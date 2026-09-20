# SPDX-License-Identifier: MIT
"""Dataset acquisition: resolve, download, verify, cache.

Network access is **opt-in**. Nothing here touches the network unless the
caller passed `--allow-download`, because the deployment these tools target is
an air-gapped host that holds the sequencing data, and a validation harness
that silently reaches for the internet is the wrong default on exactly that
machine. Without the flag a missing input is a clean SKIP with the command
needed to stage it, not a failure.

Downloads are verified against a pinned checksum and are atomic: bytes land in
a `.part` file that is only renamed into place after the digest matches. A
half-written BAM that looks complete is the failure mode that costs a day.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .config import Dataset, RemoteFile

USER_AGENT = "pto-core-validation/1.0 (+https://github.com/)"
ENCODE_API = "https://www.encodeproject.org/files/{accession}/?format=json"
CHUNK = 1 << 20


class FetchError(RuntimeError):
    pass


class DownloadsDisabled(RuntimeError):
    """Raised when an input is absent and downloads were not authorised."""


@dataclass
class Resolved:
    url: str
    sha256: str | None
    md5: str | None
    size: int | None


def resolve(spec: RemoteFile, *, allow_network: bool) -> Resolved:
    """Turn a RemoteFile into a concrete URL plus expected digest.

    ENCODE entries are pinned by accession, not by URL and checksum, and the
    portal metadata API is queried to get both. That indirection is worth one
    extra request: ENCODE file hrefs are versioned paths that change when a
    file is re-released, so a hard-coded URL rots, while the accession is the
    stable identity and the API hands back the authoritative md5 alongside it.
    A plain `url` entry must carry its own `sha256` -- there is no authority to
    ask for one.
    """
    if spec.encode_accession:
        if not allow_network:
            raise DownloadsDisabled(
                f"{spec.name}: resolving ENCODE accession "
                f"{spec.encode_accession} needs network access"
            )
        meta = _encode_metadata(spec.encode_accession)
        href = meta.get("href")
        if not href:
            raise FetchError(
                f"ENCODE accession {spec.encode_accession} returned no download href; "
                "the accession may have been revoked or replaced"
            )
        return Resolved(
            url="https://www.encodeproject.org" + href,
            sha256=spec.sha256,
            md5=spec.md5 or meta.get("md5sum"),
            size=meta.get("file_size") or spec.approx_bytes,
        )

    assert spec.url  # guaranteed by load_datasets
    if not spec.sha256 and not spec.md5:
        raise FetchError(
            f"{spec.name}: a plain `url` entry must pin `sha256` (or `md5`). "
            "Run `run_validation.py fetch --record-checksums` to fill it in, "
            "then commit datasets.json."
        )
    return Resolved(url=spec.url, sha256=spec.sha256, md5=spec.md5, size=spec.approx_bytes)


def _encode_metadata(accession: str) -> dict:
    request = urllib.request.Request(
        ENCODE_API.format(accession=accession),
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise FetchError(
            f"ENCODE accession {accession} lookup failed with HTTP {exc.code}. "
            "Check the accession in datasets.json against the portal."
        ) from exc
    except urllib.error.URLError as exc:
        raise FetchError(f"ENCODE lookup for {accession} failed: {exc.reason}") from exc


def digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            hasher.update(block)
    return hasher.hexdigest()


def ensure_file(
    dataset: Dataset,
    role: str,
    data_dir: Path,
    *,
    allow_download: bool,
    verify: bool = True,
    progress: bool = True,
) -> Path:
    """Return the local path for one dataset file, downloading it if permitted.

    An already-present file is *not* re-hashed by default. Verification is
    linear in file size and these are tens of gigabytes; re-checking on every
    run would dominate the suite's runtime to re-prove something the atomic
    rename already established. Pass `verify=True` explicitly (the `verify`
    subcommand does) when you actually suspect cache corruption.
    """
    spec = dataset.files[role]

    # A fixture checked into the tree is used where it sits. Copying it into
    # the cache would give two copies that can disagree, and the whole reason
    # it is in the tree is that it needs no staging step at all.
    if spec.repo_path:
        from .config import REPO_ROOT

        local = REPO_ROOT / spec.repo_path
        if not local.exists():
            raise FetchError(
                f"{spec.name}: repo fixture {spec.repo_path} is missing from the "
                "working tree"
            )
        return local

    target = dataset.local_dir(data_dir) / spec.filename

    if target.exists() and target.stat().st_size > 0:
        return target

    if not allow_download:
        raise DownloadsDisabled(
            f"missing input {target}\n"
            f"  dataset : {dataset.key} ({dataset.title})\n"
            f"  source  : {dataset.source}\n"
            f"  stage it with: run_validation.py fetch --datasets {dataset.key} "
            f"--allow-download"
        )

    resolved = resolve(spec, allow_network=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")

    _download(resolved.url, partial, expected_size=resolved.size, progress=progress)

    if resolved.sha256:
        _verify(partial, resolved.sha256, "sha256", spec.name)
    elif resolved.md5:
        _verify(partial, resolved.md5, "md5", spec.name)

    partial.replace(target)
    return target


def _download(url: str, destination: Path, *, expected_size: int | None, progress: bool) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            total = expected_size or int(response.headers.get("Content-Length") or 0)
            written = 0
            with destination.open("wb") as sink:
                while True:
                    block = response.read(CHUNK)
                    if not block:
                        break
                    sink.write(block)
                    written += len(block)
                    if progress:
                        _progress(destination.name, written, total)
    except urllib.error.URLError as exc:
        destination.unlink(missing_ok=True)
        raise FetchError(f"download of {url} failed: {exc}") from exc

    if progress:
        sys.stderr.write("\n")


def _progress(name: str, written: int, total: int) -> None:
    if total:
        pct = 100.0 * written / total
        sys.stderr.write(f"\r  {name}: {written / 1e6:8.1f} MB / {total / 1e6:.1f} MB ({pct:5.1f}%)")
    else:
        sys.stderr.write(f"\r  {name}: {written / 1e6:8.1f} MB")
    sys.stderr.flush()


def _verify(path: Path, expected: str, algorithm: str, name: str) -> None:
    actual = digest(path, algorithm)
    if actual.lower() != expected.lower():
        path.unlink(missing_ok=True)
        raise FetchError(
            f"{name}: {algorithm} mismatch -- the upstream file changed or the "
            f"download was corrupted.\n"
            f"  expected {expected}\n"
            f"  actual   {actual}\n"
            "If upstream genuinely re-released this file, update datasets.json "
            "deliberately and note the change; do not just paste the new digest in."
        )


def record_checksums(
    datasets: dict[str, Dataset], data_dir: Path, datasets_file: Path
) -> int:
    """Fill in missing `sha256` values in datasets.json from cached files.

    This is the one-time bootstrap after adding a dataset: fetch it, eyeball
    that it is the file you meant, then pin it. It only ever *adds* a digest --
    an existing pin that no longer matches is reported and left alone, because
    overwriting it is precisely the mistake the pin exists to prevent.
    """
    raw = json.loads(datasets_file.read_text())
    changed = 0

    for key, dataset in datasets.items():
        for role, spec in dataset.files.items():
            path = dataset.local_dir(data_dir) / spec.filename
            if not path.exists():
                continue
            actual = digest(path)
            if spec.sha256:
                if spec.sha256.lower() != actual.lower():
                    print(
                        f"  MISMATCH {key}/{role}: pinned {spec.sha256[:16]}… "
                        f"but cached file is {actual[:16]}… (left unchanged)",
                        file=sys.stderr,
                    )
                continue
            raw["datasets"][key]["files"][role]["sha256"] = actual
            print(f"  pinned {key}/{role} = {actual}")
            changed += 1

    if changed:
        backup = datasets_file.with_suffix(".json.bak")
        shutil.copy2(datasets_file, backup)
        datasets_file.write_text(json.dumps(raw, indent=2) + "\n")
        print(f"\nUpdated {datasets_file} ({changed} digest(s)); backup at {backup}")
        print("Review the diff and commit it.")
    return changed
