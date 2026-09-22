# Deploying `cuttag_profiler`'s dashboard

This directory stands up `cuttag_profiler serve` — the embedded, single-binary
HTTP dashboard — on a machine you control, for one operator to inspect their
own CUT&Tag/CUT&RUN signal matrices. It is **not** meant for a managed,
multi-tenant deployment: a build for that shape sets `PTO_CLOUD_BUILD` (see
`CMakeLists.txt`), which compiles the server out of the binary entirely.

Two ways to run it:

- **Locally**: build `cuttag_profiler` and run `cuttag_profiler serve
  --data-root <dir>` directly. No container, no auth flags needed — it binds
  loopback by default and refuses to run confinement-disabled.
- **On your own cloud instance**, when the data has to live closer to where it
  was generated than your laptop: `provision_aws.sh` + `Dockerfile` +
  `docker-compose.yml` here do that. The rest of this file is about that path.

```bash
git clone <repo> pto-core && cd pto-core
sudo ./deploy/provision_aws.sh
```

`provision_aws.sh` is idempotent and can also be pasted as EC2 user-data. It
installs Docker, builds the image, brings the service up, and prints the exact
SSH tunnel command for your instance at the end.

## Security posture — read this before running it

**The dashboard is single-operator, not multi-tenant, and was never designed
to be reachable from the open internet.** That remains true even though it is
not literally unauthenticated:

- **It does require a token.** Since `SECURITY_HTTP_2026-08-15` finding H1c
  closed, every server mints a random 256-bit session token per run and
  requires it (as `?t=` or `Authorization: Bearer`) on every `/api/` request.
  `--allow-remote` — which this deployment always passes, because binding
  `0.0.0.0` inside the container is unavoidable if the port is to be published
  at all — additionally *refuses to start* combined with `--no-auth` or with
  confinement disabled. There is no supported way to run this image with
  authentication off.
- **But it is one shared secret for the run, not a login system.** There is no
  per-user identity, no rotation, no revocation, and no rate limiting on the
  token check beyond what the transport gives you. Anyone who can read the
  operator's terminal, shell history, or container logs has it — that is the
  design's own stated bound (see `docs/SECURITY_HTTP_2026-08-15.md`, H1c).
  That is an appropriate bar for "reachable only over an SSH tunnel I control"
  and not an appropriate bar for "reachable from the public internet."
- **`--data-root` confinement (finding H2) is closed and mandatory.** The
  binary refuses to start without it, and the entrypoint here still passes it
  explicitly so confinement is visible in the image rather than inferred from
  a default. A finding that is *not* fully closed: `M6` (arbitrary file read
  via a crafted region-file path) is closed for BED/GTF region files but only
  partly closed for BAM paths, because `sam_index_load` needs a real
  filesystem path. Keep the data root scoped to only what this run needs to
  see.

**Given that, port 8080 must never be bound to a public interface.** The
supported access pattern is an SSH tunnel:

```bash
ssh -N -L 8080:127.0.0.1:8080 ubuntu@<instance>
# then open http://127.0.0.1:8080/ on your laptop
```

`docker-compose.yml` publishes `127.0.0.1:8080:8080` — loopback on the host
side — which is the actual security boundary here, not the container's own
`PROFILER_BIND=0.0.0.0` (which only matters inside the container's network
namespace) and not `ufw`/security-group rules alone: Docker installs its DNAT
rules ahead of `ufw`'s `INPUT` chain, so a port published on every interface
bypasses a host firewall rule you might otherwise trust. `provision_aws.sh`
configures `ufw` to allow only SSH and deliberately does not open 8080.

Your EC2 security group should have exactly one inbound rule:

```
Type: SSH   Protocol: TCP   Port: 22   Source: <your.ip.address>/32
```

Do not add an inbound rule for 8080. If you genuinely need direct exposure —
multiple people need dashboard access without an SSH key each — put an
authenticating reverse proxy in front (ALB + Cognito/OIDC, or nginx with mTLS)
and restrict the source range; do not rely on the run token as internet-facing
auth.

## Getting the URL (and its token) out of a detached container

`docker compose up -d` runs the service in the background, so you will not see
the startup banner directly. The server prints the dashboard URL — including
`?t=<token>` — once, to stdout, at startup:

```bash
cd deploy
docker compose logs cuttag-profiler | grep 'dashboard ready'
```

Open that URL through the SSH tunnel above, not by publishing the container's
port. The token is stripped from the address bar by the dashboard's own JS
after the first load, so it will not end up in your browser history from that
point on — but it will remain in the container logs for the life of the run.

## What the image already hardens, so you don't have to

From `Dockerfile` and `docker-compose.yml`: the process runs as an
unprivileged user with every capability dropped, the root filesystem is
read-only with a size-capped `tmpfs` for `/tmp`, `no-new-privileges` is set,
the data volume is mounted read-only, and container resource limits (4 GiB
memory, 2 CPUs) bound the blast radius of a malformed-input OOM on top of the
application's own request-size and bin-count ceilings. None of this is a
substitute for the network boundary above — it limits what a bug in the
BAM/BED parser can do once reached, not who can reach it.
