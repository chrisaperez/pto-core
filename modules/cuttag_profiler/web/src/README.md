# Frontend source

`web/index.html` is the dashboard, and it is deliberately a single
dependency-free file rather than a compiled React/Tailwind bundle.

## Why not React

The original design called for a pre-compiled React + Tailwind bundle embedded
as byte arrays. Building it that way would have meant that anyone compiling the
project needed Node.js, a lockfile-pinned `node_modules`, and a working
registry mirror — on the same machine whose selling point is that it has no
outbound network. It also means a reviewer cannot read the shipped UI: they
would be auditing minified output.

The dashboard is one form, three canvases and a table. Hand-written HTML, CSS
and ES2020 covers that in ~700 readable lines with:

- **no build step** — `cmake --build` is the whole toolchain;
- **no third-party JavaScript** in the binary, so the CSP can be
  `default-src 'self'` with nothing to whitelist;
- **an auditable artefact** — the exact text served to the browser is the text
  in the repository.

Charts are drawn on `<canvas>` rather than with a plotting library. The
meta-profile is a filled line chart; the heatmap is rasterised at its natural
resolution (one pixel per bin per band) onto an offscreen canvas and scaled
with `drawImage`, which keeps it correct on HiDPI displays. Note that
`putImageData` is *not* usable for this: it writes raw device pixels and
ignores the context transform, so the image lands in the top-left quadrant at
half size on a 2x screen.

## If you do want a React build

Nothing in the C++ layer cares. `cmake/embed_assets.cmake` globs
`web/*.{html,css,js,svg,png,ico,woff2,json}` and skips `web/src/` entirely, so:

1. put the React sources here in `web/src/`;
2. configure the bundler to emit into `web/` (`index.html` plus hashed assets);
3. rebuild — the generated `embedded_assets.cpp` picks the new files up, and
   every one of them is served from the binary as before.

Keep the output self-contained. Any `<script src="https://…">` or webfont URL
that survives bundling will be blocked by the Content-Security-Policy set in
`src/http_server.cpp`, and would break the air-gap guarantee that justifies
this tool's existence.

## API contract

The dashboard talks to three endpoints, all documented in the main README:

| Endpoint | Method | Purpose |
| --- | --- | --- |
| `/api/health` | GET | version and core count |
| `/api/chroms?bam=…` | GET | header contigs, mapped-read total |
| `/api/profile` | POST | compute a matrix; returns profile, heatmap, stats |
| `/api/matrix?token=…` | GET | full matrix as TSV, keyed by a run token |
