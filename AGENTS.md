# Agent notes

- Use `pixi run` for Python and project tooling.
- Keep implementations minimal and direct.
- Preserve the batch presets, sorted mesh order, per-attempt seeds, and sampling distributions.
- Run `pixi run test` after code changes.
- Build Elmer with `pixi run install-elmer`; do not depend on sibling repositories.
- Batch presets launch long simulations against external mesh datasets. Do not launch full batches just to test changes.
