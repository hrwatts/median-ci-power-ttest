# Contributing

This repository is primarily a research compendium for a manuscript rather than a general-purpose software package, so the most helpful contributions are usually:

- spotting mistakes in the statistical argument or interpretation
- identifying unclear wording in the manuscript or overview docs
- reporting reproducibility problems
- suggesting cleaner numerical experiments or presentation improvements

## Before Opening A Change

- Background reading: [README.md](README.md) and [docs/overview.md](docs/overview.md).
- Suggestions that would materially change the manuscript's claims should begin with an issue or discussion before a large patch is prepared.
- Keep the repo focused. Small, well-scoped improvements are much easier to review than broad refactors.

## Reproducibility Expectations

- Preserve the ability to rebuild the paper from source.
- Do not commit transient LaTeX or Python cache files.
- Changes to numerical results should be accompanied by regenerated outputs in `results/` and `figures/`.
- Keep documentation changes aligned with the actual implementation and manuscript wording.

## Style

- Prefer clear names and small functions over clever shortcuts.
- Keep statistical claims precise and appropriately qualified.
- Prefer readability for a research audience over compactness.
