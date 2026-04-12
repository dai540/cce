# Design guide

## Why the package is small

This package was rebuilt to minimize size on disk and dependency
surface.

The design rules are:

- use only source files on `main`
- keep public functions few and explicit
- avoid large bundled data
- avoid external tutorial downloads
- make pkgdown reproducible from source only

## Output contract

Both workflows return the same components:

- curves
- effects
- diagnostics
- metadata

That common contract matters more than algorithmic breadth. The package
is optimized for lightweight reuse.

## Workflow split

The package separates two use cases:

- [`fit_cce_vs()`](https://dai540.github.io/cce/reference/fit_cce_vs.md)
  for observed `A` versus `SOC`
- [`project_soc_only()`](https://dai540.github.io/cce/reference/project_soc_only.md)
  for assumption-based planning from `SOC`

This keeps the implementation simple while still covering the main
planning questions addressed by the package.

## Tutorial data policy

Tutorials use only:

- synthetic data from
  [`cce_demo_data()`](https://dai540.github.io/cce/reference/cce_demo_data.md)
- [`survival::veteran`](https://rdrr.io/pkg/survival/man/veteran.html),
  which is already distributed with `survival`

No large external data is downloaded by the package or by the
documentation.
