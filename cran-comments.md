# CRAN submission comments — causalsim 0.1.0

## Test environments

- Local: macOS 26.1 (aarch64-apple-darwin), R 4.5.2
- win-builder: Windows Server 2022 x64, R-devel (2026-08-17 r90424 ucrt)

## R CMD check results

Local (`R CMD check --as-cran`, incl. `--run-donttest`): 0 errors | 0 warnings | 0 notes

win-builder (R-devel): 0 errors | 0 warnings | 1 note

The one NOTE is:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Chayce Reed <Chayce.Reed.HSE@dartmouth.edu>'

New submission

Possibly misspelled words in DESCRIPTION:
  DGPs (7:6)
  confounder (9:54)
```

Both parts are expected:

- "New submission" is expected, as this is the initial submission of the package.
- "DGPs" and "confounder" are not misspellings. "DGPs" is the plural of DGP
  (data generating process), used throughout the package, and "confounder" is
  standard causal-inference terminology.

## Notes to reviewer

This is an initial submission.

Examples for `causalsim_eval()` and `causalsim_grid()` are wrapped in
`\donttest{}` because they run Monte Carlo simulation loops that exceed the
5-second example time limit on CRAN check machines. Both functions are fully
covered by the test suite.
