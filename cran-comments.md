# CRAN submission comments — causalsim 0.1.0

## Test environments

- macOS (local), R 4.5.2
- GitHub Actions: ubuntu-latest, R release

## R CMD CHECK results

0 errors | 0 warnings | 0 notes

## Notes to reviewer

This is an initial submission.

Examples for `causalsim_eval()` and `causalsim_grid()` are wrapped in
`\donttest{}` because they run Monte Carlo simulation loops that exceed the
5-second example time limit on CRAN check machines. The functions are fully
covered by the test suite.
