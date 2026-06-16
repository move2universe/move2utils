# Contributing to move2utils

`move2utils` is developed on **MPCDF GitLab**
(<https://gitlab.mpcdf.mpg.de/anenvi/r-packages/move2utils>); this
GitHub repository is its public release. You are very welcome to:

- **Report bugs / request features** by opening an issue here, ideally with a
  minimal reproducible example on bundled data (`inst/extdata/` or the
  `move2` example data).
- **Propose changes** via a pull request. Because development happens on
  GitLab, accepted PRs are ported there and flow back on the next release,
  with authorship preserved.

`move2utils` is intentionally minimal in scope: utilities that operate
directly on `move2` objects. Proposals outside that scope are likely better
placed in their own package.

New functionality should come with tests and `roxygen2` docs; please keep
`R CMD check` green.
