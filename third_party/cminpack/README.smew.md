# Bundled cminpack source

SMEW vendors the double-precision hybrid nonlinear solver from
[devernay/cminpack](https://github.com/devernay/cminpack), tag `v1.3.14`.

- Upstream commit: `48c2b6ecd180ad134c626365e3092ba5dd5463a7`
- Imported: 2026-09-23
- Local modifications to copied upstream files: none
- Build policy: normal SMEW builds use these committed files; they do not
  download cminpack.

The copied source files are `hybrd.c`, `hybrd1.c`, `dogleg.c`, `dpmpar.c`,
`enorm.c`, `fdjac1.c`, `qform.c`, `qrfac.c`, `r1mpyq.c`, and `r1updt.c`.
The copied headers are `cminpack.h`, `cminpackP.h`, and `minpackP.h`.
`CopyrightMINPACK.txt`, `README.md`, and `readme.txt` are copied upstream
notices/documentation. The first nine C files are the dependency closure for
`hybrd`; `hybrd1.c` is included for its convenience wrapper. No other
cminpack algorithms or optional BLAS/LAPACK implementations are included.

Files from upstream `src/` and `include/` are placed directly in this directory
without changing their contents. Their SHA-256 hashes are in `SHA256SUMS`.
To update, select an explicit upstream tag and full commit, inspect the hybrid
dependency closure and notices, copy the same explicit file list, review the
upstream diff, and regenerate `SHA256SUMS`. Keep SMEW integration code outside
the copied files.
