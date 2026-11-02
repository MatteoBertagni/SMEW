Release packaging and publication
=================================

The **Build and publish distributions** GitHub Actions workflow is defined in
``.github/workflows/python-publish.yml``. It builds one source archive and uses
cibuildwheel to build all platform wheels from that archive.

Packaging checks
----------------

Every wheel is installed and tested outside the source tree, including
Numba-to-cminpack calls and short Python/compiled simulation comparisons. Native
extensions and license notices must be present; missing extensions fail the job.
The separate **Tests** workflow runs the full model regressions.

Source archives include the Cython and cminpack sources, the wheel test suite,
and the documentation sources and example notebook. They exclude local binaries,
Numba caches, generated documentation, and large regression baseline files.

Preview builds
--------------

Pull requests and manual workflow runs produce downloadable source and wheel
artifacts without publishing to PyPI. Use these runs to check packaging before
releasing.

Publish a release
-----------------

1. Set a new version in ``smew/_version.py`` and commit the change before
   creating the release tag. PyPI versions cannot be overwritten.
2. Publish a GitHub release for that tag. This triggers the same source and
   wheel builds used for preview runs.
3. Check the **Build and publish distributions** workflow. It uploads the
   source archive and wheels to PyPI only after every platform passes.
