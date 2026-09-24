# SPDX-License-Identifier: AGPL-3.0-only
"""Small helpers shared by SMEW modules."""


def require_backend(backend):
    """Validate execution choice and check that the native solver is installed."""
    if backend not in ("python", "compiled"):
        raise ValueError(
            f"Unknown backend {backend!r}; choose 'python' or 'compiled'."
        )
    if backend == "compiled":
        try:
            from smew._native import _minpack
        except ImportError as exc:
            raise RuntimeError(
                "The compiled backend is unavailable. Install a native wheel, "
                "build the extensions with 'make install-native', or select "
                "backend='python'."
            ) from exc
        return _minpack
    return None
