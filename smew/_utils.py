# SPDX-License-Identifier: AGPL-3.0-only
"""Small helpers shared by SMEW modules."""


def require_backend(backend):
    """Validate backend selection until compiled execution is available."""
    if backend not in ("python", "compiled"):
        raise ValueError(
            f"Unknown backend {backend!r}; choose 'python' or 'compiled'."
        )
    if backend == "compiled":
        raise RuntimeError(
            "The compiled backend has not been implemented yet. "
            "Use backend='python' for now."
        )
