"""Backward-compatible imports for :mod:`MSAI.python.ms2.diagnostics`."""

try:
    from .ms2.diagnostics import *  # noqa: F403
except ImportError:
    from ms2.diagnostics import *  # type: ignore[import-not-found] # noqa: F403
