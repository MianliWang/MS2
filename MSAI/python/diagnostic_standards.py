"""Backward-compatible imports for :mod:`MSAI.python.ms2.diagnostics`."""

try:
    from .ms2.diagnostics import *  # noqa: F401,F403
except ImportError:
    from ms2.diagnostics import *  # type: ignore # noqa: F401,F403
