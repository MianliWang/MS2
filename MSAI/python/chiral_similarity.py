"""Backward-compatible imports for :mod:`MSAI.python.ms2.similarity`.

New code should import from ``MSAI.python.ms2.similarity`` directly.
"""

try:
    from .ms2.similarity import *  # noqa: F403
except ImportError:  # Support ``python MSAI/python/<script>.py`` callers.
    from ms2.similarity import *  # type: ignore[import-not-found] # noqa: F403
