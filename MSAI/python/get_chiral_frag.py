"""Legacy CLI facade for the canonical MS2 analysis pipeline.

The implementation now lives in :mod:`MSAI.python.ms2.pipeline`.
"""

try:
    from .ms2.pipeline import *  # noqa: F403
    from .ms2.pipeline import main
except ImportError:
    from ms2.pipeline import *  # type: ignore[import-not-found] # noqa: F403
    from ms2.pipeline import main  # type: ignore[import-not-found]


if __name__ == "__main__":
    main()
