"""Legacy CLI facade for the canonical MS2 analysis pipeline.

The implementation now lives in :mod:`MSAI.python.ms2.pipeline`.
"""

try:
    from .ms2.pipeline import *  # noqa: F401,F403
    from .ms2.pipeline import main
except ImportError:
    from ms2.pipeline import *  # type: ignore # noqa: F401,F403
    from ms2.pipeline import main  # type: ignore


if __name__ == "__main__":
    main()
