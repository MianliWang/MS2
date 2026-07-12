"""Legacy CLI facade for :mod:`MSAI.python.ms2_review.report`."""

try:
    from .ms2_review.report import *  # noqa: F401,F403
    from .ms2_review.report import main
except ImportError:
    from ms2_review.report import *  # type: ignore # noqa: F401,F403
    from ms2_review.report import main  # type: ignore


if __name__ == "__main__":
    main()
