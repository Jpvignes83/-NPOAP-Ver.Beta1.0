"""Core package for the photometric analysis tool.

This package exposes the main classes and functions used throughout the
application:

- ``ImageProcessor`` for CCD reduction and calibration
- ``CalculateMaster`` to build master calibration frames
- ``process_alignment`` and ``parallel_process_images`` for image alignment
- ``PerformPlateSolving`` for astrometric solving
- ``JD_BJD_Calculator`` for astronomical time calculations
- ``process_photometry`` for basic photometric analysis
- ``validate_fits_header`` to check FITS header coordinates
"""

from .image_processor import ImageProcessor
from .calculate_master import CalculateMaster
from .align_images import process_alignment, parallel_process_images
#from .jd_bjd_calculator import JD_BJD_Calculator
from .header_manager import validate_fits_header

# Optional components that require extra dependencies
try:  # pragma: no cover - optional dependency
    from .perform_plate_solving import PerformPlateSolving
except Exception:  # ImportError or runtime error
    PerformPlateSolving = None  # type: ignore

try:  # pragma: no cover - optional dependency
    from .photometry_processing import process_photometry
except Exception:  # ImportError or runtime error
    process_photometry = None  # type: ignore

__all__ = [
    "ImageProcessor",
    "CalculateMaster",
    "process_alignment",
    "parallel_process_images",
    "PerformPlateSolving",
    "JD_BJD_Calculator",
    "process_photometry",
    "validate_fits_header",
]
