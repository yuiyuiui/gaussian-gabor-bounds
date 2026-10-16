"""Gaussian Gabor spectral intervals with explicit, non-certified error indicators."""

from .sections import GaussianPlan, Resolution, ResourceLimitError, choose_resolution

__version__ = "0.1.0"
__all__ = ["GaussianPlan", "Resolution", "ResourceLimitError", "choose_resolution"]
