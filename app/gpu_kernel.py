"""Numba engine: CUDA kernel on NVIDIA GPUs, parallel-CPU fallback otherwise."""
import numpy as np
from numba import njit, prange
from .config import GUIDE_LENGTH

try:
    from numba import cuda
    GPU_AVAILABLE = cuda.is_available()
except Exception:                       # no driver / no GPU
    cuda, GPU_AVAILABLE = None, False

L = GUIDE_LENGTH