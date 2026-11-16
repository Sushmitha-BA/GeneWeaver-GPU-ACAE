"""DNA helpers: encoding, reverse complement, mock genome generator."""
import random
import numpy as np

_LUT = np.full(256, 4, dtype=np.uint8)           # 4 = N / unknown
for _i, _b in enumerate("ACGT"):
    _LUT[ord(_b)] = _i
    _LUT[ord(_b.lower())] = _i
_DEC = "ACGTN"


def encode(seq: str) -> np.ndarray:
    """'ACGT' -> array([0,1,2,3])."""
    return _LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def decode(arr) -> str:
    return "".join(_DEC[int(i)] for i in arr)


def rc_encoded(arr: np.ndarray) -> np.ndarray:
    """Reverse complement of an encoded array (A<->T, C<->G, N stays N)."""
    return np.where(arr < 4, 3 - arr, 4).astype(np.uint8)[::-1].copy()