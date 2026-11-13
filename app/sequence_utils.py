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