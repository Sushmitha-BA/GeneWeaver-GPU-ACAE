"""Pure-Python CPU baseline (Week 1) - deliberately simple and slow."""
from . config import GUIDE_LENGTH


def scan_cpu(chunk, guide, max_mm):
    """Return [(index, mismatches)] for every site with an NGG PAM and
    <= max_mm mismatches against the guide. chunk/guide are uint8 arrays."""
    c, g = chunk.tolist(), guide.tolist()
    L = GUIDE_LENGTH
    hits = []
    for i in range(len(c) - L - 2):
        if c[i + L + 1] != 2 or c[i + L + 2] != 2:   # PAM = N G G (G = 2)
            continue
        mm = 0
        for j in range(L):
            if c[i + j] != g[j]:
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm:
            hits.append((i, mm))
    return hits