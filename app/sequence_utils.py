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


def make_mock_genome(path, length=1_000_000, guide=None, plants=20, seed=42):
    """Write a random FASTA genome and plant near-matches of `guide`
    (0-3 mismatches, followed by an NGG PAM) so a scan always finds hits."""
    rng = random.Random(seed)
    seq = [rng.choice("ACGT") for _ in range(length)]
    if guide:
        for _ in range(plants):
            site = list(guide)
            for p in rng.sample(range(len(guide)), rng.randint(0, 3)):
                site[p] = rng.choice([b for b in "ACGT" if b != site[p]])
            site += [rng.choice("ACGT"), "G", "G"]
            pos = rng.randint(0, length - len(site) - 1)
            seq[pos:pos + len(site)] = site
    with open(path, "w") as fh:
        fh.write(">mock_chr1 GeneWeaver mock genome\n")
        s = "".join(seq)
        for i in range(0, len(s), 80):
            fh.write(s[i:i + 80] + "\n")
    return path