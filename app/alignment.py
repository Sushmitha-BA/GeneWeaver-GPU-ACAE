
"""Basic CRISPR alignment with NGG PAM detection."""


def reverse_complement(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def find_matches(genome, guide, max_mismatches=3):
    genome = genome.upper()
    guide = guide.upper()

    if not guide or any(b not in "ACGT" for b in guide):
        raise ValueError("Guide must contain only A, C, G and T")

    if max_mismatches < 0:
        raise ValueError("max_mismatches cannot be negative")

    results = []
    size = len(guide)

    for i in range(len(genome) - size - 2 + 1):
        target = genome[i:i + size]
        pam = genome[i + size:i + size + 3]

        # Forward strand: NGG PAM after target
        if pam[1:] == "GG":
            mismatches = sum(a != b for a, b in zip(guide, target))

            if mismatches <= max_mismatches:
                results.append({
                    "position": i,
                    "strand": "+",
                    "target": target,
                    "pam": pam,
                    "mismatches": mismatches
                })

        # Reverse strand: CCN PAM before target
        reverse_pam = genome[i:i + 3]
        if reverse_pam[:2] == "CC":
            reverse_target = genome[i + 3:i + 3 + size]

            if len(reverse_target) == size:
                matched_target = reverse_complement(reverse_target)
                mismatches = sum(
                    a != b for a, b in zip(guide, matched_target)
                )

                if mismatches <= max_mismatches:
                    results.append({
                        "position": i + 3,
                        "strand": "-",
                        "target": matched_target,
                        "pam": reverse_complement(reverse_pam),
                        "mismatches": mismatches
                    })

    return results
