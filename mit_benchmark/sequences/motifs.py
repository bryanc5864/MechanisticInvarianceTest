"""Promoter motif definitions and utilities."""

import random
from typing import List, Tuple

# Consensus sequences for E. coli σ70 promoters
MINUS_35_CONSENSUS = "TTGACA"
MINUS_10_CONSENSUS = "TATAAT"
MINUS_10_BROKEN = "TGTAAT"  # T→G at position 2 weakens the -10 box

# Compensatory elements
UP_PROXIMAL_CONSENSUS = "AAAAAARNR"  # R=A/G, N=any base
EXTENDED_10_CONSENSUS = "TGT"  # Extended -10 element (TGn motif)

# Position weight matrices for scoring (derived from E. coli promoter compilations)
# Values are log-odds scores relative to background (0.25 for each nucleotide)
MINUS_35_PWM = {
    0: {'T': 1.5, 'A': -0.5, 'C': -1.0, 'G': -1.0},
    1: {'T': 1.8, 'A': -1.0, 'C': -1.0, 'G': -0.5},
    2: {'G': 1.5, 'A': -0.5, 'T': -0.5, 'C': -1.0},
    3: {'A': 1.8, 'T': -0.5, 'C': -1.0, 'G': -1.0},
    4: {'C': 1.5, 'A': -0.5, 'T': -0.5, 'G': -1.0},
    5: {'A': 1.5, 'T': -0.5, 'C': -0.5, 'G': -1.0},
}

MINUS_10_PWM = {
    0: {'T': 1.8, 'A': -0.5, 'C': -1.0, 'G': -1.0},
    1: {'A': 1.8, 'T': -0.5, 'C': -1.0, 'G': -1.0},
    2: {'T': 1.5, 'A': -0.5, 'C': -0.5, 'G': -0.5},
    3: {'A': 1.8, 'T': -0.5, 'C': -1.0, 'G': -1.0},
    4: {'A': 1.5, 'T': -0.5, 'C': -0.5, 'G': -1.0},
    5: {'T': 1.8, 'A': -0.5, 'C': -1.0, 'G': -1.0},
}


def expand_iupac(pattern: str) -> List[str]:
    """Expand IUPAC ambiguity codes into all possible sequences.

    Args:
        pattern: DNA sequence with IUPAC codes (R=A/G, N=A/C/G/T, etc.)

    Returns:
        List of all possible concrete sequences
    """
    iupac_codes = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],  # Purine
        'Y': ['C', 'T'],  # Pyrimidine
        'S': ['C', 'G'],  # Strong
        'W': ['A', 'T'],  # Weak
        'K': ['G', 'T'],  # Keto
        'M': ['A', 'C'],  # Amino
        'B': ['C', 'G', 'T'],  # Not A
        'D': ['A', 'G', 'T'],  # Not C
        'H': ['A', 'C', 'T'],  # Not G
        'V': ['A', 'C', 'G'],  # Not T
        'N': ['A', 'C', 'G', 'T'],  # Any
    }

    results = ['']
    for char in pattern.upper():
        if char not in iupac_codes:
            raise ValueError(f"Unknown IUPAC code: {char}")
        new_results = []
        for seq in results:
            for base in iupac_codes[char]:
                new_results.append(seq + base)
        results = new_results

    return results


def sample_iupac(pattern: str) -> str:
    """Sample a single concrete sequence from an IUPAC pattern.

    Args:
        pattern: DNA sequence with IUPAC codes

    Returns:
        A single concrete DNA sequence
    """
    iupac_codes = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['C', 'G'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T'],
    }

    result = []
    for char in pattern.upper():
        if char not in iupac_codes:
            raise ValueError(f"Unknown IUPAC code: {char}")
        result.append(random.choice(iupac_codes[char]))

    return ''.join(result)


def score_pwm(sequence: str, pwm: dict) -> float:
    """Score a sequence against a position weight matrix.

    Args:
        sequence: DNA sequence to score
        pwm: Position weight matrix (dict of position -> {base: score})

    Returns:
        Sum of log-odds scores
    """
    if len(sequence) != len(pwm):
        raise ValueError(f"Sequence length {len(sequence)} != PWM length {len(pwm)}")

    total_score = 0.0
    for i, base in enumerate(sequence.upper()):
        if base in pwm[i]:
            total_score += pwm[i][base]
        else:
            total_score += -2.0  # Penalty for unexpected base

    return total_score


def generate_random_spacer(length: int = 17) -> str:
    """Generate a random spacer sequence.

    Args:
        length: Length of spacer (default 17bp)

    Returns:
        Random DNA sequence
    """
    return ''.join(random.choice('ACGT') for _ in range(length))


def generate_random_flank(length: int) -> str:
    """Generate random flanking sequence with slight AT bias (typical of E. coli).

    Args:
        length: Length of flanking sequence

    Returns:
        Random DNA sequence with AT bias
    """
    # E. coli has ~50% AT content, promoter regions often slightly higher
    weights = [0.28, 0.22, 0.22, 0.28]  # A, C, G, T
    return ''.join(random.choices('ACGT', weights=weights, k=length))


def mutate_minus_10(sequence: str, position: int = 1) -> str:
    """Create a broken -10 box by mutating a key position.

    The canonical mutation is T→G at position 2 (0-indexed position 1),
    which significantly weakens σ70 recognition.

    Args:
        sequence: Original -10 sequence (6bp)
        position: Position to mutate (default 1, the second T)

    Returns:
        Mutated -10 sequence
    """
    if len(sequence) != 6:
        raise ValueError(f"Expected 6bp -10 box, got {len(sequence)}")

    seq_list = list(sequence)
    original = seq_list[position]

    # Mutate to a different base (prefer G for T→G canonical mutation)
    if original == 'T':
        seq_list[position] = 'G'
    elif original == 'A':
        seq_list[position] = 'C'
    else:
        # Pick a random different base
        alternatives = [b for b in 'ACGT' if b != original]
        seq_list[position] = random.choice(alternatives)

    return ''.join(seq_list)


def get_compensation_elements() -> Tuple[str, str]:
    """Get compensatory elements (UP element and extended -10).

    Returns:
        Tuple of (UP element sequence, extended -10 sequence)
    """
    # Sample a concrete UP element from the consensus
    up_element = sample_iupac(UP_PROXIMAL_CONSENSUS)
    extended_10 = EXTENDED_10_CONSENSUS

    return up_element, extended_10


def scramble_sequence(sequence: str) -> str:
    """Scramble a sequence while preserving nucleotide composition.

    Args:
        sequence: DNA sequence to scramble

    Returns:
        Scrambled sequence with same composition
    """
    seq_list = list(sequence)
    random.shuffle(seq_list)
    return ''.join(seq_list)
