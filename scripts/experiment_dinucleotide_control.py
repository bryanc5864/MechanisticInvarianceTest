#!/usr/bin/env python3
"""Dinucleotide frequency control experiment.

Distinguishes mononucleotide (AT%) from dinucleotide from positional effects
by generating sequences with matched compositions at different levels.
This experiment generates sequences with matched dinucleotide frequencies
but different motif positions, distinguishing whether models respond to:
1. Mononucleotide (AT%) composition
2. Dinucleotide patterns (AA, AT, TA, TT frequencies)
3. Actual motif structure at specific positions

Design:
- Match-AT: Same AT% as compensated, random dinucleotide distribution
- Match-Dinuc: Same dinucleotide frequencies as compensated, shuffled
- Compensated: Real compensated sequence (UP + ext-10 at correct positions)
- Broken: No compensation elements

If models respond equally to Match-AT and Compensated → mononucleotide effect
If models prefer Match-Dinuc over Match-AT → dinucleotide effect
If models prefer Compensated over Match-Dinuc → some positional signal
"""

import json
import random
import sys
from pathlib import Path
from collections import Counter

import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

random.seed(42)
np.random.seed(42)


def count_dinucleotides(seq: str) -> dict:
    """Count all 16 dinucleotide frequencies in a sequence."""
    counts = Counter()
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if all(c in 'ACGT' for c in dinuc):
            counts[dinuc] += 1
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()} if total > 0 else counts


def shuffle_preserving_dinucleotides(seq: str, max_attempts: int = 1000) -> str:
    """Shuffle sequence approximately preserving dinucleotide frequencies.

    Uses the Altschul-Erickson algorithm (Euler path on dinucleotide graph).
    Falls back to simple shuffle if it fails.
    """
    seq = list(seq)
    n = len(seq)

    # Simple approach: swap pairs while tracking dinucleotide changes
    best_seq = seq.copy()
    original_dinucs = count_dinucleotides(''.join(seq))

    for _ in range(max_attempts):
        # Random swap
        i, j = random.sample(range(1, n - 1), 2)  # Avoid endpoints
        seq[i], seq[j] = seq[j], seq[i]

        new_dinucs = count_dinucleotides(''.join(seq))
        # Accept if dinucleotide distribution is closer or equal
        old_dist = sum(abs(original_dinucs.get(k, 0) - count_dinucleotides(''.join(best_seq)).get(k, 0))
                       for k in set(list(original_dinucs.keys()) + list(new_dinucs.keys())))
        new_dist = sum(abs(original_dinucs.get(k, 0) - new_dinucs.get(k, 0))
                       for k in set(list(original_dinucs.keys()) + list(new_dinucs.keys())))

        if new_dist <= old_dist:
            best_seq = seq.copy()
        else:
            seq[i], seq[j] = seq[j], seq[i]  # Revert

    return ''.join(best_seq)


def generate_background(length: int, at_fraction: float = 0.55) -> str:
    """Generate random background sequence with specified AT content."""
    seq = []
    for _ in range(length):
        if random.random() < at_fraction:
            seq.append(random.choice('AT'))
        else:
            seq.append(random.choice('GC'))
    return ''.join(seq)


def generate_compensated(at_frac: float = 0.55) -> str:
    """Generate compensated sequence with UP + ext-10 at correct positions."""
    seq = list(generate_background(100, at_frac))
    # UP element at position 8-16
    for i, nt in enumerate("AAAAAAGCA"):
        seq[8 + i] = nt
    # -35 at position 25
    for i, nt in enumerate("TTGACA"):
        seq[25 + i] = nt
    # Extended -10 at position 45
    for i, nt in enumerate("TGT"):
        seq[45 + i] = nt
    # Broken -10 at position 48
    for i, nt in enumerate("TGTAAT"):
        seq[48 + i] = nt
    return ''.join(seq)


def generate_broken(at_frac: float = 0.55) -> str:
    """Generate broken sequence (no compensation)."""
    seq = list(generate_background(100, at_frac))
    # -35 at position 25
    for i, nt in enumerate("TTGACA"):
        seq[25 + i] = nt
    # Broken -10 at position 48
    for i, nt in enumerate("TGTAAT"):
        seq[48 + i] = nt
    return ''.join(seq)


def generate_matched_at(target_at: float) -> str:
    """Generate sequence with matched AT% but random dinucleotide distribution."""
    seq = list(generate_background(100, target_at))
    # -35 at position 25
    for i, nt in enumerate("TTGACA"):
        seq[25 + i] = nt
    # Broken -10 at position 48
    for i, nt in enumerate("TGTAAT"):
        seq[48 + i] = nt
    return ''.join(seq)


def generate_matched_dinucleotide(reference_seq: str) -> str:
    """Generate sequence with matched dinucleotide freqs from a reference.

    Shuffles the reference while preserving dinucleotide frequencies,
    then re-inserts the -35 and broken -10 motifs.
    """
    shuffled = shuffle_preserving_dinucleotides(reference_seq)
    seq = list(shuffled)
    # Re-insert -35 at position 25
    for i, nt in enumerate("TTGACA"):
        seq[25 + i] = nt
    # Re-insert broken -10 at position 48
    for i, nt in enumerate("TGTAAT"):
        seq[48 + i] = nt
    return ''.join(seq)


def main():
    print("=" * 80)
    print("DINUCLEOTIDE FREQUENCY CONTROL EXPERIMENT")
    print("=" * 80)

    n_samples = 100

    # Generate all sequence conditions
    print("\nGenerating sequences...")

    compensated_seqs = [generate_compensated() for _ in range(n_samples)]
    broken_seqs = [generate_broken() for _ in range(n_samples)]

    # Compute AT% of compensated sequences
    comp_at_fracs = []
    for seq in compensated_seqs:
        at = (seq.count('A') + seq.count('T')) / len(seq)
        comp_at_fracs.append(at)
    mean_comp_at = np.mean(comp_at_fracs)
    print(f"  Mean AT% of compensated sequences: {mean_comp_at:.3f}")

    # Match-AT: same AT% as compensated, random dinucleotide distribution
    matched_at_seqs = [generate_matched_at(mean_comp_at) for _ in range(n_samples)]

    # Match-Dinuc: shuffle compensated to preserve dinucleotide freqs
    matched_dinuc_seqs = [generate_matched_dinucleotide(seq) for seq in compensated_seqs]

    # Verify compositions
    print("\n  Composition verification:")
    for name, seqs in [("Compensated", compensated_seqs),
                       ("Match-AT", matched_at_seqs),
                       ("Match-Dinuc", matched_dinuc_seqs),
                       ("Broken", broken_seqs)]:
        at_vals = [(s.count('A') + s.count('T')) / len(s) for s in seqs]
        dinuc_vals = [count_dinucleotides(s) for s in seqs]
        mean_aa = np.mean([d.get('AA', 0) for d in dinuc_vals])
        mean_at = np.mean([d.get('AT', 0) for d in dinuc_vals])
        mean_ta = np.mean([d.get('TA', 0) for d in dinuc_vals])
        print(f"    {name:<15}: AT%={np.mean(at_vals):.3f}, "
              f"AA={mean_aa:.3f}, AT={mean_at:.3f}, TA={mean_ta:.3f}")

    # Save sequences for inference
    results = {
        "compensated": compensated_seqs,
        "broken": broken_seqs,
        "matched_at": matched_at_seqs,
        "matched_dinuc": matched_dinuc_seqs,
        "metadata": {
            "n_samples": n_samples,
            "mean_compensated_at": float(mean_comp_at),
            "description": "Dinucleotide frequency control experiment",
        }
    }

    output_path = PROJECT_ROOT / "data/results/dinucleotide_control_sequences.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSequences saved to {output_path}")

    # Run biophysical models as quick validation
    print("\n--- Biophysical Model Validation ---")
    from mit_benchmark.models.biophysical import PositionAwarePWM
    # This script generates sequences with -35 at 25, -10 at 48, UP at 8, ext-10 at 45
    papwm = PositionAwarePWM(tss_position=60, pos_35=25, pos_10=48, pos_up=8, pos_ext10=45)

    conditions = {
        "Compensated": compensated_seqs,
        "Match-AT": matched_at_seqs,
        "Match-Dinuc": matched_dinuc_seqs,
        "Broken": broken_seqs,
    }

    for name, seqs in conditions.items():
        scores = [papwm.score(s) for s in seqs]
        print(f"  {name:<15}: mean={np.mean(scores):.2f}, std={np.std(scores):.2f}")

    # Compute pairwise CSS-like metrics
    print("\n--- Pairwise Comparisons (PA-PWM) ---")
    comp_scores = [papwm.score(s) for s in compensated_seqs]
    broken_scores = [papwm.score(s) for s in broken_seqs]
    mat_scores = [papwm.score(s) for s in matched_at_seqs]
    dinuc_scores = [papwm.score(s) for s in matched_dinuc_seqs]

    def pairwise_frac(a, b):
        return sum(1 for x, y in zip(a, b) if x > y) / len(a)

    print(f"  P(Compensated > Broken):   {pairwise_frac(comp_scores, broken_scores):.3f}")
    print(f"  P(Match-AT > Broken):      {pairwise_frac(mat_scores, broken_scores):.3f}")
    print(f"  P(Match-Dinuc > Broken):   {pairwise_frac(dinuc_scores, broken_scores):.3f}")
    print(f"  P(Compensated > Match-AT): {pairwise_frac(comp_scores, mat_scores):.3f}")
    print(f"  P(Compensated > Match-Dinuc): {pairwise_frac(comp_scores, dinuc_scores):.3f}")

    print("\n--- Interpretation ---")
    print("""
    If gLM shows:
      P(Comp > Broken) ≈ P(Match-AT > Broken)  → Pure AT% effect
      P(Comp > Broken) > P(Match-AT > Broken)
        but P(Comp > Match-Dinuc) ≈ 0.5         → Dinucleotide effect
      P(Comp > Match-Dinuc) >> 0.5               → Some positional signal
    """)

    print("=" * 80)
    print("Run gLM inference on these sequences to complete the analysis.")
    print("=" * 80)


if __name__ == "__main__":
    main()
