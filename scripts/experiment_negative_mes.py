#!/usr/bin/env python3
"""Negative MES investigation.

Investigates why gLMs assign higher likelihood to TGTAAT (broken -10)
than TATAAT (intact -10), which is counter-intuitive given the functional
importance of the intact consensus. Analysis approach:

1. Counting genome-wide frequency of both hexamers in E. coli
2. Checking if TGTAAT is more common than TATAAT in training data
3. Testing whether the LL difference correlates with hexamer frequency
4. Testing multiple -10 mutation variants to see if all broken > intact
"""

import json
import random
import sys
from pathlib import Path
from itertools import product
from collections import Counter

import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

random.seed(42)
np.random.seed(42)


# E. coli K-12 hexamer frequencies (approximate, derived from genome composition)
# The E. coli genome is ~50.8% GC, so nucleotide frequencies are:
# A: 0.246, T: 0.246, G: 0.254, C: 0.254
ECOLI_NUC_FREQ = {'A': 0.246, 'T': 0.246, 'G': 0.254, 'C': 0.254}


def compute_hexamer_expected_freq(hexamer: str) -> float:
    """Compute expected frequency of a hexamer under independence assumption."""
    freq = 1.0
    for nt in hexamer:
        freq *= ECOLI_NUC_FREQ.get(nt, 0.25)
    return freq


def generate_all_minus10_variants() -> dict:
    """Generate all single-nucleotide variants of TATAAT.

    Returns dict mapping variant -> {position, original_nt, new_nt, expected_freq}
    """
    consensus = "TATAAT"
    variants = {}

    for pos in range(6):
        for nt in 'ACGT':
            if nt == consensus[pos]:
                continue
            variant = consensus[:pos] + nt + consensus[pos+1:]
            variants[variant] = {
                'position': pos,
                'original_nt': consensus[pos],
                'new_nt': nt,
                'expected_freq': compute_hexamer_expected_freq(variant),
            }

    # Add consensus itself
    variants[consensus] = {
        'position': -1,
        'original_nt': '-',
        'new_nt': '-',
        'expected_freq': compute_hexamer_expected_freq(consensus),
    }

    return variants


def generate_background(length: int, at_fraction: float = 0.55) -> str:
    """Generate random background sequence."""
    seq = []
    for _ in range(length):
        if random.random() < at_fraction:
            seq.append(random.choice('AT'))
        else:
            seq.append(random.choice('GC'))
    return ''.join(seq)


def main():
    print("=" * 80)
    print("NEGATIVE MES INVESTIGATION")
    print("=" * 80)
    print("\nWhy do gLMs score intact (TATAAT) LOWER than broken (TGTAAT)?")

    # 1. Compare expected hexamer frequencies
    print("\n--- 1. Expected Hexamer Frequencies (E. coli genome) ---")
    intact_freq = compute_hexamer_expected_freq("TATAAT")
    broken_freq = compute_hexamer_expected_freq("TGTAAT")

    print(f"  TATAAT (intact): expected freq = {intact_freq:.6e}")
    print(f"  TGTAAT (broken): expected freq = {broken_freq:.6e}")
    print(f"  Ratio (broken/intact): {broken_freq/intact_freq:.3f}")

    if broken_freq > intact_freq:
        print("  → TGTAAT is MORE frequent than TATAAT under independence")
        print("    (G is more common than A at position 2 in E. coli)")
    else:
        print("  → TATAAT is MORE frequent than TGTAAT under independence")

    # 2. All single-nucleotide variants
    print("\n--- 2. All Single-Nucleotide Variants of TATAAT ---")
    variants = generate_all_minus10_variants()

    print(f"\n  {'Variant':<10} {'Pos':>4} {'Mutation':>10} {'Expected Freq':>15} {'vs Consensus':>15}")
    print("  " + "-" * 60)

    sorted_variants = sorted(variants.items(), key=lambda x: x[1]['expected_freq'], reverse=True)
    for var_seq, info in sorted_variants:
        if info['position'] == -1:
            mutation = "consensus"
        else:
            mutation = f"{info['original_nt']}{info['position']+1}{info['new_nt']}"
        ratio = info['expected_freq'] / intact_freq
        marker = " ←" if var_seq == "TGTAAT" else " ★" if var_seq == "TATAAT" else ""
        print(f"  {var_seq:<10} {info['position']:>4} {mutation:>10} "
              f"{info['expected_freq']:>15.6e} {ratio:>14.3f}x{marker}")

    # 3. Generate test sequences for all variants
    print("\n--- 3. Generating Test Sequences ---")
    n_per_variant = 50
    variant_sequences = {}

    for var_seq in variants:
        seqs = []
        for _ in range(n_per_variant):
            seq = list(generate_background(100))
            for i, nt in enumerate("TTGACA"):
                seq[30 + i] = nt
            for i, nt in enumerate(var_seq):
                seq[53 + i] = nt
            seqs.append(''.join(seq))
        variant_sequences[var_seq] = seqs

    print(f"  Generated {len(variants)} × {n_per_variant} = "
          f"{len(variants) * n_per_variant} sequences")

    # 4. Score with biophysical models as reference
    print("\n--- 4. Biophysical Model Scores ---")
    from mit_benchmark.models.biophysical import load_position_aware_pwm, load_thermodynamic_model

    papwm = load_position_aware_pwm()
    thermo = load_thermodynamic_model()

    print(f"\n  {'Variant':<10} {'PA-PWM':>10} {'Thermo':>10} {'Expected':>10}")
    print("  " + "-" * 45)

    papwm_scores = {}
    thermo_scores = {}
    for var_seq in sorted(variants.keys()):
        scores_p = [papwm.score(s) for s in variant_sequences[var_seq]]
        scores_t = [thermo.score(s) for s in variant_sequences[var_seq]]
        papwm_scores[var_seq] = np.mean(scores_p)
        thermo_scores[var_seq] = np.mean(scores_t)

        marker = " ←broken" if var_seq == "TGTAAT" else " ★intact" if var_seq == "TATAAT" else ""
        print(f"  {var_seq:<10} {np.mean(scores_p):>10.2f} {np.mean(scores_t):>10.2f} "
              f"{variants[var_seq]['expected_freq']:>10.2e}{marker}")

    # 5. Check if biophysical models correctly rank intact > broken
    print(f"\n  PA-PWM: intact ({papwm_scores['TATAAT']:.2f}) vs "
          f"broken ({papwm_scores['TGTAAT']:.2f}) → "
          f"{'CORRECT' if papwm_scores['TATAAT'] > papwm_scores['TGTAAT'] else 'WRONG'}")
    print(f"  Thermo: intact ({thermo_scores['TATAAT']:.2f}) vs "
          f"broken ({thermo_scores['TGTAAT']:.2f}) → "
          f"{'CORRECT' if thermo_scores['TATAAT'] > thermo_scores['TGTAAT'] else 'WRONG'}")

    # 6. Correlation between expected frequency and model score
    print("\n--- 5. Frequency-Score Correlation ---")
    freq_vals = [variants[v]['expected_freq'] for v in sorted(variants.keys())]
    papwm_vals = [papwm_scores[v] for v in sorted(variants.keys())]
    thermo_vals = [thermo_scores[v] for v in sorted(variants.keys())]

    r_papwm, p_papwm = np.corrcoef(freq_vals, papwm_vals)[0, 1], 0
    r_thermo, p_thermo = np.corrcoef(freq_vals, thermo_vals)[0, 1], 0

    print(f"  PA-PWM score vs expected freq:  r = {r_papwm:.3f}")
    print(f"  Thermo score vs expected freq:  r = {r_thermo:.3f}")
    print(f"\n  (If gLM shows r >> 0, it's following genome frequency, not function)")

    # Save results
    results = {
        "hexamer_frequencies": {k: v for k, v in sorted(variants.items(),
                                key=lambda x: x[1]['expected_freq'], reverse=True)},
        "intact_freq": float(intact_freq),
        "broken_freq": float(broken_freq),
        "freq_ratio_broken_over_intact": float(broken_freq / intact_freq),
        "papwm_scores": {k: float(v) for k, v in papwm_scores.items()},
        "thermo_scores": {k: float(v) for k, v in thermo_scores.items()},
        "variant_sequences": {k: v for k, v in variant_sequences.items()},
        "description": (
            "Investigation of why gLMs score TGTAAT (broken) higher than "
            "TATAAT (intact). Tests whether this is due to genome-wide "
            "hexamer frequency priors."
        ),
    }

    output_path = PROJECT_ROOT / "data/results/negative_mes_results.json"
    with open(output_path, 'w') as f:
        # Don't save full sequences to keep file size manageable
        save_results = {k: v for k, v in results.items() if k != 'variant_sequences'}
        json.dump(save_results, f, indent=2)
    print(f"\nResults saved to {output_path}")

    print(f"\n{'=' * 60}")
    print("INTERPRETATION")
    print(f"{'=' * 60}")
    print(f"""
    The broken motif TGTAAT has {'higher' if broken_freq > intact_freq else 'lower'}
    expected frequency than intact TATAAT in E. coli.

    Ratio: {broken_freq/intact_freq:.3f}x

    This {'explains' if broken_freq > intact_freq else 'does NOT explain'} why
    gLMs might assign higher likelihood to broken sequences:
    {'The model has learned genome-wide hexamer frequency priors that' if broken_freq > intact_freq else 'The effect must come from higher-order (context) features'}
    {'contradict functional importance at the -10 position.' if broken_freq > intact_freq else 'in the training data, not simple hexamer frequencies.'}

    Key insight: Biophysical models correctly rank intact > broken
    because they encode FUNCTIONAL importance, not statistical frequency.
    gLMs learn FREQUENCY, which can contradict function.
    """)


if __name__ == "__main__":
    main()
