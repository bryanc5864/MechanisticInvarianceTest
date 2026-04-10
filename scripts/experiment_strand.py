#!/usr/bin/env python3
"""Experiment 3.3: Strand Orientation Test
Test if model knows promoter elements are strand-specific.
"""

import argparse
import json
import random
import sys
import numpy as np
import torch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

random.seed(42)
np.random.seed(42)

COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def reverse_complement(seq):
    """Get reverse complement of DNA sequence."""
    return ''.join(COMPLEMENT[b] for b in reversed(seq))


def generate_background(length, at_fraction=0.55):
    """Generate sequence with specified AT content."""
    seq = []
    for _ in range(length):
        if random.random() < at_fraction:
            seq.append(random.choice('AT'))
        else:
            seq.append(random.choice('GC'))
    return ''.join(seq)


def generate_forward_promoter():
    """Generate promoter on forward strand."""
    seq = list(generate_background(100))

    # -35 box: TTGACA at position 30
    minus_35 = "TTGACA"
    for i, base in enumerate(minus_35):
        seq[30 + i] = base

    # -10 box: TATAAT at position 53
    minus_10 = "TATAAT"
    for i, base in enumerate(minus_10):
        seq[53 + i] = base

    return ''.join(seq)


def generate_reverse_in_place():
    """Put reverse complement of motifs in same positions (wrong orientation)."""
    seq = list(generate_background(100))

    # Reverse complement of -35: TGTCAA at position 30
    minus_35_rc = reverse_complement("TTGACA")  # TGTCAA
    for i, base in enumerate(minus_35_rc):
        seq[30 + i] = base

    # Reverse complement of -10: ATTATA at position 53
    minus_10_rc = reverse_complement("TATAAT")  # ATTATA
    for i, base in enumerate(minus_10_rc):
        seq[53 + i] = base

    return ''.join(seq)


def generate_reverse_promoter():
    """Full reverse complement (correct for opposite strand)."""
    fwd = generate_forward_promoter()
    return reverse_complement(fwd)


def generate_scrambled_motifs():
    """Scramble motif sequences (same composition, wrong sequence)."""
    seq = list(generate_background(100))

    # Scrambled -35 (same composition as TTGACA)
    minus_35_scrambled = list("TTGACA")
    random.shuffle(minus_35_scrambled)
    for i, base in enumerate(minus_35_scrambled):
        seq[30 + i] = base

    # Scrambled -10 (same composition as TATAAT)
    minus_10_scrambled = list("TATAAT")
    random.shuffle(minus_10_scrambled)
    for i, base in enumerate(minus_10_scrambled):
        seq[53 + i] = base

    return ''.join(seq)


def main():
    import os
    from scripts.run_inference import get_model_wrapper

    parser = argparse.ArgumentParser(description="Strand Orientation Experiment")
    parser.add_argument("--model", type=str, default="hyenadna", help="Model name (e.g., hyenadna, evo2_1b, caduceus)")
    parser.add_argument("--gpu", type=int, default=0, help="GPU device ID")
    args = parser.parse_args()

    os.environ['CUDA_VISIBLE_DEVICES'] = str(args.gpu)
    model_name = args.model
    device = "cuda" if args.gpu >= 0 else "cpu"

    n_per_condition = 50

    conditions = {
        'forward': [generate_forward_promoter() for _ in range(n_per_condition)],
        'reverse_in_place': [generate_reverse_in_place() for _ in range(n_per_condition)],
        'full_reverse': [generate_reverse_promoter() for _ in range(n_per_condition)],
        'scrambled': [generate_scrambled_motifs() for _ in range(n_per_condition)],
    }

    # Load model
    print(f"Loading {model_name}...")
    model = get_model_wrapper(model_name, device)
    model.load_model()

    def compute_ll(sequence):
        return model.compute_log_likelihood(sequence)

    # Compute LL for all sequences
    print("Computing log-likelihoods...")
    results = {}
    for condition, seqs in conditions.items():
        print(f"  {condition}...")
        lls = [compute_ll(s) for s in seqs]
        results[condition] = {
            'mean': np.mean(lls),
            'std': np.std(lls),
            'values': lls
        }

    # Print results
    print("\n" + "="*70)
    print("EXPERIMENT 3.3: STRAND ORIENTATION RESULTS")
    print("="*70)

    print("\nMean LL by condition:")
    print("-"*60)
    print(f"{'Condition':<20} {'Mean LL':<15} {'Std':<10} {'Description'}")
    print("-"*60)

    descriptions = {
        'forward': 'Correct strand orientation',
        'reverse_in_place': 'RC motifs in same position (wrong)',
        'full_reverse': 'Full RC (correct for - strand)',
        'scrambled': 'Scrambled motifs (control)',
    }

    for condition in conditions:
        mean = results[condition]['mean']
        std = results[condition]['std']
        desc = descriptions[condition]
        print(f"{condition:<20} {mean:<15.3f} {std:<10.3f} {desc}")

    # Key tests
    print("\n" + "="*70)
    print("KEY FINDINGS")
    print("="*70)

    fwd_mean = results['forward']['mean']
    rev_mean = results['reverse_in_place']['mean']
    full_rev_mean = results['full_reverse']['mean']
    scrambled_mean = results['scrambled']['mean']

    print(f"\n1. Forward vs RC-in-place:  {fwd_mean - rev_mean:+.3f}")
    print(f"   (If strand-aware: Forward >> RC-in-place)")

    print(f"\n2. Forward vs Full-RC:      {fwd_mean - full_rev_mean:+.3f}")
    print(f"   (Both should be similar if strand-blind)")

    print(f"\n3. Forward vs Scrambled:    {fwd_mean - scrambled_mean:+.3f}")
    print(f"   (Motif recognition test)")

    # Pairwise comparisons
    from scipy import stats

    print("\n4. Statistical tests:")
    fwd_lls = results['forward']['values']
    for condition in ['reverse_in_place', 'full_reverse', 'scrambled']:
        other_lls = results[condition]['values']
        t_stat, p_val = stats.ttest_ind(fwd_lls, other_lls)
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        print(f"   Forward vs {condition}: t={t_stat:+.2f}, p={p_val:.4f} {sig}")

    # Unload model
    model.unload_model()

    # Save results
    output_path = Path(f'data/results/strand_{model_name}_results.json')
    save_results = {k: {'mean': float(v['mean']), 'std': float(v['std'])} for k, v in results.items()}
    with open(output_path, 'w') as f:
        json.dump(save_results, f, indent=2)

    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
