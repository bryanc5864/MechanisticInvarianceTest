#!/usr/bin/env python3
"""Experiment 3.2: Spacing Sensitivity
Test if model knows optimal 17bp spacing between -35 and -10.
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


def generate_background(length, at_fraction=0.55):
    """Generate sequence with specified AT content."""
    seq = []
    for _ in range(length):
        if random.random() < at_fraction:
            seq.append(random.choice('AT'))
        else:
            seq.append(random.choice('GC'))
    return ''.join(seq)


def generate_promoter_with_spacing(spacing):
    """Generate promoter with specified spacing between -35 and -10."""
    # -35 box at fixed position 30
    minus_35_start = 30
    minus_35 = "TTGACA"  # 6bp

    # -10 box position depends on spacing
    minus_10_start = minus_35_start + 6 + spacing
    minus_10 = "TATAAT"  # 6bp (intact for this test)

    # Need sequence long enough
    total_len = max(100, minus_10_start + 6 + 20)

    seq = list(generate_background(total_len))

    # Insert -35 box
    for i, base in enumerate(minus_35):
        seq[minus_35_start + i] = base

    # Insert -10 box
    for i, base in enumerate(minus_10):
        seq[minus_10_start + i] = base

    return ''.join(seq[:100])  # Truncate to 100bp


def main():
    import os
    from scripts.run_inference import get_model_wrapper

    parser = argparse.ArgumentParser(description="Spacing Sensitivity Experiment")
    parser.add_argument("--model", type=str, default="hyenadna", help="Model name (e.g., hyenadna, evo2_1b, caduceus)")
    parser.add_argument("--gpu", type=int, default=0, help="GPU device ID")
    args = parser.parse_args()

    os.environ['CUDA_VISIBLE_DEVICES'] = str(args.gpu)
    model_name = args.model
    device = "cuda" if args.gpu >= 0 else "cpu"

    # Spacings to test (optimal is 17bp)
    spacings = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    n_per_spacing = 50

    # Generate sequences
    sequences = {}
    for spacing in spacings:
        sequences[spacing] = [generate_promoter_with_spacing(spacing) for _ in range(n_per_spacing)]

    # Load model
    print(f"Loading {model_name}...")
    model = get_model_wrapper(model_name, device)
    model.load_model()

    def compute_ll(sequence):
        return model.compute_log_likelihood(sequence)

    # Compute LL for all sequences
    print("Computing log-likelihoods...")
    results = {}
    for spacing in spacings:
        print(f"  Spacing = {spacing}bp...")
        lls = [compute_ll(s) for s in sequences[spacing]]
        results[spacing] = {
            'mean': np.mean(lls),
            'std': np.std(lls),
            'values': lls
        }

    # Print results
    print("\n" + "="*70)
    print("EXPERIMENT 3.2: SPACING SENSITIVITY RESULTS")
    print("="*70)

    print("\nMean LL by spacing (optimal is 17bp):")
    print("-"*50)
    print(f"{'Spacing':<12} {'Mean LL':<15} {'Std':<10} {'vs Optimal':<12}")
    print("-"*50)

    optimal_mean = results[17]['mean']
    for spacing in spacings:
        mean = results[spacing]['mean']
        std = results[spacing]['std']
        diff = mean - optimal_mean
        marker = " <-- OPTIMAL" if spacing == 17 else ""
        print(f"{spacing:<12} {mean:<15.3f} {std:<10.3f} {diff:+.3f}{marker}")

    # Find peak
    means = [results[s]['mean'] for s in spacings]
    peak_spacing = spacings[np.argmax(means)]

    print("\n" + "="*70)
    print("KEY FINDINGS")
    print("="*70)

    print(f"\n1. Peak LL spacing:       {peak_spacing}bp")
    print(f"   (Optimal is 17bp)")

    print(f"\n2. Range of means:")
    print(f"   Max: {max(means):.3f} (at {peak_spacing}bp)")
    print(f"   Min: {min(means):.3f} (at {spacings[np.argmin(means)]}bp)")
    print(f"   Difference: {max(means) - min(means):.3f}")

    # If mechanistic, expect clear peak at 17 with decay
    # If not, expect flat or noisy

    from scipy import stats
    corr, p = stats.pearsonr(spacings, means)
    print(f"\n3. Correlation (spacing vs LL): r = {corr:.3f}, p = {p:.4f}")

    # Quadratic fit (expect peak at 17)
    coeffs = np.polyfit(spacings, means, 2)
    peak_from_fit = -coeffs[1] / (2 * coeffs[0])
    print(f"4. Peak from quadratic fit: {peak_from_fit:.1f}bp")

    # Unload model
    model.unload_model()

    # Save results
    output_path = Path(f'data/results/spacing_{model_name}_results.json')
    save_results = {str(k): {'mean': float(v['mean']), 'std': float(v['std'])} for k, v in results.items()}
    with open(output_path, 'w') as f:
        json.dump(save_results, f, indent=2)

    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
