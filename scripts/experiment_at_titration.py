#!/usr/bin/env python3
"""Experiment 2.1: AT Titration
Test if model LL correlates with background AT content independent of motifs.
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

# Sequence layout
MINUS_35_START = 30
MINUS_35_END = 36
MINUS_10_START = 53
MINUS_10_END = 59
UP_START = 15
UP_END = 24
EXT_10_START = 50
EXT_10_END = 53


def generate_background(length, at_fraction):
    """Generate sequence with specified AT content."""
    seq = []
    for _ in range(length):
        if random.random() < at_fraction:
            seq.append(random.choice('AT'))
        else:
            seq.append(random.choice('GC'))
    return ''.join(seq)


def generate_promoter_with_at_background(at_fraction, motif_config):
    """Generate 100bp promoter with specified background AT and motifs."""
    seq = list(generate_background(100, at_fraction))

    # Insert -35 box
    minus_35 = motif_config.get('minus_35', 'TTGACA')
    for i, base in enumerate(minus_35):
        seq[MINUS_35_START + i] = base

    # Insert -10 box
    minus_10 = motif_config.get('minus_10', 'TATAAT')
    for i, base in enumerate(minus_10):
        seq[MINUS_10_START + i] = base

    # Optionally insert UP element
    if motif_config.get('up_element'):
        up = motif_config['up_element']
        for i, base in enumerate(up[:9]):
            seq[UP_START + i] = base

    # Optionally insert extended -10
    if motif_config.get('extended_10'):
        ext = motif_config['extended_10']
        for i, base in enumerate(ext[:3]):
            seq[EXT_10_START + i] = base

    return ''.join(seq)


def main():
    import os
    from scripts.run_inference import get_model_wrapper

    parser = argparse.ArgumentParser(description="AT Titration Experiment")
    parser.add_argument("--model", type=str, default="hyenadna", help="Model name (e.g., hyenadna, evo2_1b, caduceus)")
    parser.add_argument("--gpu", type=int, default=0, help="GPU device ID")
    args = parser.parse_args()

    os.environ['CUDA_VISIBLE_DEVICES'] = str(args.gpu)
    model_name = args.model
    device = "cuda" if args.gpu >= 0 else "cpu"

    # AT levels to test
    at_levels = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    n_per_level = 50

    # Motif configurations
    configs = {
        'intact': {'minus_35': 'TTGACA', 'minus_10': 'TATAAT'},
        'broken': {'minus_35': 'TTGACA', 'minus_10': 'TGTAAT'},
        'compensated': {'minus_35': 'TTGACA', 'minus_10': 'TGTAAT',
                       'up_element': 'AAAAAAGAG', 'extended_10': 'TGT'},
    }

    # Generate sequences
    sequences = {}
    for config_name, config in configs.items():
        sequences[config_name] = {}
        for at in at_levels:
            sequences[config_name][at] = [
                generate_promoter_with_at_background(at, config)
                for _ in range(n_per_level)
            ]

    # Load model
    print(f"Loading {model_name}...")
    model = get_model_wrapper(model_name, device)
    model.load_model()

    def compute_ll(sequence):
        return model.compute_log_likelihood(sequence)

    # Compute LL for all sequences
    print("Computing log-likelihoods...")
    results = {}
    for config_name in configs:
        results[config_name] = {}
        for at in at_levels:
            print(f"  {config_name} AT={at}...")
            lls = [compute_ll(s) for s in sequences[config_name][at]]
            results[config_name][at] = {
                'mean': np.mean(lls),
                'std': np.std(lls),
                'values': lls
            }

    # Print results
    print("\n" + "="*70)
    print("EXPERIMENT 2.1: AT TITRATION RESULTS")
    print("="*70)

    print("\nMean LL by AT content and motif configuration:")
    print("-"*70)
    print(f"{'AT%':<8}", end="")
    for config_name in configs:
        print(f"{config_name:<20}", end="")
    print()
    print("-"*70)

    for at in at_levels:
        print(f"{at:<8.1f}", end="")
        for config_name in configs:
            mean = results[config_name][at]['mean']
            std = results[config_name][at]['std']
            print(f"{mean:>8.2f} ± {std:<8.2f}", end="")
        print()

    # Compute correlations
    print("\n" + "="*70)
    print("CORRELATIONS (LL vs AT content)")
    print("="*70)

    for config_name in configs:
        all_at = []
        all_ll = []
        for at in at_levels:
            all_at.extend([at] * n_per_level)
            all_ll.extend(results[config_name][at]['values'])

        corr = np.corrcoef(all_at, all_ll)[0, 1]
        print(f"  {config_name}: r = {corr:.3f}")

    # Key test: Does motif matter at fixed AT?
    print("\n" + "="*70)
    print("MOTIF EFFECT AT FIXED AT CONTENT")
    print("="*70)

    for at in at_levels:
        intact_mean = results['intact'][at]['mean']
        broken_mean = results['broken'][at]['mean']
        comp_mean = results['compensated'][at]['mean']

        print(f"\nAT = {at}:")
        print(f"  Intact - Broken:      {intact_mean - broken_mean:+.3f}")
        print(f"  Compensated - Broken: {comp_mean - broken_mean:+.3f}")

    # Unload model
    model.unload_model()

    # Save results
    output_path = Path(f'data/results/at_titration_{model_name}_results.json')
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Convert to serializable format
    save_results = {}
    for config_name in results:
        save_results[config_name] = {}
        for at in results[config_name]:
            save_results[config_name][str(at)] = {
                'mean': float(results[config_name][at]['mean']),
                'std': float(results[config_name][at]['std']),
            }

    with open(output_path, 'w') as f:
        json.dump(save_results, f, indent=2)

    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()
