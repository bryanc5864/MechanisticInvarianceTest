#!/usr/bin/env python3
"""Run model inference on benchmark sequences.

Usage:
    python scripts/run_inference.py \
        --sequences data/sequences/all_sequences.json \
        --models evo2_1b,dnabert2,nt_500m \
        --output data/results/ \
        --gpu 0
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np
from tqdm import tqdm

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mit_benchmark.sequences.generator import SequenceGenerator, PromoterSequence
from mit_benchmark.utils.config import RESULTS_DIR, MODEL_CONFIGS


def get_model_wrapper(model_name: str, device: str):
    """Get the appropriate model wrapper.

    Args:
        model_name: Model identifier
        device: Device for inference

    Returns:
        Model wrapper instance
    """
    from mit_benchmark.models.autoregressive import Evo2Wrapper, HyenaDNAWrapper
    from mit_benchmark.models.masked_lm import (
        DNABERT2Wrapper,
        NucleotideTransformerWrapper,
        GROVERWrapper,
        CaduceusWrapper,
    )
    from mit_benchmark.models.baselines import KmerBaseline, PWMBaseline, RandomBaseline

    model_map = {
        "evo2_1b": lambda: Evo2Wrapper(device=device),
        "dnabert2": lambda: DNABERT2Wrapper(device=device),
        "nt_500m": lambda: NucleotideTransformerWrapper(device=device),
        "hyenadna": lambda: HyenaDNAWrapper(device=device),
        "grover": lambda: GROVERWrapper(device=device),
        "caduceus": lambda: CaduceusWrapper(device=device),
        "kmer": lambda: KmerBaseline(),
        "pwm": lambda: PWMBaseline(),
        "random": lambda: RandomBaseline(),
    }

    if model_name not in model_map:
        raise ValueError(f"Unknown model: {model_name}. Available: {list(model_map.keys())}")

    return model_map[model_name]()


def run_inference_single_model(
    model_name: str,
    sequences: List[PromoterSequence],
    device: str,
    checkpoint_path: Optional[Path] = None,
) -> Dict[str, List[float]]:
    """Run inference for a single model.

    Args:
        model_name: Model identifier
        sequences: List of PromoterSequence objects
        device: Device for inference
        checkpoint_path: Path to save checkpoints

    Returns:
        Dictionary mapping class labels to lists of log-likelihoods
    """
    print(f"\nRunning inference for: {model_name}")
    print("-" * 40)

    # Get model wrapper
    model = get_model_wrapper(model_name, device)

    # Load model
    model.load_model()

    # Group sequences by class
    by_class: Dict[str, List[PromoterSequence]] = {}
    for seq in sequences:
        if seq.class_label not in by_class:
            by_class[seq.class_label] = []
        by_class[seq.class_label].append(seq)

    # Run inference
    results: Dict[str, List[float]] = {}

    for class_label in sorted(by_class.keys()):
        class_seqs = by_class[class_label]
        print(f"  Class {class_label}: {len(class_seqs)} sequences")

        # Extract sequence strings
        seq_strings = [s.sequence for s in class_seqs]

        # Compute log-likelihoods
        lls = model.compute_batch_log_likelihoods(seq_strings, batch_size=8)
        results[class_label] = lls

        # Save checkpoint
        if checkpoint_path:
            checkpoint = {
                "model": model_name,
                "completed_classes": list(results.keys()),
                "results": {k: [float(x) for x in v] for k, v in results.items()},
            }
            with open(checkpoint_path, 'w') as f:
                json.dump(checkpoint, f)

    # Unload model to free memory
    model.unload_model()

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Run model inference on MIT benchmark sequences"
    )
    parser.add_argument(
        "--sequences",
        type=Path,
        required=True,
        help="Path to sequences JSON file",
    )
    parser.add_argument(
        "--models",
        type=str,
        default="kmer,pwm,random",
        help="Comma-separated list of models to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=RESULTS_DIR,
        help="Output directory for results",
    )
    parser.add_argument(
        "--gpu",
        type=int,
        default=0,
        help="GPU device ID",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from checkpoints if available",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("MIT Benchmark Inference")
    print("=" * 60)

    # Set GPU
    if args.gpu >= 0:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
        device = "cuda"
    else:
        device = "cpu"

    print(f"Device: {device}")
    print(f"Sequences: {args.sequences}")
    print(f"Models: {args.models}")
    print(f"Output: {args.output}")
    print()

    # Load sequences
    print("Loading sequences...")
    sequences = SequenceGenerator.load_sequences(args.sequences)
    print(f"  Loaded {len(sequences)} sequences")
    print()

    # Parse models
    model_names = [m.strip() for m in args.models.split(",")]
    print(f"Models to run: {model_names}")
    print()

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Run inference for each model
    all_results = {}

    for model_name in model_names:
        checkpoint_path = args.output / f"{model_name}_checkpoint.json"
        result_path = args.output / f"{model_name}_results.json"

        # Check for existing results
        if args.resume and result_path.exists():
            print(f"Skipping {model_name} (results exist)")
            with open(result_path, 'r') as f:
                all_results[model_name] = json.load(f)
            continue

        try:
            # Run inference
            results = run_inference_single_model(
                model_name,
                sequences,
                device,
                checkpoint_path,
            )

            # Save results
            with open(result_path, 'w') as f:
                json.dump({k: [float(x) for x in v] for k, v in results.items()}, f, indent=2)

            all_results[model_name] = {k: [float(x) for x in v] for k, v in results.items()}
            print(f"  Results saved to: {result_path}")

            # Remove checkpoint
            if checkpoint_path.exists():
                checkpoint_path.unlink()

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Save combined results
    combined_path = args.output / "all_results.json"
    with open(combined_path, 'w') as f:
        json.dump(all_results, f, indent=2)

    print()
    print("=" * 60)
    print(f"Inference complete! Results saved to: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
