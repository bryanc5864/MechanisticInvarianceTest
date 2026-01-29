#!/usr/bin/env python3
"""Compute evaluation metrics from inference results.

Usage:
    python scripts/compute_metrics.py \
        --results data/results/ \
        --output data/results/metrics.json
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict

import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mit_benchmark.evaluation.metrics import compute_all_metrics, MetricsResult
from mit_benchmark.utils.config import RESULTS_DIR


def main():
    parser = argparse.ArgumentParser(
        description="Compute MIT benchmark evaluation metrics"
    )
    parser.add_argument(
        "--results",
        type=Path,
        default=RESULTS_DIR,
        help="Directory containing inference results",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output file for metrics (default: results/metrics.json)",
    )

    args = parser.parse_args()

    if args.output is None:
        args.output = args.results / "metrics.json"

    print("=" * 60)
    print("MIT Benchmark Metrics Computation")
    print("=" * 60)
    print(f"Results directory: {args.results}")
    print(f"Output: {args.output}")
    print()

    # Find all result files (exclude combined all_results.json)
    result_files = [
        f for f in args.results.glob("*_results.json")
        if f.name != "all_results.json"
    ]
    if not result_files:
        print("ERROR: No result files found!")
        print(f"Looking in: {args.results}")
        sys.exit(1)

    print(f"Found {len(result_files)} result files:")
    for f in result_files:
        print(f"  - {f.name}")
    print()

    # Compute metrics for each model
    all_metrics: Dict[str, Dict] = {}

    for result_file in result_files:
        model_name = result_file.stem.replace("_results", "")
        print(f"Computing metrics for: {model_name}")

        # Load results
        with open(result_file, 'r') as f:
            results = json.load(f)

        # Convert to numpy arrays
        predictions = {k: np.array(v) for k, v in results.items()}

        # Compute metrics
        metrics = compute_all_metrics(predictions, model_name)

        # Print summary
        print(f"  CSS: {metrics.css:.3f} (p={metrics.css_pvalue:.4f})")
        print(f"  MES Natural: {metrics.mes_natural:.2f}")
        print(f"  MES Synthetic: {metrics.mes_synthetic:.2f}")
        print(f"  CIR: {metrics.cir:.2f}")
        print(f"  CM: {metrics.cm:.2f}")
        print(f"  SCR: {metrics.scr:.3f}")
        print()

        all_metrics[model_name] = metrics.to_dict()

    # Save metrics
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(all_metrics, f, indent=2)

    print("=" * 60)
    print(f"Metrics saved to: {args.output}")
    print("=" * 60)

    # Print comparison table
    print()
    print("Model Comparison (CSS):")
    print("-" * 40)
    sorted_models = sorted(all_metrics.items(), key=lambda x: x[1]["css"], reverse=True)
    for model_name, metrics in sorted_models:
        sig = "***" if metrics["css_pvalue"] < 0.001 else "**" if metrics["css_pvalue"] < 0.01 else "*" if metrics["css_pvalue"] < 0.05 else ""
        print(f"  {model_name:20s} CSS={metrics['css']:.3f} {sig}")


if __name__ == "__main__":
    main()
