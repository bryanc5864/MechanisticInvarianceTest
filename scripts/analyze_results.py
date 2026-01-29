#!/usr/bin/env python3
"""Analyze benchmark results and generate reports/figures.

Usage:
    python scripts/analyze_results.py \
        --metrics data/results/metrics.json \
        --output figures/
"""

import argparse
import json
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mit_benchmark.evaluation.metrics import MetricsResult
from mit_benchmark.evaluation.analysis import (
    run_statistical_tests,
    compare_models,
    generate_report,
)
from mit_benchmark.utils.config import FIGURES_DIR


def main():
    parser = argparse.ArgumentParser(
        description="Analyze MIT benchmark results"
    )
    parser.add_argument(
        "--metrics",
        type=Path,
        required=True,
        help="Path to metrics JSON file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=FIGURES_DIR,
        help="Output directory for report and figures",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip plot generation",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("MIT Benchmark Analysis")
    print("=" * 60)
    print(f"Metrics file: {args.metrics}")
    print(f"Output directory: {args.output}")
    print()

    # Load metrics
    print("Loading metrics...")
    with open(args.metrics, 'r') as f:
        metrics_dict = json.load(f)

    # Convert to MetricsResult objects
    all_metrics = {}
    for model_name, m in metrics_dict.items():
        all_metrics[model_name] = MetricsResult(**m)

    print(f"  Loaded metrics for {len(all_metrics)} models")
    print()

    # Run statistical tests
    print("Running statistical tests...")
    stat_tests = run_statistical_tests(all_metrics)
    for test in stat_tests:
        sig = "✓" if test.p_value < 0.05 else "✗"
        print(f"  {test.test_name}: p={test.p_value:.4f} {sig}")
    print()

    # Create comparison table
    print("Model Comparison Table:")
    print("-" * 60)
    comparison_df = compare_models(all_metrics)
    print(comparison_df.to_string(index=False))
    print()

    # Generate report
    print("Generating report...")
    report_path = generate_report(
        all_metrics,
        args.output,
        include_plots=not args.no_plots,
    )
    print(f"  Report saved to: {report_path}")
    print()

    # Summary
    print("=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    print()
    print("Key Findings:")
    print("-" * 40)

    # Best model
    best_model = max(all_metrics.items(), key=lambda x: x[1].css)
    print(f"Best Model: {best_model[0]} (CSS={best_model[1].css:.3f})")

    # Models with significant CSS
    sig_models = [m for m, r in all_metrics.items() if r.css_pvalue < 0.05 and r.css > 0.5]
    if sig_models:
        print(f"Models with significant compensation sensitivity: {', '.join(sig_models)}")
    else:
        print("No models showed significant compensation sensitivity (p < 0.05)")

    # Average CSS
    avg_css = sum(m.css for m in all_metrics.values()) / len(all_metrics)
    print(f"Average CSS across models: {avg_css:.3f}")

    if not args.no_plots:
        print()
        print("Figures generated:")
        for fig in args.output.glob("*.png"):
            print(f"  - {fig.name}")


if __name__ == "__main__":
    main()
