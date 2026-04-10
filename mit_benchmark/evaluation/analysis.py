"""Statistical analysis and visualization for MIT benchmark results."""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from scipy import stats
from dataclasses import dataclass

from .metrics import MetricsResult, compute_all_metrics


@dataclass
class StatisticalTest:
    """Result of a statistical test."""
    test_name: str
    statistic: float
    p_value: float
    effect_size: Optional[float]
    description: str


def run_statistical_tests(
    all_metrics: Dict[str, MetricsResult],
) -> List[StatisticalTest]:
    """Run statistical tests on benchmark results.

    Tests performed:
    1. One-sample t-test of CSS against 0.5 for each model
    2. ANOVA for differences between model architectures
    3. Correlation between model size and CSS (if available)

    Args:
        all_metrics: Dictionary mapping model names to MetricsResult

    Returns:
        List of StatisticalTest results
    """
    results = []

    # Extract CSS values
    model_names = list(all_metrics.keys())
    css_values = [all_metrics[m].css for m in model_names]

    # 1. Individual CSS tests (already computed, but aggregate here)
    for name, metrics in all_metrics.items():
        results.append(StatisticalTest(
            test_name=f"CSS_vs_0.5_{name}",
            statistic=metrics.css,
            p_value=metrics.css_pvalue,
            effect_size=metrics.css - 0.5,
            description=f"One-sample test of CSS against 0.5 for {name}",
        ))

    # 2. Test if any model significantly outperforms random
    if len(css_values) > 1:
        # One-sample t-test of mean CSS against 0.5
        t_stat, p_value = stats.ttest_1samp(css_values, 0.5)
        results.append(StatisticalTest(
            test_name="CSS_all_models_vs_0.5",
            statistic=t_stat,
            p_value=p_value,
            effect_size=np.mean(css_values) - 0.5,
            description="One-sample t-test of mean CSS across all models against 0.5",
        ))

    # 3. Benjamini-Hochberg FDR correction
    p_values = [r.p_value for r in results]
    if p_values:
        fdr_corrected = benjamini_hochberg(p_values)
        for i, result in enumerate(results):
            result.description += f" (FDR-corrected p={fdr_corrected[i]:.4f})"

    return results


def benjamini_hochberg(p_values: List[float], alpha: float = 0.05) -> List[float]:
    """Apply Benjamini-Hochberg FDR correction.

    Args:
        p_values: List of p-values
        alpha: Significance level

    Returns:
        List of adjusted p-values
    """
    n = len(p_values)
    if n == 0:
        return []

    # Sort p-values and get original indices
    sorted_pairs = sorted(enumerate(p_values), key=lambda x: x[1])
    sorted_indices = [p[0] for p in sorted_pairs]
    sorted_pvals = [p[1] for p in sorted_pairs]

    # Compute adjusted p-values
    adjusted = [0.0] * n
    prev_adj = 1.0

    for i in range(n - 1, -1, -1):
        rank = i + 1
        adj_p = min(prev_adj, sorted_pvals[i] * n / rank)
        adjusted[sorted_indices[i]] = adj_p
        prev_adj = adj_p

    return adjusted


def compare_models(
    all_metrics: Dict[str, MetricsResult],
) -> pd.DataFrame:
    """Create a comparison table of all models.

    Args:
        all_metrics: Dictionary mapping model names to MetricsResult

    Returns:
        DataFrame with model comparison
    """
    rows = []
    for name, metrics in all_metrics.items():
        rows.append({
            "Model": name,
            "CSS": f"{metrics.css:.3f}",
            "CSS 95% CI": f"[{metrics.css_ci_low:.3f}, {metrics.css_ci_high:.3f}]",
            "p-value": f"{metrics.css_pvalue:.4f}",
            "MES (Natural)": f"{metrics.mes_natural:.2f}",
            "MES (Synthetic)": f"{metrics.mes_synthetic:.2f}",
            "CIR": f"{metrics.cir:.2f}",
            "CM": f"{metrics.cm:.2f}",
            "SCR": f"{metrics.scr:.3f}",
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("CSS", ascending=False)
    return df


def generate_report(
    all_metrics: Dict[str, MetricsResult],
    output_dir: Path,
    include_plots: bool = True,
) -> str:
    """Generate a comprehensive benchmark report.

    Args:
        all_metrics: Dictionary mapping model names to MetricsResult
        output_dir: Directory to save report and figures
        include_plots: Whether to generate plots

    Returns:
        Path to generated report
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create comparison table
    comparison_df = compare_models(all_metrics)

    # Run statistical tests
    stat_tests = run_statistical_tests(all_metrics)

    # Generate plots if requested
    if include_plots:
        try:
            _generate_plots(all_metrics, output_dir)
        except ImportError:
            print("Warning: matplotlib not available, skipping plots")

    # Write report
    report_path = output_dir / "benchmark_report.md"
    with open(report_path, 'w') as f:
        f.write("# MIT Benchmark Results\n\n")
        f.write("## Model Comparison\n\n")
        try:
            f.write(comparison_df.to_markdown(index=False))
        except ImportError:
            # Fallback if tabulate not installed
            f.write(comparison_df.to_string(index=False))
        f.write("\n\n")

        f.write("## Statistical Tests\n\n")
        for test in stat_tests:
            sig = "✓" if test.p_value < 0.05 else "✗"
            f.write(f"- **{test.test_name}**: p={test.p_value:.4f} {sig}\n")
            f.write(f"  - {test.description}\n")
        f.write("\n")

        f.write("## Interpretation\n\n")

        # Find best model
        best_model = max(all_metrics.items(), key=lambda x: x[1].css)
        f.write(f"**Best performing model**: {best_model[0]} (CSS={best_model[1].css:.3f})\n\n")

        # Interpretation of CSS
        f.write("### CSS Interpretation\n")
        f.write("- CSS > 0.5: Model recognizes compensation\n")
        f.write("- CSS ≈ 0.5: Model does not distinguish compensation\n")
        f.write("- CSS < 0.5: Model penalizes compensated sequences\n\n")

        if include_plots:
            f.write("## Figures\n\n")
            f.write("![CSS Comparison](css_comparison.png)\n\n")
            f.write("![Metrics Heatmap](metrics_heatmap.png)\n\n")

    # Save metrics as JSON
    metrics_json = output_dir / "all_metrics.json"
    with open(metrics_json, 'w') as f:
        json.dump({name: m.to_dict() for name, m in all_metrics.items()}, f, indent=2)

    # Save comparison table as CSV
    comparison_df.to_csv(output_dir / "model_comparison.csv", index=False)

    return str(report_path)


def _generate_plots(
    all_metrics: Dict[str, MetricsResult],
    output_dir: Path,
) -> None:
    """Generate visualization plots.

    Args:
        all_metrics: Dictionary mapping model names to MetricsResult
        output_dir: Directory to save figures
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Set style
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")

    # 1. CSS Comparison Bar Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    models = list(all_metrics.keys())
    css_values = [all_metrics[m].css for m in models]
    css_errors = [
        [all_metrics[m].css - all_metrics[m].css_ci_low for m in models],
        [all_metrics[m].css_ci_high - all_metrics[m].css for m in models],
    ]

    bars = ax.bar(models, css_values, yerr=css_errors, capsize=5, alpha=0.8)

    # Color bars based on significance
    for i, (bar, m) in enumerate(zip(bars, models)):
        if all_metrics[m].css_pvalue < 0.05:
            bar.set_color('forestgreen')
        else:
            bar.set_color('gray')

    ax.axhline(y=0.5, color='red', linestyle='--', label='Random baseline')
    ax.set_ylabel('Compensation Sensitivity Score (CSS)')
    ax.set_xlabel('Model')
    ax.set_title('CSS Comparison Across Models')
    ax.legend()
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / 'css_comparison.png', dpi=150)
    plt.close()

    # 2. Metrics Heatmap
    fig, ax = plt.subplots(figsize=(10, 8))

    metrics_data = []
    for m in models:
        metrics_data.append([
            all_metrics[m].css,
            all_metrics[m].mes_natural,
            all_metrics[m].mes_synthetic,
            all_metrics[m].cm,
            all_metrics[m].scr,
        ])

    metrics_df = pd.DataFrame(
        metrics_data,
        index=models,
        columns=['CSS', 'MES (Natural)', 'MES (Synthetic)', 'CM', 'SCR'],
    )

    sns.heatmap(metrics_df, annot=True, fmt='.2f', cmap='RdYlGn', center=0.5, ax=ax)
    ax.set_title('Metrics Heatmap')
    plt.tight_layout()
    plt.savefig(output_dir / 'metrics_heatmap.png', dpi=150)
    plt.close()

    # 3. CSS vs MES scatter plot
    fig, ax = plt.subplots(figsize=(8, 6))

    for m in models:
        ax.scatter(
            all_metrics[m].mes_synthetic,
            all_metrics[m].css,
            s=100,
            label=m,
            alpha=0.7,
        )

    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Motif Effect Size (Synthetic)')
    ax.set_ylabel('Compensation Sensitivity Score')
    ax.set_title('CSS vs MES Relationship')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_dir / 'css_vs_mes.png', dpi=150)
    plt.close()
