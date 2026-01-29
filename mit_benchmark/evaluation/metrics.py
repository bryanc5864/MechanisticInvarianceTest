"""Evaluation metrics for the MIT benchmark.

Primary metric: Compensation Sensitivity Score (CSS)
Secondary metrics: MES, CIR, CM, SCR
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
from scipy import stats


@dataclass
class MetricsResult:
    """Container for all evaluation metrics."""
    # Primary metric
    css: float  # Compensation Sensitivity Score
    css_ci_low: float  # 95% CI lower bound
    css_ci_high: float  # 95% CI upper bound
    css_pvalue: float  # P-value for CSS > 0.5

    # Secondary metrics
    mes_natural: float  # Motif Effect Size (natural intact vs broken)
    mes_synthetic: float  # Motif Effect Size (synthetic intact vs broken)
    cir: float  # Context Independence Ratio
    cm: float  # Compensation Magnitude
    scr: float  # Scramble Control Ratio

    # Bootstrap CIs for secondary metrics
    mes_natural_ci_low: float = 0.0
    mes_natural_ci_high: float = 0.0
    mes_synthetic_ci_low: float = 0.0
    mes_synthetic_ci_high: float = 0.0
    cm_ci_low: float = 0.0
    cm_ci_high: float = 0.0
    scr_ci_low: float = 0.0
    scr_ci_high: float = 0.0
    scr_pvalue: float = 1.0  # P-value for SCR > 0.5

    # Additional statistics
    n_samples: int = 0
    model_name: str = ""

    def to_dict(self) -> Dict:
        return asdict(self)


def compute_css(
    ll_broken: np.ndarray,
    ll_compensated: np.ndarray,
) -> Tuple[float, float, float, float]:
    """Compute Compensation Sensitivity Score.

    CSS measures the fraction of cases where compensated sequences
    have higher likelihood than broken sequences.

    CSS = mean(LL(compensated) > LL(broken))

    A model that understands compensation should have CSS > 0.5.
    A model that only sees surface-level motif damage should have CSS ≈ 0.5.

    Args:
        ll_broken: Log-likelihoods for broken sequences (class D)
        ll_compensated: Log-likelihoods for compensated sequences (class E)

    Returns:
        Tuple of (CSS, CI_low, CI_high, p_value)
    """
    ll_broken = np.asarray(ll_broken)
    ll_compensated = np.asarray(ll_compensated)

    # Pairwise comparisons
    # For matched pairs, compare directly
    # For unmatched, use all pairs
    n = min(len(ll_broken), len(ll_compensated))
    comparisons = ll_compensated[:n] > ll_broken[:n]

    css = np.mean(comparisons)

    # Bootstrap confidence interval
    n_bootstrap = 1000
    bootstrap_css = []
    for _ in range(n_bootstrap):
        idx = np.random.randint(0, n, size=n)
        bootstrap_css.append(np.mean(comparisons[idx]))

    ci_low = np.percentile(bootstrap_css, 2.5)
    ci_high = np.percentile(bootstrap_css, 97.5)

    # One-sample t-test against 0.5
    t_stat, p_value = stats.ttest_1samp(comparisons.astype(float), 0.5)
    # One-tailed test for CSS > 0.5
    if css > 0.5:
        p_value = p_value / 2
    else:
        p_value = 1 - p_value / 2

    return css, ci_low, ci_high, p_value


def _compute_mes_raw(ll_intact: np.ndarray, ll_broken: np.ndarray) -> float:
    """Compute Cohen's d effect size (internal helper)."""
    if len(ll_intact) == 0 or len(ll_broken) == 0:
        return 0.0

    mean_intact = np.mean(ll_intact)
    mean_broken = np.mean(ll_broken)

    n1, n2 = len(ll_intact), len(ll_broken)
    var1 = np.var(ll_intact, ddof=1) if n1 > 1 else 0.0
    var2 = np.var(ll_broken, ddof=1) if n2 > 1 else 0.0

    if n1 + n2 <= 2:
        return 0.0

    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

    if pooled_std == 0 or np.isnan(pooled_std):
        if mean_intact > mean_broken:
            return 10.0
        elif mean_intact < mean_broken:
            return -10.0
        return 0.0

    result = (mean_intact - mean_broken) / pooled_std
    return np.clip(result, -10.0, 10.0)


def compute_mes(
    ll_intact: np.ndarray,
    ll_broken: np.ndarray,
    return_ci: bool = False,
    n_bootstrap: int = 1000,
) -> float:
    """Compute Motif Effect Size.

    MES measures how much the model distinguishes intact from broken
    promoters using Cohen's d.

    MES = (mean(LL_intact) - mean(LL_broken)) / pooled_std

    A higher MES indicates the model recognizes the importance of
    the intact -10 box.

    Args:
        ll_intact: Log-likelihoods for intact sequences
        ll_broken: Log-likelihoods for broken sequences
        return_ci: If True, return (MES, ci_low, ci_high) tuple
        n_bootstrap: Number of bootstrap samples for CI

    Returns:
        Cohen's d effect size, or (d, ci_low, ci_high) if return_ci=True
    """
    ll_intact = np.asarray(ll_intact)
    ll_broken = np.asarray(ll_broken)

    mes = _compute_mes_raw(ll_intact, ll_broken)

    if not return_ci:
        return mes

    # Bootstrap CI
    bootstrap_mes = []
    n1, n2 = len(ll_intact), len(ll_broken)
    for _ in range(n_bootstrap):
        idx1 = np.random.randint(0, n1, size=n1)
        idx2 = np.random.randint(0, n2, size=n2)
        bootstrap_mes.append(_compute_mes_raw(ll_intact[idx1], ll_broken[idx2]))

    ci_low = np.percentile(bootstrap_mes, 2.5)
    ci_high = np.percentile(bootstrap_mes, 97.5)

    return mes, ci_low, ci_high


def compute_cir(
    mes_synthetic: float,
    mes_natural: float,
) -> float:
    """Compute Context Independence Ratio.

    CIR measures whether the model's motif sensitivity transfers
    from synthetic to natural contexts.

    CIR = MES_synthetic / MES_natural

    CIR ≈ 1 indicates consistent motif sensitivity.
    CIR >> 1 indicates overfitting to synthetic patterns.
    CIR << 1 indicates poor transfer to synthetic contexts.

    Args:
        mes_synthetic: MES for synthetic sequences (C vs D)
        mes_natural: MES for natural sequences (A vs B)

    Returns:
        Context Independence Ratio
    """
    if mes_natural == 0:
        return float('inf') if mes_synthetic > 0 else 1.0

    return mes_synthetic / mes_natural


def compute_cm(
    ll_broken: np.ndarray,
    ll_compensated: np.ndarray,
    ll_intact: np.ndarray,
) -> float:
    """Compute Compensation Magnitude.

    CM measures how much of the likelihood loss from breaking the -10
    is recovered by adding compensatory elements.

    CM = (mean(LL_compensated) - mean(LL_broken)) / (mean(LL_intact) - mean(LL_broken))

    CM = 0: No recovery (compensation not recognized)
    CM = 1: Full recovery (compensation fully compensates)
    CM > 1: Over-recovery (compensation adds more than it restores)

    Args:
        ll_broken: Log-likelihoods for broken sequences
        ll_compensated: Log-likelihoods for compensated sequences
        ll_intact: Log-likelihoods for intact sequences

    Returns:
        Compensation magnitude (fraction of recovery)
    """
    ll_broken = np.asarray(ll_broken)
    ll_compensated = np.asarray(ll_compensated)
    ll_intact = np.asarray(ll_intact)

    mean_broken = np.mean(ll_broken)
    mean_compensated = np.mean(ll_compensated)
    mean_intact = np.mean(ll_intact)

    denominator = mean_intact - mean_broken
    if denominator == 0:
        return 0.0

    return (mean_compensated - mean_broken) / denominator


def compute_scr(
    ll_compensated: np.ndarray,
    ll_scrambled: np.ndarray,
    return_ci: bool = False,
    n_bootstrap: int = 1000,
):
    """Compute Scramble Control Ratio.

    SCR tests whether the model's response to compensation is due to
    specific motif structure vs. general sequence composition.

    SCR = mean(LL_compensated > LL_scrambled)

    SCR >> 0.5 indicates the model recognizes motif structure.
    SCR ≈ 0.5 indicates the model only sees composition.

    Args:
        ll_compensated: Log-likelihoods for compensated sequences (E)
        ll_scrambled: Log-likelihoods for scrambled sequences (H)
        return_ci: If True, return (SCR, ci_low, ci_high, p_value) tuple
        n_bootstrap: Number of bootstrap samples for CI

    Returns:
        Scramble control ratio, or (SCR, ci_low, ci_high, p_value) if return_ci=True
    """
    ll_compensated = np.asarray(ll_compensated)
    ll_scrambled = np.asarray(ll_scrambled)

    n = min(len(ll_compensated), len(ll_scrambled))
    comparisons = ll_compensated[:n] > ll_scrambled[:n]
    scr = np.mean(comparisons)

    if not return_ci:
        return scr

    # Bootstrap CI
    bootstrap_scr = []
    for _ in range(n_bootstrap):
        idx = np.random.randint(0, n, size=n)
        bootstrap_scr.append(np.mean(comparisons[idx]))

    ci_low = np.percentile(bootstrap_scr, 2.5)
    ci_high = np.percentile(bootstrap_scr, 97.5)

    # P-value (one-tailed test for SCR > 0.5)
    t_stat, p_value = stats.ttest_1samp(comparisons.astype(float), 0.5)
    if scr > 0.5:
        p_value = p_value / 2
    else:
        p_value = 1 - p_value / 2

    return scr, ci_low, ci_high, p_value


def compute_all_metrics(
    predictions: Dict[str, np.ndarray],
    model_name: str,
) -> MetricsResult:
    """Compute all evaluation metrics.

    Args:
        predictions: Dictionary mapping class labels to log-likelihood arrays
            Expected keys: A, B, C, D, E, F, G, H
        model_name: Name of the model

    Returns:
        MetricsResult with all metrics
    """
    # Extract predictions by class
    ll_a = np.asarray(predictions.get("A", []))  # Natural intact
    ll_b = np.asarray(predictions.get("B", []))  # Natural broken
    ll_c = np.asarray(predictions.get("C", []))  # Synthetic intact
    ll_d = np.asarray(predictions.get("D", []))  # Synthetic broken
    ll_e = np.asarray(predictions.get("E", []))  # Synthetic compensated
    ll_f = np.asarray(predictions.get("F", []))  # Over-compensated
    ll_g = np.asarray(predictions.get("G", []))  # Natural compensated
    ll_h = np.asarray(predictions.get("H", []))  # Scrambled

    # Primary metric: CSS (broken vs compensated)
    css, css_ci_low, css_ci_high, css_pvalue = compute_css(ll_d, ll_e)

    # Secondary metrics with bootstrap CIs
    mes_natural, mes_nat_ci_low, mes_nat_ci_high = compute_mes(
        ll_a, ll_b, return_ci=True
    )
    mes_synthetic, mes_syn_ci_low, mes_syn_ci_high = compute_mes(
        ll_c, ll_d, return_ci=True
    )
    cir = compute_cir(mes_synthetic, mes_natural)
    cm = compute_cm(ll_d, ll_e, ll_c)
    scr, scr_ci_low, scr_ci_high, scr_pvalue = compute_scr(
        ll_e, ll_h, return_ci=True
    )

    # Bootstrap CI for CM
    n_boot = 1000
    cm_boots = []
    n_d, n_e, n_c = len(ll_d), len(ll_e), len(ll_c)
    if n_d > 0 and n_e > 0 and n_c > 0:
        for _ in range(n_boot):
            d_idx = np.random.randint(0, n_d, size=n_d)
            e_idx = np.random.randint(0, n_e, size=n_e)
            c_idx = np.random.randint(0, n_c, size=n_c)
            denom = np.mean(ll_c[c_idx]) - np.mean(ll_d[d_idx])
            if denom != 0:
                cm_boots.append(
                    (np.mean(ll_e[e_idx]) - np.mean(ll_d[d_idx])) / denom
                )
            else:
                cm_boots.append(0.0)
        cm_ci_low = np.percentile(cm_boots, 2.5)
        cm_ci_high = np.percentile(cm_boots, 97.5)
    else:
        cm_ci_low, cm_ci_high = 0.0, 0.0

    n_samples = sum(len(predictions.get(k, [])) for k in "ABCDEFGH")

    return MetricsResult(
        css=css,
        css_ci_low=css_ci_low,
        css_ci_high=css_ci_high,
        css_pvalue=css_pvalue,
        mes_natural=mes_natural,
        mes_synthetic=mes_synthetic,
        cir=cir,
        cm=cm,
        scr=scr,
        mes_natural_ci_low=mes_nat_ci_low,
        mes_natural_ci_high=mes_nat_ci_high,
        mes_synthetic_ci_low=mes_syn_ci_low,
        mes_synthetic_ci_high=mes_syn_ci_high,
        cm_ci_low=cm_ci_low,
        cm_ci_high=cm_ci_high,
        scr_ci_low=scr_ci_low,
        scr_ci_high=scr_ci_high,
        scr_pvalue=scr_pvalue,
        n_samples=n_samples,
        model_name=model_name,
    )
