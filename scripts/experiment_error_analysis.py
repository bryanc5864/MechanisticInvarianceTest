#!/usr/bin/env python3
"""Error analysis experiment.

Analyzes patterns in model successes and failures to understand what
sequence features predict correct classification. Key questions:
1. Which compensated sequences does HyenaDNA score correctly (> broken)?
2. Are failures correlated with GC content, specific background patterns,
   or motif variant properties?
3. Do all models fail on the same sequences (systematic) or different ones (random)?
"""

import json
import sys
from pathlib import Path
from collections import Counter

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def compute_sequence_properties(seq: str) -> dict:
    """Compute compositional properties of a sequence."""
    n = len(seq)
    at_content = (seq.count('A') + seq.count('T')) / n
    gc_content = (seq.count('G') + seq.count('C')) / n

    # Dinucleotide frequencies
    dinucs = Counter()
    for i in range(n - 1):
        dinucs[seq[i:i+2]] += 1

    # AT runs (longest consecutive A/T stretch)
    max_at_run = 0
    current_run = 0
    for c in seq:
        if c in 'AT':
            current_run += 1
            max_at_run = max(max_at_run, current_run)
        else:
            current_run = 0

    # GC runs
    max_gc_run = 0
    current_run = 0
    for c in seq:
        if c in 'GC':
            current_run += 1
            max_gc_run = max(max_gc_run, current_run)
        else:
            current_run = 0

    # Regional AT content (upstream, middle, downstream)
    upstream_at = (seq[:30].count('A') + seq[:30].count('T')) / 30
    middle_at = (seq[30:60].count('A') + seq[30:60].count('T')) / 30
    downstream_at = (seq[60:].count('A') + seq[60:].count('T')) / 40

    total_dinucs = sum(dinucs.values())

    return {
        'at_content': at_content,
        'gc_content': gc_content,
        'max_at_run': max_at_run,
        'max_gc_run': max_gc_run,
        'upstream_at': upstream_at,
        'middle_at': middle_at,
        'downstream_at': downstream_at,
        'aa_freq': dinucs.get('AA', 0) / total_dinucs if total_dinucs else 0,
        'tt_freq': dinucs.get('TT', 0) / total_dinucs if total_dinucs else 0,
        'at_freq': dinucs.get('AT', 0) / total_dinucs if total_dinucs else 0,
        'ta_freq': dinucs.get('TA', 0) / total_dinucs if total_dinucs else 0,
        'cg_freq': dinucs.get('CG', 0) / total_dinucs if total_dinucs else 0,
    }


def main():
    print("=" * 80)
    print("ERROR ANALYSIS EXPERIMENT")
    print("=" * 80)

    # Load sequences
    seq_path = PROJECT_ROOT / "data/sequences/all_sequences.json"
    with open(seq_path) as f:
        sequences = json.load(f)

    if isinstance(sequences, list):
        sequences = {s['id']: s for s in sequences}

    # Load model results
    results_dir = PROJECT_ROOT / "data/results"
    model_scores = {}

    # Patterns to exclude from model results
    exclude_patterns = [
        'all_results', 'metrics', 'biophysical_comparison',
        'at_titration', 'positional_sweep', 'spacing_', 'strand_',
        'dinucleotide_control', 'negative_mes', 'error_analysis',
        'mpra_library',
    ]

    for result_file in results_dir.glob("*_results.json"):
        skip = False
        for pattern in exclude_patterns:
            if pattern in result_file.name:
                skip = True
                break
        if skip:
            continue
        model_name = result_file.stem.replace('_results', '')
        with open(result_file) as f:
            raw_scores = json.load(f)

        # Handle both formats:
        # Old format: per-sequence keys like {"A_000": -5.3, "A_001": -5.1, ...}
        # New format: per-class arrays like {"A": [-5.3, -5.1, ...], "B": [...], ...}
        first_key = next(iter(raw_scores))
        if isinstance(raw_scores[first_key], list):
            # New format: convert to per-sequence keys
            converted = {}
            for class_label, values in raw_scores.items():
                for idx, val in enumerate(values):
                    seq_id = f"{class_label}_{idx:03d}"
                    converted[seq_id] = val
            model_scores[model_name] = converted
        else:
            model_scores[model_name] = raw_scores

    if not model_scores:
        print("No model results found. Run inference first.")
        return

    print(f"\nModels found: {list(model_scores.keys())}")

    # Get D (broken) and E (compensated) sequences
    d_ids = sorted([k for k in sequences if sequences[k]['class_label'] == 'D'])
    e_ids = sorted([k for k in sequences if sequences[k]['class_label'] == 'E'])
    h_ids = sorted([k for k in sequences if sequences[k]['class_label'] == 'H'])

    n = min(len(d_ids), len(e_ids))
    d_ids = d_ids[:n]
    e_ids = e_ids[:n]

    analysis_results = {}

    for model_name, scores in model_scores.items():
        print(f"\n{'=' * 60}")
        print(f"MODEL: {model_name}")
        print(f"{'=' * 60}")

        # Check which pairs the model gets right (E > D)
        correct = []
        incorrect = []
        for d_id, e_id in zip(d_ids, e_ids):
            if d_id not in scores or e_id not in scores:
                continue
            d_score = scores[d_id]
            e_score = scores[e_id]

            d_seq = sequences[d_id]['sequence']
            e_seq = sequences[e_id]['sequence']

            d_props = compute_sequence_properties(d_seq)
            e_props = compute_sequence_properties(e_seq)

            pair_info = {
                'd_id': d_id, 'e_id': e_id,
                'd_score': d_score, 'e_score': e_score,
                'diff': e_score - d_score,
                'd_props': d_props, 'e_props': e_props,
                'at_diff': e_props['at_content'] - d_props['at_content'],
            }

            if e_score > d_score:
                correct.append(pair_info)
            else:
                incorrect.append(pair_info)

        total = len(correct) + len(incorrect)
        css = len(correct) / total if total > 0 else 0.5

        print(f"\n  CSS: {css:.3f} ({len(correct)}/{total} correct)")

        if not correct or not incorrect:
            print("  Insufficient data for error analysis.")
            continue

        # Analyze differences between correct and incorrect pairs
        print(f"\n  --- Property Comparison: Correct vs Incorrect ---")

        properties_to_check = [
            ('at_diff', 'AT% difference (E - D)'),
            ('e_props.at_content', 'AT% of compensated seq'),
            ('d_props.at_content', 'AT% of broken seq'),
            ('e_props.max_at_run', 'Max AT run in compensated'),
            ('e_props.upstream_at', 'Upstream AT% in compensated'),
            ('e_props.aa_freq', 'AA dinuc freq in compensated'),
        ]

        model_analysis = {
            'css': css,
            'n_correct': len(correct),
            'n_incorrect': len(incorrect),
            'property_tests': {},
        }

        for prop_key, prop_name in properties_to_check:
            if '.' in prop_key:
                parts = prop_key.split('.')
                correct_vals = [p[parts[0]][parts[1]] for p in correct]
                incorrect_vals = [p[parts[0]][parts[1]] for p in incorrect]
            else:
                correct_vals = [p[prop_key] for p in correct]
                incorrect_vals = [p[prop_key] for p in incorrect]

            mean_correct = np.mean(correct_vals)
            mean_incorrect = np.mean(incorrect_vals)

            # t-test for difference
            if len(correct_vals) > 1 and len(incorrect_vals) > 1:
                t_stat, p_val = stats.ttest_ind(correct_vals, incorrect_vals)
            else:
                t_stat, p_val = 0.0, 1.0

            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""

            print(f"  {prop_name:<35}: correct={mean_correct:.4f}, "
                  f"incorrect={mean_incorrect:.4f}, p={p_val:.4f} {sig}")

            model_analysis['property_tests'][prop_key] = {
                'name': prop_name,
                'mean_correct': float(mean_correct),
                'mean_incorrect': float(mean_incorrect),
                'p_value': float(p_val),
                't_stat': float(t_stat),
            }

        # Score difference distribution
        correct_diffs = [p['diff'] for p in correct]
        incorrect_diffs = [p['diff'] for p in incorrect]

        print(f"\n  --- Score Difference Distribution ---")
        print(f"  Correct pairs:   mean diff = {np.mean(correct_diffs):+.3f}, "
              f"std = {np.std(correct_diffs):.3f}")
        print(f"  Incorrect pairs: mean diff = {np.mean(incorrect_diffs):+.3f}, "
              f"std = {np.std(incorrect_diffs):.3f}")

        model_analysis['correct_diff_mean'] = float(np.mean(correct_diffs))
        model_analysis['incorrect_diff_mean'] = float(np.mean(incorrect_diffs))

        analysis_results[model_name] = model_analysis

    # Cross-model agreement analysis
    if len(model_scores) > 1:
        print(f"\n{'=' * 60}")
        print("CROSS-MODEL AGREEMENT")
        print(f"{'=' * 60}")

        model_names = list(model_scores.keys())
        # For each pair, check if correct/incorrect across models
        for i, m1 in enumerate(model_names):
            for j, m2 in enumerate(model_names):
                if i >= j:
                    continue
                both_correct = 0
                both_incorrect = 0
                disagree = 0
                total = 0
                for d_id, e_id in zip(d_ids, e_ids):
                    if (d_id not in model_scores[m1] or e_id not in model_scores[m1] or
                            d_id not in model_scores[m2] or e_id not in model_scores[m2]):
                        continue
                    m1_correct = model_scores[m1][e_id] > model_scores[m1][d_id]
                    m2_correct = model_scores[m2][e_id] > model_scores[m2][d_id]
                    total += 1
                    if m1_correct and m2_correct:
                        both_correct += 1
                    elif not m1_correct and not m2_correct:
                        both_incorrect += 1
                    else:
                        disagree += 1

                if total > 0:
                    agreement = (both_correct + both_incorrect) / total
                    print(f"  {m1} vs {m2}: agreement={agreement:.3f} "
                          f"(both correct={both_correct}, both wrong={both_incorrect}, "
                          f"disagree={disagree})")

        analysis_results['cross_model'] = {
            'description': 'See terminal output for cross-model agreement'
        }

    # Save results
    output_path = results_dir / "error_analysis_results.json"
    with open(output_path, 'w') as f:
        json.dump(analysis_results, f, indent=2)
    print(f"\nResults saved to {output_path}")

    print(f"\n{'=' * 60}")
    print("INTERPRETATION")
    print(f"{'=' * 60}")
    print("""
    If correct pairs have significantly higher AT% difference → AT-driven
    If correct pairs have longer AT runs → Run-length heuristic
    If no property differs significantly → Errors are stochastic
    If models agree on failures → Systematic failure mode
    If models disagree → Random/model-specific failures
    """)


if __name__ == "__main__":
    main()
