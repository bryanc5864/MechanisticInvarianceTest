# Mechanistic Invariance Test for Genomic Language Models

A benchmark for evaluating whether genomic language models (gLMs) encode mechanistic understanding of bacterial transcriptional regulation, or merely learn shallow statistical associations.

## Abstract

We introduce the **Mechanistic Invariance Test (MIT)**, a benchmark that probes whether gLMs understand regulatory compensation in E. coli σ70 promoters. The core principle: promoter function can be maintained through compensatory regulatory elements even when primary motifs are damaged. A model with true mechanistic understanding should recognize that *broken + compensated ≈ intact*.

We evaluate 5 gLMs spanning autoregressive (HyenaDNA, Evo2-1B), masked language model (NT-500M, GROVER), and bidirectional SSM (Caduceus) architectures. Our key finding: while HyenaDNA shows statistically significant compensation sensitivity (CSS=0.630, p=0.004), extended experiments reveal this signal is driven entirely by AT content correlation (r=0.78-0.96), not mechanistic understanding. All tested models fail diagnostic tests for positional encoding, spacing sensitivity, and strand orientation.

In contrast, biophysical baselines with ~100 parameters of explicit biological knowledge achieve CSS=1.000 and correctly encode all mechanistic constraints. This demonstrates that the MIT benchmark is solvable, and that current gLMs have not learned equivalent principles despite orders of magnitude more parameters.

## Key Results

| Model | CSS | SCR | Significant? | Interpretation |
|-------|-----|-----|--------------|----------------|
| PA-PWM (biophysical) | 1.000 | 0.980 | - | Solves MIT perfectly |
| RPA-PWM (relative) | 1.000 | 0.920 | - | Relative constraints suffice |
| HyenaDNA | 0.630 | 0.480 | Yes (p=0.004) | AT-driven, not mechanistic |
| Evo2-1B | 0.600 | 0.460 | No (FDR p=0.09) | AT-driven |
| NT-500M | 0.540 | 0.400 | No | Near random |
| GROVER | 0.520 | 0.520 | No | Near random |
| Caduceus | 0.490 | 0.420 | No | At chance |

**CSS** = Compensation Sensitivity Score (P(compensated > broken))
**SCR** = Scramble Control Ratio (P(structured > scrambled))

See **[RESULTS.md](RESULTS.md)** for comprehensive analysis including extended experiments on AT titration, positional sweeps, spacing sensitivity, and strand orientation.

## Biological Background

In E. coli σ70 promoters:
- The **-10 box** (consensus: TATAAT) is critical for σ70 recognition
- A **broken -10 box** (e.g., TGTAAT with T→G mutation) reduces promoter strength >100-fold
- **Compensatory elements** can restore function:
  - **UP element**: AT-rich region upstream of -35 (contacts α subunit)
  - **Extended -10**: TGT triplet immediately upstream of -10
  - **Optimal spacing**: 17±1 bp between -35 and -10

A model encoding mechanistic understanding should:
1. Recognize that compensated sequences score higher than broken sequences (CSS > 0.5)
2. Distinguish structured compensation from scrambled composition (SCR > 0.5)
3. Encode correct positional constraints (UP upstream of -35, spacing = 17bp)
4. Be strand-aware (promoter elements are directional)

## Sequence Classes

| Class | Name | N | Description |
|-------|------|---|-------------|
| A | Natural Intact | 100 | Real promoters with strong -10 box |
| B | Natural Broken | 100 | Real promoters with mutated -10, no compensation |
| C | Synthetic Intact | 100 | Consensus -35 (TTGACA) and -10 (TATAAT) |
| D | Synthetic Broken | 100 | Consensus -35, broken -10 (TGTAAT) |
| E | Synthetic Compensated | 100 | Broken -10 + UP element + extended -10 |
| F | Over-Compensated | 50 | Broken -10 + all compensatory elements |
| G | Natural Compensated | 50 | Real promoters with compensation |
| H | Scrambled Control | 50 | Same composition as E, scrambled motifs |

**Total: 650 sequences**

## Installation

### Requirements
- Python 3.9+
- PyTorch 2.0+ with CUDA support
- NVIDIA GPU with 8GB+ memory (for gLM inference)

### Setup

```bash
# Clone repository
git clone https://github.com/bryanc5864/MechanisticInvarianceTest.git
cd DOT-ICL

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install package
pip install -e .

# Or install dependencies directly
pip install -r requirements.txt
```

### Model-Specific Dependencies

**Evo2-1B** (requires separate environment due to transformer-engine):
```bash
pip install evo2 transformer-engine flash-attn
```

**Caduceus** (requires mamba-ssm):
```bash
pip install mamba-ssm==2.2.4 causal-conv1d==1.5.0.post8
```

## Quick Start

### Full Reproducibility (Single Command)

```bash
# Run complete experiment pipeline
python scripts/run_all_experiments.py --gpu=1

# Results saved to:
# - data/results/metrics.json (all metrics)
# - data/results/*.json (per-model scores)
# - figures/*.png (visualizations)
# - logs/<timestamp>/experiment.log (full log)
```

### Step-by-Step Reproduction

```bash
# 1. Generate benchmark sequences (650 sequences, 8 classes)
python scripts/generate_sequences.py --output data/sequences/

# 2. Run baseline models (CPU, fast)
python scripts/run_inference.py \
    --sequences data/sequences/all_sequences.json \
    --models kmer,pwm,random \
    --output data/results/

# 3. Run gLM inference (GPU required)
python scripts/run_inference.py \
    --sequences data/sequences/all_sequences.json \
    --models hyenadna,grover,nt_500m \
    --output data/results/ \
    --gpu 0

# 4. Compute metrics
python scripts/compute_metrics.py \
    --results data/results/ \
    --output data/results/metrics.json

# 5. Generate figures and analysis
python scripts/analyze_results.py \
    --metrics data/results/metrics.json \
    --output figures/
```

### Extended Experiments

```bash
# Biophysical model comparison (PA-PWM, RPA-PWM, Thermo)
python scripts/run_biophysical_comparison.py

# RPA-PWM: Fair baseline with relative constraints only
python scripts/rpa_pwm.py

# Diagnostic experiments
python scripts/experiment_at_titration.py --model hyenadna
python scripts/experiment_positional_sweep.py --model hyenadna
python scripts/experiment_spacing.py --model hyenadna
python scripts/experiment_strand.py --model hyenadna

# Extended validation experiments
python scripts/experiment_dinucleotide_control.py
python scripts/experiment_error_analysis.py
python scripts/experiment_negative_mes.py
python scripts/experiment_mpra_validation.py
```

## Models Evaluated

| Model | Type | Architecture | Scoring Method |
|-------|------|--------------|----------------|
| HyenaDNA-medium | Autoregressive | Hyena | Log-likelihood |
| Evo2-1B | Autoregressive | StripedHyena | Log-likelihood |
| NT-500M | Masked LM | Transformer | Pseudo-log-likelihood |
| GROVER | Masked LM | Transformer | Mean embedding |
| Caduceus | Bidirectional | Mamba SSM | Pseudo-log-likelihood |
| PA-PWM | Biophysical | PWM | Position-aware scoring |
| RPA-PWM | Biophysical | PWM | Relative position scoring |
| k-mer | Statistical | - | TF-IDF |

## Metrics

| Metric | Definition | Interpretation |
|--------|------------|----------------|
| **CSS** | P(LL_compensated > LL_broken) | >0.5 = recognizes compensation |
| **SCR** | P(LL_structured > LL_scrambled) | >0.5 = recognizes structure |
| **MES** | Cohen's d (intact vs broken) | Effect size of motif discrimination |
| **CM** | (LL_comp - LL_broken) / (LL_intact - LL_broken) | Fraction of recovery |

## Project Structure

```
DOT-ICL/
├── mit_benchmark/
│   ├── sequences/
│   │   ├── generator.py      # Sequence generation for all 8 classes
│   │   ├── motifs.py         # σ70 promoter motif definitions
│   │   └── natural.py        # Natural promoter data from RegulonDB
│   ├── models/
│   │   ├── base.py           # Abstract model interface
│   │   ├── autoregressive.py # HyenaDNA, Evo2 wrappers
│   │   ├── masked_lm.py      # NT, GROVER, Caduceus wrappers
│   │   ├── biophysical.py    # PA-PWM, Thermodynamic models
│   │   └── baselines.py      # k-mer, PWM, random baselines
│   └── evaluation/
│       ├── metrics.py        # CSS, MES, SCR, CM implementations
│       └── analysis.py       # Statistical tests, visualization
├── scripts/
│   ├── run_all_experiments.py    # Full pipeline
│   ├── run_inference.py          # Model inference
│   ├── run_biophysical_comparison.py
│   ├── rpa_pwm.py                # Fair baseline
│   └── experiment_*.py           # Extended experiments
├── data/
│   ├── sequences/            # Generated benchmark sequences
│   └── results/              # Model predictions and metrics
├── figures/                  # Generated visualizations
├── RESULTS.md               # Comprehensive results documentation
└── requirements.txt
```

## Citation

```bibtex
@article{mit_benchmark_2026,
  title={Mechanistic Invariance Test for Genomic Language Models},
  author={Cheng, Bryan and Zhang, Jasper},
  year={2026}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.
