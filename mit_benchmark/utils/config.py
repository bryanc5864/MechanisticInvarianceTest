"""Configuration and constants for MIT benchmark."""

import os
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
SEQUENCES_DIR = DATA_DIR / "sequences"
RESULTS_DIR = DATA_DIR / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"

# Ensure directories exist
for dir_path in [SEQUENCES_DIR, RESULTS_DIR, FIGURES_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# HuggingFace configuration
HF_TOKEN = os.environ.get("HF_TOKEN", None)

# Sequence generation constants
SEQUENCE_LENGTH = 100  # Total sequence length in bp

# Promoter element positions (0-indexed)
# Standard layout: UP element (15-24), -35 box (30-35), spacer (36-52), -10 box (53-58)
UP_ELEMENT_START = 15
UP_ELEMENT_END = 24
MINUS_35_START = 30
MINUS_35_END = 36  # exclusive, 6bp
SPACER_START = 36
SPACER_END = 53  # exclusive, 17bp spacer
MINUS_10_START = 53
MINUS_10_END = 59  # exclusive, 6bp
EXTENDED_10_START = 50  # TGT just upstream of -10
EXTENDED_10_END = 53

# Consensus sequences
MINUS_35_CONSENSUS = "TTGACA"
MINUS_10_CONSENSUS = "TATAAT"
MINUS_10_BROKEN = "TGTAAT"  # T→G at position 2 (creates weaker -10)
UP_PROXIMAL_CONSENSUS = "AAAAAARNR"  # R=A/G, N=any (positions -46 to -38)
EXTENDED_10_CONSENSUS = "TGT"

# Spacer length (canonical is 17bp between -35 and -10)
SPACER_LENGTH = 17

# Number of sequences per class
CLASS_SIZES = {
    "A": 100,  # Natural intact promoters
    "B": 100,  # Natural broken promoters
    "C": 100,  # Synthetic intact promoters
    "D": 100,  # Synthetic broken promoters
    "E": 100,  # Synthetic compensated promoters
    "F": 50,   # Synthetic over-compensated promoters
    "G": 50,   # Natural compensated promoters
    "H": 50,   # Scrambled compensation control
}

# Model configurations
MODEL_CONFIGS = {
    "evo2_1b": {
        "type": "autoregressive",
        "model_id": "evo2_1b_base",
        "description": "Evo2 1B parameter autoregressive model",
    },
    "dnabert2": {
        "type": "masked_lm",
        "model_id": "zhihan1996/DNABERT-2-117M",
        "trust_remote_code": True,
        "description": "DNABERT-2 117M masked language model",
    },
    "nt_500m": {
        "type": "masked_lm",
        "model_id": "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
        "description": "Nucleotide Transformer 500M",
    },
    "hyenadna": {
        "type": "autoregressive",
        "model_id": "LongSafari/hyenadna-medium-160k-seqlen-hf",
        "trust_remote_code": True,
        "description": "HyenaDNA medium with 160k context",
    },
    "grover": {
        "type": "masked_lm",
        "model_id": "PoetschLab/GROVER",
        "description": "GROVER genomic language model",
    },
    "caduceus": {
        "type": "masked_lm",
        "model_id": "kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
        "trust_remote_code": True,
        "description": "Caduceus bidirectional DNA model",
    },
}

# Nucleotide alphabet
NUCLEOTIDES = "ACGT"

# Random seed for reproducibility
RANDOM_SEED = 42
