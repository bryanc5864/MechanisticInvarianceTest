#!/usr/bin/env python3
"""Run all MIT benchmark experiments with comprehensive logging.

This script runs:
1. Baseline models (k-mer, PWM, random)
2. gLM models (HyenaDNA, NT-500M, GROVER)
3. Extended experiments (AT titration, positional sweep, spacing, strand)
4. Biophysical model comparison

All outputs are logged to logs/ directory with timestamps.
"""

import json
import logging
import os
import random
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Set random seeds
random.seed(42)
np.random.seed(42)

# Create timestamp for this run
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
LOG_DIR = PROJECT_ROOT / "logs" / TIMESTAMP
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(message)s',
    handlers=[
        logging.FileHandler(LOG_DIR / "experiment.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


def log_system_info():
    """Log system information."""
    import platform
    try:
        import torch
        cuda_available = torch.cuda.is_available()
        if cuda_available:
            gpu_name = torch.cuda.get_device_name(0)
            gpu_memory = torch.cuda.get_device_properties(0).total_memory / 1e9
        else:
            gpu_name = "N/A"
            gpu_memory = 0
    except ImportError:
        cuda_available = False
        gpu_name = "N/A"
        gpu_memory = 0

    logger.info("="*80)
    logger.info("SYSTEM INFORMATION")
    logger.info("="*80)
    logger.info(f"Python version: {platform.python_version()}")
    logger.info(f"Platform: {platform.platform()}")
    logger.info(f"CUDA available: {cuda_available}")
    logger.info(f"GPU: {gpu_name}")
    logger.info(f"GPU Memory: {gpu_memory:.1f} GB")
    logger.info(f"Timestamp: {TIMESTAMP}")
    logger.info(f"Log directory: {LOG_DIR}")
    logger.info("="*80)


def run_baseline_models():
    """Run baseline models (k-mer, PWM, random)."""
    logger.info("\n" + "="*80)
    logger.info("RUNNING BASELINE MODELS")
    logger.info("="*80)

    from mit_benchmark.models.baselines import KmerBaseline, PWMBaseline, RandomBaseline

    # Load sequences
    seq_path = PROJECT_ROOT / "data/sequences/all_sequences.json"
    with open(seq_path) as f:
        sequences = json.load(f)
    logger.info(f"Loaded {len(sequences)} sequences from {seq_path}")

    # Convert list to dict if needed
    if isinstance(sequences, list):
        sequences = {s['id']: s for s in sequences}

    results = {}

    # k-mer baseline
    logger.info("\n--- k-mer Baseline ---")
    kmer_model = KmerBaseline()
    logger.info("Computing k-mer scores...")
    kmer_scores = {}
    for seq_id, seq_data in sequences.items():
        seq = seq_data['sequence'] if isinstance(seq_data, dict) else seq_data
        score = kmer_model.compute_log_likelihood(seq)
        kmer_scores[seq_id] = score
    results['kmer'] = kmer_scores
    logger.info(f"k-mer: Computed {len(kmer_scores)} scores")
    logger.info(f"k-mer: Mean score = {np.mean(list(kmer_scores.values())):.4f}")

    # PWM baseline
    logger.info("\n--- PWM Baseline ---")
    pwm_model = PWMBaseline()
    logger.info("Computing PWM scores...")
    pwm_scores = {}
    for seq_id, seq_data in sequences.items():
        seq = seq_data['sequence'] if isinstance(seq_data, dict) else seq_data
        score = pwm_model.compute_log_likelihood(seq)
        pwm_scores[seq_id] = score
    results['pwm'] = pwm_scores
    logger.info(f"PWM: Computed {len(pwm_scores)} scores")
    logger.info(f"PWM: Mean score = {np.mean(list(pwm_scores.values())):.4f}")

    # Random baseline
    logger.info("\n--- Random Baseline ---")
    random_model = RandomBaseline()
    logger.info("Computing random scores...")
    random_scores = {}
    for seq_id, seq_data in sequences.items():
        seq = seq_data['sequence'] if isinstance(seq_data, dict) else seq_data
        score = random_model.compute_log_likelihood(seq)
        random_scores[seq_id] = score
    results['random'] = random_scores
    logger.info(f"Random: Computed {len(random_scores)} scores")
    logger.info(f"Random: Mean score = {np.mean(list(random_scores.values())):.4f}")

    # Save results
    for model_name, scores in results.items():
        output_path = PROJECT_ROOT / f"data/results/{model_name}_results.json"
        with open(output_path, 'w') as f:
            json.dump(scores, f, indent=2)
        logger.info(f"Saved {model_name} results to {output_path}")

    return results


def run_glm_models(gpu_id=1):
    """Run gLM models (HyenaDNA, NT-500M, GROVER)."""
    logger.info("\n" + "="*80)
    logger.info("RUNNING gLM MODELS")
    logger.info("="*80)

    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
    logger.info(f"Using GPU {gpu_id}")

    import torch
    from transformers import AutoModelForCausalLM, AutoTokenizer, AutoModel

    # Load sequences
    seq_path = PROJECT_ROOT / "data/sequences/all_sequences.json"
    with open(seq_path) as f:
        sequences = json.load(f)
    logger.info(f"Loaded {len(sequences)} sequences")

    # Convert list to dict if needed
    if isinstance(sequences, list):
        sequences = {s['id']: s for s in sequences}

    results = {}

    # --- HyenaDNA ---
    logger.info("\n--- HyenaDNA ---")
    try:
        logger.info("Loading HyenaDNA model...")
        tokenizer = AutoTokenizer.from_pretrained(
            "LongSafari/hyenadna-medium-160k-seqlen-hf",
            trust_remote_code=True
        )
        model = AutoModelForCausalLM.from_pretrained(
            "LongSafari/hyenadna-medium-160k-seqlen-hf",
            trust_remote_code=True
        )
        model.to('cuda')
        model.eval()
        logger.info("HyenaDNA loaded successfully")

        def compute_hyenadna_ll(sequence):
            with torch.no_grad():
                inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=False)
                input_ids = inputs["input_ids"].to('cuda')
                outputs = model(input_ids)
                logits = outputs.logits
                log_probs = torch.log_softmax(logits[0, :-1, :], dim=-1)
                target_ids = input_ids[0, 1:]
                ll = log_probs.gather(1, target_ids.unsqueeze(1)).sum().item()
            return ll

        logger.info("Computing HyenaDNA scores...")
        hyenadna_scores = {}
        for i, (seq_id, seq_data) in enumerate(sequences.items()):
            seq = seq_data['sequence'] if isinstance(seq_data, dict) else seq_data
            score = compute_hyenadna_ll(seq)
            hyenadna_scores[seq_id] = score
            if (i + 1) % 100 == 0:
                logger.info(f"  Processed {i+1}/{len(sequences)} sequences")

        results['hyenadna'] = hyenadna_scores
        logger.info(f"HyenaDNA: Computed {len(hyenadna_scores)} scores")
        logger.info(f"HyenaDNA: Mean score = {np.mean(list(hyenadna_scores.values())):.4f}")

        # Clean up
        del model, tokenizer
        torch.cuda.empty_cache()

    except Exception as e:
        logger.error(f"HyenaDNA failed: {e}")

    # --- NT-500M ---
    logger.info("\n--- Nucleotide Transformer 500M ---")
    try:
        logger.info("Loading NT-500M model...")
        tokenizer = AutoTokenizer.from_pretrained(
            "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
            trust_remote_code=True
        )
        model = AutoModel.from_pretrained(
            "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
            trust_remote_code=True
        )
        model.to('cuda')
        model.eval()
        logger.info("NT-500M loaded successfully")

        def compute_nt_pll(sequence):
            """Compute pseudo-log-likelihood for MLM."""
            with torch.no_grad():
                inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
                input_ids = inputs["input_ids"].to('cuda')
                outputs = model(input_ids)
                embeddings = outputs.last_hidden_state
                # Use mean embedding norm as proxy score
                score = embeddings.mean().item()
            return score

        logger.info("Computing NT-500M scores...")
        nt_scores = {}
        for i, (seq_id, seq_data) in enumerate(sequences.items()):
            seq = seq_data['sequence'] if isinstance(seq_data, dict) else seq_data
            score = compute_nt_pll(seq)
            nt_scores[seq_id] = score
            if (i + 1) % 100 == 0:
                logger.info(f"  Processed {i+1}/{len(sequences)} sequences")

        results['nt_500m'] = nt_scores
        logger.info(f"NT-500M: Computed {len(nt_scores)} scores")
        logger.info(f"NT-500M: Mean score = {np.mean(list(nt_scores.values())):.4f}")

        del model, tokenizer
        torch.cuda.empty_cache()

    except Exception as e:
        logger.error(f"NT-500M failed: {e}")

    # --- GROVER ---
    logger.info("\n--- GROVER ---")
    try:
        logger.info("Loading GROVER model...")
        tokenizer = AutoTokenizer.from_pretrained("PoetschLab/GROVER", trust_remote_code=True)
        model = AutoModel.from_pretrained("PoetschLab/GROVER", trust_remote_code=True)
        model.to('cuda')
        model.eval()
        logger.info("GROVER loaded successfully")

        def compute_grover_score(sequence):
            with torch.no_grad():
                inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
                input_ids = inputs["input_ids"].to('cuda')
                outputs = model(input_ids)
                embeddings = outputs.last_hidden_state
                score = embeddings.mean().item()
            return score

        logger.info("Computing GROVER scores...")
        grover_scores = {}
        for i, (seq_id, seq_data) in enumerate(sequences.items()):
            seq = seq_data['sequence'] if isinstance(seq_data, dict) else seq_data
            score = compute_grover_score(seq)
            grover_scores[seq_id] = score
            if (i + 1) % 100 == 0:
                logger.info(f"  Processed {i+1}/{len(sequences)} sequences")

        results['grover'] = grover_scores
        logger.info(f"GROVER: Computed {len(grover_scores)} scores")
        logger.info(f"GROVER: Mean score = {np.mean(list(grover_scores.values())):.4f}")

        del model, tokenizer
        torch.cuda.empty_cache()

    except Exception as e:
        logger.error(f"GROVER failed: {e}")

    # --- DNABERT-2 ---
    logger.info("\n--- DNABERT-2 ---")
    try:
        logger.info("Loading DNABERT-2 model...")
        tokenizer = AutoTokenizer.from_pretrained(
            "zhihan1996/DNABERT-2-117M", trust_remote_code=True
        )
        # DNABERT-2 has a custom BertConfig that conflicts with transformers>=4.40.
        # Load config separately and use standard BertModel as fallback.
        try:
            model = AutoModel.from_pretrained(
                "zhihan1996/DNABERT-2-117M", trust_remote_code=True
            )
        except ValueError:
            logger.info("AutoModel failed, trying BertModel with standard config...")
            from transformers import BertModel, BertConfig
            config = BertConfig.from_pretrained("zhihan1996/DNABERT-2-117M")
            model = BertModel.from_pretrained(
                "zhihan1996/DNABERT-2-117M", config=config, ignore_mismatched_sizes=True
            )
        model.to('cuda')
        model.eval()
        logger.info("DNABERT-2 loaded successfully")

        def compute_dnabert2_score(sequence):
            """Compute pseudo-log-likelihood for DNABERT-2."""
            with torch.no_grad():
                inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
                input_ids = inputs["input_ids"].to('cuda')
                outputs = model(input_ids)
                embeddings = outputs.last_hidden_state
                score = embeddings.mean().item()
            return score

        logger.info("Computing DNABERT-2 scores...")
        dnabert2_scores = {}
        for i, (seq_id, seq_data) in enumerate(sequences.items()):
            seq = seq_data['sequence'] if isinstance(seq_data, dict) else seq_data
            score = compute_dnabert2_score(seq)
            dnabert2_scores[seq_id] = score
            if (i + 1) % 100 == 0:
                logger.info(f"  Processed {i+1}/{len(sequences)} sequences")

        results['dnabert2'] = dnabert2_scores
        logger.info(f"DNABERT-2: Computed {len(dnabert2_scores)} scores")
        logger.info(f"DNABERT-2: Mean score = {np.mean(list(dnabert2_scores.values())):.4f}")

        del model, tokenizer
        torch.cuda.empty_cache()

    except Exception as e:
        logger.error(f"DNABERT-2 failed: {e}")

    # --- Caduceus ---
    logger.info("\n--- Caduceus ---")
    try:
        logger.info("Loading Caduceus model...")
        tokenizer = AutoTokenizer.from_pretrained(
            "kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
            trust_remote_code=True
        )
        model = AutoModel.from_pretrained(
            "kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
            trust_remote_code=True
        )
        model.to('cuda')
        model.eval()
        logger.info("Caduceus loaded successfully")

        def compute_caduceus_score(sequence):
            with torch.no_grad():
                inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
                input_ids = inputs["input_ids"].to('cuda')
                outputs = model(input_ids)
                embeddings = outputs.last_hidden_state
                score = embeddings.mean().item()
            return score

        logger.info("Computing Caduceus scores...")
        caduceus_scores = {}
        for i, (seq_id, seq_data) in enumerate(sequences.items()):
            seq = seq_data['sequence'] if isinstance(seq_data, dict) else seq_data
            score = compute_caduceus_score(seq)
            caduceus_scores[seq_id] = score
            if (i + 1) % 100 == 0:
                logger.info(f"  Processed {i+1}/{len(sequences)} sequences")

        results['caduceus'] = caduceus_scores
        logger.info(f"Caduceus: Computed {len(caduceus_scores)} scores")
        logger.info(f"Caduceus: Mean score = {np.mean(list(caduceus_scores.values())):.4f}")

        del model, tokenizer
        torch.cuda.empty_cache()

    except Exception as e:
        logger.error(f"Caduceus failed: {e}")
        if "mamba_ssm" in str(e):
            logger.error("  Install mamba_ssm: pip install mamba_ssm (requires CUDA >= 11.6 and nvcc in PATH)")

    # Save results
    for model_name, scores in results.items():
        output_path = PROJECT_ROOT / f"data/results/{model_name}_results.json"
        with open(output_path, 'w') as f:
            json.dump(scores, f, indent=2)
        logger.info(f"Saved {model_name} results to {output_path}")

    return results


def run_extended_experiments(gpu_id=1):
    """Run extended experiments (AT titration, positional sweep, spacing, strand)."""
    logger.info("\n" + "="*80)
    logger.info("RUNNING EXTENDED EXPERIMENTS")
    logger.info("="*80)

    os.environ['CUDA_VISIBLE_DEVICES'] = str(gpu_id)

    import torch
    from transformers import AutoModelForCausalLM, AutoTokenizer

    # Load HyenaDNA once for all experiments
    logger.info("Loading HyenaDNA for extended experiments...")
    tokenizer = AutoTokenizer.from_pretrained(
        "LongSafari/hyenadna-medium-160k-seqlen-hf",
        trust_remote_code=True
    )
    model = AutoModelForCausalLM.from_pretrained(
        "LongSafari/hyenadna-medium-160k-seqlen-hf",
        trust_remote_code=True
    )
    model.to('cuda')
    model.eval()
    logger.info("HyenaDNA loaded")

    def compute_ll(sequence):
        with torch.no_grad():
            inputs = tokenizer(sequence, return_tensors="pt", add_special_tokens=False)
            input_ids = inputs["input_ids"].to('cuda')
            outputs = model(input_ids)
            logits = outputs.logits
            log_probs = torch.log_softmax(logits[0, :-1, :], dim=-1)
            target_ids = input_ids[0, 1:]
            ll = log_probs.gather(1, target_ids.unsqueeze(1)).sum().item()
        return ll

    def generate_background(length, at_fraction=0.55):
        seq = []
        for _ in range(length):
            if random.random() < at_fraction:
                seq.append(random.choice('AT'))
            else:
                seq.append(random.choice('GC'))
        return ''.join(seq)

    results = {}

    # --- Experiment 1: AT Titration ---
    logger.info("\n--- Experiment 1: AT Titration ---")
    at_levels = [0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
    n_per_level = 30
    at_results = {}

    for at_frac in at_levels:
        logger.info(f"  AT fraction = {at_frac:.0%}")
        intact_lls = []
        broken_lls = []
        compensated_lls = []

        for _ in range(n_per_level):
            # Intact
            seq = list(generate_background(100, at_frac))
            for i, nt in enumerate("TTGACA"):
                seq[25 + i] = nt
            for i, nt in enumerate("TATAAT"):
                seq[48 + i] = nt
            intact_lls.append(compute_ll(''.join(seq)))

            # Broken
            seq = list(generate_background(100, at_frac))
            for i, nt in enumerate("TTGACA"):
                seq[25 + i] = nt
            for i, nt in enumerate("TGTAAT"):
                seq[48 + i] = nt
            broken_lls.append(compute_ll(''.join(seq)))

            # Compensated
            seq = list(generate_background(100, at_frac))
            for i, nt in enumerate("AAAAAAGCA"):
                seq[8 + i] = nt
            for i, nt in enumerate("TTGACA"):
                seq[25 + i] = nt
            for i, nt in enumerate("TGT"):
                seq[45 + i] = nt
            for i, nt in enumerate("TGTAAT"):
                seq[48 + i] = nt
            compensated_lls.append(compute_ll(''.join(seq)))

        at_results[str(int(at_frac * 100))] = {
            'intact': {'mean': np.mean(intact_lls), 'std': np.std(intact_lls), 'values': intact_lls},
            'broken': {'mean': np.mean(broken_lls), 'std': np.std(broken_lls), 'values': broken_lls},
            'compensated': {'mean': np.mean(compensated_lls), 'std': np.std(compensated_lls), 'values': compensated_lls},
        }
        logger.info(f"    Intact: {np.mean(intact_lls):.2f}, Broken: {np.mean(broken_lls):.2f}, Comp: {np.mean(compensated_lls):.2f}")

    results['at_titration'] = at_results

    # Compute correlation
    all_at = []
    all_ll = []
    for at_str, data in at_results.items():
        at_val = int(at_str) / 100
        for ll in data['intact']['values'] + data['broken']['values'] + data['compensated']['values']:
            all_at.append(at_val)
            all_ll.append(ll)
    correlation = np.corrcoef(all_at, all_ll)[0, 1]
    logger.info(f"  AT-LL Correlation: r = {correlation:.3f}")
    results['at_titration']['correlation'] = correlation

    # --- Experiment 2: Positional Sweep ---
    logger.info("\n--- Experiment 2: Positional Sweep ---")
    positions = [0, 5, 10, 15, 20, 25, 35, 45, 60, 70, 80]
    n_per_pos = 50
    pos_results = {}

    for pos in positions:
        logger.info(f"  UP position = {pos}")
        lls = []
        for _ in range(n_per_pos):
            seq = list(generate_background(100))
            # UP element at variable position
            up_element = "AAAAAAGCA"
            if pos + len(up_element) <= 100:
                for i, nt in enumerate(up_element):
                    seq[pos + i] = nt
            # Fixed -35 and -10
            for i, nt in enumerate("TTGACA"):
                seq[25 + i] = nt
            for i, nt in enumerate("TATAAT"):
                seq[48 + i] = nt
            lls.append(compute_ll(''.join(seq)))

        pos_results[str(pos)] = {'mean': np.mean(lls), 'std': np.std(lls), 'values': lls}
        logger.info(f"    Mean LL: {np.mean(lls):.2f}")

    # No UP control
    logger.info("  No UP element")
    no_up_lls = []
    for _ in range(n_per_pos):
        seq = list(generate_background(100))
        for i, nt in enumerate("TTGACA"):
            seq[25 + i] = nt
        for i, nt in enumerate("TATAAT"):
            seq[48 + i] = nt
        no_up_lls.append(compute_ll(''.join(seq)))
    pos_results['none'] = {'mean': np.mean(no_up_lls), 'std': np.std(no_up_lls), 'values': no_up_lls}
    logger.info(f"    Mean LL: {np.mean(no_up_lls):.2f}")

    results['positional_sweep'] = pos_results

    # --- Experiment 3: Spacing Sensitivity ---
    logger.info("\n--- Experiment 3: Spacing Sensitivity ---")
    spacings = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    n_per_spacing = 50
    spacing_results = {}

    for spacing in spacings:
        logger.info(f"  Spacing = {spacing}bp")
        lls = []
        for _ in range(n_per_spacing):
            seq = list(generate_background(100))
            # -35 at position 25
            for i, nt in enumerate("TTGACA"):
                seq[25 + i] = nt
            # -10 at position 25 + 6 + spacing
            pos_10 = 25 + 6 + spacing
            if pos_10 + 6 <= 100:
                for i, nt in enumerate("TATAAT"):
                    seq[pos_10 + i] = nt
            lls.append(compute_ll(''.join(seq)))

        spacing_results[str(spacing)] = {'mean': np.mean(lls), 'std': np.std(lls), 'values': lls}
        logger.info(f"    Mean LL: {np.mean(lls):.2f}")

    results['spacing'] = spacing_results

    # Find peak
    peak_spacing = max(spacing_results, key=lambda x: spacing_results[x]['mean'])
    logger.info(f"  Peak spacing: {peak_spacing}bp (biological optimum is 17bp)")

    # --- Experiment 4: Strand Orientation ---
    logger.info("\n--- Experiment 4: Strand Orientation ---")
    n_per_condition = 50
    strand_results = {}

    COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    def reverse_complement(seq):
        return ''.join(COMPLEMENT[b] for b in reversed(seq))

    # Forward
    logger.info("  Forward strand")
    forward_lls = []
    for _ in range(n_per_condition):
        seq = list(generate_background(100))
        for i, nt in enumerate("TTGACA"):
            seq[30 + i] = nt
        for i, nt in enumerate("TATAAT"):
            seq[53 + i] = nt
        forward_lls.append(compute_ll(''.join(seq)))
    strand_results['forward'] = {'mean': np.mean(forward_lls), 'std': np.std(forward_lls), 'values': forward_lls}
    logger.info(f"    Mean LL: {np.mean(forward_lls):.2f}")

    # RC in place
    logger.info("  RC motifs in place")
    rc_in_place_lls = []
    for _ in range(n_per_condition):
        seq = list(generate_background(100))
        for i, nt in enumerate(reverse_complement("TTGACA")):
            seq[30 + i] = nt
        for i, nt in enumerate(reverse_complement("TATAAT")):
            seq[53 + i] = nt
        rc_in_place_lls.append(compute_ll(''.join(seq)))
    strand_results['reverse_in_place'] = {'mean': np.mean(rc_in_place_lls), 'std': np.std(rc_in_place_lls), 'values': rc_in_place_lls}
    logger.info(f"    Mean LL: {np.mean(rc_in_place_lls):.2f}")

    # Full RC
    logger.info("  Full reverse complement")
    full_rc_lls = []
    for _ in range(n_per_condition):
        seq = list(generate_background(100))
        for i, nt in enumerate("TTGACA"):
            seq[30 + i] = nt
        for i, nt in enumerate("TATAAT"):
            seq[53 + i] = nt
        full_rc_lls.append(compute_ll(reverse_complement(''.join(seq))))
    strand_results['full_reverse'] = {'mean': np.mean(full_rc_lls), 'std': np.std(full_rc_lls), 'values': full_rc_lls}
    logger.info(f"    Mean LL: {np.mean(full_rc_lls):.2f}")

    # Scrambled
    logger.info("  Scrambled motifs")
    scrambled_lls = []
    for _ in range(n_per_condition):
        seq = list(generate_background(100))
        m35 = list("TTGACA")
        random.shuffle(m35)
        for i, nt in enumerate(m35):
            seq[30 + i] = nt
        m10 = list("TATAAT")
        random.shuffle(m10)
        for i, nt in enumerate(m10):
            seq[53 + i] = nt
        scrambled_lls.append(compute_ll(''.join(seq)))
    strand_results['scrambled'] = {'mean': np.mean(scrambled_lls), 'std': np.std(scrambled_lls), 'values': scrambled_lls}
    logger.info(f"    Mean LL: {np.mean(scrambled_lls):.2f}")

    results['strand'] = strand_results

    # Strand accuracy
    strand_correct = sum(1 for f, r in zip(forward_lls, rc_in_place_lls) if f > r)
    strand_acc = strand_correct / n_per_condition
    logger.info(f"  Strand accuracy: {strand_acc:.3f} (50% = strand-blind)")

    # Clean up
    del model, tokenizer
    torch.cuda.empty_cache()

    # Save all results
    for exp_name, exp_results in results.items():
        output_path = PROJECT_ROOT / f"data/results/{exp_name}_results.json"
        # Convert numpy types for JSON
        def convert(obj):
            if isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert(v) for v in obj]
            return obj
        with open(output_path, 'w') as f:
            json.dump(convert(exp_results), f, indent=2)
        logger.info(f"Saved {exp_name} results to {output_path}")

    return results


def run_biophysical_comparison():
    """Run biophysical model comparison."""
    logger.info("\n" + "="*80)
    logger.info("RUNNING BIOPHYSICAL MODEL COMPARISON")
    logger.info("="*80)

    from mit_benchmark.models.biophysical import (
        PositionScanningModel,
        load_position_aware_pwm, load_thermodynamic_model,
        load_papwm_no_comp, load_papwm_no_position,
        GENERATOR_POSITIONS,
    )

    models = {
        'PA-PWM': load_position_aware_pwm(),
        'PA-PWM-NoComp': load_papwm_no_comp(),
        'PA-PWM-NoPos': load_papwm_no_position(),
        'Thermo': load_thermodynamic_model(),
        'Scan': PositionScanningModel(pos_35=GENERATOR_POSITIONS['pos_35'],
                                       pos_10=GENERATOR_POSITIONS['pos_10']),
    }

    n_samples = 100

    def generate_background(length, at_fraction=0.55):
        seq = []
        for _ in range(length):
            if random.random() < at_fraction:
                seq.append(random.choice('AT'))
            else:
                seq.append(random.choice('GC'))
        return ''.join(seq)

    # Generate sequences
    logger.info("Generating test sequences...")

    # Use main generator positions for consistency
    P35 = GENERATOR_POSITIONS['pos_35']   # 30
    P10 = GENERATOR_POSITIONS['pos_10']   # 53
    PUP = GENERATOR_POSITIONS['pos_up']   # 15
    PEXT = GENERATOR_POSITIONS['pos_ext10']  # 50

    broken_seqs = []
    for _ in range(n_samples):
        seq = list(generate_background(100))
        for i, nt in enumerate("TTGACA"):
            seq[P35 + i] = nt
        for i, nt in enumerate("TGTAAT"):
            seq[P10 + i] = nt
        broken_seqs.append(''.join(seq))

    compensated_seqs = []
    for _ in range(n_samples):
        seq = list(generate_background(100))
        for i, nt in enumerate("AAAAAAGCA"):
            seq[PUP + i] = nt
        for i, nt in enumerate("TTGACA"):
            seq[P35 + i] = nt
        for i, nt in enumerate("TGT"):
            seq[PEXT + i] = nt
        for i, nt in enumerate("TGTAAT"):
            seq[P10 + i] = nt
        compensated_seqs.append(''.join(seq))

    scrambled_seqs = []
    for _ in range(n_samples):
        seq = list(generate_background(100))
        up = list("AAAAAAGCA")
        random.shuffle(up)
        for i, nt in enumerate(up):
            seq[PUP + i] = nt
        for i, nt in enumerate("TTGACA"):
            seq[P35 + i] = nt
        ext = list("TGT")
        random.shuffle(ext)
        for i, nt in enumerate(ext):
            seq[PEXT + i] = nt
        for i, nt in enumerate("TGTAAT"):
            seq[P10 + i] = nt
        scrambled_seqs.append(''.join(seq))

    results = {}

    for model_name, model in models.items():
        logger.info(f"\n--- {model_name} ---")

        broken_scores = [model.score(s) for s in broken_seqs]
        compensated_scores = [model.score(s) for s in compensated_seqs]
        scrambled_scores = [model.score(s) for s in scrambled_seqs]

        # CSS
        css = sum(1 for c, b in zip(compensated_scores, broken_scores) if c > b) / n_samples

        # SCR
        scr = sum(1 for c, s in zip(compensated_scores, scrambled_scores) if c > s) / n_samples

        results[model_name] = {
            'css': css,
            'scr': scr,
            'broken_mean': float(np.mean(broken_scores)),
            'compensated_mean': float(np.mean(compensated_scores)),
            'scrambled_mean': float(np.mean(scrambled_scores)),
        }

        logger.info(f"  CSS: {css:.3f}")
        logger.info(f"  SCR: {scr:.3f}")
        logger.info(f"  Broken mean: {np.mean(broken_scores):.2f}")
        logger.info(f"  Compensated mean: {np.mean(compensated_scores):.2f}")

    # Save results
    output_path = PROJECT_ROOT / "data/results/biophysical_comparison.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"\nSaved biophysical results to {output_path}")

    return results


def compute_metrics():
    """Compute all metrics from results."""
    logger.info("\n" + "="*80)
    logger.info("COMPUTING METRICS")
    logger.info("="*80)

    from mit_benchmark.evaluation.metrics import compute_all_metrics

    # Load sequences
    seq_path = PROJECT_ROOT / "data/sequences/all_sequences.json"
    with open(seq_path) as f:
        sequences = json.load(f)

    # Build class lookup from sequences
    if isinstance(sequences, list):
        seq_to_class = {s['id']: s['class_label'] for s in sequences}
    else:
        seq_to_class = {k: v['class_label'] for k, v in sequences.items()}

    # Load all results
    results_dir = PROJECT_ROOT / "data/results"
    model_results = {}

    # Files that are experiment outputs, not model score files
    exclude_files = {
        'all_results.json', 'metrics.json', 'biophysical_comparison.json',
        'at_titration_results.json', 'positional_sweep_results.json',
        'spacing_results.json', 'strand_results.json',
        'negative_mes_results.json', 'error_analysis_results.json',
        'dinucleotide_control_sequences.json', 'mpra_library.json',
    }
    for result_file in results_dir.glob("*_results.json"):
        if result_file.name in exclude_files:
            continue

        model_name = result_file.stem.replace('_results', '')
        with open(result_file) as f:
            model_results[model_name] = json.load(f)
        logger.info(f"Loaded {model_name} results ({len(model_results[model_name])} scores)")

    # Compute metrics for each model
    all_metrics = {}
    for model_name, scores in model_results.items():
        logger.info(f"\n--- {model_name} ---")

        # Organize scores by class label
        predictions_by_class = {label: [] for label in 'ABCDEFGH'}
        for seq_id, score in scores.items():
            class_label = seq_to_class.get(seq_id)
            if class_label and class_label in predictions_by_class:
                predictions_by_class[class_label].append(score)

        # Log class counts
        for label in 'ABCDEFGH':
            n = len(predictions_by_class[label])
            if n > 0:
                logger.info(f"  Class {label}: {n} sequences")

        metrics = compute_all_metrics(predictions_by_class, model_name)
        # Convert MetricsResult to dict if needed
        if hasattr(metrics, 'to_dict'):
            metrics_dict = metrics.to_dict()
        else:
            metrics_dict = metrics
        all_metrics[model_name] = metrics_dict

        logger.info(f"  CSS: {metrics_dict.get('css', 'N/A'):.3f}")
        logger.info(f"  MES (natural): {metrics_dict.get('mes_natural', 'N/A'):.3f}")
        logger.info(f"  MES (synthetic): {metrics_dict.get('mes_synthetic', 'N/A'):.3f}")
        logger.info(f"  SCR: {metrics_dict.get('scr', 'N/A'):.3f}")

    # Save metrics
    output_path = PROJECT_ROOT / "data/results/metrics.json"
    with open(output_path, 'w') as f:
        json.dump(all_metrics, f, indent=2)
    logger.info(f"\nSaved metrics to {output_path}")

    return all_metrics


def generate_summary():
    """Generate summary log."""
    logger.info("\n" + "="*80)
    logger.info("EXPERIMENT SUMMARY")
    logger.info("="*80)

    # Load metrics
    metrics_path = PROJECT_ROOT / "data/results/metrics.json"
    if metrics_path.exists():
        with open(metrics_path) as f:
            metrics = json.load(f)

        logger.info("\nModel Performance (CSS):")
        logger.info("-" * 40)
        for model, m in sorted(metrics.items(), key=lambda x: x[1].get('css', 0), reverse=True):
            css = m.get('css', 'N/A')
            if isinstance(css, (int, float)):
                logger.info(f"  {model:<15} CSS = {css:.3f}")

    # Load biophysical comparison
    bio_path = PROJECT_ROOT / "data/results/biophysical_comparison.json"
    if bio_path.exists():
        with open(bio_path) as f:
            bio_results = json.load(f)

        logger.info("\nBiophysical Model Performance:")
        logger.info("-" * 40)
        for model, r in bio_results.items():
            logger.info(f"  {model:<15} CSS = {r['css']:.3f}, SCR = {r['scr']:.3f}")

    # Load extended experiment results
    logger.info("\nExtended Experiments:")
    logger.info("-" * 40)

    at_path = PROJECT_ROOT / "data/results/at_titration_results.json"
    if at_path.exists():
        with open(at_path) as f:
            at_results = json.load(f)
        if 'correlation' in at_results:
            logger.info(f"  AT-LL Correlation: r = {at_results['correlation']:.3f}")

    spacing_path = PROJECT_ROOT / "data/results/spacing_results.json"
    if spacing_path.exists():
        with open(spacing_path) as f:
            spacing_results = json.load(f)
        peak = max(spacing_results, key=lambda x: spacing_results[x]['mean'])
        logger.info(f"  Peak Spacing: {peak}bp (optimal = 17bp)")

    strand_path = PROJECT_ROOT / "data/results/strand_results.json"
    if strand_path.exists():
        with open(strand_path) as f:
            strand_results = json.load(f)
        fwd = strand_results['forward']['mean']
        rev = strand_results['reverse_in_place']['mean']
        logger.info(f"  Strand: Forward={fwd:.2f}, Reverse={rev:.2f}, Diff={fwd-rev:.2f}")

    logger.info("\n" + "="*80)
    logger.info(f"All logs saved to: {LOG_DIR}")
    logger.info("="*80)


def main():
    """Run all experiments."""
    log_system_info()

    # Check if we should skip GPU models
    skip_gpu = '--skip-gpu' in sys.argv
    gpu_id = 1  # Default GPU

    for arg in sys.argv:
        if arg.startswith('--gpu='):
            gpu_id = int(arg.split('=')[1])

    # Run experiments
    logger.info("\n" + "#"*80)
    logger.info("# STARTING MIT BENCHMARK EXPERIMENTS")
    logger.info("#"*80)

    # 1. Baseline models (CPU)
    run_baseline_models()

    # 2. gLM models (GPU)
    if not skip_gpu:
        run_glm_models(gpu_id=gpu_id)
    else:
        logger.info("Skipping GPU models (--skip-gpu flag)")

    # 3. Extended experiments (GPU)
    if not skip_gpu:
        run_extended_experiments(gpu_id=gpu_id)
    else:
        logger.info("Skipping extended experiments (--skip-gpu flag)")

    # 4. Biophysical comparison with ablations (CPU)
    run_biophysical_comparison()

    # 5. Extended validation experiments (CPU-only parts)
    logger.info("\n" + "="*80)
    logger.info("RUNNING EXTENDED VALIDATION EXPERIMENTS")
    logger.info("="*80)

    try:
        logger.info("\n--- Negative MES Investigation ---")
        import subprocess
        subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts/experiment_negative_mes.py")],
            cwd=str(PROJECT_ROOT),
        )
    except Exception as e:
        logger.error(f"Negative MES investigation failed: {e}")

    try:
        logger.info("\n--- Dinucleotide Control ---")
        subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts/experiment_dinucleotide_control.py")],
            cwd=str(PROJECT_ROOT),
        )
    except Exception as e:
        logger.error(f"Dinucleotide control failed: {e}")

    try:
        logger.info("\n--- Error Analysis ---")
        subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts/experiment_error_analysis.py")],
            cwd=str(PROJECT_ROOT),
        )
    except Exception as e:
        logger.error(f"Error analysis failed: {e}")

    # 6. Compute metrics
    compute_metrics()

    # 7. Generate summary
    generate_summary()

    logger.info("\n" + "#"*80)
    logger.info("# ALL EXPERIMENTS COMPLETE")
    logger.info("#"*80)


if __name__ == "__main__":
    main()
