"""Masked language model wrappers using pseudo-log-likelihood."""

import torch
import numpy as np
from typing import List, Optional
from tqdm import tqdm

from .base import BaseGLM


def compute_pseudo_log_likelihood(
    model,
    tokenizer,
    sequence: str,
    device: str = "cuda",
    mask_token_id: Optional[int] = None,
) -> float:
    """Compute pseudo-log-likelihood for masked language models.

    PLL(x) = sum_i log P(x_i | x_{-i})

    This masks each position one at a time and sums the log probabilities
    of the original tokens given the masked context.

    Args:
        model: HuggingFace model
        tokenizer: HuggingFace tokenizer
        sequence: DNA sequence
        device: Device for inference
        mask_token_id: Token ID for mask token (auto-detected if None)

    Returns:
        Pseudo-log-likelihood score
    """
    # Tokenize
    inputs = tokenizer(
        sequence,
        return_tensors="pt",
        add_special_tokens=True,
    )
    input_ids = inputs["input_ids"].to(device)
    attention_mask = inputs.get("attention_mask", torch.ones_like(input_ids)).to(device)

    # Get mask token ID
    if mask_token_id is None:
        mask_token_id = tokenizer.mask_token_id
        if mask_token_id is None:
            # Try common alternatives
            if hasattr(tokenizer, 'mask_token'):
                mask_token_id = tokenizer.convert_tokens_to_ids(tokenizer.mask_token)
            else:
                raise ValueError("Could not find mask token ID")

    # Get special token positions to skip
    special_ids = set()
    if tokenizer.cls_token_id is not None:
        special_ids.add(tokenizer.cls_token_id)
    if tokenizer.sep_token_id is not None:
        special_ids.add(tokenizer.sep_token_id)
    if tokenizer.pad_token_id is not None:
        special_ids.add(tokenizer.pad_token_id)

    total_ll = 0.0
    seq_len = input_ids.shape[1]

    # Batch multiple masked positions for efficiency
    # We'll mask one position at a time but can batch across sequences
    for pos in range(seq_len):
        original_token = input_ids[0, pos].item()

        # Skip special tokens
        if original_token in special_ids:
            continue

        # Create masked input
        masked_ids = input_ids.clone()
        masked_ids[0, pos] = mask_token_id

        # Forward pass
        with torch.no_grad():
            outputs = model(masked_ids, attention_mask=attention_mask)
            logits = outputs.logits

        # Get log probability of original token
        log_probs = torch.log_softmax(logits[0, pos, :], dim=-1)
        total_ll += log_probs[original_token].item()

    return total_ll


class DNABERT2Wrapper(BaseGLM):
    """Wrapper for DNABERT-2 masked language model.

    DNABERT-2 uses BPE tokenization and is trained with masked LM objective.
    """

    def __init__(
        self,
        model_name: str = "zhihan1996/DNABERT-2-117M",
        device: str = "cuda",
    ):
        """Initialize DNABERT-2 wrapper.

        Args:
            model_name: HuggingFace model ID
            device: Device for inference
        """
        super().__init__(model_name=model_name, device=device)

    def load_model(self) -> None:
        """Load DNABERT-2 model."""
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        print(f"Loading DNABERT-2 model: {self.model_name}")

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.model_name,
            trust_remote_code=True,
        )
        self.model = AutoModelForMaskedLM.from_pretrained(
            self.model_name,
            trust_remote_code=True,
            torch_dtype=torch.float16,
        )
        self.model.to(self.device)
        self.model.eval()
        self._loaded = True
        print(f"DNABERT-2 model loaded on {self.device}")

    def compute_log_likelihood(self, sequence: str) -> float:
        """Compute pseudo-log-likelihood for a DNA sequence.

        Args:
            sequence: DNA sequence

        Returns:
            Pseudo-log-likelihood score
        """
        if not self._loaded:
            self.load_model()

        return compute_pseudo_log_likelihood(
            self.model,
            self.tokenizer,
            sequence.upper(),
            self.device,
        )

    def compute_batch_log_likelihoods(
        self,
        sequences: List[str],
        batch_size: int = 8,
    ) -> List[float]:
        """Compute PLLs for multiple sequences."""
        results = []
        for seq in tqdm(sequences, desc="DNABERT-2 inference"):
            results.append(self.compute_log_likelihood(seq))
        return results


class NucleotideTransformerWrapper(BaseGLM):
    """Wrapper for Nucleotide Transformer (InstaDeep).

    NT uses 6-mer tokenization and is trained on multi-species genomes.
    """

    def __init__(
        self,
        model_name: str = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species",
        device: str = "cuda",
    ):
        """Initialize NT wrapper.

        Args:
            model_name: HuggingFace model ID
            device: Device for inference
        """
        super().__init__(model_name=model_name, device=device)

    def load_model(self) -> None:
        """Load Nucleotide Transformer model."""
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        print(f"Loading Nucleotide Transformer: {self.model_name}")

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.model_name,
            trust_remote_code=True,
        )
        self.model = AutoModelForMaskedLM.from_pretrained(
            self.model_name,
            trust_remote_code=True,
        )
        self.model.to(self.device)
        self.model.eval()
        self._loaded = True
        print(f"Nucleotide Transformer loaded on {self.device}")

    def compute_log_likelihood(self, sequence: str) -> float:
        """Compute pseudo-log-likelihood.

        Args:
            sequence: DNA sequence

        Returns:
            PLL score
        """
        if not self._loaded:
            self.load_model()

        return compute_pseudo_log_likelihood(
            self.model,
            self.tokenizer,
            sequence.upper(),
            self.device,
        )

    def compute_batch_log_likelihoods(
        self,
        sequences: List[str],
        batch_size: int = 8,
    ) -> List[float]:
        """Compute PLLs for multiple sequences."""
        results = []
        for seq in tqdm(sequences, desc="NT inference"):
            results.append(self.compute_log_likelihood(seq))
        return results


class GROVERWrapper(BaseGLM):
    """Wrapper for GROVER genomic language model.

    GROVER is trained on microbial genomes with masked LM objective.
    """

    def __init__(
        self,
        model_name: str = "PoetschLab/GROVER",
        device: str = "cuda",
    ):
        """Initialize GROVER wrapper.

        Args:
            model_name: HuggingFace model ID
            device: Device for inference
        """
        super().__init__(model_name=model_name, device=device)

    def load_model(self) -> None:
        """Load GROVER model."""
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        print(f"Loading GROVER model: {self.model_name}")

        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self.model = AutoModelForMaskedLM.from_pretrained(
            self.model_name,
            torch_dtype=torch.float16,
        )
        self.model.to(self.device)
        self.model.eval()
        self._loaded = True
        print(f"GROVER model loaded on {self.device}")

    def compute_log_likelihood(self, sequence: str) -> float:
        """Compute pseudo-log-likelihood.

        Args:
            sequence: DNA sequence

        Returns:
            PLL score
        """
        if not self._loaded:
            self.load_model()

        return compute_pseudo_log_likelihood(
            self.model,
            self.tokenizer,
            sequence.upper(),
            self.device,
        )

    def compute_batch_log_likelihoods(
        self,
        sequences: List[str],
        batch_size: int = 8,
    ) -> List[float]:
        """Compute PLLs for multiple sequences."""
        results = []
        for seq in tqdm(sequences, desc="GROVER inference"):
            results.append(self.compute_log_likelihood(seq))
        return results


class CaduceusWrapper(BaseGLM):
    """Wrapper for Caduceus bidirectional DNA model.

    Caduceus uses Mamba architecture with bidirectional capability.
    """

    def __init__(
        self,
        model_name: str = "kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16",
        device: str = "cuda",
    ):
        """Initialize Caduceus wrapper.

        Args:
            model_name: HuggingFace model ID
            device: Device for inference
        """
        super().__init__(model_name=model_name, device=device)

    def load_model(self) -> None:
        """Load Caduceus model."""
        from transformers import AutoModelForMaskedLM, AutoTokenizer

        print(f"Loading Caduceus model: {self.model_name}")

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.model_name,
            trust_remote_code=True,
        )
        self.model = AutoModelForMaskedLM.from_pretrained(
            self.model_name,
            trust_remote_code=True,
            torch_dtype=torch.float16,
        )
        self.model.to(self.device)
        self.model.eval()
        self._loaded = True
        print(f"Caduceus model loaded on {self.device}")

    def compute_log_likelihood(self, sequence: str) -> float:
        """Compute pseudo-log-likelihood.

        Caduceus is bidirectional, so PLL is appropriate.
        Caduceus (Mamba-based) does not use attention_mask, so we
        compute PLL directly without passing it.

        Args:
            sequence: DNA sequence

        Returns:
            PLL score
        """
        if not self._loaded:
            self.load_model()

        sequence = sequence.upper()

        # Tokenize
        inputs = self.tokenizer(
            sequence,
            return_tensors="pt",
            add_special_tokens=True,
        )
        input_ids = inputs["input_ids"].to(self.device)

        # Get mask token ID
        mask_token_id = self.tokenizer.mask_token_id
        if mask_token_id is None:
            if hasattr(self.tokenizer, 'mask_token'):
                mask_token_id = self.tokenizer.convert_tokens_to_ids(self.tokenizer.mask_token)
            else:
                raise ValueError("Could not find mask token ID")

        # Get special token positions to skip
        special_ids = set()
        if self.tokenizer.cls_token_id is not None:
            special_ids.add(self.tokenizer.cls_token_id)
        if self.tokenizer.sep_token_id is not None:
            special_ids.add(self.tokenizer.sep_token_id)
        if self.tokenizer.pad_token_id is not None:
            special_ids.add(self.tokenizer.pad_token_id)

        total_ll = 0.0
        seq_len = input_ids.shape[1]

        for pos in range(seq_len):
            original_token = input_ids[0, pos].item()

            if original_token in special_ids:
                continue

            masked_ids = input_ids.clone()
            masked_ids[0, pos] = mask_token_id

            with torch.no_grad():
                outputs = self.model(masked_ids)
                logits = outputs.logits

            log_probs = torch.log_softmax(logits[0, pos, :], dim=-1)
            total_ll += log_probs[original_token].item()

        return total_ll

    def compute_batch_log_likelihoods(
        self,
        sequences: List[str],
        batch_size: int = 8,
    ) -> List[float]:
        """Compute PLLs for multiple sequences."""
        results = []
        for seq in tqdm(sequences, desc="Caduceus inference"):
            results.append(self.compute_log_likelihood(seq))
        return results
