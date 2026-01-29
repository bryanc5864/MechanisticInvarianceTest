"""Baseline models for comparison (k-mer, PWM, random)."""

import numpy as np
from collections import Counter
from typing import List, Dict, Optional
from tqdm import tqdm

from .base import BaseGLM
from ..sequences.motifs import MINUS_35_PWM, MINUS_10_PWM, score_pwm


class KmerBaseline(BaseGLM):
    """K-mer frequency baseline model.

    Scores sequences based on k-mer frequencies from E. coli genome.
    Uses TF-IDF style scoring.
    """

    def __init__(
        self,
        k: int = 6,
        model_name: str = "kmer_baseline",
        device: str = "cpu",
    ):
        """Initialize k-mer baseline.

        Args:
            k: K-mer size (default 6 for hexamers)
            model_name: Model identifier
            device: Device (always CPU for this model)
        """
        super().__init__(model_name=model_name, device="cpu")
        self.k = k
        self.kmer_freqs: Dict[str, float] = {}
        self.total_kmers = 0

    def load_model(self) -> None:
        """Initialize k-mer frequency table.

        For simplicity, we use uniform background frequencies.
        In practice, these would be computed from E. coli genome.
        """
        print(f"Initializing {self.k}-mer baseline")

        # Initialize with pseudocounts (Laplace smoothing)
        # There are 4^k possible k-mers
        n_kmers = 4 ** self.k
        self.background_prob = 1.0 / n_kmers

        # Generate expected frequencies based on E. coli GC content (~50%)
        # This is simplified - real implementation would use genome-derived frequencies
        self._generate_background_freqs()

        self._loaded = True
        print(f"{self.k}-mer baseline initialized")

    def _generate_background_freqs(self) -> None:
        """Generate background k-mer frequencies.

        Uses simplified model based on nucleotide frequencies.
        Real implementation would count from E. coli genome.
        """
        # E. coli nucleotide frequencies (approximate)
        nuc_freq = {'A': 0.246, 'C': 0.254, 'G': 0.254, 'T': 0.246}

        def generate_kmers(k: int, prefix: str = "") -> None:
            if k == 0:
                # Compute frequency assuming independence
                freq = 1.0
                for c in prefix:
                    freq *= nuc_freq[c]
                self.kmer_freqs[prefix] = freq
            else:
                for nuc in "ACGT":
                    generate_kmers(k - 1, prefix + nuc)

        generate_kmers(self.k)
        self.total_kmers = sum(self.kmer_freqs.values())

    def compute_log_likelihood(self, sequence: str) -> float:
        """Compute log-likelihood based on k-mer frequencies.

        Args:
            sequence: DNA sequence

        Returns:
            Log-likelihood score
        """
        if not self._loaded:
            self.load_model()

        sequence = sequence.upper()
        total_ll = 0.0

        # Count k-mers in sequence
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            if all(c in "ACGT" for c in kmer):
                freq = self.kmer_freqs.get(kmer, self.background_prob)
                total_ll += np.log(freq)

        return total_ll

    def compute_batch_log_likelihoods(
        self,
        sequences: List[str],
        batch_size: int = 8,
    ) -> List[float]:
        """Compute log-likelihoods for multiple sequences."""
        return [self.compute_log_likelihood(seq) for seq in sequences]


class PWMBaseline(BaseGLM):
    """Position Weight Matrix baseline model.

    Scores sequences by scanning for -35 and -10 promoter elements
    using position weight matrices.
    """

    def __init__(
        self,
        model_name: str = "pwm_baseline",
        device: str = "cpu",
    ):
        """Initialize PWM baseline.

        Args:
            model_name: Model identifier
            device: Device (always CPU)
        """
        super().__init__(model_name=model_name, device="cpu")

        # Positions to scan (based on standard promoter layout)
        self.minus_35_start = 30
        self.minus_10_start = 53

    def load_model(self) -> None:
        """Load PWM models."""
        print("Initializing PWM baseline")
        self.minus_35_pwm = MINUS_35_PWM
        self.minus_10_pwm = MINUS_10_PWM
        self._loaded = True
        print("PWM baseline initialized")

    def compute_log_likelihood(self, sequence: str) -> float:
        """Score sequence using PWM for promoter elements.

        Args:
            sequence: DNA sequence (expected 100bp with standard layout)

        Returns:
            Combined PWM score
        """
        if not self._loaded:
            self.load_model()

        sequence = sequence.upper()

        # Extract -35 and -10 regions
        minus_35 = sequence[self.minus_35_start:self.minus_35_start + 6]
        minus_10 = sequence[self.minus_10_start:self.minus_10_start + 6]

        # Score using PWMs
        score_35 = score_pwm(minus_35, self.minus_35_pwm)
        score_10 = score_pwm(minus_10, self.minus_10_pwm)

        # Combine scores (equal weighting)
        return score_35 + score_10

    def compute_batch_log_likelihoods(
        self,
        sequences: List[str],
        batch_size: int = 8,
    ) -> List[float]:
        """Compute scores for multiple sequences."""
        return [self.compute_log_likelihood(seq) for seq in sequences]


class RandomBaseline(BaseGLM):
    """Random baseline model.

    Returns random scores for sanity checking.
    Expected CSS should be ~0.5 for this baseline.
    """

    def __init__(
        self,
        seed: int = 42,
        model_name: str = "random_baseline",
        device: str = "cpu",
    ):
        """Initialize random baseline.

        Args:
            seed: Random seed for reproducibility
            model_name: Model identifier
            device: Device (always CPU)
        """
        super().__init__(model_name=model_name, device="cpu")
        self.rng = np.random.RandomState(seed)

    def load_model(self) -> None:
        """No model to load."""
        self._loaded = True

    def compute_log_likelihood(self, sequence: str) -> float:
        """Return random score.

        Args:
            sequence: DNA sequence (ignored)

        Returns:
            Random score from standard normal
        """
        return self.rng.randn()

    def compute_batch_log_likelihoods(
        self,
        sequences: List[str],
        batch_size: int = 8,
    ) -> List[float]:
        """Return random scores for all sequences."""
        return [self.compute_log_likelihood(seq) for seq in sequences]
