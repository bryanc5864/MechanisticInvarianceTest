"""Abstract base class for genomic language model wrappers."""

from abc import ABC, abstractmethod
from typing import List, Optional, Dict, Any
import numpy as np


class BaseGLM(ABC):
    """Abstract base class for genomic language models.

    All model wrappers should inherit from this class and implement
    the required methods for computing sequence log-likelihoods.
    """

    def __init__(self, model_name: str, device: str = "cuda"):
        """Initialize the model wrapper.

        Args:
            model_name: Name/identifier for the model
            device: Device to run inference on ("cuda" or "cpu")
        """
        self.model_name = model_name
        self.device = device
        self.model = None
        self.tokenizer = None
        self._loaded = False

    @abstractmethod
    def load_model(self) -> None:
        """Load the model and tokenizer to the specified device.

        This should be called before computing likelihoods.
        Models should be loaded in eval mode with gradients disabled.
        """
        pass

    @abstractmethod
    def compute_log_likelihood(self, sequence: str) -> float:
        """Compute the log-likelihood score for a DNA sequence.

        Args:
            sequence: DNA sequence string (uppercase A, C, G, T)

        Returns:
            Log-likelihood score (higher = more likely under the model)
        """
        pass

    def compute_batch_log_likelihoods(
        self,
        sequences: List[str],
        batch_size: int = 8,
    ) -> List[float]:
        """Compute log-likelihoods for a batch of sequences.

        Default implementation processes sequences one at a time.
        Subclasses may override for more efficient batched inference.

        Args:
            sequences: List of DNA sequences
            batch_size: Batch size for processing (may be ignored)

        Returns:
            List of log-likelihood scores
        """
        return [self.compute_log_likelihood(seq) for seq in sequences]

    def unload_model(self) -> None:
        """Unload the model from memory.

        Useful for freeing GPU memory when switching between models.
        """
        if self.model is not None:
            del self.model
            self.model = None
        if self.tokenizer is not None:
            del self.tokenizer
            self.tokenizer = None
        self._loaded = False

        # Force CUDA memory cleanup
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except ImportError:
            pass

    @property
    def is_loaded(self) -> bool:
        """Check if the model is currently loaded."""
        return self._loaded

    def get_model_info(self) -> Dict[str, Any]:
        """Get information about the model.

        Returns:
            Dictionary with model metadata
        """
        return {
            "name": self.model_name,
            "device": self.device,
            "loaded": self._loaded,
            "type": self.__class__.__name__,
        }

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(model_name='{self.model_name}', device='{self.device}')"
