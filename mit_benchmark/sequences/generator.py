"""Sequence generation for all 8 classes in the MIT benchmark."""

import json
import random
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Dict, Optional

from .motifs import (
    MINUS_35_CONSENSUS,
    MINUS_10_CONSENSUS,
    MINUS_10_BROKEN,
    UP_PROXIMAL_CONSENSUS,
    EXTENDED_10_CONSENSUS,
    sample_iupac,
    scramble_sequence,
)
from .natural import NaturalPromoterFetcher


@dataclass
class PromoterSequence:
    """A promoter sequence with metadata."""
    id: str
    sequence: str  # 100bp sequence
    class_label: str  # A, B, C, D, E, F, G, H
    class_name: str  # Human-readable class name
    minus_35: str
    minus_10: str
    has_up_element: bool
    has_extended_10: bool
    is_natural: bool
    is_broken: bool
    is_compensated: bool
    metadata: Dict

    def to_dict(self) -> Dict:
        return asdict(self)


class SequenceGenerator:
    """Generate promoter sequences for all 8 benchmark classes.

    Classes:
    - A: Natural intact promoters (strong -10, no extra elements needed)
    - B: Natural broken promoters (weak -10, no compensation)
    - C: Synthetic intact promoters (consensus elements)
    - D: Synthetic broken promoters (mutated -10)
    - E: Synthetic compensated (broken -10 + UP + extended -10)
    - F: Synthetic over-compensated (broken -10 + UP + extended -10 + strong -35)
    - G: Natural compensated (based on real promoters with compensation)
    - H: Scrambled compensation control (same composition, no structure)
    """

    # Sequence layout positions (0-indexed)
    UP_START = 15
    UP_END = 24
    MINUS_35_START = 30
    MINUS_35_END = 36
    SPACER_START = 36
    SPACER_END = 53
    EXTENDED_10_START = 50
    EXTENDED_10_END = 53
    MINUS_10_START = 53
    MINUS_10_END = 59
    SEQUENCE_LENGTH = 100

    def __init__(self, seed: int = 42):
        """Initialize the generator.

        Args:
            seed: Random seed for reproducibility
        """
        self.seed = seed
        self.rng = random.Random(seed)
        self.natural_fetcher = NaturalPromoterFetcher(seed=seed)

    def _random_seq(self, length: int, at_bias: bool = True) -> str:
        """Generate random DNA sequence.

        Args:
            length: Sequence length
            at_bias: Whether to use AT-rich bias (typical of E. coli)

        Returns:
            Random DNA sequence
        """
        if at_bias:
            weights = [0.28, 0.22, 0.22, 0.28]  # A, C, G, T
        else:
            weights = [0.25, 0.25, 0.25, 0.25]
        return ''.join(self.rng.choices('ACGT', weights=weights, k=length))

    def _build_sequence(
        self,
        minus_35: str,
        minus_10: str,
        up_element: Optional[str] = None,
        extended_10: Optional[str] = None,
    ) -> str:
        """Build a 100bp promoter sequence.

        Args:
            minus_35: -35 box sequence (6bp)
            minus_10: -10 box sequence (6bp)
            up_element: UP element sequence (9bp) or None
            extended_10: Extended -10 sequence (3bp, typically TGT) or None

        Returns:
            100bp promoter sequence
        """
        seq = list(self._random_seq(self.SEQUENCE_LENGTH))

        # Place UP element (positions 15-23)
        if up_element:
            for i, base in enumerate(up_element[:9]):
                seq[self.UP_START + i] = base

        # Place -35 box (positions 30-35)
        for i, base in enumerate(minus_35):
            seq[self.MINUS_35_START + i] = base

        # Place extended -10 (positions 50-52, overlaps spacer)
        if extended_10:
            for i, base in enumerate(extended_10[:3]):
                seq[self.EXTENDED_10_START + i] = base

        # Place -10 box (positions 53-58)
        for i, base in enumerate(minus_10):
            seq[self.MINUS_10_START + i] = base

        return ''.join(seq)

    def generate_class_a(self, n: int = 100) -> List[PromoterSequence]:
        """Generate Class A: Natural intact promoters.

        These have strong -10 boxes (close to TATAAT consensus).

        Args:
            n: Number of sequences to generate

        Returns:
            List of PromoterSequence objects
        """
        natural = self.natural_fetcher.get_intact_promoters(n)
        results = []

        for i, p in enumerate(natural):
            results.append(PromoterSequence(
                id=f"A_{i:03d}",
                sequence=p.sequence,
                class_label="A",
                class_name="Natural Intact",
                minus_35=p.minus_35,
                minus_10=p.minus_10,
                has_up_element=p.has_up_element,
                has_extended_10=p.has_extended_10,
                is_natural=True,
                is_broken=False,
                is_compensated=False,
                metadata={"source": p.source, "gene": p.gene, "name": p.name},
            ))

        return results

    def generate_class_b(self, n: int = 100) -> List[PromoterSequence]:
        """Generate Class B: Natural broken promoters.

        These have weak -10 boxes (T→G mutation at position 2) and no compensation.

        Args:
            n: Number of sequences to generate

        Returns:
            List of PromoterSequence objects
        """
        natural = self.natural_fetcher.get_broken_promoters(n)
        results = []

        for i, p in enumerate(natural):
            results.append(PromoterSequence(
                id=f"B_{i:03d}",
                sequence=p.sequence,
                class_label="B",
                class_name="Natural Broken",
                minus_35=p.minus_35,
                minus_10=p.minus_10,
                has_up_element=False,
                has_extended_10=False,
                is_natural=True,
                is_broken=True,
                is_compensated=False,
                metadata={"source": p.source, "gene": p.gene, "name": p.name},
            ))

        return results

    def generate_class_c(self, n: int = 100) -> List[PromoterSequence]:
        """Generate Class C: Synthetic intact promoters.

        These use consensus -35 (TTGACA) and -10 (TATAAT) elements.

        Args:
            n: Number of sequences to generate

        Returns:
            List of PromoterSequence objects
        """
        results = []

        for i in range(n):
            sequence = self._build_sequence(
                minus_35=MINUS_35_CONSENSUS,
                minus_10=MINUS_10_CONSENSUS,
                up_element=None,
                extended_10=None,
            )

            results.append(PromoterSequence(
                id=f"C_{i:03d}",
                sequence=sequence,
                class_label="C",
                class_name="Synthetic Intact",
                minus_35=MINUS_35_CONSENSUS,
                minus_10=MINUS_10_CONSENSUS,
                has_up_element=False,
                has_extended_10=False,
                is_natural=False,
                is_broken=False,
                is_compensated=False,
                metadata={"design": "consensus"},
            ))

        return results

    def generate_class_d(self, n: int = 100) -> List[PromoterSequence]:
        """Generate Class D: Synthetic broken promoters.

        These use consensus -35 but broken -10 (TGTAAT).

        Args:
            n: Number of sequences to generate

        Returns:
            List of PromoterSequence objects
        """
        results = []

        for i in range(n):
            sequence = self._build_sequence(
                minus_35=MINUS_35_CONSENSUS,
                minus_10=MINUS_10_BROKEN,
                up_element=None,
                extended_10=None,
            )

            results.append(PromoterSequence(
                id=f"D_{i:03d}",
                sequence=sequence,
                class_label="D",
                class_name="Synthetic Broken",
                minus_35=MINUS_35_CONSENSUS,
                minus_10=MINUS_10_BROKEN,
                has_up_element=False,
                has_extended_10=False,
                is_natural=False,
                is_broken=True,
                is_compensated=False,
                metadata={"design": "broken_-10"},
            ))

        return results

    def generate_class_e(self, n: int = 100) -> List[PromoterSequence]:
        """Generate Class E: Synthetic compensated promoters.

        These have broken -10 but compensatory UP element and extended -10.

        Args:
            n: Number of sequences to generate

        Returns:
            List of PromoterSequence objects
        """
        results = []

        for i in range(n):
            # Sample UP element from consensus pattern
            up_element = sample_iupac(UP_PROXIMAL_CONSENSUS)

            sequence = self._build_sequence(
                minus_35=MINUS_35_CONSENSUS,
                minus_10=MINUS_10_BROKEN,
                up_element=up_element,
                extended_10=EXTENDED_10_CONSENSUS,
            )

            results.append(PromoterSequence(
                id=f"E_{i:03d}",
                sequence=sequence,
                class_label="E",
                class_name="Synthetic Compensated",
                minus_35=MINUS_35_CONSENSUS,
                minus_10=MINUS_10_BROKEN,
                has_up_element=True,
                has_extended_10=True,
                is_natural=False,
                is_broken=True,
                is_compensated=True,
                metadata={"design": "compensated", "up_element": up_element},
            ))

        return results

    def generate_class_f(self, n: int = 50) -> List[PromoterSequence]:
        """Generate Class F: Synthetic over-compensated promoters.

        These have broken -10 with multiple compensatory elements:
        - Strong consensus -35 (TTGACA)
        - UP element
        - Extended -10

        Args:
            n: Number of sequences to generate

        Returns:
            List of PromoterSequence objects
        """
        results = []

        # Use stronger -35 variants
        strong_35_variants = ["TTGACA", "TTGACT", "TTGATA"]

        for i in range(n):
            up_element = sample_iupac(UP_PROXIMAL_CONSENSUS)
            minus_35 = self.rng.choice(strong_35_variants)

            sequence = self._build_sequence(
                minus_35=minus_35,
                minus_10=MINUS_10_BROKEN,
                up_element=up_element,
                extended_10=EXTENDED_10_CONSENSUS,
            )

            results.append(PromoterSequence(
                id=f"F_{i:03d}",
                sequence=sequence,
                class_label="F",
                class_name="Synthetic Over-Compensated",
                minus_35=minus_35,
                minus_10=MINUS_10_BROKEN,
                has_up_element=True,
                has_extended_10=True,
                is_natural=False,
                is_broken=True,
                is_compensated=True,
                metadata={
                    "design": "over_compensated",
                    "up_element": up_element,
                    "strong_35": True,
                },
            ))

        return results

    def generate_class_g(self, n: int = 50) -> List[PromoterSequence]:
        """Generate Class G: Natural compensated promoters.

        These are natural promoters with weak -10 but compensatory elements.

        Args:
            n: Number of sequences to generate

        Returns:
            List of PromoterSequence objects
        """
        natural = self.natural_fetcher.get_compensated_promoters(n)
        results = []

        for i, p in enumerate(natural):
            results.append(PromoterSequence(
                id=f"G_{i:03d}",
                sequence=p.sequence,
                class_label="G",
                class_name="Natural Compensated",
                minus_35=p.minus_35,
                minus_10=p.minus_10,
                has_up_element=True,
                has_extended_10=True,
                is_natural=True,
                is_broken=True,
                is_compensated=True,
                metadata={"source": p.source, "gene": p.gene, "name": p.name},
            ))

        return results

    def generate_class_h(self, n: int = 50) -> List[PromoterSequence]:
        """Generate Class H: Scrambled compensation control.

        These have the same nucleotide composition as compensated promoters
        but with scrambled motifs (destroying functional structure).

        Args:
            n: Number of sequences to generate

        Returns:
            List of PromoterSequence objects
        """
        # First generate compensated sequences, then scramble the motif regions
        compensated = self.generate_class_e(n)
        results = []

        for i, p in enumerate(compensated):
            # Scramble the UP element region
            seq = list(p.sequence)

            # Scramble UP element (15-23)
            up_region = seq[self.UP_START:self.UP_END]
            self.rng.shuffle(up_region)
            seq[self.UP_START:self.UP_END] = up_region

            # Scramble extended -10 (50-52)
            ext_region = seq[self.EXTENDED_10_START:self.EXTENDED_10_END]
            self.rng.shuffle(ext_region)
            seq[self.EXTENDED_10_START:self.EXTENDED_10_END] = ext_region

            scrambled_seq = ''.join(seq)

            results.append(PromoterSequence(
                id=f"H_{i:03d}",
                sequence=scrambled_seq,
                class_label="H",
                class_name="Scrambled Control",
                minus_35=p.minus_35,
                minus_10=p.minus_10,
                has_up_element=False,  # Scrambled, so not functional
                has_extended_10=False,
                is_natural=False,
                is_broken=True,
                is_compensated=False,  # Scrambled removes compensation
                metadata={
                    "design": "scrambled_control",
                    "original_class": "E",
                },
            ))

        return results

    def generate_all_classes(self) -> Dict[str, List[PromoterSequence]]:
        """Generate sequences for all 8 classes.

        Returns:
            Dictionary mapping class label to list of sequences
        """
        return {
            "A": self.generate_class_a(100),
            "B": self.generate_class_b(100),
            "C": self.generate_class_c(100),
            "D": self.generate_class_d(100),
            "E": self.generate_class_e(100),
            "F": self.generate_class_f(50),
            "G": self.generate_class_g(50),
            "H": self.generate_class_h(50),
        }

    def save_sequences(
        self,
        sequences: Dict[str, List[PromoterSequence]],
        output_path: Path,
        format: str = "json",
    ) -> None:
        """Save sequences to file.

        Args:
            sequences: Dictionary of class label to sequences
            output_path: Path to output file
            format: Output format ("json" or "fasta")
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format == "json":
            # Flatten to list with all metadata
            all_seqs = []
            for class_label, seq_list in sequences.items():
                for seq in seq_list:
                    all_seqs.append(seq.to_dict())

            with open(output_path, 'w') as f:
                json.dump(all_seqs, f, indent=2)

        elif format == "fasta":
            with open(output_path, 'w') as f:
                for class_label, seq_list in sequences.items():
                    for seq in seq_list:
                        f.write(f">{seq.id}|{seq.class_name}|broken={seq.is_broken}|comp={seq.is_compensated}\n")
                        f.write(f"{seq.sequence}\n")

        else:
            raise ValueError(f"Unknown format: {format}")

    @staticmethod
    def load_sequences(input_path: Path) -> List[PromoterSequence]:
        """Load sequences from JSON file.

        Args:
            input_path: Path to input JSON file

        Returns:
            List of PromoterSequence objects
        """
        with open(input_path, 'r') as f:
            data = json.load(f)

        return [PromoterSequence(**d) for d in data]
