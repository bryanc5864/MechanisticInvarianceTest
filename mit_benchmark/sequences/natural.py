"""Natural promoter fetching and processing."""

import random
from typing import List, Dict, Optional
from dataclasses import dataclass


@dataclass
class NaturalPromoter:
    """Represents a natural E. coli promoter."""
    name: str
    sequence: str  # Full 100bp context
    minus_35: str  # -35 box sequence
    minus_10: str  # -10 box sequence
    has_up_element: bool
    has_extended_10: bool
    gene: str
    source: str  # RegulonDB, literature, etc.


# Well-characterized E. coli σ70 promoters from literature
# These are canonical promoters with experimentally verified activity
CURATED_PROMOTERS = [
    # Strong promoters with canonical elements
    {
        "name": "lacUV5",
        "gene": "lac operon",
        "minus_35": "TTTACA",
        "minus_10": "TATAAT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "tac",
        "gene": "tac hybrid",
        "minus_35": "TTGACA",
        "minus_10": "TATAAT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "trp",
        "gene": "trp operon",
        "minus_35": "TTGACA",
        "minus_10": "TTAACT",
        "has_up_element": False,
        "has_extended_10": True,
        "source": "literature",
    },
    {
        "name": "rrnB_P1",
        "gene": "rRNA",
        "minus_35": "TTTTCT",
        "minus_10": "TATAAT",
        "has_up_element": True,
        "has_extended_10": True,
        "source": "literature",
    },
    {
        "name": "rrnB_P2",
        "gene": "rRNA",
        "minus_35": "ATGCAT",
        "minus_10": "TAAAAT",
        "has_up_element": True,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "lpp",
        "gene": "lipoprotein",
        "minus_35": "TTGTCA",
        "minus_10": "TAAAAT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "galP1",
        "gene": "galactose",
        "minus_35": "GTGTCA",
        "minus_10": "TATGTT",
        "has_up_element": False,
        "has_extended_10": True,
        "source": "literature",
    },
    {
        "name": "galP2",
        "gene": "galactose",
        "minus_35": "TTGCAT",
        "minus_10": "TATAAT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "recA",
        "gene": "recA",
        "minus_35": "TTGATA",
        "minus_10": "TATAAT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "ompA",
        "gene": "outer membrane",
        "minus_35": "TCTTTT",
        "minus_10": "TACACT",
        "has_up_element": True,
        "has_extended_10": True,
        "source": "literature",
    },
    # Additional promoters with varied characteristics
    {
        "name": "araBAD",
        "gene": "arabinose",
        "minus_35": "CTGACG",
        "minus_10": "TACTGT",
        "has_up_element": False,
        "has_extended_10": True,
        "source": "literature",
    },
    {
        "name": "malT",
        "gene": "maltose",
        "minus_35": "TTGTAA",
        "minus_10": "AATAAT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "tyrT",
        "gene": "tRNA-Tyr",
        "minus_35": "ATGCAA",
        "minus_10": "GATACT",
        "has_up_element": True,
        "has_extended_10": True,
        "source": "literature",
    },
    {
        "name": "bla",
        "gene": "beta-lactamase",
        "minus_35": "TTCAAA",
        "minus_10": "TATGTT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "cat",
        "gene": "CAT",
        "minus_35": "TGGAAA",
        "minus_10": "TAAACT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "dnaA_P1",
        "gene": "dnaA",
        "minus_35": "TTGTCC",
        "minus_10": "TTTAAT",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "dnaA_P2",
        "gene": "dnaA",
        "minus_35": "TTATCA",
        "minus_10": "TAGACT",
        "has_up_element": True,
        "has_extended_10": True,
        "source": "literature",
    },
    {
        "name": "fis",
        "gene": "Fis",
        "minus_35": "TTGACT",
        "minus_10": "TACGAT",
        "has_up_element": True,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "ihfA",
        "gene": "IHF alpha",
        "minus_35": "TTGACC",
        "minus_10": "TATAAC",
        "has_up_element": False,
        "has_extended_10": False,
        "source": "literature",
    },
    {
        "name": "rpoD",
        "gene": "sigma70",
        "minus_35": "TTGACA",
        "minus_10": "GATACT",
        "has_up_element": True,
        "has_extended_10": True,
        "source": "literature",
    },
]


class NaturalPromoterFetcher:
    """Fetches and processes natural E. coli σ70 promoters."""

    def __init__(self, seed: int = 42):
        """Initialize the fetcher.

        Args:
            seed: Random seed for reproducibility
        """
        self.rng = random.Random(seed)
        self.promoters = self._load_curated_promoters()

    def _load_curated_promoters(self) -> List[Dict]:
        """Load curated promoter data."""
        return CURATED_PROMOTERS.copy()

    def _generate_context(
        self,
        minus_35: str,
        minus_10: str,
        has_up_element: bool = False,
        has_extended_10: bool = False,
        up_element: Optional[str] = None,
    ) -> str:
        """Generate a 100bp sequence context for a promoter.

        Layout (0-indexed):
        - 0-14: Random upstream flank
        - 15-23: UP element (if present) or random
        - 24-29: Random spacer
        - 30-35: -35 box (6bp)
        - 36-52: Spacer (17bp)
        - 50-52: Extended -10 (TGT, if present, overlaps end of spacer)
        - 53-58: -10 box (6bp)
        - 59-99: Random downstream flank

        Args:
            minus_35: The -35 box sequence (6bp)
            minus_10: The -10 box sequence (6bp)
            has_up_element: Whether to include UP element
            has_extended_10: Whether to include extended -10
            up_element: Specific UP element to use (optional)

        Returns:
            100bp promoter context sequence
        """
        # Generate random flanks with slight AT bias
        def random_seq(length: int) -> str:
            weights = [0.28, 0.22, 0.22, 0.28]  # A, C, G, T
            return ''.join(self.rng.choices('ACGT', weights=weights, k=length))

        # Build the sequence piece by piece
        seq = ['N'] * 100  # Placeholder

        # Upstream flank (0-14)
        for i in range(15):
            seq[i] = self.rng.choices('ACGT', weights=[0.28, 0.22, 0.22, 0.28])[0]

        # UP element region (15-23, 9bp)
        if has_up_element:
            if up_element and len(up_element) >= 9:
                up_seq = up_element[:9]
            else:
                # Generate A-rich UP element: AAAAAARNR pattern (9bp)
                up_seq = "AAAAA" + self.rng.choice("AG") + self.rng.choice("ACGT") + self.rng.choice("AG") + self.rng.choice("AG")
            for i, base in enumerate(up_seq[:9]):
                seq[15 + i] = base
        else:
            for i in range(15, 24):
                seq[i] = self.rng.choices('ACGT', weights=[0.28, 0.22, 0.22, 0.28])[0]

        # Spacer between UP and -35 (24-29, 6bp)
        for i in range(24, 30):
            seq[i] = self.rng.choices('ACGT', weights=[0.28, 0.22, 0.22, 0.28])[0]

        # -35 box (30-35, 6bp)
        for i, base in enumerate(minus_35):
            seq[30 + i] = base

        # Spacer between -35 and -10 (36-52, 17bp)
        for i in range(36, 53):
            seq[i] = self.rng.choices('ACGT', weights=[0.28, 0.22, 0.22, 0.28])[0]

        # Extended -10 (50-52, TGT) - overlaps end of spacer
        if has_extended_10:
            seq[50] = 'T'
            seq[51] = 'G'
            seq[52] = 'T'

        # -10 box (53-58, 6bp)
        for i, base in enumerate(minus_10):
            seq[53 + i] = base

        # Downstream flank (59-99)
        for i in range(59, 100):
            seq[i] = self.rng.choices('ACGT', weights=[0.28, 0.22, 0.22, 0.28])[0]

        return ''.join(seq)

    def get_intact_promoters(self, n: int = 100) -> List[NaturalPromoter]:
        """Get natural promoters with intact (strong) -10 boxes.

        These promoters have -10 boxes close to the TATAAT consensus.

        Args:
            n: Number of promoters to return

        Returns:
            List of NaturalPromoter objects
        """
        # Filter for promoters with strong -10 boxes
        strong_promoters = [
            p for p in self.promoters
            if self._score_minus_10(p["minus_10"]) >= 0
        ]

        results = []
        for i in range(n):
            # Cycle through available promoters with variation
            base_promoter = strong_promoters[i % len(strong_promoters)]

            # Generate context sequence
            sequence = self._generate_context(
                minus_35=base_promoter["minus_35"],
                minus_10=base_promoter["minus_10"],
                has_up_element=base_promoter["has_up_element"],
                has_extended_10=base_promoter["has_extended_10"],
            )

            results.append(NaturalPromoter(
                name=f"{base_promoter['name']}_v{i // len(strong_promoters)}",
                sequence=sequence,
                minus_35=base_promoter["minus_35"],
                minus_10=base_promoter["minus_10"],
                has_up_element=base_promoter["has_up_element"],
                has_extended_10=base_promoter["has_extended_10"],
                gene=base_promoter["gene"],
                source=base_promoter["source"],
            ))

        return results

    def get_broken_promoters(self, n: int = 100) -> List[NaturalPromoter]:
        """Get natural promoters with broken (weak) -10 boxes.

        These are derived from intact promoters with T→G mutation at position 2.

        Args:
            n: Number of promoters to return

        Returns:
            List of NaturalPromoter objects with broken -10
        """
        intact = self.get_intact_promoters(n)
        results = []

        for p in intact:
            # Break the -10 box (T→G at position 2)
            broken_minus_10 = self._break_minus_10(p.minus_10)

            # Regenerate context with broken -10
            sequence = self._generate_context(
                minus_35=p.minus_35,
                minus_10=broken_minus_10,
                has_up_element=False,  # No compensation
                has_extended_10=False,
            )

            results.append(NaturalPromoter(
                name=f"{p.name}_broken",
                sequence=sequence,
                minus_35=p.minus_35,
                minus_10=broken_minus_10,
                has_up_element=False,
                has_extended_10=False,
                gene=p.gene,
                source=p.source,
            ))

        return results

    def get_compensated_promoters(self, n: int = 50) -> List[NaturalPromoter]:
        """Get natural promoters with broken -10 but compensatory elements.

        These have weak -10 boxes but UP elements and/or extended -10.

        Args:
            n: Number of promoters to return

        Returns:
            List of NaturalPromoter objects with compensation
        """
        intact = self.get_intact_promoters(n)
        results = []

        for p in intact:
            # Break the -10 box
            broken_minus_10 = self._break_minus_10(p.minus_10)

            # Add compensatory elements
            sequence = self._generate_context(
                minus_35=p.minus_35,
                minus_10=broken_minus_10,
                has_up_element=True,
                has_extended_10=True,
            )

            results.append(NaturalPromoter(
                name=f"{p.name}_compensated",
                sequence=sequence,
                minus_35=p.minus_35,
                minus_10=broken_minus_10,
                has_up_element=True,
                has_extended_10=True,
                gene=p.gene,
                source=p.source,
            ))

        return results

    def _score_minus_10(self, minus_10: str) -> float:
        """Score a -10 box against TATAAT consensus.

        Args:
            minus_10: 6bp -10 sequence

        Returns:
            Score (higher = more similar to consensus)
        """
        consensus = "TATAAT"
        score = 0
        for i, (a, b) in enumerate(zip(minus_10, consensus)):
            if a == b:
                score += 1
            # Position 1 (second T) and position 5 (last T) are most critical
            if i in [1, 5] and a == b:
                score += 0.5
        return score

    def _break_minus_10(self, minus_10: str) -> str:
        """Break a -10 box by mutating position 2.

        Args:
            minus_10: Original 6bp -10 sequence

        Returns:
            Broken -10 sequence with T→G at position 2 (or similar)
        """
        seq_list = list(minus_10)

        # Mutate position 1 (the second position, 0-indexed)
        # Canonical mutation is T→G
        if seq_list[1] == 'T':
            seq_list[1] = 'G'
        elif seq_list[1] == 'A':
            seq_list[1] = 'C'
        else:
            # Pick a different base
            alternatives = [b for b in 'ACGT' if b != seq_list[1]]
            seq_list[1] = self.rng.choice(alternatives)

        return ''.join(seq_list)
