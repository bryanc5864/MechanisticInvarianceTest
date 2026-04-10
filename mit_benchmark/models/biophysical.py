#!/usr/bin/env python3
"""Biophysical models with explicit positional encoding.

These models encode mechanistic knowledge about σ70 promoter structure:
- Position-Aware PWM (PA-PWM): PWM scoring at expected positions
- Thermodynamic Model (Thermo): Free energy-based binding model
- Position-Scanning Model (Scan): Find best motifs, penalize wrong positions
"""

import numpy as np
from typing import Dict, Tuple, Optional


# Position Weight Matrices for σ70 promoter elements
# Values are log-odds scores (positive = favored, negative = disfavored)

PWM_MINUS_35 = [
    # Position 0: T strongly favored (TTGACA consensus)
    {'T': 1.5, 'A': -0.5, 'C': -1.5, 'G': -1.5},
    # Position 1: T strongly favored
    {'T': 1.5, 'A': -0.5, 'C': -1.5, 'G': -1.5},
    # Position 2: G strongly favored
    {'G': 1.5, 'A': -0.5, 'T': -1.0, 'C': -1.5},
    # Position 3: A strongly favored
    {'A': 1.5, 'T': -0.3, 'G': -1.5, 'C': -1.5},
    # Position 4: C strongly favored
    {'C': 1.5, 'A': -0.5, 'T': -1.0, 'G': -1.5},
    # Position 5: A strongly favored
    {'A': 1.5, 'T': -0.3, 'G': -1.5, 'C': -1.5},
]

PWM_MINUS_10 = [
    # Position 0: T strongly favored (TATAAT consensus)
    {'T': 1.5, 'A': -0.3, 'C': -2.0, 'G': -2.0},
    # Position 1: A strongly favored - CRITICAL position for -10 function
    {'A': 2.0, 'T': -0.5, 'G': -2.5, 'C': -2.5},
    # Position 2: T strongly favored
    {'T': 1.5, 'A': -0.5, 'C': -1.5, 'G': -2.0},
    # Position 3: A strongly favored
    {'A': 1.5, 'T': -0.3, 'G': -1.5, 'C': -2.0},
    # Position 4: A strongly favored
    {'A': 1.5, 'T': -0.3, 'G': -1.5, 'C': -1.5},
    # Position 5: T strongly favored
    {'T': 1.5, 'A': -0.3, 'C': -1.5, 'G': -2.0},
]

# UP element PWM (AT-rich, ~20bp upstream of -35)
PWM_UP_PROXIMAL = [
    {'A': 1.0, 'T': 0.5, 'G': -1.0, 'C': -1.0},  # A preferred
    {'A': 1.0, 'T': 0.5, 'G': -1.0, 'C': -1.0},
    {'A': 1.0, 'T': 0.5, 'G': -1.0, 'C': -1.0},
    {'A': 1.0, 'T': 0.5, 'G': -1.0, 'C': -1.0},
    {'A': 1.0, 'T': 0.5, 'G': -1.0, 'C': -1.0},
    {'A': 0.5, 'T': 0.5, 'G': -0.5, 'C': -0.5},  # RNR pattern
    {'A': 0.5, 'G': 0.5, 'T': -0.5, 'C': -0.5},
    {'A': 0.5, 'T': 0.5, 'G': -0.5, 'C': -0.5},
    {'A': 0.5, 'G': 0.5, 'T': -0.5, 'C': -0.5},
]

# Extended -10 (TGn immediately upstream of -10)
PWM_EXTENDED_10 = [
    {'T': 1.5, 'A': -0.5, 'C': -1.5, 'G': -1.5},
    {'G': 1.5, 'A': -1.0, 'T': -1.0, 'C': -1.5},
    {'A': 0.3, 'T': 0.3, 'G': 0.3, 'C': 0.3},  # n = any
]


def score_pwm(sequence: str, pwm: list, position: int) -> float:
    """Score a sequence region against a PWM.

    Args:
        sequence: DNA sequence
        pwm: List of position-specific scoring dicts
        position: Start position in sequence

    Returns:
        PWM score (sum of position scores)
    """
    if position < 0 or position + len(pwm) > len(sequence):
        return -10.0  # Out of bounds penalty

    score = 0.0
    for i, pos_scores in enumerate(pwm):
        nt = sequence[position + i]
        score += pos_scores.get(nt, -2.0)
    return score


def reverse_complement(seq: str) -> str:
    """Get reverse complement of DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(nt, 'N') for nt in reversed(seq))


class PositionAwarePWM:
    """
    Biophysical model with explicit positional encoding.

    Score = PWM(-35) + PWM(-10) + PWM(UP) + spacing_penalty + position_penalty

    This model encodes the key mechanistic principles:
    1. -35 and -10 boxes must be at specific positions relative to TSS
    2. Optimal spacing between -35 and -10 is 17bp
    3. UP element enhances transcription when at correct position
    4. Extended -10 (TGn) can compensate for weak -10
    """

    def __init__(self, tss_position: int = 60, pos_35: int = None,
                 pos_10: int = None, pos_up: int = None, pos_ext10: int = None):
        """
        Args:
            tss_position: Position of transcription start site in 100bp sequence
            pos_35: Explicit -35 position (overrides TSS-based calculation)
            pos_10: Explicit -10 position (overrides TSS-based calculation)
            pos_up: Explicit UP element position (overrides TSS-based calculation)
            pos_ext10: Explicit extended -10 position (overrides TSS-based calculation)
        """
        self.tss = tss_position

        # Use explicit positions if provided, otherwise compute from TSS
        self.pos_35 = pos_35 if pos_35 is not None else self.tss - 35
        self.pos_10 = pos_10 if pos_10 is not None else self.tss - 10
        self.pos_UP = pos_up if pos_up is not None else self.tss - 52
        self.pos_ext10 = pos_ext10 if pos_ext10 is not None else self.pos_10 - 3

        self.optimal_spacing = 17

    def score(self, sequence: str) -> float:
        """Score a promoter sequence.

        Args:
            sequence: 100bp DNA sequence

        Returns:
            Total score (higher = more promoter-like)
        """
        score = 0.0

        # 1. Score -35 box at expected position
        score += score_pwm(sequence, PWM_MINUS_35, self.pos_35)

        # 2. Score -10 box at expected position
        score_10 = score_pwm(sequence, PWM_MINUS_10, self.pos_10)
        score += score_10

        # 3. Score UP element if present (lower weight)
        up_score = score_pwm(sequence, PWM_UP_PROXIMAL, self.pos_UP)
        score += up_score * 0.5

        # 4. Spacing penalty
        actual_spacing = self.pos_10 - self.pos_35 - 6
        spacing_penalty = -abs(actual_spacing - self.optimal_spacing) * 0.5
        score += spacing_penalty

        # 5. Extended -10 bonus (especially valuable if -10 is weak)
        ext10_score = score_pwm(sequence, PWM_EXTENDED_10, self.pos_ext10)
        if ext10_score > 1.0:  # Good extended -10
            # More valuable if -10 is weak
            if score_10 < 5.0:  # Weak -10
                score += ext10_score * 1.5  # Compensation bonus
            else:
                score += ext10_score * 0.5

        return score

    def compute_log_likelihood(self, sequence: str) -> float:
        """Alias for score() to match gLM interface."""
        return self.score(sequence)


class ThermodynamicModel:
    """
    Thermodynamic model of σ70 binding based on free energy.

    ΔG_total = ΔG_35 + ΔG_10 + ΔG_spacing + ΔG_UP + ΔG_ext10

    More negative ΔG = stronger binding = better promoter

    Based on:
    - Brewster et al. (2012) Cell
    - Kinney et al. (2010) PNAS
    """

    # Binding energies - using positive scores for consensus (higher = better binding)
    # This makes the model more intuitive: higher score = better promoter
    ENERGY_MATRIX_35 = [
        # TTGACA consensus - positive for consensus, negative for non-consensus
        {'T': 2.0, 'A': 0.0, 'C': -1.0, 'G': -1.0},
        {'T': 2.0, 'A': 0.0, 'C': -1.0, 'G': -1.0},
        {'G': 2.0, 'A': 0.5, 'T': -0.5, 'C': -1.0},
        {'A': 2.0, 'T': 0.5, 'G': -1.0, 'C': -1.0},
        {'C': 2.0, 'A': 0.5, 'T': -0.5, 'G': -1.0},
        {'A': 2.0, 'T': 0.5, 'G': -1.0, 'C': -1.0},
    ]

    ENERGY_MATRIX_10 = [
        # TATAAT consensus
        {'T': 2.0, 'A': 0.5, 'C': -1.5, 'G': -2.0},
        {'A': 3.0, 'T': 0.0, 'G': -3.0, 'C': -3.0},  # Critical position - A strongly favored
        {'T': 2.0, 'A': 0.5, 'C': -1.0, 'G': -1.5},
        {'A': 2.0, 'T': 0.5, 'G': -1.0, 'C': -1.5},
        {'A': 2.0, 'T': 0.5, 'G': -1.0, 'C': -1.0},
        {'T': 2.0, 'A': 0.5, 'C': -1.0, 'G': -1.5},
    ]

    # Spacing score (17bp optimal, positive for optimal)
    SPACING_ENERGY = {
        12: -2.0, 13: -1.5, 14: -1.0, 15: -0.3, 16: -0.1,
        17: 1.0,  # Optimal - bonus
        18: -0.1, 19: -0.5, 20: -1.0, 21: -1.5, 22: -2.0,
    }

    def __init__(self, tss_position: int = 60, pos_35: int = None,
                 pos_10: int = None, pos_up: int = None, pos_ext10: int = None):
        self.tss = tss_position
        self.pos_35 = pos_35 if pos_35 is not None else self.tss - 35
        self.pos_10 = pos_10 if pos_10 is not None else self.tss - 10
        self.pos_UP = pos_up if pos_up is not None else self.tss - 52
        self.pos_ext10 = pos_ext10 if pos_ext10 is not None else self.pos_10 - 3

    def compute_dG(self, sequence: str) -> float:
        """Compute total binding free energy.

        Args:
            sequence: 100bp DNA sequence

        Returns:
            ΔG in arbitrary units (more negative = stronger binding)
        """
        # ΔG from -35 box
        dG_35 = self._score_element(sequence, self.pos_35, self.ENERGY_MATRIX_35)

        # ΔG from -10 box
        dG_10 = self._score_element(sequence, self.pos_10, self.ENERGY_MATRIX_10)

        # Spacing penalty
        spacing = self.pos_10 - self.pos_35 - 6
        dG_spacing = self.SPACING_ENERGY.get(spacing, -4.0)

        # UP element contribution (only if AT-rich at correct position)
        dG_UP = self._score_UP(sequence, self.pos_UP)

        # Extended -10 contribution
        dG_ext10 = self._score_extended_10(sequence, self.pos_ext10)

        # Compensation: extended -10 is more valuable when -10 is weak
        # Weak -10 scores low (broken TGTAAT has A->G at critical position)
        if dG_10 < 8.0 and dG_ext10 > 1.0:  # Weak -10, good ext-10
            dG_ext10 *= 2.0  # Amplify compensation effect

        return dG_35 + dG_10 + dG_spacing + dG_UP + dG_ext10

    def _score_element(self, seq: str, pos: int, matrix: list) -> float:
        """Score an element using energy matrix."""
        if pos < 0 or pos + len(matrix) > len(seq):
            return -5.0  # Out of bounds penalty

        score = 0.0
        for i, pos_scores in enumerate(matrix):
            nt = seq[pos + i]
            score += pos_scores.get(nt, -2.0)
        return score

    def _score_UP(self, seq: str, pos: int) -> float:
        """Score UP element based on AT content at correct position."""
        if pos < 0:
            return 0.0

        end_pos = min(pos + 20, len(seq))
        up_region = seq[pos:end_pos]
        if len(up_region) == 0:
            return 0.0

        at_content = (up_region.count('A') + up_region.count('T')) / len(up_region)

        # Favorable only if AT-rich (>70%) - positive scores
        if at_content > 0.80:
            return 2.5  # Strong favorable
        elif at_content > 0.70:
            return 1.5  # Moderate favorable
        elif at_content > 0.60:
            return 0.5  # Weak favorable
        return 0.0

    def _score_extended_10(self, seq: str, pos: int) -> float:
        """Score extended -10 (TGn) at correct position."""
        if pos < 0 or pos + 2 > len(seq):
            return 0.0

        motif = seq[pos:pos + 2]
        if motif == "TG":
            return 2.0  # Favorable
        elif motif[0] == 'T':
            return 0.5  # Partial
        return 0.0

    def score(self, sequence: str) -> float:
        """Score sequence (higher = better promoter)."""
        return self.compute_dG(sequence)

    def compute_log_likelihood(self, sequence: str) -> float:
        """Alias for score() to match gLM interface."""
        return self.score(sequence)


class PositionScanningModel:
    """
    Position-scanning model that finds best motif locations,
    then penalizes deviations from expected positions.

    This is more flexible than PA-PWM but still position-aware.

    Score = motif_quality + position_penalty + arrangement_penalty
    """

    def __init__(self, tss_position: int = 60, pos_35: int = None,
                 pos_10: int = None):
        self.tss = tss_position
        self.expected_35 = pos_35 if pos_35 is not None else self.tss - 35
        self.expected_10 = pos_10 if pos_10 is not None else self.tss - 10
        self.optimal_spacing = 17

    def scan_motif(self, sequence: str, pwm: list) -> Tuple[int, float]:
        """Find best position for a motif.

        Args:
            sequence: DNA sequence
            pwm: Position weight matrix

        Returns:
            (best_position, best_score)
        """
        best_pos = 0
        best_score = -float('inf')

        for pos in range(len(sequence) - len(pwm) + 1):
            score = score_pwm(sequence, pwm, pos)
            if score > best_score:
                best_score = score
                best_pos = pos

        return best_pos, best_score

    def score(self, sequence: str) -> float:
        """Score a promoter sequence.

        Args:
            sequence: 100bp DNA sequence

        Returns:
            Total score
        """
        # Find best -35 position
        pos_35, score_35 = self.scan_motif(sequence, PWM_MINUS_35)

        # Find best -10 position
        pos_10, score_10 = self.scan_motif(sequence, PWM_MINUS_10)

        # Base score from motif quality
        score = score_35 + score_10

        # Position penalty (deviation from expected)
        dev_35 = abs(pos_35 - self.expected_35)
        dev_10 = abs(pos_10 - self.expected_10)
        position_penalty = -(dev_35 + dev_10) * 0.3
        score += position_penalty

        # Arrangement penalty
        actual_spacing = pos_10 - pos_35 - 6
        spacing_penalty = -abs(actual_spacing - self.optimal_spacing) * 0.4
        score += spacing_penalty

        # Order penalty (-35 must come before -10)
        if pos_35 >= pos_10:
            score -= 5.0

        # Extended -10 bonus (scan for it near -10)
        ext10_start = max(0, pos_10 - 5)
        ext10_end = pos_10
        best_ext10 = -float('inf')
        for pos in range(ext10_start, ext10_end):
            ext_score = score_pwm(sequence, PWM_EXTENDED_10, pos)
            best_ext10 = max(best_ext10, ext_score)

        if best_ext10 > 1.0:
            # Compensation bonus if -10 is weak
            # A weak -10 (e.g., TGTAAT) scores around 3-5 vs consensus TATAAT at 9+
            if score_10 < 7.0:
                score += best_ext10 * 1.5  # Strong compensation
            else:
                score += best_ext10 * 0.3

        # UP element bonus (scan upstream of -35)
        up_start = max(0, pos_35 - 25)
        up_end = pos_35 - 5
        if up_end > up_start:
            up_region = sequence[up_start:up_end]
            at_content = (up_region.count('A') + up_region.count('T')) / len(up_region)
            if at_content > 0.75:
                score += 2.0 * (at_content - 0.5)

        return score

    def compute_log_likelihood(self, sequence: str) -> float:
        """Alias for score() to match gLM interface."""
        return self.score(sequence)


class PositionAwarePWM_NoComp(PositionAwarePWM):
    """Ablation: PA-PWM without compensation logic.

    Scores -35 and -10 at expected positions with spacing penalty,
    but does NOT score UP element or extended -10.

    If CSS ≈ 0.5 for this model, it confirms that compensation-specific
    logic is required—not just positional encoding of core elements.
    """

    def score(self, sequence: str) -> float:
        score = 0.0

        # 1. Score -35 box at expected position
        score += score_pwm(sequence, PWM_MINUS_35, self.pos_35)

        # 2. Score -10 box at expected position
        score += score_pwm(sequence, PWM_MINUS_10, self.pos_10)

        # 3. Spacing penalty
        actual_spacing = self.pos_10 - self.pos_35 - 6
        spacing_penalty = -abs(actual_spacing - self.optimal_spacing) * 0.5
        score += spacing_penalty

        # NO UP element scoring
        # NO extended -10 scoring
        # NO compensation bonus

        return score


class PositionAwarePWM_NoPosition(PositionAwarePWM):
    """Ablation: PA-PWM without positional encoding.

    Scans the entire sequence for the best -35, -10, UP, and ext-10
    matches anywhere, with no position constraints.

    If CSS is high for this model, positional encoding is not needed.
    If CSS ≈ 0.5, positional encoding is what makes PA-PWM work.
    """

    def score(self, sequence: str) -> float:
        # Scan for best -35 anywhere
        best_35 = max(
            score_pwm(sequence, PWM_MINUS_35, p)
            for p in range(len(sequence) - 6 + 1)
        )
        # Scan for best -10 anywhere
        best_10 = max(
            score_pwm(sequence, PWM_MINUS_10, p)
            for p in range(len(sequence) - 6 + 1)
        )
        score = best_35 + best_10

        # Scan for best UP element anywhere
        best_up = max(
            score_pwm(sequence, PWM_UP_PROXIMAL, p)
            for p in range(len(sequence) - len(PWM_UP_PROXIMAL) + 1)
        )
        score += best_up * 0.5

        # Scan for best extended -10 anywhere
        best_ext10 = max(
            score_pwm(sequence, PWM_EXTENDED_10, p)
            for p in range(len(sequence) - len(PWM_EXTENDED_10) + 1)
        )
        if best_ext10 > 1.0 and best_10 < 5.0:
            score += best_ext10 * 1.5
        elif best_ext10 > 1.0:
            score += best_ext10 * 0.5

        return score


# Generator-aligned positions (matching SequenceGenerator layout)
GENERATOR_POSITIONS = {
    'pos_35': 30,   # SequenceGenerator.MINUS_35_START
    'pos_10': 53,   # SequenceGenerator.MINUS_10_START
    'pos_up': 15,   # SequenceGenerator.UP_START
    'pos_ext10': 50, # SequenceGenerator.EXTENDED_10_START
}


# Convenience functions for loading models
def load_position_aware_pwm(tss: int = 60, align_to_generator: bool = True) -> PositionAwarePWM:
    """Load Position-Aware PWM model."""
    if align_to_generator:
        return PositionAwarePWM(tss_position=tss, **GENERATOR_POSITIONS)
    return PositionAwarePWM(tss_position=tss)

def load_thermodynamic_model(tss: int = 60, align_to_generator: bool = True) -> ThermodynamicModel:
    """Load Thermodynamic model."""
    if align_to_generator:
        return ThermodynamicModel(tss_position=tss, **GENERATOR_POSITIONS)
    return ThermodynamicModel(tss_position=tss)

def load_position_scanning_model(tss: int = 60) -> PositionScanningModel:
    """Load Position-Scanning model."""
    return PositionScanningModel(tss_position=tss)

def load_papwm_no_comp(tss: int = 60, align_to_generator: bool = True) -> PositionAwarePWM_NoComp:
    """Load PA-PWM ablation without compensation logic."""
    if align_to_generator:
        return PositionAwarePWM_NoComp(tss_position=tss, **GENERATOR_POSITIONS)
    return PositionAwarePWM_NoComp(tss_position=tss)

def load_papwm_no_position(tss: int = 60) -> PositionAwarePWM_NoPosition:
    """Load PA-PWM ablation without positional encoding."""
    return PositionAwarePWM_NoPosition(tss_position=tss)
