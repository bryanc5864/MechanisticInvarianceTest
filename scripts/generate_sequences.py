#!/usr/bin/env python3
"""Generate benchmark sequences for all 8 classes.

Usage:
    python scripts/generate_sequences.py --output data/sequences/
"""

import argparse
import json
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from mit_benchmark.sequences.generator import SequenceGenerator
from mit_benchmark.utils.config import SEQUENCES_DIR, CLASS_SIZES, RANDOM_SEED


def main():
    parser = argparse.ArgumentParser(
        description="Generate MIT benchmark sequences"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=SEQUENCES_DIR,
        help="Output directory for sequences",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=RANDOM_SEED,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--format",
        choices=["json", "fasta", "both"],
        default="both",
        help="Output format",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("MIT Benchmark Sequence Generation")
    print("=" * 60)
    print(f"Output directory: {args.output}")
    print(f"Random seed: {args.seed}")
    print(f"Output format: {args.format}")
    print()

    # Create generator
    generator = SequenceGenerator(seed=args.seed)

    # Generate all classes
    print("Generating sequences for all classes...")
    all_sequences = generator.generate_all_classes()

    # Print summary
    print()
    print("Class Summary:")
    print("-" * 40)
    total = 0
    for class_label, sequences in all_sequences.items():
        n = len(sequences)
        total += n
        if sequences:
            example = sequences[0]
            print(f"  Class {class_label} ({example.class_name}): {n} sequences")
    print("-" * 40)
    print(f"  Total: {total} sequences")
    print()

    # Validate sequences
    print("Validating sequences...")
    for class_label, sequences in all_sequences.items():
        for seq in sequences:
            assert len(seq.sequence) == 100, f"Sequence {seq.id} has wrong length"
            assert all(c in "ACGT" for c in seq.sequence), f"Invalid nucleotides in {seq.id}"
    print("  All sequences validated ✓")
    print()

    # Save sequences
    args.output.mkdir(parents=True, exist_ok=True)

    if args.format in ["json", "both"]:
        json_path = args.output / "all_sequences.json"
        generator.save_sequences(all_sequences, json_path, format="json")
        print(f"Saved JSON to: {json_path}")

    if args.format in ["fasta", "both"]:
        fasta_path = args.output / "all_sequences.fasta"
        generator.save_sequences(all_sequences, fasta_path, format="fasta")
        print(f"Saved FASTA to: {fasta_path}")

    # Save class-specific files
    for class_label, sequences in all_sequences.items():
        class_path = args.output / f"class_{class_label}.json"
        generator.save_sequences({class_label: sequences}, class_path, format="json")

    print()
    print("=" * 60)
    print("Sequence generation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
