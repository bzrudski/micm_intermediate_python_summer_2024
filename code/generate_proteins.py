import enum
import os.path
import random
from typing import Optional


class AminoAcid(enum.Flag):
    """Flag class representing amino acids."""

    F = enum.auto()
    L = enum.auto()
    I = enum.auto()
    M = enum.auto()
    V = enum.auto()
    P = enum.auto()
    A = enum.auto()
    W = enum.auto()
    G = enum.auto()
    S = enum.auto()
    T = enum.auto()
    Y = enum.auto()
    Q = enum.auto()
    N = enum.auto()
    C = enum.auto()
    D = enum.auto()
    E = enum.auto()
    H = enum.auto()
    K = enum.auto()
    R = enum.auto()
    NON_POLAR = F | L | I | M | V | P | A | W | G
    POLAR = S | T | Y | Q | N | C
    ACIDIC = D | E
    BASIC = H | K | R


AMINO_ACID_PROPERTIES = {
    "NON_POLAR": ["F", "L", "I", "M", "V", "P", "A", "W", "G"],
    "POLAR": ["S", "T", "Y", "Q", "N", "C"],
    "ACIDIC": ["D", "E"],
    "BASIC": ["H", "K", "R"],
}


class PeptideSequence:
    """Peptide sequence."""

    amino_acids: list[AminoAcid]

    def __init__(self, amino_acids: list[AminoAcid]) -> None:
        self.amino_acids = amino_acids.copy()

    def __str__(self) -> str:
        self_string = ""
        for aa in self.amino_acids:
            self_string += aa.name

        return self_string

    def __len__(self) -> int:
        return len(self.amino_acids)


def generate_random_peptide_sequence(length: int = 100) -> PeptideSequence:
    """Generate a random peptide sequence of a given length."""
    amino_acids = random.choices(list(AminoAcid), k=length)

    return PeptideSequence(amino_acids)


def mutate_peptide_sequence(
    seq: PeptideSequence, n: int, weights: Optional[dict[AminoAcid, float]] = None
) -> PeptideSequence:
    """Mutate a peptide sequence with weights."""
    new_amino_acid_list = seq.amino_acids.copy()

    # Randomly select the positions to change
    seq_length = len(new_amino_acid_list)
    for _ in range(n):
        # Select position to mutate
        i = random.randint(0, seq_length - 1)

        if weights is None:
            weights = {aa: 1 for aa in AminoAcid}

        # Get the new single amino acid or amino acid category
        new_amino_acid_category = random.choices(
            list(weights.keys()), weights=list(weights.values())
        )[0]

        # Randomly select a new amino acid from the category.
        # If a single amino acid, it will automatically be chosen.
        new_amino_acid = random.choice(list(new_amino_acid_category))

        new_amino_acid_list[i] = new_amino_acid

    return PeptideSequence(new_amino_acid_list)


class Species(enum.Enum):
    """Species for sequences."""

    HUMAN = "HUMAN"
    CHIMP = "CHIMP"
    MOUSE = "MOUSE"
    PIG = "PIG"
    YEAST = "YEAST"


class DiseaseState(enum.Enum):
    """Disease state."""

    HEALTHY = "HEALTHY"
    DISEASE = "DISEASE"


def generate_single_species_dataset(
    length: int,
    number_of_individuals: int,
    number_of_mutations: int,
    weights: dict[AminoAcid, float],
    species: Species = Species.HUMAN,
    disease_state: DiseaseState = DiseaseState.HEALTHY,
    output_folder: str = ".",
    initial_sequence: Optional[PeptideSequence] = None,
) -> PeptideSequence:
    """Generate a sample dataset."""

    # Create the base sequence
    base_sequence = initial_sequence or generate_random_peptide_sequence(length)

    # Randomly mutate the base
    random_sequences = [
        mutate_peptide_sequence(base_sequence, number_of_mutations, weights)
        for _ in range(number_of_individuals)
    ]

    # Write the mutated sequences to file:
    with open(
        os.path.join(output_folder, f"{species.value}_{disease_state.value}.fasta"), "w"
    ) as f:
        for i, seq in enumerate(random_sequences):
            f.write(f"> SEQ -- {species.value}_{disease_state.value}_{i}\n")
            f.write(str(seq))
            f.write("\n")

    return base_sequence


def generate_datasets_for_workshop():
    """Generate the artificial peptide sequences."""

    params = [
        (
            Species.HUMAN,
            DiseaseState.HEALTHY,
            {AminoAcid.ACIDIC: 2, AminoAcid.NON_POLAR: 1},
        ),
        (
            Species.HUMAN,
            DiseaseState.DISEASE,
            {AminoAcid.ACIDIC: 1, AminoAcid.NON_POLAR: 4},
        ),
        (
            Species.CHIMP,
            DiseaseState.HEALTHY,
            {AminoAcid.ACIDIC: 2, AminoAcid.NON_POLAR: 1, AminoAcid.POLAR: 1},
        ),
        (
            Species.CHIMP,
            DiseaseState.DISEASE,
            {AminoAcid.ACIDIC: 1, AminoAcid.NON_POLAR: 1, AminoAcid.POLAR: 5},
        ),
        (
            Species.PIG,
            DiseaseState.HEALTHY,
            {AminoAcid.ACIDIC: 1, AminoAcid.BASIC: 3, AminoAcid.POLAR: 1},
        ),
        (
            Species.PIG,
            DiseaseState.DISEASE,
            {AminoAcid.ACIDIC: 4, AminoAcid.BASIC: 3, AminoAcid.POLAR: 1},
        ),
    ]

    initial_sequence = None
    for species, disease_state, weights in params:
        initial_sequence = generate_single_species_dataset(
            length=100, number_of_individuals=50, number_of_mutations=20,
            weights=weights, species=species, disease_state=disease_state,
            output_folder="./datasets", initial_sequence=initial_sequence
        )


if __name__ == "__main__":
    generate_datasets_for_workshop()