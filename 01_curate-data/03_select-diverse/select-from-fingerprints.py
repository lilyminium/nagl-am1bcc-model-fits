from collections import defaultdict, Counter
import typing

import click
import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem

DEFAULT_ELEMENTS = ("S", "F", "Cl", "Br", "I", "P", "O", "N", "C", "H")

class AtomFingerprint(typing.NamedTuple):
    """A class to store the fingerprint of an atom."""
    element: str
    fingerprint: str


class MoleculeAtomFingerprints(typing.NamedTuple):
    """A class to store the atom fingerprints of a molecule."""

    smiles: str
    fingerprints: tuple[AtomFingerprint, ...]

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        include_elements=DEFAULT_ELEMENTS,
    ):
        from openff.toolkit import Molecule
        
        offmol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
        rdmol = offmol.to_rdkit()
        atom_fingerprints = set()

        for rdatom in rdmol.GetAtoms():
            symbol = rdatom.GetSymbol()
            # skip including fingerprint if not in include_elements
            if symbol not in include_elements:
                continue
            fingerprint: str = AllChem.GetHashedAtomPairFingerprintAsBitVect(
                rdmol,
                maxLength=2,
                nBits=2048,
                nBitsPerEntry=4,
                fromAtoms=[rdatom.GetIdx()],
            ).ToBase64()
            atom_fingerprint = AtomFingerprint(symbol, fingerprint)
            atom_fingerprints.add(atom_fingerprint)
        return cls(smiles, tuple(sorted(atom_fingerprints)))
            


def select_least_represented_molecule(
    molecule_fingerprints: set[MoleculeAtomFingerprints],
    atom_fingerprint_counts: Counter,
) -> MoleculeAtomFingerprints:
    """Select the molecule with the least represented atom fingerprint."""
    already_represented = set(atom_fingerprint_counts.keys())
    n_unrepresented = Counter({
        molecule_fingerprint: len(
            set(molecule_fingerprint.fingerprints) - already_represented
        )
        for molecule_fingerprint in molecule_fingerprints
    })
    return n_unrepresented.most_common(1)[0][0]


@click.command()
@click.option(
    "--input",
    "smiles_datasets",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    required=True,
    multiple=True,
    help=(
        "The path to the input dataset to select from."
        "This must have a `mapped_smiles` column with mapped SMILES."
    ),
)
@click.option(
    "--output",
    "output_dataset",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    required=True,
    help="The path to save the selected dataset to.",
)
@click.option(
    "--exclude",
    "exclude_datasets",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    multiple=True,
    default=[],
)
@click.option(
    "--element-order",
    "element_order",
    type=str,
    multiple=True,
    default=DEFAULT_ELEMENTS,
    show_default=True,
    help="The order to select elements in.",
)
@click.option(
    "--n-min-molecules",
    "n_min_molecules",
    type=int,
    default=4,
    show_default=True,
    help="The minimum number of molecules to select for each atom environment.",
)
def select(
    smiles_datasets: list[str],
    output_dataset: str,
    exclude_datasets: list[str] = list(),
    element_order: tuple[str, ...] = DEFAULT_ELEMENTS,
    n_min_molecules: int = 4,
    column_name: str = "mapped_smiles"
):
    import datasets
    import pyarrow as pa
    import pyarrow.compute as pc

    print("Elements:", element_order)

    # load smiles from datasets
    all_smiles = []
    for dataset_path in smiles_datasets:
        print(f"Loading from {dataset_path}")
        dataset = datasets.Dataset.load_from_disk(dataset_path)
        all_smiles.extend(dataset.data[column_name].to_pylist())

    print(f"Loaded {len(all_smiles)} original SMILES")

    exclude_smiles = []
    for dataset_path in exclude_datasets:
        print(f"Loading exclusions from {dataset_path}")
        dataset = datasets.Dataset.load_from_disk(dataset_path)
        exclude_smiles.extend(dataset.data[column_name].to_pylist())
    
    print(f"Loaded {len(exclude_smiles)} exclusions")

    all_smiles = [
        x for x in all_smiles
        if x not in exclude_smiles
    ]
    
    print(f"Selecting from {len(all_smiles)} SMILES")

    # convert smiles to atom fingerprints
    all_molecule_fingerprints = []
    for smiles in tqdm.tqdm(
        all_smiles,
        desc="Converting SMILES to atom fingerprints",
    ):
        try:
            all_molecule_fingerprints.append(
                MoleculeAtomFingerprints.from_smiles(smiles, element_order)
            )
        except BaseException as e:
            print(f"Error converting {smiles}: {e}")
            continue

    # atom fingerprints -> molecules
    atom_fingerprints_to_molecules = defaultdict(set)
    for molecule_fingerprints in tqdm.tqdm(
        all_molecule_fingerprints,
        desc="Grouping molecules by atom fingerprints",
    ):
        for atom_fingerprint in molecule_fingerprints.fingerprints:
            atom_fingerprints_to_molecules[atom_fingerprint].add(
                molecule_fingerprints
            )

    # start selecting output molecules
    selected_molecule_fingerprints = set()

    # select molecules with rare environments first
    selected_atom_environments = []
    for atom_fingerprint, molecules in tqdm.tqdm(
        atom_fingerprints_to_molecules.items(),
        desc="Selecting molecules with rare environments",
    ):
        if len(molecules) <= n_min_molecules:
            selected_molecule_fingerprints |= molecules
            selected_atom_environments.append(atom_fingerprint)

    # count number of molecules for each atom environment
    atom_fingerprint_counts = Counter()
    for molecule_fingerprints in selected_molecule_fingerprints:
        for atom_fingerprint in molecule_fingerprints.fingerprints:
            atom_fingerprint_counts[atom_fingerprint] += 1


    # select molecules with most underrepresented environments
    # for each atom fingerprint environment,
    # in specified element order
    remaining_atom_environments = (
        set(atom_fingerprints_to_molecules)
        - set(selected_atom_environments)
    )
    element_order = list(element_order)
    ordered_atom_fingerprints = sorted(
        # atom_fingerprint_counts.keys(),
        remaining_atom_environments,
        key=lambda x: (element_order.index(x[0]), x[1]),
    )
    for atom_fingerprint in tqdm.tqdm(ordered_atom_fingerprints):
        molecule_fingerprints = (
            atom_fingerprints_to_molecules[atom_fingerprint]
            - selected_molecule_fingerprints
        )
        n_remaining = n_min_molecules - atom_fingerprint_counts[atom_fingerprint]
        # if len(molecule_fingerprints) <= n_min_molecules:
        if len(molecule_fingerprints) <= n_remaining:
            selected_molecule_fingerprints |= molecule_fingerprints
            for molecule_fingerprint in molecule_fingerprints:
                for atom_fingerprint in molecule_fingerprint.fingerprints:
                    atom_fingerprint_counts[atom_fingerprint] += 1
            continue

        # for _ in range(n_min_molecules):
        for _ in range(n_remaining):
            molecule_fingerprint = select_least_represented_molecule(
                molecule_fingerprints,
                atom_fingerprint_counts,
            )
            selected_molecule_fingerprints.add(molecule_fingerprint)
            for atom_fingerprint in molecule_fingerprint.fingerprints:
                atom_fingerprint_counts[atom_fingerprint] += 1
            molecule_fingerprints -= {molecule_fingerprint}

    # convert selected molecules to dataset
    selected_smiles = [
        {"mapped_smiles": molecule_fingerprint.smiles}
        for molecule_fingerprint in selected_molecule_fingerprints
    ]

    print(f"Selected {len(selected_smiles)} molecules from {len(all_smiles)} original.")
    
    table = pa.Table.from_pylist(selected_smiles)
    dataset = datasets.Dataset(datasets.table.InMemoryTable(table))
    dataset.save_to_disk(output_dataset)


if __name__ == "__main__":
    select()
