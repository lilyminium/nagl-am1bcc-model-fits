import click
import tqdm
from typing import List

@click.command()
@click.option(
    "--input-file",
    "input_files",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    multiple=True,
    help="File containing molecules (MoleculeStore SQLITE database).",
)
@click.option(
    "--output-file",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    help="File to write selected molecules to, one on each line (*.smi).",
)
@click.option(
    "--n-min-molecules",
    type=int,
    default=5,
    help="Minimum number of molecules required to select an atom environment.",
)
@click.option(
    "--element-order",
    type=str,
    multiple=True,
    default=["S", "F", "Cl", "Br", "I", "P", "O", "N", "C"],
    help="Order of elements to use when selecting atom environments.",
)
def select_molecules(
    input_files,
    output_file,
    n_min_molecules: int = 5,
    element_order: List[str] = ["S", "F", "Cl", "Br", "I", "P", "O", "N", "C"]
):
    from openff.nagl.storage._store import MoleculeStore
    from openff.nagl._app.partitioner import FingerprintCollection

    all_smiles = set()
    original = 0

    for input_file in tqdm.tqdm(input_files, desc="loading smiles"):
        store = MoleculeStore(input_file)
        smiles = store.get_smiles()
        original += len(smiles)
        all_smiles |= set(smiles)
    

    fingerprints = FingerprintCollection.from_smiles(
        all_smiles,
        exclude_elements=tuple()
    )
    selected = fingerprints.select_atom_environments(
        n_min_molecules=n_min_molecules,
        element_order=element_order
    )

    with open(output_file, "w") as f:
        f.write("\n".join(selected))
    
    print(f"Selected {len(selected)} molecules from {len(all_smiles)} ({original} total, including duplicates)")



if __name__ == "__main__":
    select_molecules()