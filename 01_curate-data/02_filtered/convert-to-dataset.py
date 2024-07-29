import functools
import multiprocessing

import click
import datasets
import datasets.table
import pyarrow
import tqdm

def create_entry(smiles: str, dataset_name: str,):
    from openff.toolkit import Molecule

    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    mapped_smiles = molecule.to_smiles(
        isomeric=True,
        explicit_hydrogens=True,
        mapped=True,
    )
    inchi = molecule.to_inchi()

    return {
        "smiles": smiles,
        "mapped_smiles": mapped_smiles,
        "inchi": inchi,
        "dataset": dataset_name,
    }


@click.command()
@click.option(
    "--input",
    "input_file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
    help=(
        "The path to the input file containing the SMILES to convert. "
        "Each line should contain a single SMILES."
    )
)
@click.option(
    "--dataset",
    "dataset_name",
    type=str,
    required=True,
    help="The name of the dataset that the SMILES are from."
)
@click.option(
    "--output",
    "output_path",
    type=click.Path(exists=False, dir_okay=True, file_okay=True),
    required=True,
    help="The path to the dataset file or directory to save the dataset to."
)
@click.option(
    "--n-processes",
    "n_processes",
    type=int,
    default=1,
    show_default=True,
    help="The number of processes to use to convert the SMILES."
)
def convert(
    input_file: str,
    dataset_name: str,
    output_path: str,
    n_processes: int = 1,
):
    schema = pyarrow.schema(
        [
            ("smiles", pyarrow.string()),
            ("mapped_smiles", pyarrow.string()),
            ("inchi", pyarrow.string()),
            ("dataset", pyarrow.string())
        ]
    )

    with open(input_file, "r") as file:
        smiles = [line.strip() for line in file.readlines()]

    creator = functools.partial(create_entry, dataset_name=dataset_name)

    with multiprocessing.Pool(n_processes) as pool:
        all_entries = list(
            tqdm.tqdm(
                pool.imap_unordered(creator, smiles),
                total=len(smiles),
            )
        )
    
    table = pyarrow.Table.from_pylist(
        all_entries,
        schema=schema,
    )
    dataset = datasets.Dataset(datasets.table.InMemoryTable(table))
    dataset.set_format("torch")
    dataset.save_to_disk(output_path)


if __name__ == "__main__":
    convert()