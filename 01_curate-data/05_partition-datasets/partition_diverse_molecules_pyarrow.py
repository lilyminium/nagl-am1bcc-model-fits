import pathlib
from typing import Literal, List

import click
import tqdm
from click_option_group import optgroup


def load_smiles(input_files):
    import pyarrow.dataset as ds

    all_smiles = set()
    for input_file in tqdm.tqdm(input_files, desc="loading smiles"):
        dataset = ds.dataset(input_file, format="parquet")
        smiles = dataset.to_table(columns=["mapped_smiles"]).to_pydict()["mapped_smiles"]
        all_smiles |= set(smiles)

    return all_smiles

def calculate_fingerprint(smiles: str):
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect

    mol = Chem.MolFromSmiles(smiles)
    fp = GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
    return (smiles, fp)

def batch_calculate_fingerprint(batch_smiles: list[str]):
    rows = []
    for smi in tqdm.tqdm(batch_smiles, desc="Calculating fingerprints"):
        error = None
        try:
            result = calculate_fingerprint(smi)
        except BaseException as e:
            error = str(e)
        rows.append((result, error))
    return rows

@click.command()
@click.option(
    "--input",
    "input_files",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    multiple=True,
    help="File containing molecules (Pyarrow database).",
)
@click.option(
    "--output-base-name",
    type=str,
    help="Base name for output directories.",
)
@click.option(
    "--training-fraction",
    type=float,
    default=0.8,
    help="Fraction of molecules to use for training.",
)
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-workers",
    help="The number of workers to distribute the labelling across. Use -1 to request "
    "one worker per batch.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.option(
    "--worker-type",
    help="The type of worker to distribute the labelling across.",
    type=click.Choice(["lsf", "local", "slurm"]),
    default="local",
    show_default=True,
)
@optgroup.option(
    "--batch-size",
    help="The number of molecules to processes at once on a particular worker.",
    type=int,
    default=500,
    show_default=True,
)
@optgroup.group("LSF configuration", help="Options to configure LSF workers.")
@optgroup.option(
    "--memory",
    help="The amount of memory (GB) to request per LSF queue worker.",
    type=int,
    default=3,
    show_default=True,
)
@optgroup.option(
    "--walltime",
    help="The maximum wall-clock hours to request per LSF queue worker.",
    type=int,
    default=2,
    show_default=True,
)
@optgroup.option(
    "--queue",
    help="The LSF queue to submit workers to.",
    type=str,
    default="cpuqueue",
    show_default=True,
)
@optgroup.option(
    "--conda-environment",
    help="The conda environment that LSF workers should run using.",
    type=str,
)
def partition_smiles(
    input_files,
    output_base_name,
    training_fraction: float = 0.8,
    worker_type: Literal["lsf", "local"] = "local",
    queue: str = "cpuqueue",
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from rdkit import Chem
    from openff.nagl.utils._parallelization import batch_distributed
    from dask import distributed
    from rdkit.SimDivFilters import MaxMinPicker
    import pyarrow.dataset as ds
    import pyarrow.compute as pc
    import pyarrow as pa

    all_smiles = sorted(load_smiles(input_files))
    print(f"Loaded {len(all_smiles)} smiles")

    smiles_fingerprints = []

    with batch_distributed(
        all_smiles,
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    ) as batcher:
        futures = list(batcher(batch_calculate_fingerprint))
        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures),
            desc="Saving rows",
        ):
            for result, error in future.result():
                if error is not None:
                    print(error)
                else:
                    smiles_fingerprints.append(result)
    
    fingerprints = [x[1] for x in smiles_fingerprints]
    n_fingerprints = len(fingerprints)

    n_training = int(n_fingerprints * training_fraction)

    picker = MaxMinPicker()
    selected_indices = [
        i
        for i in tqdm.tqdm(
            picker.LazyBitVectorPick(fingerprints, n_fingerprints, n_training, seed=42),
            desc="Selecting training molecules",
        )
    ]

    training_smiles = []
    validation_smiles = []
    for i, (smiles, _) in enumerate(smiles_fingerprints):
        if i in selected_indices:
            training_smiles.append(smiles)
        else:
            validation_smiles.append(smiles)

    training_smiles = set(training_smiles)
    validation_smiles = set(validation_smiles)

    training_directory = pathlib.Path(f"{output_base_name}_training")
    training_directory.mkdir(exist_ok=True, parents=True)
    validation_directory = pathlib.Path(f"{output_base_name}_validation")
    validation_directory.mkdir(exist_ok=True, parents=True)

    for input_file in tqdm.tqdm(input_files):
        dataset = ds.dataset(input_file, format="parquet")
        
        training_dataset = dataset.filter(
            pc.field("mapped_smiles").isin(training_smiles)
        )
        input_name = pathlib.Path(input_file).stem
        training_dataset_path = training_directory / input_name
        ds.write_dataset(training_dataset, training_dataset_path, format="parquet")

        validation_dataset = dataset.filter(
            pc.field("mapped_smiles").isin(validation_smiles)
        )
        validation_dataset_path = validation_directory / input_name
        ds.write_dataset(validation_dataset, validation_dataset_path, format="parquet")



if __name__ == "__main__":
    partition_smiles()
