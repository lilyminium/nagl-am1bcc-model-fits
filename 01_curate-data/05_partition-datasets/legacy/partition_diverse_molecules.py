from typing import Literal, List

import click
import tqdm
from click_option_group import optgroup


def load_smiles(input_files):
    from openff.nagl.storage._store import MoleculeStore

    all_smiles = set()
    for input_file in tqdm.tqdm(input_files, desc="loading smiles"):
        if input_file.endswith(".sqlite"):
            store = MoleculeStore(input_file)
            smiles = store.get_smiles()
        else:
            with open(input_file, "r") as f:
                smiles = [x.strip() for x in f.readlines()]
        
        all_smiles |= set(smiles)

    return all_smiles

def calculate_fingerprint(smiles: str):
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect

    mol = Chem.MolFromSmiles(smiles)
    fp = GetMorganFingerprintAsBitVect(mol, 3, nBits=2048)
    return (smiles, fp)

@click.command()
@click.option(
    "--input-file",
    "input_files",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    multiple=True,
    help="File containing molecules (MoleculeStore SQLITE database).",
)
@click.option(
    "--output-base-name",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    help="Base name for output files.",
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
    from dask import distributed
    from rdkit.SimDivFilters import MaxMinPicker
    from openff.nagl._app.distributed import Manager
    from openff.nagl._cli.utils import as_batch_function_with_captured_errors

    all_smiles = sorted(load_smiles(input_files))
    print(f"Loaded {len(all_smiles)} smiles")

    manager = Manager(
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    )
    manager.set_entries(all_smiles)

    batch_func = as_batch_function_with_captured_errors(
        calculate_fingerprint,
        desc="Calculating fingerprints",
    )

    smiles_fingerprints = []
    with manager:
        futures = manager.submit_to_client(batch_func)
        for future in tqdm.tqdm(
                distributed.as_completed(futures),
                total=manager.n_batches,
                desc="Waiting for fingerprints",
            ):
                for result, error in future.result():
                    if error is not None:
                        print(error)
                    else:
                        smiles_fingerprints.append(result)
            
                future.release()
    
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

    training_file = f"{output_base_name}_training.smi"
    validation_file = f"{output_base_name}_validation.smi"

    with open(training_file, "w") as f:
        f.write("\n".join(training_smiles))

    print(f"Training set: {len(training_smiles)} molecules ({training_file})")

    with open(validation_file, "w") as f:
        f.write("\n".join(validation_smiles))

    print(f"Validation set: {len(validation_smiles)} molecules ({validation_file})")


if __name__ == "__main__":
    partition_smiles()
