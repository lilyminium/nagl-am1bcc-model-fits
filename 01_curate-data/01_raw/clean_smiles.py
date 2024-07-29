import typing
import click
import tqdm
from click_option_group import optgroup

def single_clean(smiles: str):
    from openff.toolkit.topology import Molecule
    elements = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}

    split = smiles.split(".")
    if len(split) > 1:
        smiles = max(split, key=len)

    try:
        offmol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    except BaseException as e:
        return None
    mol_elements = set([atom.symbol for atom in offmol.atoms])
    if not mol_elements.issubset(elements):
        return None
    try:
        return offmol.to_smiles()
    except BaseException as e:
        return None


def batch_clean(smiles: tuple[str, ...]) -> list[str]:
    valid = []
    for smi in tqdm.tqdm(smiles, desc="cleaning"):
        cleaned = single_clean(smi)
        if cleaned is not None:
            valid.append(cleaned)
    return valid


@click.command()
@click.option(
    "--input-file",
    help="The path to the input molecules: smiles file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--output-file",
    help="The path to save the filtered molecules to. This should be a smiles file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
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
def clean(
    input_file: str,
    output_file: str,
    worker_type: typing.Literal["lsf", "local"] = "local",
    queue: str = "cpuqueue",
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from openff.nagl.utils._parallelization import batch_distributed
    from dask import distributed

    print("starting to filter")

    with open(input_file, "r") as f:
        all_smiles = [x.strip() for x in tqdm.tqdm(f.readlines(), desc="Reading input")]
    
    clean_smiles = []

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
        futures = list(batcher(batch_clean))
        for i, future in tqdm.tqdm(
            enumerate(distributed.as_completed(futures, raise_errors=False)),
            total=len(futures),
            desc="Cleaning SMILES batches",
        ):
            batch = future.result()
            clean_smiles += batch
    
    print(f"Filtered {len(all_smiles)} SMILES to {len(clean_smiles)} SMILES")
    with open(output_file, "w") as f:
        f.write("\n".join(clean_smiles))


if __name__ == "__main__":
    clean()
