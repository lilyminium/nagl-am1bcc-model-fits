import os
import pathlib
import typing

import click
from click_option_group import optgroup
import tqdm
import torch

from openff.toolkit import Molecule

def label_molecule(smiles: str) -> torch.tensor:
    from openff.units import unit
    from openff.toolkit.utils import ToolkitRegistry, RDKitToolkitWrapper, AmberToolsToolkitWrapper

    registry = ToolkitRegistry(toolkit_precedence=[
        RDKitToolkitWrapper(),
        AmberToolsToolkitWrapper()
    ])

    offmol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
    offmol.assign_partial_charges("am1bcc", toolkit_registry=registry)
    charges = torch.tensor(
        offmol.partial_charges.m_as(unit.elementary_charge)
    ).flatten().tolist()

    return {
        "mapped_smiles": smiles,
        "am1bcc_charges_ambertools": charges,
    }
    # return charges


def batch_label(
    all_smiles: str,
):
    entries = []
    errors = []
    for smiles in tqdm.tqdm(all_smiles):
        try:
            entries.append(label_molecule(smiles))
        except BaseException as e:
            errors.append((smiles, e))
            # smiles_to_charges[smiles] = torch.tensor([])
    return entries, errors

@click.command()
@click.option(
    "--input",
    "input_dataset_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    required=True,
    help=(
        "The path to the input dataset to label."
        "This must have a `mapped_smiles` column with mapped SMILES."
    ),
)
@click.option(
    "--output",
    "output_dataset_path",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    required=True,
    help="The path to save the labelled dataset to.",
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
def label_dataset(
    input_dataset_path: str,
    output_dataset_path: str,
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
    
    import pyarrow as pa
    import pyarrow.parquet as pq
    import pyarrow.dataset as ds
    import pyarrow.compute as pc
    import datasets

    # really ensure not OpenEye
    os.environ["OE_LICENSE"] = ""

    print(f"input_dataset_path: {input_dataset_path}")

    dataset = datasets.Dataset.load_from_disk(input_dataset_path)
    all_smiles = sorted(
        dataset.data["mapped_smiles"].to_pylist(),
        key=len,
        reverse=True
    )
    print(f"Loaded {len(all_smiles)} SMILES from {input_dataset_path}")

    exclude_smiles = []
    output_directory = pathlib.Path(output_dataset_path)
    output_directory.mkdir(exist_ok=True, parents=True)

    error_directory = output_directory.parent / (str(output_directory.name) + "_errors")
    error_directory.mkdir(exist_ok=True, parents=True)
    start_index = 0

    try:
        print(f"Querying {output_directory}")
        existing = ds.dataset(output_directory)
        if existing.count_rows():
            exclude_smiles.extend(
                existing.to_table(
                    columns=["mapped_smiles"]
                ).to_pydict()["mapped_smiles"]
            )
            start_index = len(existing.files)
    except BaseException as e:
        print(e)
    print(f"Loaded {len(exclude_smiles)} SMILES already calculated")

    error_start_index = 0

    try:
        print(f"Querying {error_directory}")
        existing2 = ds.dataset(error_directory)
        if existing2.count_rows():
            exclude_smiles.extend(
                existing2.to_table(
                    columns=["mapped_smiles"]
                ).to_pydict()["mapped_smiles"]
            )
            error_start_index = len(existing2.files)
    except BaseException as e:
        print(e)


    print(f"Starting from start index {start_index}")
    if exclude_smiles:
        print(f"Excluding {len(exclude_smiles)} SMILES")

        exclude_smiles = set(exclude_smiles)
        all_smiles = sorted(
            set(all_smiles) - exclude_smiles,
            key=len,
            reverse=True
        )

    print(f"Labelling {len(all_smiles)} molecules...")


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
        futures = list(batcher(batch_label))
        for i, future in tqdm.tqdm(
            enumerate(
                distributed.as_completed(futures, raise_errors=False),
                start=error_start_index
            ),
            total=len(futures),
            desc="Updating charges",
        ):
            batch_entries, batch_errors = future.result()
            batch_table = pa.Table.from_pylist(batch_entries)
            print(batch_table.schema)

            table_path = output_directory / f"batch-{start_index:04d}.parquet"
            if len(batch_entries):
                pq.write_table(batch_table, table_path)
                print(f"Wrote {len(batch_entries)} to {table_path}")
                start_index += 1

            error_entries = [
                {
                    "mapped_smiles": smiles,
                    "error": str(error)
                }
                for smiles, error in batch_errors
            ]
            error_table = pa.Table.from_pylist(error_entries)
            error_path = error_directory / f"batch-{i:04d}.parquet"
            pq.write_table(error_table, error_path)


if __name__ == "__main__":
    label_dataset()