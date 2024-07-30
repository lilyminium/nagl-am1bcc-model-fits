import os
import functools
import pathlib
import yaml
import io
import gc
import pickle
import glob
import torch

import typing
import tqdm
import click
from click_option_group import optgroup
import yaml
import pytorch_lightning as pl
import torch

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.dataset as ds

from openff.nagl.toolkits.openff import capture_toolkit_warnings

SCHEMA = pa.schema([
    pa.field("pickled", pa.binary()),
])

def batch_convert_original(
    batch_number,
    source=None,
    atom_features=[],
    bond_features=[],
    columns=[],
):
    from openff.nagl.nn._dataset import DGLMoleculeDatasetEntry
    from openff.nagl.toolkits.openff import capture_toolkit_warnings

    dataset = ds.dataset(source, format="parquet")
    all_files = sorted(dataset.files)

    with capture_toolkit_warnings():
        all_rows = []
        for num in batch_number:
            table = ds.dataset(all_files[num], format="parquet").to_table()
            rows = table.to_pylist()
            all_rows.extend([
                {k: v for k, v in x.items() if k in columns}
                for x in rows
            ])
        rows_ = []
        for row in tqdm.tqdm(all_rows, desc="Converting rows"):
            try:
                entry = DGLMoleculeDatasetEntry._from_unfeaturized_pyarrow_row(row, atom_features, bond_features)
            except BaseException as e:
                print(e)
                continue
            f = io.BytesIO()
            pickle.dump(entry, f)
            rows_.append(f.getvalue())

        batch_with_features = pa.RecordBatch.from_arrays(
            [rows_],
            schema=SCHEMA
        )
        return batch_with_features
    



@click.command()
@click.option(
    "--model-config-file",
    "model_config_files",
    help="The path to a YAML configuration file for the model.",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    default=tuple(),
    multiple=True,
)
@click.option(
    "--data-config-file",
    "data_config_files",
    help="The path to a YAML configuration file for the data.",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    default=tuple(),
    multiple=True,
)
@click.option(
    "--data-cache-directory",
    help="Path to cached data",
    type=click.Path(exists=False, file_okay=False, dir_okay=True),
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
def train_model(
    data_cache_directory: str,
    model_config_files: typing.Tuple[str, ...] = tuple(),
    data_config_files: typing.Tuple[str, ...] = tuple(),
    learning_rate: float = 0.001,
    n_processes: int = 4,
    worker_type: str = "local",
    queue: str = "cpuqueue",
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from openff.nagl.config.training import TrainingConfig
    from openff.nagl.training.training import DGLMoleculeDataModule
    from openff.nagl.nn._dataset import _get_hashed_arrow_dataset_path, _LazyDGLMoleculeDataset
    from openff.nagl.utils._parallelization import batch_distributed

    from dask import distributed
    print("starting")


    model_config_kwargs = {}
    for model_file in model_config_files:
        with open(model_file, "r") as f:
            model_config_kwargs.update(yaml.load(f, Loader=yaml.FullLoader))
    
    data_config_kwargs = {}
    for data_file in data_config_files:
        with open(data_file, "r") as f:
            data_config_kwargs.update(yaml.load(f, Loader=yaml.FullLoader))

    for v in data_config_kwargs.values():
        v["cache_directory"] = os.path.abspath(data_cache_directory)

    optimizer_config_kwargs = {
        "optimizer": "Adam",
        "learning_rate": learning_rate
    }

    training_config = TrainingConfig(
        model=model_config_kwargs,
        data=data_config_kwargs,
        optimizer=optimizer_config_kwargs,
    )

    data_module = DGLMoleculeDataModule(
        training_config,
        n_processes=n_processes
    )

    atom_features = training_config.model.atom_features
    bond_features = training_config.model.bond_features

    for dataset in ("training", "validation", "test"):
        print(f"---- {dataset} ----")
        data_config = getattr(data_module.config.data, dataset)
        cachedir = pathlib.Path(data_config.cache_directory)
        if data_config.sources:
            columns = data_config.get_required_target_columns()
            if "mapped_smiles" not in columns:
                columns = ["mapped_smiles"] + columns

            for source in data_config.sources:
                file_path = _get_hashed_arrow_dataset_path(
                    source,
                    atom_features,
                    bond_features,
                    columns,
                ).with_suffix(".arrow")
                output_file_path = cachedir / file_path
                if output_file_path.exists():
                    print(f"Skipping {source} as it has already been processed.")
                    continue

                print(f"Processing {source} to {output_file_path}.")

                dataset = ds.dataset(source, format="parquet")
                n_files = len(dataset.files)

                batch_func = functools.partial(
                    batch_convert_original,
                    source=source,
                    atom_features=training_config.model.atom_features,
                    bond_features=training_config.model.bond_features,
                    columns=columns,
                )

                with pa.OSFile(str(output_file_path), "wb") as sink:
                    with pa.ipc.new_file(sink, SCHEMA) as writer:
                        with batch_distributed(
                            # batch_func,
                            list(range(n_files)),
                            n_workers=n_workers,
                            batch_size=batch_size,
                            worker_type=worker_type,
                            queue=queue,
                            conda_environment=conda_environment,
                            memory=memory,
                            walltime=walltime,
                        ) as batcher:

                            futures = list(
                                batcher(
                                    batch_convert_original,
                                    source=source,
                                    atom_features=training_config.model.atom_features,
                                    bond_features=training_config.model.bond_features,
                                    columns=columns,
                                )
                            )
                            for future in tqdm.tqdm(
                                distributed.as_completed(futures, raise_errors=False),
                                total=len(futures),
                                desc="Saving rows",
                            ):
                                batch = future.result()
                                writer.write_batch(batch)
            


if __name__ == "__main__":
    print("before calling train")
    train_model()
