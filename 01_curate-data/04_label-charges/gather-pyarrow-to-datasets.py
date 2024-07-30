import pathlib
import typing

import click
from click_option_group import optgroup
import tqdm
import torch


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
    "--charges",
    "charge_dataset_path",
    type=click.Path(exists=True, dir_okay=True, file_okay=True),
    required=True,
    help=(
        "The path to the charge dataset."
    ),
)
@click.option(
    "--output",
    "output_dataset_path",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    required=True,
    help="The path to save the labelled dataset to.",
)
@click.option(
    "--column",
    "column_name",
    type=str,
    default="am1bccelf10_charges_openeye",
)
def label_dataset(
    input_dataset_path: str,
    charge_dataset_path: str,
    output_dataset_path: str,
    column_name: str = "am1bccelf10_charges_openeye",
):
    from openff.nagl.utils._parallelization import batch_distributed
    from dask import distributed
    
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.compute as pc
    import datasets

    dataset = datasets.Dataset.load_from_disk(input_dataset_path)
    all_smiles = sorted(
        dataset.data["mapped_smiles"].to_pylist(),
        key=len,
        reverse=True
    )

    print(f"Labelling {len(all_smiles)} molecules...")

    charge_datasets = ds.dataset(charge_dataset_path)
    print(f"Loaded {charge_datasets.count_rows()} charges...")
    all_charges = {
        row["mapped_smiles"]: row[column_name]
        for row in charge_datasets.to_table().to_pylist()
    }
    print(f"Loaded {len(all_charges)} unique charges...")

    expression = pc.field("mapped_smiles").isin(all_charges)
    rows = dataset.data.to_pylist()
    seen_smiles = set()
    unique_rows = []
    for row in tqdm.tqdm(rows):
        if row["mapped_smiles"] in all_charges:
            if row["mapped_smiles"] in seen_smiles:
                continue
            unique_rows.append(row)
            seen_smiles.add(row["mapped_smiles"])

    table = pa.Table.from_pylist(unique_rows)
    print(table.schema)
    new_all_smiles = table.to_pydict()["mapped_smiles"]

    new_column = [
        all_charges[smiles]
        for smiles in new_all_smiles
    ]
    field = pa.field(column_name, pa.list_(pa.float64()))
    table = table.append_column(field, [new_column])

    # huggingface_table = datasets.table.InMemoryTable(table)
    new_dataset = datasets.Dataset(table)
    new_dataset.set_format("torch")
    new_dataset.save_to_disk(output_dataset_path)

    print(f"Done labelling {len(all_smiles)} molecules.")
    print(f"Saved to {output_dataset_path}.")

if __name__ == "__main__":
    label_dataset()
