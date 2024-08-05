import pathlib
from typing import Tuple, List, Literal

import click
import tqdm
from click_option_group import optgroup


def generate_conformers_and_label(smiles, ambertools=False, openeye=True):
    import numpy as np
    from openff.units import unit
    from openff.toolkit import Molecule
    from openff.toolkit.utils import OpenEyeToolkitWrapper, AmberToolsToolkitWrapper
    from openff.nagl.toolkits.openff import capture_toolkit_warnings
    from openff.recharge.grids import GridGenerator, MSKGridSettings

    with capture_toolkit_warnings():
        offmol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

    print(smiles)

    try:
        offmol.generate_conformers(
            n_conformers=1000,
            rms_cutoff=0.05 * unit.angstrom,
            make_carboxylic_acids_cis=True,
            toolkit_registry=OpenEyeToolkitWrapper(),
        )
    except:
        raise ValueError(smiles)
    offmol.apply_elf_conformer_selection(
        limit=10,
        toolkit_registry=OpenEyeToolkitWrapper(),
    )

    if openeye:
        openeye_charges = []
    if ambertools:
        ambertools_charges = []

    errors = []

    for conf in tqdm.tqdm(offmol.conformers, desc="Computing charges"):
        if ambertools:
            try:
                offmol.assign_partial_charges(
                    "am1bcc",
                    use_conformers=[conf],
                    toolkit_registry=AmberToolsToolkitWrapper(),
                )
                ambertools_charges.append(
                    offmol.partial_charges.m_as(unit.elementary_charge)
                )
            except BaseException as e:
                errors.append(str(e))

        if openeye:
            offmol.assign_partial_charges(
                "am1bcc",
                use_conformers=[conf],
                toolkit_registry=OpenEyeToolkitWrapper(),
            )
            openeye_charges.append(
                offmol.partial_charges.m_as(unit.elementary_charge)
            )

    if openeye:
        openeye_charge = np.array(openeye_charges).mean(axis=0)
    if ambertools:
        if len(ambertools_charges):
            ambertools_charge = np.array(ambertools_charges).mean(axis=0)
        else:
            raise ValueError(f"No AmberTools conformers for {smiles}")

    all_conformers = []
    all_am1bcc_esps_openeye = []
    all_am1bcc_esps_ambertools = []
    all_invs = []
    all_n_esps = []
    all_am1bcc_dipoles_openeye = []
    all_am1bcc_dipoles_ambertools = []

    settings = MSKGridSettings()

    for conf in offmol.conformers:
        coordinates = conf.m_as(unit.angstrom)
        grid = GridGenerator.generate(offmol, conf, settings)

        displacement = grid[:, None, :] - conf[None, :, :]
        distance = (displacement ** 2).sum(axis=-1) ** 0.5
        distance = distance.m_as(unit.bohr)
        inv_distance = 1 / distance

        all_conformers.extend(coordinates.flatten())
        all_n_esps.append(len(inv_distance))
        all_invs.extend(inv_distance.flatten())
        
        if openeye:
            am1bcc_esp_openeye = inv_distance @ openeye_charge
            am1bcc_dipole_openeye = openeye_charge @ coordinates
            all_am1bcc_esps_openeye.extend(am1bcc_esp_openeye)
            all_am1bcc_dipoles_openeye.extend(am1bcc_dipole_openeye)

        if ambertools:
            am1bcc_esp_ambertools = inv_distance @ ambertools_charge
            am1bcc_dipole_ambertools = ambertools_charge @ coordinates
            all_am1bcc_esps_ambertools.extend(am1bcc_esp_ambertools)
            all_am1bcc_dipoles_ambertools.extend(am1bcc_dipole_ambertools)

    mapped_smiles = offmol.to_smiles(mapped=True)

    row = [
        smiles,
        mapped_smiles,
        np.array(all_conformers),
        len(offmol.conformers),
        np.array(all_n_esps),
        np.array(all_invs),
    ]

    if openeye:
        row += [
            np.array(openeye_charge),
            np.array(all_am1bcc_esps_openeye),
            np.array(all_am1bcc_dipoles_openeye),
        ]

    if ambertools:
        row += [
            np.array(ambertools_charge),
            np.array(all_am1bcc_esps_ambertools),
            np.array(all_am1bcc_dipoles_ambertools),
        ]

    if len(errors):
        error = f"{smiles}\n" + "\n".join(errors)
    else:
        error = ""

    return row, error

def batch_label(smiles_batch, ambertools=False, openeye=True):
    import pyarrow as pa

    schema_fields = [
        pa.field("smiles", pa.string()),
        pa.field("mapped_smiles", pa.string()),
        pa.field("conformers", pa.list_(pa.float32())),
        pa.field("n_conformers", pa.int64()),
        pa.field("esp_lengths", pa.list_(pa.int64())),
        pa.field("esp_grid_inverse_distances", pa.list_(pa.float32())),        
    ]

    if openeye:
        schema_fields += [
            pa.field("am1bcc_charges", pa.list_(pa.float32())),
            pa.field("am1bcc_esps", pa.list_(pa.float32())),
            pa.field("am1bcc_dipoles", pa.list_(pa.float32())),
        ]
    if ambertools:
        schema_fields += [
            pa.field("am1bcc_charges_ambertools", pa.list_(pa.float32())),
            pa.field("am1bcc_esps_ambertools", pa.list_(pa.float32())),
            pa.field("am1bcc_dipoles_ambertools", pa.list_(pa.float32())),
        ]

    schema = pa.schema(schema_fields)

    rows = []
    error = ""
    for smi in tqdm.tqdm(smiles_batch, desc="Labelling SMILES"):
        try:
            output_row, error_ = generate_conformers_and_label(
                smi,
                ambertools=ambertools,
                openeye=openeye,
            )
            error += error_
        except Exception as e:
            error += f"{smi}\n" + str(e)
        else:
            rows.append(output_row)

    
    batch = pa.RecordBatch.from_arrays(
        list(zip(*rows)),
        schema=schema
    )
    
    return batch, error

@click.command()
@click.option(
    "--dataset-file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="File containing FreeSolv rows.",
)
@click.option(
    "--output-path",
    help="The path to the output database to save the labelled molecules in.",
    type=click.Path(exists=False, file_okay=True, dir_okay=True),
    required=True,
)
@click.option(
    "--ambertools/--no-ambertools",
    help="Whether to use AmberTools to compute partial charges.",
    default=False,
    show_default=True,
    is_flag=True,
)
@click.option(
    "--openeye/--no-openeye",
    help="Whether to use OpenEye to compute partial charges.",
    default=True,
    show_default=True,
    is_flag=True,
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

def generate_all(
    dataset_file: str,
    output_path: str,
    ambertools: bool = False,
    openeye: bool = True,
    worker_type: Literal["lsf", "local"] = "local",
    queue: str = "cpuqueue",
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from openff.nagl.utils._parallelization import batch_distributed
    from dask import distributed
    import pyarrow.dataset as ds
    import pyarrow as pa

    with open(dataset_file, "r") as file:
        all_smiles = [x.strip() for x in file.readlines()]


    output_path = pathlib.Path(output_path)
    counter = 0

    if output_path.exists():
        existing = ds.dataset(output_path)
        existing_smiles = existing.to_table(columns=["smiles"]).to_pydict()["smiles"]
        all_smiles = [x for x in all_smiles if x not in existing_smiles]

        n_batches = int(sorted(output_path.glob("batch-*"))[-1].stem.split("-")[-1])
        counter = n_batches + 1

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
        futures = list(
            batcher(batch_label, ambertools=ambertools, openeye=openeye)
        )
        for i, future in tqdm.tqdm(
            enumerate(distributed.as_completed(futures, raise_errors=False)),
            total=len(futures),
            desc="Saving rows",
        ):
            batch, error = future.result()
            if error:
                print(error)
            batch_path = output_path / f"batch-{counter + i:05d}"
            ds.write_dataset(batch, batch_path, format="parquet")
            future.release()


if __name__ == "__main__":
    generate_all()