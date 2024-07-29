import functools
from typing import Tuple, List, Literal

import click
import tqdm
from click_option_group import optgroup


def generate_conformers_and_label(
    smiles: str,
    partial_charge_methods: Tuple[str, ...] = tuple(),
    bond_order_methods: Tuple[str, ...] = tuple()
):
    from openff.units import unit
    from openff.toolkit import Molecule
    from openff.toolkit.utils import OpenEyeToolkitWrapper
    from openff.nagl.storage.record import MoleculeRecord
    from openff.toolkit.utils.exceptions import UndefinedStereochemistryError
    from openff.nagl.toolkits.openff import capture_toolkit_warnings

    with capture_toolkit_warnings():
        offmol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

        # try:
        #     offmol = Molecule.from_mapped_smiles(smiles)
        # except UndefinedStereochemistryError:
        #     offmol = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)
        #     stereo = offmol.enumerate_stereoisomers(
        #         undefined_only=True,
        #         max_isomers=1,
        #         rationalise=True
        #     )
        #     if not stereo:
        #         stereo = [offmol]
        #     offmol = Molecule.from_mapped_smiles(stereo[0].to_smiles(mapped=True), allow_undefined_stereo=True)

    offmol.generate_conformers(
        n_conformers=1000,
        rms_cutoff=0.05 * unit.angstrom,
        make_carboxylic_acids_cis=True,
        toolkit_registry=OpenEyeToolkitWrapper(),
    )
    offmol.apply_elf_conformer_selection(
        limit=10,
        toolkit_registry=OpenEyeToolkitWrapper(),
    )

    record = MoleculeRecord.from_openff(
        offmol,
        partial_charge_methods=partial_charge_methods,
        bond_order_methods=bond_order_methods,
        generate_conformers=False
    )
    return record


@click.command()
@click.option(
    "--smiles-file",
    "smiles_files",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="File containing SMILES strings, one on each line.",
    multiple=True
)
@click.option(
    "--output-file",
    help="The path to the SQLite database (.sqlite) to save the labelled molecules in.",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--partial-charge-method",
    "partial_charge_methods",
    help="The partial charge methods to compute",
    multiple=True,
    default=tuple(),
    show_default=True,
)
@click.option(
    "--deduplicate",
    help="Whether to remove duplicate molecules from the input SMILES files.",
    is_flag=True,
    default=False,
    show_default=True,
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
    smiles_files: List[str],
    output_file: str,
    partial_charge_methods: Tuple[str, ...],
    deduplicate: bool = False,
    worker_type: Literal["lsf", "local"] = "local",
    queue: str = "cpuqueue",
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from openff.nagl.storage._store import MoleculeStore
    from openff.nagl._app.distributed import Manager
    from dask import distributed
    from openff.nagl._cli.utils import as_batch_function_with_captured_errors
    from openff.toolkit import Molecule

    all_smiles = set()
    for smiles_file in tqdm.tqdm(smiles_files, desc="loading molecules"):
        with open(smiles_file, "r") as f:
            smiles = [x.strip() for x in f.readlines()]
            all_smiles |= set(smiles)

    print(f"Loaded {len(all_smiles)} unique SMILES.")

    if deduplicate:
        canonical_smiles = set()
        for smiles in tqdm.tqdm(all_smiles, desc="canonicalizing SMILES"):
            canonical_smiles.add(
                Molecule.from_smiles(
                    smiles,
                    allow_undefined_stereo=True,
                ).to_smiles()
            )

        all_smiles = canonical_smiles
        print(f"Canonicalized to {len(all_smiles)} unique SMILES.")

    manager = Manager(
        batch_size=batch_size,
        worker_type=worker_type,
        queue=queue,
        conda_environment=conda_environment,
        memory=memory,
        walltime=walltime,
        n_workers=n_workers,
    )
    manager.set_entries(sorted(all_smiles, key=len, reverse=True))

    single_func = functools.partial(
        generate_conformers_and_label,
        partial_charge_methods=partial_charge_methods,
    )

    batch_func = as_batch_function_with_captured_errors(
        single_func,
        desc="Labelling records"
    )

    store = MoleculeStore(output_file)

    n_results = 0

    error_file = output_file.replace(".sqlite", "_errors.dat")
    with open(error_file, "w") as f:
        with manager:
            futures = manager.submit_to_client(batch_func)

            for future in tqdm.tqdm(
                distributed.as_completed(futures),
                total=manager.n_batches,
                desc="Saving records",
            ):
                batch_results = []
                for result, error in future.result():
                    if error:
                        f.write(error)
                        f.write("\n")
                    else:
                        try:
                            store.store(result)
                            n_results += 1
                        except BaseException as e:
                            f.write(error)
                            f.write("\n")

                        # batch_results.append(result)
                # n_results += len(batch_results)
                # store.store(batch_results)

                future.release()

    print(f"Saved {n_results} records to {output_file}.")


if __name__ == "__main__":
    generate_all()