import collections
import typing
import multiprocessing

import click
from click_option_group import optgroup

import tqdm
from openff.toolkit import Molecule, ForceField
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import pandas as pd

from openff.toolkit.utils import RDKitToolkitWrapper

RDKIT_WRAPPER = RDKitToolkitWrapper()

def load_valid_smiles(smi: str):
    try:
        mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
    except Exception as e:
        return f"Failed to parse {smi}: {e}"
    else:
        return mol.to_rdkit()
    

def _flip_r_and_s(stereo: str):
    if stereo == "R":
        return "S"
    elif stereo == "S":
        return "R"
    else:
        raise ValueError(f"Unknown stereochemistry {stereo}")


def are_enantiomers(mol1, mol2):
    # assert they map onto each other and that stereoisomer enumeration
    # keep same atom order
    assert mol1.is_isomorphic_with(
        mol2,
        atom_stereochemistry_matching=False,
        bond_stereochemistry_matching=False,
        strip_pyrimidal_n_atom_stereo=True
    )

    # check if they are enantiomers
    # i.e. all stereocenters are inverted
    rdmol1 = mol1.to_rdkit()
    rdmol2 = mol2.to_rdkit()

    stereocenters1 = Chem.FindMolChiralCenters(rdmol1, includeUnassigned=False)
    stereocenters2 = Chem.FindMolChiralCenters(rdmol2, includeUnassigned=False)

    flip_1_stereocenters = [
        (_, _flip_r_and_s(stereo))
        for _, stereo in stereocenters1
    ]
    return flip_1_stereocenters == stereocenters2



def handle_single_molecule(smi: str):
    # if multiple fragments, keep largest only
    smi = max(smi.split("."), key=len)
    # load molecule to ensure it's valid
    try:
        mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
    except Exception as e:
        return None, f"Failed to parse {smi}: {e}"
    
    if mol.n_atoms > 90:
        return None, f"Too many atoms ({mol.n_atoms})"

    # check force field can assign parameters
    ff = ForceField("openff-2.1.0.offxml")
    try:
        ff.create_interchange(mol.to_topology())
    except Exception as e:
        return [], f"Could not assign parameters to {smi}: {e}"
    
    # count stereocenters
    rdmol = mol.to_rdkit()
    Chem.AssignStereochemistry(rdmol, force=True, cleanIt=True)
    n_stereocenters = len(Chem.FindMolChiralCenters(rdmol, includeUnassigned=True))
    # n_stereocenters = rdMolDescriptors.CalcNumAtomStereoCenters(rdmol)
    if n_stereocenters <= 1:
        return [], f"Molecule has <=1 stereocenters, {smi}"
    
    # enumerate stereoisomers
    # doesn't include self
    mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)

    try:
        isomers_mol = mol.enumerate_stereoisomers(
            undefined_only=False,
            max_isomers=100,
            rationalise=True,
            # toolkit_registry=RDKIT_WRAPPER
        )
    except Exception as e:
        # conversion back to mol can still fail on stereochemistry
        return [], f"Failed to enumerate stereoisomers, {e}, {smi}"
    
    isomers = []
    # if all stereocenters defined in input mol, include
    try:
        isomer = Molecule.from_mapped_smiles(mol.to_smiles(mapped=True))
    except Exception as e:
        pass
    else:
        isomers.append(isomer)

    # check if all stereocenters have been defined to Toolkit satisfaction
    for isomer_mol in isomers_mol:
        try:
            isomer = Molecule.from_mapped_smiles(isomer_mol.to_smiles(mapped=True))
        except Exception as e:
            continue
        else:
            isomers.append(isomer)

    # check enantiomers
    if not len(isomers):
        return [], f"No stereoisomers found, {smi}"
    
    diastereomers = [isomers.pop(0)]
    for isomer in isomers:
        if not are_enantiomers(diastereomers[0], isomer):
            diastereomers.append(isomer)
    if len(diastereomers) <= 1:
        return [], f"No diastereomers found, {smi}"

    diastereomer_smiles = [
        mol.to_smiles(mapped=True)
        for mol in diastereomers
    ]
    return (mol.n_atoms, n_stereocenters, diastereomer_smiles), None
    

def batch_handle_molecules(smiles: list[str]):
    entries = []
    errors = []
    for smi in tqdm.tqdm(smiles):
        entry, error = handle_single_molecule(smi)
        if error:
            errors.append(error)
        else:
            entries.append(entry)
    return entries, errors


@click.command()
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
def main(
    nprocs: int = 16,
    verbose: bool = True,
    worker_type: typing.Literal["lsf", "local"] = "local",
    queue: str = "cpuqueue",
    conda_environment: str = "openff-nagl",
    memory: int = 4,  # GB
    walltime: int = 32,  # hours
    batch_size: int = 300,
    n_workers: int = -1,
):
    from openff.nagl.utils._parallelization import batch_distributed, as_batch_function
    from dask import distributed
    from dask.distributed import get_client

    # load chembl; pass through Toolkit first
    input_file = "input/all-chembl.smi"

    with open(input_file, "r") as f:
        all_smiles = [line.strip() for line in f]
    n_smiles = len(all_smiles)

    all_smiles = sorted(all_smiles, key=len)[:100000]

    all_outputs = []
    all_errors = []
    threshold = 5000
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
        futures = list(batcher(
            batch_handle_molecules,
        ))

        for future in tqdm.tqdm(
            distributed.as_completed(futures, raise_errors=False),
            total=len(futures)
        ):
            output, error = future.result()
            all_outputs.extend(output)
            all_errors.extend(error)

            if len(all_outputs) >= threshold:
                break
                client = get_client()
                client.cancel(futures)
    
    by_n_atoms = collections.defaultdict(list)
    for diastereomers in all_outputs:
        n_atoms, n_stereo, smiles = diastereomers
        # check if this molecule has already been seen
        for existing in by_n_atoms[n_atoms]:
            mol1 = Molecule.from_smiles(existing[-1][0])
            mol2 = Molecule.from_smiles(smiles[0])
            if mol1.is_isomorphic_with(
                mol2,
                atom_stereochemistry_matching=False,
                bond_stereochemistry_matching=False,
                strip_pyrimidal_n_atom_stereo=True
            ):
                break
        else:
            by_n_atoms[n_atoms].append((n_stereo, smiles))

    # outputs = [
    #     handle_single_molecule(smi)
    #     for smi in tqdm.tqdm(all_smiles[:110])
    # ]

    # # with multiprocessing.Pool(nprocs) as pool:
    # #     outputs = list(
    # #         tqdm.tqdm(
    # #             pool.imap(handle_single_molecule, all_smiles[:10000]),
    # #             total=n_smiles
    # #     ))

    # by_n_atoms = collections.defaultdict(list)
    # for diastereomers, error in outputs:
    #     if error:
    #         continue
    #     n_atoms, n_stereo, smiles = diastereomers
    #     # check if this molecule has already been seen
    #     for existing in by_n_atoms[n_atoms]:
    #         mol1 = Molecule.from_smiles(existing[-1][0])
    #         mol2 = Molecule.from_smiles(smiles[0])
    #         if mol1.is_isomorphic_with(
    #             mol2,
    #             atom_stereochemistry_matching=False,
    #             bond_stereochemistry_matching=False,
    #             strip_pyrimidal_n_atom_stereo=True
    #         ):
    #             break
    #     else:
    #         by_n_atoms[n_atoms].append((n_stereo, smiles))

    if verbose:
        for error in all_errors:
            print(error)

    all_entries = []
    mol_id = 1
    for n_atoms in sorted(by_n_atoms.keys()):
        for n_stereo, diastereomers in by_n_atoms[n_atoms]:
            for j, isomer in enumerate(diastereomers, 1):
                entry = {
                    "smiles": isomer,
                    "n_atoms": n_atoms,
                    "n_stereocenters": n_stereo,
                    "mol_id": mol_id,
                    "diastereomer_id": j
                }
                all_entries.append(entry)
            mol_id += 1

    df = pd.DataFrame(all_entries)
    df.to_csv("input/chembl-diastereomers.csv", index=False)

    with open("output/chembl-diastereomers.smi", "w") as f:
        f.write("\n".join(df.smiles.values))

    

if __name__ == "__main__":
    main()
