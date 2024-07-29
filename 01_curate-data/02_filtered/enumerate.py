import functools
import multiprocessing
from typing import Set

import click
from click_option_group import optgroup
import tqdm

def _enumerate_tautomers(
    smiles: str,
    enumerate_tautomers: bool,
    max_tautomers: int,
    enumerate_protomers: bool,
    max_protomers: int,
) -> Set[str]:

    from openff.toolkit.utils.exceptions import RadicalsNotSupportedError


    found_forms = {smiles}

    from openff.toolkit.topology import Molecule
    from openff.toolkit.utils import (
        OpenEyeToolkitWrapper,
        RDKitToolkitWrapper,
        ToolkitRegistry,
    )

    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

    if enumerate_tautomers:
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[RDKitToolkitWrapper, OpenEyeToolkitWrapper],
            exception_if_unavailable=False,
        )

        found_forms.update(
            tautomer.to_smiles()
            for tautomer in molecule.enumerate_tautomers(
                max_states=max_tautomers,
                toolkit_registry=toolkit_registry
            )
        )

    if enumerate_protomers:  # pragma: no cover
        from openeye import oechem, oequacpac

        oe_molecule: oechem.OEMol = molecule.to_openeye()
        i = 0
        for oe_protomer in oequacpac.OEGetReasonableProtomers(oe_molecule):
            try:
                offmol = Molecule.from_openeye(oe_protomer, allow_undefined_stereo=True)
            except RadicalsNotSupportedError:
                continue
            found_forms.add(offmol.to_smiles())
            i += 1
            if i >= max_protomers:
                break

    all_forms = []
    for smi in found_forms:
        try:
            all_forms.append((smi, smiles_to_inchi_key(smi)))
        except:
            pass

    return all_forms
    # return found_forms



def smiles_to_inchi_key(smiles: str) -> str:
    from openff.toolkit.topology import Molecule

    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    return molecule.to_inchikey(fixed_hydrogens=True)



@click.command(
    # "enumerate",
    short_help="Enumerate all reasonable tautomers of a molecule set.",
    help="Enumerates all reasonable tautomers (as determine by the OpenEye toolkit) "
    "of a specified set of molecules.",
)
@click.option(
    "--input-file",
    help="The path to the input molecules. This should either be an SDF or a GZipped "
    "SDF file.",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--output-file",
    help="The path to save the enumerated molecules to. This should either be an SDF or "
    "a GZipped SDF file.",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--tautomers/--no-tautomers",
    "enumerate_tautomers",
    help="Whether to enumerate possible tautomers or not.",
    default=False,
    show_default=True,
)
@click.option(
    "--max-tautomers",
    help="The maximum number of tautomers to generate per input molecule.",
    type=int,
    default=16,
    show_default=True,
)
@click.option(
    "--protomers/--no-protomers",
    "enumerate_protomers",
    help="Whether to enumerate the possible protontation states or not. "
    "(requires oequacpac)",
    default=False,
    show_default=True,
)
@click.option(
    "--max-protomers",
    help="The maximum number of protontation states to generate per input molecule.",
    type=int,
    default=16,
    show_default=True,
)
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-processes",
    help="The number of processes to parallelize the enumeration over.",
    type=int,
    default=1,
    show_default=True,
)
def enumerate_cli(
    input_file: str,
    output_file: str,
    enumerate_tautomers: bool,
    max_tautomers: int,
    enumerate_protomers: bool,
    max_protomers: int,
    n_processes: int,
):

    print(
        f" - Enumerating"
        f"{' tautomers' if enumerate_tautomers else ''}"
        f"{'/' if enumerate_protomers and enumerate_tautomers else ''}"
        f"{' protomers' if enumerate_protomers else ''}"
    )

    with open(input_file, "r") as f:
        all_smiles = [x.strip() for x in tqdm.tqdm(f.readlines(), desc="Reading input")]

    generator = functools.partial(
        _enumerate_tautomers,
        enumerate_tautomers=enumerate_tautomers,
        max_tautomers=max_tautomers,
        enumerate_protomers=enumerate_protomers,
        max_protomers=max_protomers
    )

    unique_molecules = set()
    enumerated_smiles = []
    with multiprocessing.Pool(processes=n_processes) as pool:
        for meric_smiles in tqdm.tqdm(
            pool.imap(generator, all_smiles),
            total=len(all_smiles),
            desc="enumerating smiles"
        ):
            for smi, inchi_key in meric_smiles:
                # inchi_key = smiles_to_inchi_key(smi)
                if inchi_key not in unique_molecules:
                    enumerated_smiles.append(smi)
                    unique_molecules.add(inchi_key)

    with open(output_file, "w") as f:
        f.write("\n".join(enumerated_smiles))

    print(f"Wrote {len(enumerated_smiles)} *meric SMILES from {len(all_smiles)} original to {output_file}")

if __name__ == "__main__":
    enumerate_cli()