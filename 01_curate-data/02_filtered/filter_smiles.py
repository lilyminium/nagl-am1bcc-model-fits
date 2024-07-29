
import functools
import logging
from typing import TYPE_CHECKING, Iterable, List, Tuple, Union
import multiprocessing


import click
from click_option_group import optgroup

from openff.units.elements import MASSES, SYMBOLS

if TYPE_CHECKING:
    from openff.toolkit.topology.molecule import Molecule, unit

logger = logging.getLogger(__name__)

INV_SYMBOLS = {v: k for k, v in SYMBOLS.items()}

def get_atomic_number(el: Union[int, str]) -> int:
    if isinstance(el, int):
        return el
    return INV_SYMBOLS[el]


def apply_filter(
    molecule: "Molecule",
    allowed_elements: Tuple[int],
    min_mass: "unit.Quantity",
    max_mass: "unit.Quantity",
    n_rotatable_bonds: int,
) -> bool:
    mass = sum(MASSES[atom.atomic_number] for atom in molecule.atoms)

    return (
        all(atom.atomic_number in allowed_elements for atom in molecule.atoms)
        and mass > min_mass
        and mass < max_mass
        and len(molecule.find_rotatable_bonds()) <= n_rotatable_bonds
    )


def split_and_apply_filter(
    molecule: "Molecule",
    allowed_elements: Tuple[int],
    min_mass: "unit.Quantity",
    max_mass: "unit.Quantity",
    n_rotatable_bonds: int,
    only_retain_largest: bool = True,
    as_smiles: bool = False,
    mapped_smiles: bool = False,
):
    from openff.toolkit.topology import Molecule

    if isinstance(molecule, str):
        if only_retain_largest:
            split = molecule.split(".")
            molecule = max(split, key=len)
        molecule = Molecule.from_smiles(molecule, allow_undefined_stereo=True)
    try:
        if only_retain_largest:
            split_smiles = molecule.to_smiles().split(".")
            if len(split_smiles) > 1:
                largest = max(split_smiles, key=len)
                molecule = Molecule.from_smiles(largest, allow_undefined_stereo=True)
                logger.debug(f"Keeping '{largest}' from '{split_smiles}'")
        valid = apply_filter(
            molecule,
            allowed_elements=allowed_elements,
            min_mass=min_mass,
            max_mass=max_mass,
            n_rotatable_bonds=n_rotatable_bonds,
        )
        if valid:
            if as_smiles:
                return molecule.to_smiles(mapped=mapped_smiles)
            else:
                return molecule
    except Exception as e:
        logger.warning(f"Failed to process molecule {molecule}, {e}")


def filter_molecules(
    molecules: Iterable["Molecule"],
    only_retain_largest: bool = True,
    allowed_elements: Tuple[Union[str, int], ...] = (
        "H",
        "C",
        "N",
        "O",
        "F",
        "P",
        "S",
        "Cl",
        "Br",
        "I",
    ),
    min_mass: "unit.Quantity" = 250,
    max_mass: "unit.Quantity" = 350,
    n_rotatable_bonds: int = 7,
    n_processes: int = 1,
    as_smiles: bool = False,
) -> Iterable["Molecule"]:
    import tqdm
    from openff.toolkit.topology.molecule import unit

    from openff.nagl.utils.openff import capture_toolkit_warnings

    allowed_elements = [get_atomic_number(x) for x in allowed_elements]

    if not isinstance(min_mass, unit.Quantity):
        min_mass = min_mass * unit.amu
    if not isinstance(max_mass, unit.Quantity):
        max_mass = max_mass * unit.amu

    with capture_toolkit_warnings():
        filterer = functools.partial(
            split_and_apply_filter,
            only_retain_largest=only_retain_largest,
            allowed_elements=allowed_elements,
            n_rotatable_bonds=n_rotatable_bonds,
            min_mass=min_mass,
            max_mass=max_mass,
            as_smiles=as_smiles,
        )
        with multiprocessing.Pool(processes=n_processes) as pool:
            for molecule in tqdm.tqdm(
                pool.imap(filterer, molecules), desc="filtering molecules"
            ):
                if molecule is not None:
                    yield molecule



@click.command(
    short_help="Filter undesirable chemistries and counter-ions.",
    help="Filters a set of molecules based on the criteria specified by:\n\n"
    "    [1] Bleiziffer, Patrick, Kay Schaller, and Sereina Riniker. 'Machine learning "
    "of partial charges derived from high-quality quantum-mechanical calculations.' "
    "JCIM 58.3 (2018): 579-590.\n\nIn particular molecules are only retained if they "
    "have a weight between 250 and 350 g/mol, have less than seven rotatable bonds and "
    "are composed of only H, C, N, O, F, P, S, Cl, Br, and I.\n\nThis script will also "
    "optionally remove any counter-ions by retaining only the largest molecule if "
    "multiple components are present.",
)
@click.option(
    "--input-file",
    help="The path to the input molecules: SDF, zipped SDF, or smiles file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--output-file",
    help="The path to save the filtered molecules to. This should be an SDF or smiles file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--only-retain-largest",
    is_flag=True,
    help="If specified counter ions (and molecules) will be removed.",
    default=False,
    show_default=True,
)
@click.option(
    "--min-mass",
    type=float,
    help="Minimum mass (g/mol)",
    default=250,
    show_default=True,
)
@click.option(
    "--max-mass",
    type=float,
    help="Maximum mass (g/mol)",
    default=250,
    show_default=True,
)
@click.option(
    "--n-rotatable-bonds",
    type=int,
    help="Number of rotatable bonds",
    default=7,
    show_default=True,
)
@click.option(
    "--element",
    "allowed_elements",
    type=str,
    multiple=True,
    help="Allowed elements",
    default=("H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"),
    show_default=True,
)
@click.option(
    "--exclude-smiles",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    multiple=True,
)
@optgroup.group("Parallelization configuration")
@optgroup.option(
    "--n-processes",
    help="The number of processes to parallelize the filtering over.",
    type=int,
    default=1,
    show_default=True,
)
def filter_cli(
    input_file: str,
    output_file: str,
    allowed_elements: Tuple[str, ...],
    only_retain_largest: bool,
    exclude_smiles: Tuple[str, ...] = tuple(),
    min_mass: float = 250,
    max_mass: float = 350,
    n_rotatable_bonds: int = 7,
    n_processes: int = 1,
):
    from collections import defaultdict
    import tqdm
    from openff.toolkit.topology import Molecule
    from openff.nagl.toolkits.openff import capture_toolkit_warnings

    print("loading excluded smiles")
    all_excluded_smiles = []
    for filename in exclude_smiles:
        with open(filename, "r") as f:
            read = [x.strip() for x in f.readlines()]
            all_excluded_smiles.extend([x for x in read if x])
    print(f"loaded {len(all_excluded_smiles)} excluded smiles")

    with capture_toolkit_warnings():
        offmols = [
            Molecule.from_smiles(x, allow_undefined_stereo=True)
            for x in tqdm.tqdm(all_excluded_smiles, desc="Loading excluded SMILES")
        ]
    all_offmols = defaultdict(list)
    for offmol in offmols:
        all_offmols[offmol.n_atoms].append(offmol)

    print("starting to filter")

    with open(input_file, "r") as f:
        smiles = [x.strip() for x in tqdm.tqdm(f.readlines(), desc="Reading input")]

    not_radicals = []
    with capture_toolkit_warnings():
        for smi in tqdm.tqdm(smiles, desc="filtering out radicals"):
            try:
                offmol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
                offmol = Molecule.from_openeye(offmol.to_openeye(), allow_undefined_stereo=True)
                offmol = Molecule.from_rdkit(offmol.to_rdkit(), allow_undefined_stereo=True)
            except:
                pass
            else:
                not_radicals.append(offmol.to_smiles())

    counter = 0
    with open(output_file, "w") as f:
        for smi in filter_molecules(
            not_radicals,
            only_retain_largest=only_retain_largest,
            allowed_elements=allowed_elements,
            min_mass=min_mass,
            max_mass=max_mass,
            n_rotatable_bonds=n_rotatable_bonds,
            n_processes=n_processes,
            as_smiles=True
        ):
            offmol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
            candidates = all_offmols[offmol.n_atoms]
            if any([x.is_isomorphic_with(offmol) for x in candidates]):
                continue
            counter += 1
            f.write(f"{offmol.to_smiles()}\n")
    
    print(f"Filtered {len(smiles)} SMILES to {counter} SMILES")


if __name__ == "__main__":
    filter_cli()