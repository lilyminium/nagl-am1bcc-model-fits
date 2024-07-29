import click
import collections
import itertools
import tqdm

def main():
    """
    This script enumerates the entire possible space of molecules <= 3 heavy atoms,
    with elements supported by NAGL. Charge states ranging -4 to +4 are also included.
    It generates all possible combinations of atoms and bonds between them.

    Radicals and invalid molecules are filtered out by using the OpenFF Toolkit
    to parse the resulting smiles.
    """
    from openff.toolkit import Molecule

    elements = ["H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"]

    atoms = list(elements)
    for charge_states in ["-4", "-3", "-2", "-1", "+1", "+2", "+3", "+4"]:
        for atom in elements:
            atoms.append(f"[{atom}{charge_states}]")
    print(f"Generated {len(atoms)} atoms")
    
    small_small_molecules = set(atoms)

    bond_types = ["-", "=", "#", ":"]
    # two elements should be unique
    for bond_type in tqdm.tqdm(bond_types):
        for el1, el2 in itertools.product(atoms, repeat=2):
            small_small_molecules.add(f"{el1}{bond_type}{el2}")
    
    # three elements may not be unique so sort by key
    molecules_by_key = collections.defaultdict(set)
    for el1, el2, el3 in tqdm.tqdm(itertools.product(atoms, repeat=3)):
        for bt1, bt2 in itertools.product(bond_types, repeat=2):
            # key is all bonds, sorted
            bond1_elements = sorted([el1, el2])
            bond1 = (bond1_elements[0], bt1, bond1_elements[1])
            bond2_elements = sorted([el2, el3])
            bond2 = (bond2_elements[0], bt2, bond2_elements[1])

            key = tuple(sorted([bond1, bond2]))
            forward = f"{el1}{bt1}{el2}{bt2}{el3}"
            molecules_by_key[key].add(forward)
            reverse = f"{el3}{bt2}{el2}{bt1}{el1}"
            if reverse != forward and reverse in molecules_by_key[key]:
                molecules_by_key[key].remove(reverse)
            for bt3 in bond_types:
                bond3_elements = sorted([el1, el3])
                bond3 = (bond3_elements[0], bt3, bond3_elements[1])
                key = tuple(sorted([bond1, bond2, bond3]))
                forward = f"{el1}1{bt1}{el2}{bt2}{el3}{bt3}1"
                reverse = f"{el3}1{bt3}{el2}{bt2}{el1}{bt1}1"
                molecules_by_key[key].add(forward)
                if reverse != forward and reverse in molecules_by_key[key]:
                    molecules_by_key[key].remove(reverse)

    # keep only unique molecules from molecules_by_key
    for candidates in tqdm.tqdm(
        molecules_by_key.values(),
        desc="Checking unique molecules",
        total=len(molecules_by_key)
    ):
        if len(candidates) == 1:
            small_small_molecules.update(candidates)
        else:
            candidates = list(candidates)
            valid_molecules = []
            for candidate in candidates:
                try:
                    mol = Molecule.from_smiles(candidate, allow_undefined_stereo=True)
                    valid_molecules.append((candidate, mol))
                except Exception as e:
                    pass
            
            if valid_molecules:
                unique = [valid_molecules.pop(0)]

                for candidate, mol in valid_molecules:
                    if not any(
                        [
                            mol.is_isomorphic_with(other)
                            for _, other in unique
                        ]
                    ):
                        unique.append((candidate, mol))
                    del mol
            
                small_small_molecules.update([smi for smi, _ in unique])
            



    mapped_smiles = []
    valid_smiles = []
    errors = []
    small_small_molecules = sorted(small_small_molecules, key=len, reverse=True)
    for smi in tqdm.tqdm(small_small_molecules, desc="Checking molecules"):
        try:
            mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)            
        except Exception as e:
            errors.append((smi, e))
            continue
        
        mapped_smiles.append(mol.to_smiles(mapped=True))
        valid_smiles.append(smi)

    with open("input/small-small-molecules.smi", "w") as f:
        f.write("\n".join(valid_smiles))
    print(f"Saved {len(valid_smiles)} small-small molecules to input/small-small-molecules.smi")
    with open("output/mapped-small-small-molecules.smi", "w") as f:
        f.write("\n".join(mapped_smiles))
    print(f"Saved {len(mapped_smiles)} small-small molecules to output/mapped-small-small-molecules.smi")

    with open("small-small-errors.txt", "w") as f:
        for smi, e in errors:
            f.write(f"{smi} **** {e}\n")

if __name__ == "__main__":
    main()
