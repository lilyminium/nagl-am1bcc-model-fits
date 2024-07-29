
import click
import tqdm

@click.command()
@click.option(
    "--input-file",
    help="Input FASTA file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--output-smiles-file",
    help="Output file to write normal SMILES to",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--output-mapped-smiles-file",
    help="Output file to write mapped SMILES to",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
)
def main(input_file, output_smiles_file, output_mapped_smiles_file):
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions
    from openff.toolkit.topology import Molecule

    with open(input_file, "r") as file:
        all_codes = [x.strip() for x in file.readlines()]


    unmapped_smiles = []
    mapped_smiles = []
    acetyl_rxn = rdChemReactions.ReactionFromSmarts(
        "[C:4](=[O:5])-[C:6](-[C:7])-[N:1](-[H:2])-[H:3]>>"
        "[C:4](=[O:5])-[C:6](-[C:7])-[N:1](-[H:2])-[C:3](=O)C"
    )
    nme_rxn = rdChemReactions.ReactionFromSmarts(
        "[N:1][C:2]([C:3])[C:4](=[O:5])O>>"
        "[N:1]-[C:2](-[C:3])-[C:4](=[O:5])-NC"
    )
    arg_rxn = rdChemReactions.ReactionFromSmarts(
        "[C:8]-[N+0:1]-[C:2](-[N:3](-[H:4])[H:5])=[NH1:6]-[H:7]>>"
        "[C:8]-[N+0:1]-[C:2](-[N:3](-[H:4])[H:5])=[NH2+1:6]"
    )
    for codes in tqdm.tqdm(all_codes):
        try:
            rdmol = Chem.AddHs(Chem.MolFromFASTA(codes))
            rdmol = acetyl_rxn.RunReactants((rdmol,))[0][0]
            rdmol = nme_rxn.RunReactants((rdmol,))[0][0]
            rdmol = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(rdmol)))
            n_R = len([x for x in codes if x == "R"])
            for i in range(n_R):
                rdmol = arg_rxn.RunReactants((rdmol,))[0][0]
                rdmol = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(rdmol)))
            rdmol = Chem.MolFromSmiles(Chem.MolToSmiles(rdmol))
            unmapped_smiles.append(Chem.MolToSmiles(rdmol, isomericSmiles=True))
            offmol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)
            mapped_smiles.append(offmol.to_smiles(mapped=True))
        except Exception as e:
            print(f"Error processing {codes}: {e}")

    with open(output_mapped_smiles_file, "w") as file:
        file.write("\n".join(mapped_smiles))

    with open(output_smiles_file, "w") as file:
        file.write("\n".join(unmapped_smiles))


if __name__ == "__main__":
    main()
