import itertools

import click
import tqdm

from rdkit import Chem
from rdkit.Chem import rdChemReactions

@click.command()
@click.option(
    "--output-smiles-file",
    help="SMILES file output",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--output-fasta-file",
    help="FASTA file output",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--max-length",
    type=int,
    default=4,
    show_default=True
)
def main(output_smiles_file, output_fasta_file, max_length):
    natural_codes = [
        # # Positively charged
        "R",
        "H",
        "K",
        # Negatively charged
        "D",
        "E",
        # Polar uncharged
        "S",
        "T",
        "N",
        "Q",
        # # Special cases
        "G",
        "P",
        "C",
        # # Hydrophobic
        "A",
        "V",
        "I",
        "L",
        "M",
        # # Hydrophobic
        "F",
        "Y",
        "W",
    ]

    all_codes = []

    # Define the set of 'reactions' that will set the correct protonation state of each
    # residue.
    for i in range(1, max_length + 1):
        for codes in itertools.combinations_with_replacement(natural_codes, r=i):
            if codes not in all_codes:
                all_codes.append(codes)

    acetyl_rxn = rdChemReactions.ReactionFromSmarts(
        "[C:4](=[O:5])-[C:6]-[N:1](-[H:2])-[H:3]>>"
        "[C:4](=[O:5])-[C:6]-[N:1](-[H:2])-[C:3](=O)C"
    )
    nme_rxn = rdChemReactions.ReactionFromSmarts(
        "[N:1][C:2][C:4](=[O:5])O>>"
        "[N:1]-[C:2]-[C:4](=[O:5])-NC"
    )
    arg_rxn = rdChemReactions.ReactionFromSmarts(
        "[C:8]-[N+0:1]-[C:2](-[N:3](-[H:4])[H:5])=[NH1:6]-[H:7]>>"
        "[C:8]-[N+0:1]-[C:2](-[N:3](-[H:4])[H:5])=[NH2+1:6]"
    )

    smiles = []
    for codes in tqdm.tqdm(all_codes):
        fasta = "".join(codes)
        rdmol = Chem.MolFromFASTA(fasta)
        smiles.append(Chem.MolToSmiles(rdmol))

        # run reactions
        if len(codes) > 2:
            if codes[0] != "P":
                rdmol = Chem.AddHs(rdmol)
                rdmol = acetyl_rxn.RunReactants((rdmol,))[0][0]
            if codes[-1] != "P":
                rdmol = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(rdmol)))
                rdmol = nme_rxn.RunReactants((rdmol,))[0][0]
            n_R = len([x for x in codes if x == "R"])
            for i in range(n_R):
                rdmol = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(rdmol)))
                rdmol = arg_rxn.RunReactants((rdmol,))[0][0]
                
            rdmol = Chem.MolFromSmiles(Chem.MolToSmiles(rdmol))
            smiles.append(Chem.MolToSmiles(rdmol))


    with open(output_smiles_file, "w") as file:
        file.write("\n".join(smiles))

    with open(output_fasta_file, "w") as file:
        texts = ["".join(x) for x in all_codes]
        file.write("\n".join(texts))

    print(f"Generated {len(all_codes)} peptides")

if __name__ == "__main__":
    main()
