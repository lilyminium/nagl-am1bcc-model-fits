"""
This script generates peptide chains 5 residues long, with and without post-translational modifications.
The following modifications are attempted:

* lipidation
* glycosylation
* methylation
* phosphorylation
* acetylation
* nitrosylation
* ubiquitination

"""

import itertools

import tqdm
import pandas as pd
from openff.toolkit import Molecule

from rdkit import Chem
from rdkit.Chem import rdChemReactions

PATTERNS = dict(
    # lipidation
    n_myristoylation = "[N:1](-[H:2])-[H:3]>>[N:1](-[H:2])C(=O)CCCCCCCCCCCCC",
    s_palmitoylation = "[S:1]-[H:2]>>[S:1]-C(=O)CCCCCCCCCCCCCCC",

    # phosphorylation
    phosphorylation_1 = "[O:1]-[H:2]>>[O:1]-P(-[O-])(=O)-[O-]",
    phosphorylation_2 = "[O:1]-[H:2]>>[O:1]-P(-[O-])(=O)-O-P(-[O-])(=O)-[O-]",

    # nitrosylation
    nitrosylation = "[S:1]-[H:2]>>[S:1]-N(=O)",

    # acetylation
    acetylation = "[N+1:1](-[H:2])-[H:3]>>[N:1](-[H:2])-C(=O)C",

    # glycosylation
    n_glycosylation_glcnac = "[C:1](=[O:2])-[N:3](-[H:4])-[H:5]>>[C:1](=[O:2])-[N:3](-[H:4])-C1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H](NC(=O)C)1",
    o_glycosylation_galnac = "[O:1]-[H:2]>>[O:1]-C1O[C@@H](CO)[C@@H](O)[C@@H](O)[C@H](NC(=O)C)1",

    # methylation
    n3_methylation_1 = "[N+1:1]-[H:2]>>[N+1:1]-C",
    n3_methylation_2 = "[N+1:1](-[H:2])-[H:3]>>[N+1:1](C)C",
    n3_methylation_3 = "[N+1:1](-[H:2])(-[H:3])-[H:4]>>[N+1:1](C)(C)C",

    arg_methylation_1 = "[C:1](=[N+1:2](-[H:3])-[H:4])-[N:5](-[H:6])-[H:7]>>[C:1](=[N+1:2](-[H:3])-[H:4])-[N:5](-[H:6])C",
    arg_methylation_2a = "[C:1](=[N+1:2](-[H:3])-[H:4])-[N:5](-[H:6])-[H:7]>>[C:1](=[N+1:2](-[H:3])-[H:4])-[N:5](C)C",
    arg_methylation_2b = "[C:1](=[N+1:2](-[H:3])-[H:4])-[N:5](-[H:6])-[H:7]>>[C:1](=[N+1:2](-[H:3])C)-[N:5](-[H:6])C",
)

PTM_REACTIONS = {
    k: rdChemReactions.ReactionFromSmarts(v)
    for k, v in PATTERNS.items()
}


def rdmol_to_smiles(rdkit_molecule):
    offmol = Molecule.from_smiles(
        Chem.MolToSmiles(rdkit_molecule),
        allow_undefined_stereo=True
    )
    return offmol.to_smiles()

def main():
    cap_ace = "[N:10]-[C:8](=[O:9])[C:7]-[N:11]-[C:1](=[O:2])-[C:3]-[N:4](-[H:5])-[H:6]>>[N:10]-[C:8](=[O:9])[C:7]-[N:11]-[C:1](=[O:2])-[C:3]-[N:4](-[H:5])C(=O)C"
    cap_nme = "[N:10]-[C:9]-[C:7](=[O:8])-[N:6]-[C:5]-[C:1](=[O:2])-[O:3]-[H:4]>>[N:10]-[C:9]-[C:7](=[O:8])-[N:6]-[C:5]-[C:1](=[O:2])-NC"
    cap_ace_rxn = rdChemReactions.ReactionFromSmarts(cap_ace)
    cap_nme_rxn = rdChemReactions.ReactionFromSmarts(cap_nme)

    valid_aas = "ACDEFGHIKLMNPQRSTVWY"

    smiles_data = {
        "smiles": [],
        "modification": [],
        "fasta": [],
    }
    for triplet in tqdm.tqdm(itertools.combinations_with_replacement(valid_aas, r=3)):
        peptide_aas = f"GV{''.join(triplet)}VG"
        peptide = Chem.AddHs(Chem.MolFromFASTA(peptide_aas))
        
        ace_peptide = cap_ace_rxn.RunReactants((peptide,))
        capped_peptide = cap_nme_rxn.RunReactants(ace_peptide[0])[0][0]
        smiles_data["smiles"].append(rdmol_to_smiles(capped_peptide))
        smiles_data["fasta"].append(peptide_aas)
        smiles_data["modification"].append("plain")
        
        for reaction_name, reaction in PTM_REACTIONS.items():
            products = reaction.RunReactants((capped_peptide,))
            if not products:
                continue
            product = products[0][0]
            smiles_data["smiles"].append(rdmol_to_smiles(product))
            smiles_data["fasta"].append(peptide_aas)
            smiles_data["modification"].append(reaction_name)
    
    df = pd.DataFrame(smiles_data)
    df.to_csv("input/peptides_post-translational_modifications.csv")

    with open("output/peptides_post-translational_modifications.smi", "w") as f:
        f.write("\n".join(df.smiles.values))
