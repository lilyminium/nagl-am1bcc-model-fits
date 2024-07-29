import click
import tqdm


@click.command()
@click.option(
    "--input-file",
    help="The path to a file for molecules (sqlite, smiles, sdf)",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
)
@click.option(
    "--output-file",
    help="The path to the output library charge collection file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
    required=True
)
def extract_smiles(
    input_file: str,
    output_file: str
):
    from openeye import oechem
    from openff.toolkit.topology import Molecule

    all_smiles = set()

    counter = 0
    with oechem.oemolistream(input_file) as ifs:
        for oemol in tqdm.tqdm(ifs.GetOEMols(), desc="extracting smiles"):
            offmol = Molecule.from_openeye(oemol, allow_undefined_stereo=True)
            all_smiles.add(offmol.to_smiles())
            counter += 1

    with open(output_file, "w") as f:
        f.write("\n".join(sorted(all_smiles)))

    print(f"Wrote {len(all_smiles)} SMILES from {counter} input OEMols to {output_file}")



if __name__ == "__main__":
    extract_smiles()