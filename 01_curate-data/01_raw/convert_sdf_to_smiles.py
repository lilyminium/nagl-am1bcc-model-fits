import click
import tqdm


@click.command()
@click.option(
    "--sdf-file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    help="File containing molecules (SDF)",
)
@click.option(
    "--output-file",
    help="The path to the SMILES file (.smi) to save the SMILES strings in.",
)
def convert_to_smiles(
    sdf_file: str,
    output_file: str,
):
    from openff.toolkit import Molecule

    with open(output_file, "w") as f:
        from openeye import oechem

        stream = oechem.oemolistream()
        stream.open(sdf_file)
        is_sdf = stream.GetFormat() == oechem.OEFormat_SDF
        for oemol in tqdm.tqdm(stream.GetOEMols()):
            if is_sdf and hasattr(oemol, "GetConfIter"):
                for conf in oemol.GetConfIter():
                    confmol = conf.GetMCMol()

            else:
                confmol = oemol
        
            try:
                molecule = Molecule.from_openeye(confmol, allow_undefined_stereo=True)
            except BaseException as e:
                print(e)
            else:
                f.write(f"{molecule.to_smiles()}\n")


if __name__ == "__main__":
    convert_to_smiles()