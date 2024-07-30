import click

@click.command()
@click.option(
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
    required=True,
)
@click.option(
    "--output",
    "output_path",
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    required=True,
)
def export(
    input_path: str,
    output_path: str,
):
    """
    This script adds a version number to the GNNModel.
    It also adds a ChemicalDomain to describe the allowed elements and forbidden patterns.
    """
    import torch
    from openff.nagl.domains import ChemicalDomain

    data = torch.load(input_path)
    data["hyperparameters"]["config"]["version"] = "0.1"
    domain = ChemicalDomain(
        allowed_elements=[1, 6, 7, 8, 9, 15, 16, 17, 35, 53],
        forbidden_patterns=[
            '[#15:1]#[*:2]',
            '[#15:1]-[#53:2]',
            '[#15:1]:[*:2]',
            '[#15:1]=[!#6&!#7&!#8&!#16:2]',
            '[#16:1]#[*:2]',
            '[#16:1]-[#35,#53:2]',
            '[#16:1]=[!#15&!#6&!#7&!#8:2]',
            '[#17:1]#,:,=[*:2]',
            '[#17:1]-[!#15&!#16!#6&!#7&!#8:2]',
            '[#1:1]#,:,=[*:2]',
            '[#1:1]-[#1:2]',
            '[#35:1]#,:,=[*:2]',
            '[#35:1]-[!#15&!#6&!#7&!#8:2]',
            '[#53:1]#,:,=[*:2]',
            '[#53:1]-[!#53&!#6&!#7&!#8:2]',
            '[#7:1]#[!#7&!#6:2]',
            '[#8:1]#[*:2]',
            '[#9:1]#,:,=[*:2]',
            '[#9:1]-[!#15&!#16!#6&!#7&!#8:2]'
        ]
    )

    data["hyperparameters"]["chemical_domain"] = domain.dict()
    torch.save(data, output_path)


if  __name__ == "__main__":
    export()