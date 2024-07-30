# Labelling charges (legacy)

These are the scripts that were actually used to generate data for training and validation. Two versions are provided for different versions of NAGL.

In both:

* 1000 conformers with a 0.05 A RMS cutoff are generated for each molecule. Carboxylic acids are made cis, and the OpenEye toolkit is used.
* A maximum 10 ELF conformers are selected using the OpenEye toolkit.
* Each conformer is labelled using the "am1bcc" partial charge method with the AmberTools and/or OpenEye toolkits. 


``generate_conformers_and_label.py`` uses an older (v0.2.x) version of NAGL that still stores MoleculeRecords in a MoleculeStore database.

``generate_conformers_and_label-pyarrow_both.py`` uses a newer (v0.3.x+) version that stores data in a PyArrow table.

Environments are provided for both, where `openff-nagl-refactor` is the environment used for the older (v0.2.x) version of NAGL, and `openff-nagl-pyarrow` is for the newer (v0.3.x+) version of NAGL. While there are minor differences between the two environments, the main relevant software packages are the same:

* OpenEye == 2022.2.1
* RDKit == 2022.03.05
* AmberTools == 22.0


These are provided for clarity only; if reproducing this work it is highly recommended to use an input file to generate a new environment appropriate to your machine.