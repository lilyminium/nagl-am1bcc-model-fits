# Labelling charges (legacy)

These are the scripts that were actually used to generate data for training and validation. Two versions are provided for different versions of NAGL.

In both:

* 1000 conformers with a 0.05 A RMS cutoff are generated for each molecule. Carboxylic acids are made cis, and the OpenEye toolkit is used.
* A maximum 10 ELF conformers are selected using the OpenEye toolkit.
* Each conformer is labelled using the "am1bcc" partial charge method with the AmberTools and/or OpenEye toolkits. 


``generate_conformers_and_label.py`` uses an older (v0.2.x) version of NAGL that still stores MoleculeRecords in a MoleculeStore database.

``generate_conformers_and_label-pyarrow_both.py`` uses a newer (v0.3.x+) version that stores data in a PyArrow table.

Environments are provided for both.