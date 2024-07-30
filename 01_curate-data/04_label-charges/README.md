# Labelling charges

The scripts in here label SMILES with charges.

The scripts in the top-level directory label SMILES with either the "am1bccelf10" charge method (OpenEye) or "am1bcc" (AmberTools), via the OpenFF toolkit. These scripts were, e.g., used to generate data for the lookup table in the rc3 release candidate. 

The scripts in the legacy/ directory were used to generate data for training/validation/testing. They manually enumerate conformers and generate "elf10" charges by averaging AM1-BCC charges for both.


A full conda environment is provided in `openff-nagl-test.yaml` for full clarity. For reproduction, it is recommended to generate your own environments with an input file as in the top level.