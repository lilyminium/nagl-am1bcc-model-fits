# Labelling charges

The scripts in here label SMILES with charges.

The scripts in the top-level directory label SMILES with either the "am1bccelf10" charge method (OpenEye) or "am1bcc" (AmberTools). These scripts were used to generate data for the lookup table in the rc3 release candidate. 

The scripts in the legacy/ directory were used to generate data for training/validation. They manually enumerate conformers and generate "elf10" charges by averaging AM1-BCC charges for both.
