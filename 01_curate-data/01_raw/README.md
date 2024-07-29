# Obtaining raw data

This directory contains helper scripts for downloading, cleaning, or otherwise initially processing raw data.

## Scripts

* generate-short-peptide-chains.py: generates peptide chains from amino acid codes ranging up to 4 
* generate-peptide-chains.py: generates peptide chains from amino acid codes ranging 5-30 residues in length
* cap-peptides.py: caps chains with an acetyl and n-methyl cap
* generate-small-small-molecules.py: generates small molecules with up to 3 heavy atoms


## input

The files in the ``input/`` directory contain raw datasets.
Not all these are included for licensing reasons.

* ChEMBL_eps_78.sdf.gz, from Bleiziffer et al 2018
* ZINC_eps_78.sdf.gz,
* enamine-10240.sdf.gz, from enamine
* enamine-50240.sdf.gz, from enamine
* fda.smi
* pdb.smi
* all-chembl.smi


## output

This directory contains files with a smiles pattern on each line.

