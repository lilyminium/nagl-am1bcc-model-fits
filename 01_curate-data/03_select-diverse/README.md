# Selecting diverse molecules

The code here is similar to the selection process of diverse molecules from [Bleiziffer et al. 2018](https://doi.org/10.1021/acs.jcim.7b00663).

The code:

* assigns an AtomPairFingerprint to each atom with length 2, 2048 bits, and 4 bits per entry
* firstly include all molecules that contain a "rare" atom environment, where a rare atom environment is one that is present in <= N_MIN_MOLECULES
* for each of the remaining "non-rare" atom environments, roughly ordered by element ("S", "F", "Cl", "Br", "I", "P", "O", "N", "C", "H"):
    * the molecule with the highest number of not-yet-represented atom environments is added to the the final dataset, iteratively, until each atom environment has at least N_MIN_MOLECULES molecules


The code in the `legacy/` directory uses an earlier version (v0.2.x) of NAGL to do this selection. The implementation is more convoluted but follows the same algorithm.
