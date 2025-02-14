# molden2qmc
molden format file to Casino gwfn.data file format converter.

This python script is used to create qwfn.data file which used for input in the [CASINO](https://vallico.net/casinoqmc/what-is-casino/) program.
The gwfn.data file contains molecular geometry and molecular orbitals expanded in a Gaussian basis set.

In general any file that confirms [MOLDEN](https://www.theochem.ru.nl/molden/molden_format.html) specification can be used,

but different quantum chemistry code not always fully comply with the required specifications,
in such cases, the script corrects these discrepancies.

The script supports the Molden file generated by the following programs:

0. — [TURBOMOLE](http://www.turbomole.com/)
1. — [PSI4](http://www.psicode.org/)
2. — [CFOUR](http://www.cfour.de/)
3. — [ORCA](https://orcaforum.cec.mpg.de/)
4. — [DALTON](http://daltonprogram.org/)
5. — [MOLPRO](https://www.molpro.net/)
6. — [NWCHEM](https://nwchemgit.github.io/)
7. — [QCHEM](http://www.q-chem.com/)


## Specific features
Only **CFOUR 2.1** version supports proper ordering of g-orbitals in MOLDEN file so recommended for use.

Only in case of **NWCHEM** 6.8.1 version able to generate [MOLDEN](https://nwchemgit.github.io/Properties.html#moldenfile) file without a patch.
**molden_norm none** is valid option to generate molden2qmc compatible format.
I want to pay attention to that only spherical-harmonic (5 d, 7 f, 9 g, ...) angular functions supported.

In case of **QChem** also supported only spherical-harmonic angular functions. You should checkout
PURECART variable which handle the angular form of the basis functions.

Latest versions tested:
- TURBOMOLE - V6.6
- Psi4      - 1.1 release
- CFOUR     - 2.1
- ORCA      - Version 4.2.0
- Dalton    - 2018.2
- Molpro    - Version 2012.1
- NWChem    - 6.8.1
- Q-Chem    - 4.4
