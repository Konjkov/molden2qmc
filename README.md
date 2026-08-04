# molden2qmc

Converter from the [MOLDEN](https://www.theochem.ru.nl/molden/molden_format.html) file format
to the `gwfn.data` format used by the [CASINO](https://vallico.net/casinoqmc/what-is-casino/)
quantum Monte Carlo program, and to the [QMCPACK](https://qmcpack.org/) HDF5 format.

The `gwfn.data` file contains the molecular geometry and the molecular orbitals expanded in a
Gaussian basis set. In general, any file conforming to the MOLDEN specification can be used, but
different quantum chemistry codes do not always fully comply with the specification — in such
cases the script corrects the discrepancies, which is why the code that produced the file must
be given explicitly on the command line.

## Installation

```shell
pip install .
```

This installs the `molden2qmc` command. HDF5 output for QMCPACK is optional and needs an extra
dependency:

```shell
pip install '.[qmcpack]'
```

The scripts can also be run directly from a checkout without installing, e.g. `./molden2qmc.py`.

## Usage

```shell
molden2qmc <code> <input_file> [output_file]
```

`<code>` is the number of the quantum chemistry code that produced the MOLDEN file:

| code | program                                          |
|------|--------------------------------------------------|
| 0    | [TURBOMOLE](http://www.turbomole.com/)           |
| 1    | [PSI4](http://www.psicode.org/) or PySCF         |
| 2    | [CFOUR](http://www.cfour.de/) 2.0beta            |
| 3    | [ORCA](https://orcaforum.cec.mpg.de/) 3.X – 5.X  |
| 4    | [DALTON](http://daltonprogram.org/) 2016         |
| 5    | [MOLPRO](https://www.molpro.net/), MOLCAS or openMOLCAS |
| 6    | [NWCHEM](https://nwchemgit.github.io/)           |
| 7    | [QCHEM](http://www.q-chem.com/) 4.X              |

`output_file` defaults to `gwfn.data` (or `Mol.orbs.h5` with `--qmcpack`).

Options:

- `--pseudoatoms none|all|<numbers>` — the presence of a pseudopotential cannot be detected from
  the MOLDEN file, so the atoms it was used for must be listed: `none`, `all`, or a whitespace
  separated list of atom numbers starting from 1.
- `--qmcpack` — write an HDF5 file in the QMCPACK format instead of `gwfn.data`.
- `--cell '1 0 0; 0 1 0; 0 0 1'` — set the lattice vectors of a periodic system, which a MOLDEN
  file does not contain.

Example:

```shell
molden2qmc 3 orca.molden.input gwfn.data --pseudoatoms none
```

### make_geminal_casl.py

Generates the `GEMINAL` section of a CASINO `parameters.casl` file from a `gwfn.data` file in the
given run directory (default: the current one):

```shell
./make_geminal_casl.py [rundir] --channels p:2 d:2 --mirror
```

Without `--channels` an HF-like section is generated. See the docstring at the top of the script
for the meaning of shells, anchors and seeds.

### multideterminant.py

Work in progress, and deliberately not part of the installed package: there is no standard for a
multideterminant extension of the MOLDEN format analogous to MOLDEN itself, so the script has to
parse the native output of each code and remains incomplete. It is kept in the repository and can
be run directly from a checkout.

## Specific features

Only **CFOUR 2.1** supports proper ordering of g-orbitals in the MOLDEN file, so it is the
recommended version.

Only **NWCHEM** 6.8.1 is able to generate a [MOLDEN](https://nwchemgit.github.io/Properties.html#moldenfile)
file without a patch. `molden_norm none` is a valid option for generating a molden2qmc compatible
file. Note that only spherical-harmonic (5 d, 7 f, 9 g, …) angular functions are supported.

For **QChem** only spherical-harmonic angular functions are supported as well; check the
`PURECART` variable, which controls the angular form of the basis functions.

For **ORCA**, code `3` covers versions 3.X through 5.X. Newer versions (6.X) have not been tested;
they may well work, since nothing is known to have changed in the MOLDEN output, but that has not
been verified.

## Tested versions

| program   | version       |
|-----------|---------------|
| TURBOMOLE | V6.6          |
| PSI4      | 1.1 release   |
| CFOUR     | 2.1           |
| ORCA      | 5.X           |
| DALTON    | 2018.2        |
| MOLPRO    | 2012.1        |
| NWCHEM    | 6.8.1         |
| Q-CHEM    | 4.4           |

## Tests

```shell
python -m unittest test_molden2qmc.py
```
