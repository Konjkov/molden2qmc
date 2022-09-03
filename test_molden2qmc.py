#!/usr/bin/env python3

__author__ = 'Vladimir Konjkov'
__version__ = '4.0.3'

import unittest
import filecmp
import numpy as np
import molden2qmc

molden2qmc.__version__ = __version__

np.set_printoptions(threshold=np.inf, suppress=True, linewidth=10000)


def mo_matrix(m, col=0, skip=4, spin='Alpha'):
    """
    Convert MO-coefficients to numpy array, remove first four MO witch is 1s of N atoms.

    :param m: Molden converter class i.e. Orca, CFour, Turbomole, ...
    :param col: column which sign set to '+'
    :return: MO-coefficients without first four 1s-orbitals witch close to degenerate.
    """
    mo = np.empty((m.nbasis_functions(), m.nbasis_functions()), 'd')
    for n, orbital in enumerate(m.mo_matrix):
        m = 0
        if orbital['SPIN'] != spin:
            continue
        for ao in orbital['MO']:
            mo[n, m:m+len(ao['DATA'])] = ao['DATA']
            m += len(ao['DATA'])
    return (mo.T*np.sign(mo[:, col])).T[skip:, :]


class test_Turbomole(unittest.TestCase):
    base_dir = 'test/N4/TURBOMOLE/'
    molden_file = 'molden.input'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            turbomole = molden2qmc.Turbomole(f)
        turbomole.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/SVP_Turbomole/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(turbomole), mo_matrix(orca), atol=0.0001))

    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            turbomole = molden2qmc.Turbomole(f)
        turbomole.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/cc-pVTZ_Turbomole/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(turbomole), mo_matrix(orca), atol=0.002))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            turbomole = molden2qmc.Turbomole(f)
        turbomole.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/cc-pVQZ_Turbomole/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(turbomole), mo_matrix(orca), atol=0.001))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            turbomole = molden2qmc.Turbomole(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_PSI4(unittest.TestCase):
    base_dir = 'test/N4/PSI4/'
    molden_file = 'N4.n4.molden'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            psi4 = molden2qmc.PSI4(f)
        psi4.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA4/RHF/SVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(psi4), mo_matrix(orca), atol=0.0001))

    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            psi4 = molden2qmc.PSI4(f)
        psi4.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/cc-pVTZ/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(psi4), mo_matrix(orca), atol=0.001))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            psi4 = molden2qmc.PSI4(f)
        psi4.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/cc-pVQZ/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(psi4), mo_matrix(orca), atol=0.001))
    

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            psi4 = molden2qmc.PSI4(f)
        psi4.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        #with open('test/N4/ORCA4/UHF/SVP/N4.molden.input', "r") as f:
        #    orca = molden2qmc.Orca(f)
        #self.assertTrue(np.allclose(mo_matrix(psi4), mo_matrix(orca), atol=0.0001))


class test_CFour(unittest.TestCase):
    base_dir = 'test/N4/CFOUR/'
    molden_file = 'MOLDEN'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            cfour = molden2qmc.CFour(f)
        cfour.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/SVP_CFOUR/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(cfour), mo_matrix(orca), atol=0.0001))

    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            cfour = molden2qmc.CFour(f)
        cfour.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/cc-pVTZ_CFOUR/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(cfour), mo_matrix(orca), atol=0.001))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            cfour = molden2qmc.CFour(f)
        cfour.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/cc-pVQZ_CFOUR/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(cfour), mo_matrix(orca), atol=0.001))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.CFour(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_Orca(unittest.TestCase):
    base_dir = 'test/N4/ORCA/'
    molden_file = 'N4.molden.input'

    def test_RHF_SVP_Dalton(self):
        test_dir = 'RHF/SVP_Dalton/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_SVP_CFOUR(self):
        test_dir = 'RHF/SVP_CFOUR/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_TZVP_Dalton(self):
        test_dir = 'RHF/TZVP_Dalton/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_QZVP_Dalton(self):
        test_dir = 'RHF/QZVP_Dalton/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVDZ(self):
        test_dir = 'RHF/cc-pVDZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVTZ_Turbomole(self):
        test_dir = 'RHF/cc-pVTZ_Turbomole/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVQZ_Turbomole(self):
        test_dir = 'RHF/cc-pVQZ_Turbomole/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Orca(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_Dalton(unittest.TestCase):
    base_dir = 'test/N4/DALTON/'
    molden_file = 'molden.inp'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            dalton = molden2qmc.Dalton(f)
        dalton.gwfn()
        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/SVP_Dalton/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(dalton), mo_matrix(orca), atol=0.001))

    def test_RHF_TZVP(self):
        test_dir = 'RHF/TZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            dalton = molden2qmc.Dalton(f)
        dalton.gwfn()
        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/TZVP_Dalton/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(dalton), mo_matrix(orca), atol=0.0005))

    def test_RHF_QZVP(self):
        test_dir = 'RHF/QZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            dalton = molden2qmc.Dalton(f)
        dalton.gwfn()
        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/QZVP_Dalton/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(dalton, col=2), mo_matrix(orca, col=2), atol=0.001))

#    def test_ROHF_SVP(self):
#        test_dir = 'ROHF/SVP/'
#        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
#            dalton = molden2qmc.Dalton(f)
#        dalton.gwfn()
#        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
#        with open('test/N4/ORCA4/ROHF/SVP/N4.molden.input', "r") as f:
#            orca = molden2qmc.Orca(f)
#        self.assertTrue(np.allclose(mo_matrix(dalton), mo_matrix(orca), atol=0.0001))


class test_Molpro(unittest.TestCase):
    base_dir = 'test/N4/MOLPRO/'
    molden_file = 'n4.molden'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molpro = molden2qmc.Molpro(f)
        molpro.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/SVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(molpro), mo_matrix(orca), atol=0.001))

    def test_RHF_TZVP(self):
        test_dir = 'RHF/TZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molpro = molden2qmc.Molpro(f)
        molpro.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/TZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(molpro), mo_matrix(orca), atol=0.01))

    def test_RHF_QZVP(self):
        test_dir = 'RHF/QZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molpro = molden2qmc.Molpro(f)
        molpro.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/QZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(molpro), mo_matrix(orca), atol=0.01))

    @unittest.skip("Not implemented")
    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molden2qmc.Molpro(f).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_NwChem(unittest.TestCase):
    base_dir = 'test/N4/NWCHEM/'
    molden_file = 'N4.molden'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f)
        nwchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/SVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(nwchem), mo_matrix(orca), atol=0.001))

    def test_RHF_TZVP(self):
        test_dir = 'RHF/TZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f)
        nwchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/TZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(nwchem), mo_matrix(orca), atol=0.001))

    def test_RHF_QZVP(self):
        test_dir = 'RHF/QZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f)
        nwchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/QZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(nwchem), mo_matrix(orca), atol=0.001))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f)
        nwchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        #with open('test/N4/ORCA4/UHF/SVP/N4.molden.input', "r") as f:
        #    orca = molden2qmc.Orca(f)
        #self.assertTrue(np.allclose(mo_matrix(nwchem), mo_matrix(orca), atol=0.001))

    @unittest.skip("Cartesian basis not supported in NWChem")
    def test_RHF_SVP_cart(self):
        test_dir = 'RHF/SVP_cart/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f)
        nwchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/SVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(nwchem), mo_matrix(orca), atol=0.001))

    @unittest.skip("Cartesian basis not supported in NWChem")
    def test_RHF_TZVP_cart(self):
        test_dir = 'RHF/TZVP_cart/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f)
        nwchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/TZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(nwchem), mo_matrix(orca), atol=0.001))

    @unittest.skip("Cartesian basis not supported in NWChem")
    def test_RHF_QZVP_cart(self):
        test_dir = 'RHF/QZVP_cart/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f)
        nwchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/QZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(nwchem), mo_matrix(orca), atol=0.001))


class test_QChem(unittest.TestCase):
    base_dir = 'test/N4/QCHEM/'
    molden_file = 'N4.molden'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f)
        qchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/SVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(qchem), mo_matrix(orca), atol=0.001))

    def test_RHF_TZVP(self):
        test_dir = 'RHF/TZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f)
        qchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/TZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(qchem), mo_matrix(orca), atol=0.001))

    def test_RHF_QZVP(self):
        test_dir = 'RHF/QZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f)
        qchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/QZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(qchem), mo_matrix(orca), atol=0.001))

    @unittest.skip("This basis is segmented contracted in QChem")
    def test_RHF_cc_pVDZ(self):
        test_dir = 'RHF/cc-pVDZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f)
        qchem.gwfn()
        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/cc-pVDZ/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(qchem), mo_matrix(orca), atol=0.001))

    @unittest.skip("This basis is segmented contracted in QChem")
    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f)
        qchem.gwfn()
        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/cc-pVTZ/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(qchem), mo_matrix(orca), atol=0.001))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f)
        qchem.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/cc-pVQZ/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(qchem), mo_matrix(orca), atol=0.001))

    def test_UHF_cc_pVDZ(self):
        test_dir = 'UHF/cc-pVDZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f)
        qchem.gwfn()


class test_Orca4(unittest.TestCase):
    base_dir = 'test/N4/ORCA4/'
    molden_file = 'N4.molden.input'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            orca4 = molden2qmc.Orca(f)
        orca4.gwfn()
        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/SVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(orca4), mo_matrix(orca), atol=0.001))

    def test_RHF_TZVP(self):
        test_dir = 'RHF/TZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            orca4 = molden2qmc.Orca(f)
        orca4.gwfn()
        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/TZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(orca4), mo_matrix(orca), atol=0.001))

    def test_RHF_QZVP(self):
        test_dir = 'RHF/QZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            orca4 = molden2qmc.Orca(f)
        orca4.gwfn()
        #self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        with open('test/N4/ORCA/RHF/QZVP/N4.molden.input', "r") as f:
            orca = molden2qmc.Orca(f)
        self.assertTrue(np.allclose(mo_matrix(orca4), mo_matrix(orca), atol=0.001))

import os

class test_QMCPACK(unittest.TestCase):
    base_dir = 'test/N4/QMCPACK/'
    molden_file = 'N4.molden'

    def test_MOLPRO_RHF_SVP(self):
        test_dir = 'MOLPRO/RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molpro = molden2qmc.Molpro(f,qmcpack_normalization=True)
        molpro.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5')
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_ORCA4_UHF_SVP(self):
        test_dir = 'ORCA4/UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            orca = molden2qmc.Orca(f,qmcpack_normalization=True)
        orca.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5')
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_ORCA4_RHF_SVP(self):
        test_dir = 'ORCA4/RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            orca = molden2qmc.Orca(f,qmcpack_normalization=True)
        orca.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5')
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_CFOUR_UHF_SVP(self):
        test_dir = 'CFOUR/UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            cfour = molden2qmc.CFour(f,qmcpack_normalization=True)
        cfour.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5')
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_MOLPRO_RHF_TZVP(self):
        test_dir = 'MOLPRO/RHF/TZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            molpro = molden2qmc.Molpro(f,qmcpack_normalization=True)
        molpro.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_NWCHEM_RHF_SVP(self):
        test_dir = 'NWCHEM/RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f,qmcpack_normalization=True)
        nwchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_NWCHEM_UHF_SVP(self):
        test_dir = 'NWCHEM/UHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f,qmcpack_normalization=True)
        nwchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_NWCHEM_RHF_TZVP(self):
        test_dir = 'NWCHEM/RHF/TZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            nwchem = molden2qmc.NwChem(f,qmcpack_normalization=True)
        nwchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_PSI4_RHF_ccpVTZ(self):
        test_dir = 'PSI4/RHF/cc-pVTZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            psi4 = molden2qmc.PSI4(f,qmcpack_normalization=True)
        psi4.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_PSI4_RHF_SVP(self):
        test_dir = 'PSI4/RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            psi4 = molden2qmc.PSI4(f,qmcpack_normalization=True)
        psi4.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_QCHEM_RHF_ccpVDZ(self):
        test_dir = 'QCHEM/RHF/cc-pVDZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f,qmcpack_normalization=True)
        qchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_QCHEM_RHF_ccpVTZ(self):
        test_dir = 'QCHEM/RHF/cc-pVTZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f,qmcpack_normalization=True)
        qchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_QCHEM_RHF_SVP(self):
        test_dir = 'QCHEM/RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f,qmcpack_normalization=True)
        qchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_QCHEM_RHF_SVP(self):
        test_dir = 'QCHEM/UHF/cc-pVDZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f,qmcpack_normalization=True)
        qchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_QCHEM_RHF_TZVP(self):
        test_dir = 'QCHEM/RHF/TZVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            qchem = molden2qmc.QChem(f,qmcpack_normalization=True)
        qchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_QCHEM_d_ane(self):
        test_dir = '../../QCHEM/d-ane/'
        ane_file='methane.molden'
        with open(self.base_dir + test_dir + ane_file, "r") as f:
            qchem = molden2qmc.QChem(f,qmcpack_normalization=True)
        qchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_QCHEM_f_ane(self):
        test_dir = '../../QCHEM/f-ane/'
        ane_file='methane.molden'
        with open(self.base_dir + test_dir + ane_file, "r") as f:
            qchem = molden2qmc.QChem(f,qmcpack_normalization=True)
        qchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_QCHEM_g_ane(self):
        test_dir = '../../QCHEM/g-ane/'
        ane_file='methane.molden'
        with open(self.base_dir + test_dir + ane_file, "r") as f:
            qchem = molden2qmc.QChem(f,qmcpack_normalization=True)
        qchem.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5')
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_TURBOMOLE_RHF_SVP(self):
        test_dir = 'TURBOMOLE/RHF/SVP/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            turbomole = molden2qmc.Turbomole(f,qmcpack_normalization=True)
        turbomole.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

    def test_TURBOMOLE_RHF_ccpVTZ(self):
        test_dir = 'TURBOMOLE/RHF/cc-pVTZ/'
        with open(self.base_dir + test_dir + self.molden_file, "r") as f:
            turbomole = molden2qmc.Turbomole(f,qmcpack_normalization=True)
        turbomole.qmcpack()
        res=os.system('h5diff -d 0.000001 Mol.orbs.h5 ' + self.base_dir + test_dir + 'Mol.orbs.h5') 
        os.remove("Mol.orbs.h5")
        self.assertTrue(res == 0)

if __name__ == '__main__':
    """python -m unittest test_molden2qmc.test_QChem"""
    unittest.main()
