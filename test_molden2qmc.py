#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = 'Vladimir Konjkov'
__version__ = '2.6.0'

import unittest
import filecmp
import numpy as np
import molden2qmc
#from shutil import copyfile
#        copyfile('gwfn.data', self.base_dir + test_dir + 'gwfn.data')

molden2qmc.__version__ = __version__


def mo_matrix(m, col=0):
    """
    Convert MO-coefficients to numpy array, remove first four MO witch is 1s of N atoms.

    :param m: Molden converter class i.e. Orca, CFour, Turbomole, ...
    :param col: column which sign set to '+'
    :return: MO-coefficients with out first four 1s-orbitals witch close to degenerate.
    """
    mo = np.empty((m.nbasis_functions(), m.nbasis_functions()), 'd')
    for n, orbital1 in enumerate(m.mo_matrix):
        m = 0
        for ao in orbital1['MO']:
            mo[n,m:m+len(ao['DATA'])] = ao['DATA']
            m += len(ao['DATA'])
    return (mo.T*np.sign(mo[:, col])).T[4:, :]


class test_Turbomole(unittest.TestCase):
    base_dir = 'test/N4/TURBOMOLE/'
    molden_file = 'molden.input'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        turbomole = molden2qmc.Turbomole(open(self.base_dir + test_dir + self.molden_file, "r"))
        turbomole.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/SVP_Turbomole/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(turbomole), mo_matrix(orca), atol=0.0001))

    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ/'
        turbomole = molden2qmc.Turbomole(open(self.base_dir + test_dir + self.molden_file, "r"))
        turbomole.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/cc-pVTZ_Turbomole/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(turbomole), mo_matrix(orca), atol=0.002))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ/'
        turbomole = molden2qmc.Turbomole(open(self.base_dir + test_dir + self.molden_file, "r"))
        turbomole.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/cc-pVQZ_Turbomole/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(turbomole), mo_matrix(orca), atol=0.001))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        molden2qmc.Turbomole(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_PSI4(unittest.TestCase):
    base_dir = 'test/N4/PSI4/'
    molden_file = 'N4.out.n4.molden'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        psi4 = molden2qmc.PSI4(open(self.base_dir + test_dir + self.molden_file, "r"))
        psi4.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/SVP_Psi4/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(psi4), mo_matrix(orca), atol=0.0001))

    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ/'
        psi4 = molden2qmc.PSI4(open(self.base_dir + test_dir + self.molden_file, "r"))
        psi4.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/cc-pVTZ_Psi4/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(psi4), mo_matrix(orca), atol=0.001))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ/'
        psi4 = molden2qmc.PSI4(open(self.base_dir + test_dir + self.molden_file, "r"))
        psi4.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/cc-pVQZ_Psi4/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(psi4), mo_matrix(orca), atol=0.001))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        molden2qmc.PSI4(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_CFour(unittest.TestCase):
    base_dir = 'test/N4/CFOUR/'
    molden_file = 'MOLDEN'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP_patch/'
        cfour = molden2qmc.CFour(open(self.base_dir + test_dir + self.molden_file, "r"))
        cfour.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/SVP_CFOUR/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(cfour), mo_matrix(orca), atol=0.0001))

    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ_patch/'
        cfour = molden2qmc.CFour(open(self.base_dir + test_dir + self.molden_file, "r"))
        cfour.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/cc-pVTZ_CFOUR/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(cfour), mo_matrix(orca), atol=0.001))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ_patch/'
        cfour = molden2qmc.CFour(open(self.base_dir + test_dir + self.molden_file, "r"))
        cfour.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/cc-pVQZ_CFOUR/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(cfour), mo_matrix(orca), atol=0.001))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        molden2qmc.CFour(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_Orca(unittest.TestCase):
    base_dir = 'test/N4/ORCA/'
    molden_file = 'N4.molden.input'

    def test_RHF_SVP_Molpro(self):
        test_dir = 'RHF/SVP_Molpro/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_SVP_Dalton(self):
        test_dir = 'RHF/SVP_Dalton/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_SVP_PSI4(self):
        test_dir = 'RHF/SVP_Psi4/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_SVP_CFOUR(self):
        test_dir = 'RHF/SVP_CFOUR/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_TZVP_Molpro(self):
        test_dir = 'RHF/TZVP_Molpro/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_TZVP_Dalton(self):
        test_dir = 'RHF/TZVP_Dalton/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_QZVP_Molpro(self):
        test_dir = 'RHF/QZVP_Molpro/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_QZVP_Dalton(self):
        test_dir = 'RHF/QZVP_Dalton/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVTZ(self):
        test_dir = 'RHF/cc-pVTZ_Turbomole/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVTZ_PSI4(self):
        test_dir = 'RHF/cc-pVTZ_Psi4/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVQZ(self):
        test_dir = 'RHF/cc-pVQZ_Turbomole/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_RHF_cc_pVQZ_PSI4(self):
        test_dir = 'RHF/cc-pVQZ_Psi4/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        molden2qmc.Orca(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_Dalton(unittest.TestCase):
    base_dir = 'test/N4/DALTON/'
    molden_file = 'molden.inp'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        dalton = molden2qmc.Dalton(open(self.base_dir + test_dir + self.molden_file, "r"))
        dalton.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/SVP_Dalton/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(dalton), mo_matrix(orca), atol=0.001))

    def test_RHF_TZVP(self):
        test_dir = 'RHF/TZVP/'
        dalton = molden2qmc.Dalton(open(self.base_dir + test_dir + self.molden_file, "r"))
        dalton.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/TZVP_Dalton/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(dalton), mo_matrix(orca), atol=0.001))

    def test_RHF_QZVP(self):
        test_dir = 'RHF/QZVP/'
        dalton = molden2qmc.Dalton(open(self.base_dir + test_dir + self.molden_file, "r"))
        dalton.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/QZVP_Dalton/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(dalton, col=2), mo_matrix(orca, col=2), atol=0.002))

    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        molden2qmc.Dalton(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


class test_Molpro(unittest.TestCase):
    base_dir = 'test/N4/MOLPRO/'
    molden_file = 'n4.molden'

    def test_RHF_SVP(self):
        test_dir = 'RHF/SVP/'
        molpro = molden2qmc.Molpro(open(self.base_dir + test_dir + self.molden_file, "r"))
        molpro.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/SVP_Molpro/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(molpro), mo_matrix(orca), atol=0.001))

    def test_RHF_TZVP(self):
        test_dir = 'RHF/TZVP/'
        molpro = molden2qmc.Molpro(open(self.base_dir + test_dir + self.molden_file, "r"))
        molpro.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/TZVP_Molpro/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(molpro), mo_matrix(orca), atol=0.01))

    def test_RHF_QZVP(self):
        test_dir = 'RHF/QZVP/'
        molpro = molden2qmc.Molpro(open(self.base_dir + test_dir + self.molden_file, "r"))
        molpro.gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))
        orca = molden2qmc.Orca(open('test/N4/ORCA/RHF/QZVP_Molpro/N4.molden.input', "r"))
        self.assertTrue(np.allclose(mo_matrix(molpro), mo_matrix(orca), atol=0.01))

    @unittest.skip("Not implemented")
    def test_UHF_SVP(self):
        test_dir = 'UHF/SVP/'
        molden2qmc.Molpro(open(self.base_dir + test_dir + self.molden_file, "r")).gwfn()
        self.assertTrue(filecmp.cmp(self.base_dir + test_dir + 'gwfn.data', 'gwfn.data'))


if __name__ == '__main__':
    unittest.main()
    #suite = unittest.TestLoader().loadTestsFromTestCase(test_PSI4)
    #unittest.TextTestRunner(verbosity=2).run(suite)
