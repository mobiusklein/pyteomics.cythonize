from os import path
import pyteomics
import pickle
#pyteomics.__path__ = [path.abspath(path.join(path.dirname(__file__), path.pardir, 'pyteomics'))]
import unittest
import random
from pyteomics import auxiliary, mass
from pyteomics import cmass, cparser
cmass.nist_mass = mass.nist_mass
import gzip



class MassTest(unittest.TestCase):
    def setUp(self):
        self.mass_data = {
            'A' : {0: (1.0, 1.0),
                   1: (1.0, 0.5),
                   2: (2.0, 0.5)},
            'B' : {0: (2.0, 1.0),
                   2: (2.0, 0.5),
                   3: (3.0, 0.5)},
            'C' : {0: (3.0, 1.0),
                   3: (3.0, 0.5),
                   4: (4.0, 0.5)},
            'D' : {0: (4.0, 1.0),
                   4: (4.0, 0.5),
                   5: (5.0, 0.5)},
            'E' : {0: (5.0, 1.0),
                   5: (5.0, 0.5),
                   6: (6.0, 0.5)},
            'F' : {0: (6.0, 1.0),
                   6: (6.0, 0.7),
                   7: (7.0, 0.3)},
            'H+': {0: (5.0, 1.0),
                   5: (5.0, 1.0)},
            }

        self.mass_H = cmass.nist_mass['H'][0][0]
        self.mass_O = cmass.nist_mass['O'][0][0]
        self.test_aa_mass = {'X': 1.0, 'Y': 2.0, 'Z': 3.0}
        self.random_peptides = [
            ''.join([random.choice('XYZ') for i in range(20)])
            for i in range(10)]

        self.aa_comp = {'X':   cmass.Composition({'A': 1},
                                                mass_data=self.mass_data),
                        'Y':   cmass.Composition({'B': 1},
                                                mass_data=self.mass_data),
                        'Z':   cmass.Composition({'C': 1},
                                                mass_data=self.mass_data),
                        'H-':  cmass.Composition({'D': 1},
                                                mass_data=self.mass_data),
                        '-OH': cmass.Composition({'E': 1},
                                                mass_data=self.mass_data),
                        }

        self.ion_comp = {'M': cmass.Composition({},
                                               mass_data=self.mass_data),
                         'a': cmass.Composition({'A': -1},
                                               mass_data=self.mass_data)}
        self.mods = {'a': cmass.Composition(A=1),
                'b': cmass.Composition(B=1)}
        self.d = {atom: 1 for atom in 'ABCDE'}

    def test_fast_mass(self):
        for pep in self.random_peptides:
            self.assertAlmostEqual(
                cmass.fast_mass(pep, aa_mass=self.test_aa_mass),
                sum(pep.count(aa) * m
                     for aa, m in self.test_aa_mass.items())
                + self.mass_H * 2.0 + self.mass_O)

    def test_fast_mass2(self):
        for pep in self.random_peptides:
            self.assertAlmostEqual(
                cmass.fast_mass2(pep, aa_mass=self.test_aa_mass),
                sum(pep.count(aa) * m
                     for aa, m in self.test_aa_mass.items())
                + self.mass_H * 2.0 + self.mass_O)


    def test_Composition_dict(self):
        # Test Composition from a dict.
        self.assertEqual(
            cmass.Composition(self.d, mass_data=self.mass_data), self.d)

    def test_Composition_formula(self):
        # Test Composition from a formula.
        self.assertEqual(self.d,
            cmass.Composition(formula='ABCDE',
                         mass_data={atom: {0: (1.0, 1.0)} for atom in 'ABCDE'}))

    def test_Composition_seq(self):
        # Test Composition from a sequence.
        self.assertEqual(self.d,
            cmass.Composition(sequence='XYZ', aa_comp=self.aa_comp))

    def test_Composition_pseq(self):
        # Test Composition from a parsed sequence.
        self.assertEqual(
            cmass.Composition(parsed_sequence=['X', 'Y', 'Z'],
                             aa_comp=self.aa_comp),
            {atom: 1 for atom in 'ABC'})

    def test_Composition_sseq(self):
        # Test Composition from a split sequence.
        self.assertEqual(
            cmass.Composition(split_sequence=[('X',), ('Y',), ('Z',)],
                             aa_comp=self.aa_comp),
            {atom: 1 for atom in 'ABC'})

    def test_Composition_sum(self):
        # Test sum of Composition objects.
        self.assertEqual(
            cmass.Composition(sequence='XXY', aa_comp=self.aa_comp)
            + cmass.Composition(sequence='YZZ', aa_comp=self.aa_comp),
            {atom: 2 for atom in 'ABCDE'})

    def test_Composition_sub(self):
        # Test subtraction of Composition objects
        self.assertEqual({}
            - cmass.Composition(sequence='XYZ', aa_comp=self.aa_comp),
            {atom: -1 for atom in 'ABCDE'})

    def test_Composition_mul(self):
        # Test multiplication of Composition by integers
        self.assertEqual(
                2 * cmass.Composition(sequence='XYZ', aa_comp=self.aa_comp),
                {atom: 2 for atom in 'ABCDE'})
        self.assertEqual(
                cmass.Composition(sequence='XYZ', aa_comp=self.aa_comp) * 2,
                {atom: 2 for atom in 'ABCDE'})

    def test_Composition_positional(self):
        # Test creation from positional args
        ac = self.aa_comp.copy()
        ac.update(self.mods)
        self.assertEqual(cmass.Composition('aXbYZ', aa_comp=ac),
                {'A': 2, 'B': 2, 'C': 1, 'D': 1, 'E': 1})
        self.assertEqual(cmass.Composition('AB2C3', mass_data=self.mass_data),
                {'A': 1, 'B': 2, 'C': 3})

    def test_calculate_mass(self):
        # Calculate mass by a formula.
        self.assertEqual(
            cmass.calculate_mass(formula='ABCDE', mass_data=self.mass_data),
            sum([self.mass_data[atom][0][0] for atom in 'ABCDE']))

        # Calculate mass by a sequence.
        self.assertEqual(
            cmass.calculate_mass(sequence='XYZ',
                                aa_comp=self.aa_comp,
                                mass_data=self.mass_data),
            sum([self.mass_data[atom][0][0] for atom in 'ABCDE']))

        # Calculate mass by a parsed sequence.
        self.assertEqual(
            cmass.calculate_mass(parsed_sequence=['H-','X','Y','Z','-OH'],
                                aa_comp=self.aa_comp,
                                mass_data=self.mass_data),
            sum([self.mass_data[atom][0][0] for atom in 'ABCDE']))

        # Calculate average mass by a formula.
        self.assertEqual(
            cmass.calculate_mass(formula='ABCDE',
                                average=True,
                                mass_data=self.mass_data),
            sum([self.mass_data[atom][isotope][0]
                 * self.mass_data[atom][isotope][1]
                 for atom in 'ABCDE'
                 for isotope in self.mass_data[atom] if isotope != 0]))

        # Calculate m/z of an ion.
        for charge in [1,2,3]:
            self.assertEqual(
                cmass.calculate_mass(formula='ABCDE',
                                    ion_type='M',
                                    charge=charge,
                                    mass_data=self.mass_data),
                cmass.calculate_mass(formula='ABCDE'+'H+%d' % (charge,),
                                    mass_data=self.mass_data))

            self.assertEqual(
                cmass.calculate_mass(formula='ABCDE',
                                    ion_type='M',
                                    charge=charge,
                                    mass_data=self.mass_data),
                (cmass.calculate_mass(formula='ABCDE',
                                    mass_data=self.mass_data)
                 + self.mass_data['H+'][0][0] * charge
                 ) / charge)

            self.assertRaises(
                auxiliary.PyteomicsError,
                cmass.calculate_mass,
                **{'formula': 'ABCDEH+%d' % charge,
                   'ion_type': 'M',
                   'charge': charge,
                   'mass_data': self.mass_data})

        # Sanity check.
        for pep in self.random_peptides:
            self.assertEqual(cmass.calculate_mass(
                sequence=pep, aa_comp=self.aa_comp, mass_data=self.mass_data,
                ion_comp=self.ion_comp),
                cmass.calculate_mass(
                    parsed_sequence=cparser.parse(
                        pep, labels=['X', 'Y', 'Z'], show_unmodified_termini=True),
                    aa_comp=self.aa_comp, mass_data=self.mass_data,
                    ion_comp=self.ion_comp))

    def test_composition_objects_are_pickleable(self):
        dict_ = cmass.Composition(self.d, mass_data=self.mass_data)
        formula = cmass.Composition(formula='ABCDE',
                         mass_data={atom: {0: (1.0, 1.0)} for atom in 'ABCDE'})
        sequence = cmass.Composition(sequence='XYZ', aa_comp=self.aa_comp)
        parsed_sequence = cmass.Composition(parsed_sequence=['X', 'Y', 'Z'],
                             aa_comp=self.aa_comp)
        split_sequence = cmass.Composition(split_sequence=[('X',), ('Y',), ('Z',)],
                             aa_comp=self.aa_comp)

        self.assertEqual(dict_, pickle.loads(pickle.dumps(dict_)))
        self.assertEqual(formula, pickle.loads(pickle.dumps(formula)))
        self.assertEqual(sequence, pickle.loads(pickle.dumps(sequence)))
        self.assertEqual(parsed_sequence, pickle.loads(pickle.dumps(parsed_sequence)))
        self.assertEqual(split_sequence, pickle.loads(pickle.dumps(split_sequence)))


if __name__ == '__main__':
    unittest.main()
