import numpy as np
from unittest import TestCase
from zlib import crc32

from pycqed.measurement.randomized_benchmarking.clifford_group import(
    clifford_lookuptable, clifford_group_single_qubit)

import pycqed.measurement.randomized_benchmarking.randomized_benchmarking \
    as rb

from pycqed.measurement.randomized_benchmarking.clifford_decompositions \
    import(gate_decomposition)

from pycqed.measurement.randomized_benchmarking import \
    two_qubit_clifford_group as tqc
from pycqed.measurement.randomized_benchmarking.generate_clifford_hash_tables import construct_clifford_lookuptable

class TestLookuptable(TestCase):
    def test_unique_mapping(self):
        for row in clifford_lookuptable:
            self.assertFalse(len(row) > len(set(row)))

    def test_sum_of_rows(self):
        expected_sum = np.sum(range(len(clifford_group_single_qubit)))
        for row in clifford_lookuptable:
            self.assertEqual(np.sum(row), expected_sum)

    def test_element_index_in_group(self):
        for row in clifford_lookuptable:
            for el in row:
                self.assertTrue(el < len(clifford_group_single_qubit))


class TestCalculateNetClifford(TestCase):
    def test_identity_does_nothing(self):
        id_seq = np.zeros(5)
        net_cl = rb.calculate_net_clifford(id_seq)
        self.assertEqual(net_cl, 0)

        for i in range(len(clifford_group_single_qubit)):
            id_seq[3] = i
            net_cl = rb.calculate_net_clifford(id_seq)
            self.assertEqual(net_cl, i)

    def test_pauli_squared_is_ID(self):
        for cl in [0, 3, 6, 9, 12]:  # 12 is Hadamard
            net_cl = rb.calculate_net_clifford([cl, cl])
            self.assertEqual(net_cl, 0)


class TestRecoveryClifford(TestCase):
    def testInversionRandomSequence(self):
        random_cliffords = np.random.randint(0, len(clifford_group_single_qubit), 100)
        net_cl = rb.calculate_net_clifford(random_cliffords)

        for des_cl in range(len(clifford_group_single_qubit)):
            rec_cliff = rb.calculate_recovery_clifford(net_cl, des_cl)
            comb_seq = np.append(random_cliffords, rec_cliff)

            comb_net_cl_simple = rb.calculate_net_clifford([net_cl, rec_cliff])
            comb_net_cl = rb.calculate_net_clifford(comb_seq)

            self.assertEqual(comb_net_cl, des_cl)
            self.assertEqual(comb_net_cl_simple, des_cl)


class TestRB_sequence(TestCase):
    def test_net_cliff(self):
        for i in range(len(clifford_group_single_qubit)):
            rb_seq = rb.randomized_benchmarking_sequence(500, desired_net_cl=i)
            net_cliff = rb.calculate_net_clifford(rb_seq)
            self.assertEqual(net_cliff, i)

    def test_seed_reproduces(self):
        rb_seq_a = rb.randomized_benchmarking_sequence(500, seed=5)
        rb_seq_b = rb.randomized_benchmarking_sequence(500, seed=None)
        rb_seq_c = rb.randomized_benchmarking_sequence(500, seed=5)
        rb_seq_d = rb.randomized_benchmarking_sequence(500, seed=None)
        self.assertTrue((rb_seq_a == rb_seq_c).all())
        self.assertTrue((rb_seq_a != rb_seq_b).any)
        self.assertTrue((rb_seq_c != rb_seq_b).any)
        self.assertTrue((rb_seq_b != rb_seq_d).any)


class TestGateDecomposition(TestCase):
    def test_unique_elements(self):
        for gate in gate_decomposition:
            self.assertEqual(gate_decomposition.count(gate), 1)

    def test_average_number_of_gates(self):
        from itertools import chain
        avg_nr_gates = len(list(chain(*gate_decomposition)))/24
        self.assertEqual(avg_nr_gates, 1.875)


######################################################################
# Two qubit clifford group below
######################################################################

class TestHashedLookuptables(TestCase):
    def test_single_qubit_hashtable_constructed(self):
        hash_table = construct_clifford_lookuptable(tqc.SingleQubitClifford,
                                                    np.arange(24))

        for i in range(24):
            Cl = tqc.SingleQubitClifford(i)
            target_hash = crc32(Cl.pauli_transfer_matrix.tobytes())
            table_idx = hash_table.index(target_hash)
            self.assertEqual(table_idx, i)

    def test_single_qubit_hashtable_file(self):
        hash_table = tqc.get_single_qubit_clifford_hash_table()

        for i in range(24):
            Cl = tqc.SingleQubitClifford(i)
            target_hash = crc32(Cl.pauli_transfer_matrix.tobytes())
            table_idx = hash_table.index(target_hash)
            self.assertEqual(table_idx, i)



