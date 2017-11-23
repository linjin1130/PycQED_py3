import os
import unittest
import numpy as np

try:
    from pycqed.measurement.openql_experiments import single_qubit_oql as sqo
    from pycqed.measurement.openql_experiments.generate_CCL_cfg import  \
        generate_config
    from openql import openql as ql

    rootDir = os.path.dirname(os.path.realpath(__file__))
    curdir = os.path.dirname(__file__)
    config_fn = os.path.join(curdir, 'test_cfg_CCL.json')

    output_dir = os.path.join(curdir, 'test_output')
    ql.set_output_dir(output_dir)

    class Test_configuration_files(unittest.TestCase):

        def test_openQL_config_valid(self):
            test_config_fn = os.path.join(curdir, 'test_gen_cfg_CCL.json')
            generate_config(filename=test_config_fn,
                            mw_pulse_duration=20, ro_duration=300,
                            init_duration=200000)
            # If this compiles we conclude that the generated config is valid
            sqo.AllXY(qubit_idx=0, platf_cfg=test_config_fn)

    class Test_single_qubit_seqs_CCL(unittest.TestCase):

        def test_vsm_timing_cal_sequence(self):
            sqo.vsm_timing_cal_sequence(qubit_idx=1,
                                        platf_cfg=config_fn)

        def test_CW_RO_seq(self):
            sqo.CW_RO_sequence(qubit_idx=0, platf_cfg=config_fn)

        def test_pulsed_spec_seq(self):
            sqo.pulsed_spec_seq(qubit_idx=0, spec_pulse_length=80e-9,
                                platf_cfg=config_fn)

        def test_allxy(self):
            # Only test if it compiles
            sqo.AllXY(qubit_idx=0, platf_cfg=config_fn)

        def test_T1(self):
            # Only test if it compiles
            sqo.T1(times=np.arange(0, 1e-6, 20e-9),
                   qubit_idx=0, platf_cfg=config_fn)

        def test_Ramsey(self):
            # Only test if it compiles
            sqo.Ramsey(times=np.arange(0, 1e-6, 20e-9),
                       qubit_idx=0, platf_cfg=config_fn)

        def test_echo(self):
            # Only test if it compiles
            sqo.echo(times=np.arange(0, 2e-6, 40e-9),
                     qubit_idx=0, platf_cfg=config_fn)

        def test_single_elt_on(self):
            sqo.single_elt_on(qubit_idx=0, platf_cfg=config_fn)

        def test_flipping(self):
            number_of_flips = np.arange(10)
            sqo.flipping(qubit_idx=0, number_of_flips=number_of_flips,
                         platf_cfg=config_fn)

        def test_butterfly(self):
            sqo.butterfly(qubit_idx=0, initialize=True,
                          platf_cfg=config_fn)

        def test_off_on(self):
            sqo.off_on(0, pulse_comb='off', platf_cfg=config_fn)
            sqo.off_on(0, pulse_comb='on', platf_cfg=config_fn)
            sqo.off_on(0, pulse_comb='off_on', platf_cfg=config_fn)


except ImportError as e:
    class TestMissingDependency(unittest.TestCase):

        @unittest.skip('Missing dependency - ' + str(e))
        def test_fail(self):
            pass

if __name__ == '__main__':
    unittest.main()