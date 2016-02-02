# Note to Xiang, remove the imports that are not used :)
import time
import numpy as np
import sys
# import visa
import unittest
import logging

qcpath = 'D:\GitHubRepos\Qcodes'
if qcpath not in sys.path:
    sys.path.append(qcpath)

from qcodes.instrument.visa import VisaInstrument
from qcodes.utils import validators as vals

# cython drivers for encoding and decoding
import pyximport
pyximport.install(setup_args={"script_args": ["--compiler=msvc"],
                              "include_dirs": np.get_include()},
                  reload_support=True)

from ._controlbox import codec as c
from ._controlbox import Assembler
from . import QuTech_ControlBoxdriver as qcb
from ._controlbox import defHeaders_CBox_v3 as defHeaders # File containing bytestring commands

'''
@author: Xiang Fu
The driver for ControlBox version 3. This is inherited from the driver of
ControlBox version 2.
'''

class QuTech_ControlBox_v3(qcb.QuTech_ControlBox):
    def __init__(self, *args, **kwargs):
        super(QuTech_ControlBox_v3, self).__init__(*args, **kwargs)

        # self.add_parameter('acquisition_mode',
        #                    set_cmd=self._do_set_acquisition_mode,
        #                    get_cmd=self._do_get_acquisition_mode,
        #                    vals=vals.Anything())

        self.add_parameter('core_state',
                           set_cmd=self._do_set_core_state,
                           get_cmd=self._do_get_core_state,
                           vals=vals.Anything())
        self.add_parameter('trigger_source',
                           set_cmd=self._do_set_trigger_source,
                           get_cmd=self._do_get_trigger_source,
                           vals=vals.Anything())

        self.set_master_controller_working_state(0, 0, 0)

    def run_test_suite(self):
            from importlib import reload  # Useful for testing
            from ._controlbox import test_suite
            reload(test_suite)
            # pass the CBox to the module so it can be used in the tests
            self.c = c  # make the codec callable from the testsuite
            test_suite.CBox = self
            suite = unittest.TestLoader().loadTestsFromTestCase(test_suite.CBox_tests)
            unittest.TextTestRunner(verbosity=2).run(suite)

    def _do_get_firmware_version(self):
        logging.warning("This function is obselete, now try function get_master_controller_params.")

    def get_master_controller_params(self):
        message = self.create_message(defHeaders.ReadVersion)
        (stat, mesg) = self.serial_write(message)
        if stat:
            param_msg = self.serial_read()

        # Decoding the message
        # version number, first 3 bytes
        v_list = []
        for byte in bytes(param_msg):
            v_list.append(int(byte) - 128)  # -128 is to remove the highest bit 1 for data_bytes
            #print(byte, ": ", "{0:{fill}8b}".format(byte, fill='0'))

        version = str(v_list[0])+'.'+str(v_list[1]) + \
            '.'+str(v_list[2])

        self._acquisition_mode = v_list[3]
        self._core_state = v_list[4]
        self._trigger_source = v_list[5]
        self.adc_offset = (v_list[6] << 4) + v_list[7]
        self.signal_delay= (v_list[8] << 4) + v_list[9]
        self.integration_length = (v_list[10] << 7) + v_list[11]
        a11 = (v_list[12] << 7) + v_list[13]
        a12 = (v_list[14] << 7) + v_list[15]
        a21 = (v_list[16] << 7) + v_list[17]
        a22 = (v_list[18] << 7) + v_list[19]
        self.lin_trans_coeffs =  np.array([a11, a12, a21, a22])
        self.threshold0 = (v_list[20] << 21) + (v_list[21] << 14) + (v_list[22] << 7) + v_list[23];
        self.threshold1 = (v_list[24] << 21) + (v_list[25] << 14) + (v_list[26] << 7) + v_list[27];
        self.log_length = (v_list[28] << 7) + v_list[29];
        self.nr_samples = (v_list[30] << 7) + v_list[31];
        self.avg_size = v_list[32]
        self.nr_averages = 2**self.avg_size

        self.snapshot()

        return version

    def _do_set_log_length(self, length):
        '''
        set the number of measurements of the log in test mode.

        @param length : the number of measurements range (1, 8000)
        @return stat : 0 if the upload succeeded and 1 if the upload failed.
        '''
        v = self.get_master_controller_params()
        if v.startswith('v2.15'):
            n_bytes = 3
        else:
            # logging.warning('Version != 2.15 using old protocol for log length')
            n_bytes = 2

        cmd = defHeaders.UpdateLoggerMaxCounterHeader
        data_bytes = c.encode_byte(length-1, 7,
                                   expected_number_of_bytes=n_bytes)
        # Changed from version 2.15 onwards
        message = c.create_message(cmd, data_bytes)
        (stat, mesg) = self.serial_write(message)
        if stat:
            self._log_length = length
        else:
            raise Exception('Failed to set log_length')
        return (stat, message)

    def _do_set_core_state(self, core_state):
        if self.get('acquisition_mode') is not None:
            tmp_acquisition_mode = self.get('acquisition_mode')
            print('_do_set_core_state\got_acquisition_mode: ', tmp_acquisition_mode)
        else:
            tmp_acquisition_mode = 0   # idle state

        if self.get('trigger_source') is not None:
            tmp_trigger_source = self.get('trigger_source')
            print('_do_set_core_state\got_trigger_source: ', tmp_trigger_source)
        else:
            tmp_trigger_source = 0      # internal trigger

        self.set_master_controller_working_state(core_state, tmp_acquisition_mode,
                                  tmp_trigger_source)

    def _do_set_acquisition_mode(self, acquisition_mode):
        if not('core_state' in self.parameters):
            logging.warning('Core state needs to be fixed for acq')
            return

        if self.get('core_state') is not None:
            tmp_core_state = self.get('core_state')
            print('_do_set_acquisition_mode\got_core_state: ', tmp_core_state)
        else:
            tmp_core_state = 0   # idle state

        if not('trigger_source' in self.parameters):
            logging.warning('MasterController trigger source needs to be fixed for acq')
            return

        if self.get('trigger_source') is not None:
            tmp_trigger_source = self.get('trigger_source')
            print('_do_set_acquisition_mode\got_trigger_source: ', tmp_trigger_source)
        else:
            tmp_trigger_source = 0      # internal trigger

        self.set_master_controller_working_state(tmp_core_state, acquisition_mode,
                                  tmp_trigger_source)

    def _do_set_trigger_source(self, trigger_source):
        if self.get('core_state') is not None:
            tmp_core_state = self.get('core_state')
            print('_do_set_trigger_source\got_core_state: ', tmp_core_state)
        else:
            tmp_core_state = 0   # idle state

        if self.get('acquisition_mode') is not None:
            tmp_acquisition_mode = self.get('acquisition_mode')
            print('_do_set_trigger_source\got_acquisition_mode: ', tmp_acquisition_mode)
        else:
            tmp_acquisition_mode = 0   # idle state

        self.set_master_controller_working_state(tmp_core_state, tmp_acquisition_mode,
                                  trigger_source)

    def _do_get_core_state(self):
        return self._core_state

    def _do_get_trigger_source(self):
        return self._trigger_source

    def _do_get_acquisition_mode(self):
        return self._acquisition_mode

    def _set_touch_n_go_parameters(self):

        '''
        Touch 'n Go is only valid for ControlBox version 2, so this function is
        obselete.
        '''
        print("Warning: Touch 'n Go is only valid for ControlBox version 2. omit")

    def set_master_controller_working_state(self, core_state, acquisition_mode, trigger_source):
        '''
        @param core_states: activate the core or disable it,
                        0 = idle,
                        1 = active
        @param acquisition_modes: the data collection mode in MasterController,
                         0 = idle (default),
                         1 = integration logging mode,
                         2 = integration averaging mode,
                         3 = input averaging mode,
                         4 = integration streaming mode
        @param trigger_sources: trigger source of the MasterController,
                        0 = internal trigger (default),
                        1 = external trigger,
                        2 = mixed trigger
        @return stat : True if the upload succeeded and False if the upload failed
        '''
        core_state = str(core_state)
        acquisition_mode = str(acquisition_mode)
        trigger_source = str(trigger_source)

        acquisition_mode_int = None
        for i in range(len(defHeaders.acquisition_modes)):
            if acquisition_mode.upper() in defHeaders.acquisition_modes[i].upper():
                acquisition_mode_int = i
                break
        if acquisition_mode_int is None:
            raise KeyError('acquisition_mode %s not recognized')
        if acquisition_mode_int == 1 and self._log_length > 8000:
            logging.warning('Log length can be max 8000 in int. log. mode')

        core_state_int = None
        for i in range(len(defHeaders.core_states)):
            if core_state.upper() in defHeaders.core_states[i].upper():
                core_state_int = i
                break
        if core_state_int is None:
            raise KeyError('core_state %s not recognized')

        trigger_source_int = None
        for i in range(len(defHeaders.trigger_sources)):
            if trigger_source.upper() in defHeaders.trigger_sources[i].upper():
                trigger_source_int = i
                break
        if trigger_source_int is None:
            raise KeyError('trigger_state %s not recognized')

        # Here the actual acquisition_mode is set
        cmd = defHeaders.UpdateMCWorkingState
        data_bytes = c.encode_byte(core_state_int, 7, expected_number_of_bytes=1)
        data_bytes += c.encode_byte(acquisition_mode_int, 7, expected_number_of_bytes=1)
        data_bytes += c.encode_byte(trigger_source_int, 7, expected_number_of_bytes=1)
        message = c.create_message(cmd, data_bytes)

        (stat, mesg) = self.serial_write(message)
        if stat:
            self._acquisition_mode = defHeaders.acquisition_modes[acquisition_mode_int]
            self._core_state = defHeaders.core_states[core_state_int]
            self._trigger_source = defHeaders.trigger_sources[trigger_source_int]
        else:
            raise Exception('Failed to set acquisition_mode')

        print('acquisition_mode:', self.get('acquisition_mode'))  # ensure updating of the value
        print('core_state:', self.get('core_state'))
        print('trigger_source:', self.get('trigger_source'))

        return (stat, message)

    def load_instructions(self, asm_file):
        '''
        set the weights of the integregration

        @param instructions : the instructions, an array of 32-bit instructions
        @return stat : 0 if the upload succeeded and 1 if the upload failed.
        '''

        asm = Assembler.Assembler(asm_file)

        instructions = asm.convert_to_instructions()

        if instructions == False :
            print("Error: the assembly file is of errors.")
            return False

        # for instruction in instructions:
        #     print(format(instruction, 'x').zfill(8))

        # Check the instruction list length
        if len(instructions) == 0 :
            raise ValueError("The instruction list is empty.")

        cmd = defHeaders.LoadInstructionsHeader
        data_bytes = bytearray()

        instr_length = 32
        data_bytes.extend(c.encode_array(
                        self.convert_arrary_to_signed(instructions, instr_length),
                        data_bits_per_byte = 4,
                        bytes_per_value = 8))

        message = c.create_message(cmd, bytes(data_bytes))
        (stat, mesg) = self.serial_write(message)

        if not stat:
            raise Exception('Failed to load instructions')

        return (stat, mesg)


    def set_conditional_tape(self, awg_nr, tape_nr, tape):
        '''
        NOTE: The ControlBox does not support timing tape till now.

        set the conditional tape content for an awg

        @param awg : the awg of the dac, (0,1,2).
        @param tape_nr : the number of the tape, integer ranging (0~6)
        @param tape : the array of entries, with a maximum number of entries 512.
            Every entry is an integer has the following structure:
                |WaitingTime (9bits) | PUlse number (3 bits) | EndofSegment marker (1bit)|
            WaitingTime: The waiting time before the end of last pulse or trigger, in ns.
            Pulse number: 0~7, indicating which pulse to be output
            EndofSegment marker: 1 if the entry is the last entry of the tape, otherwise 0.
        @return stat : 0 if the upload succeeded and 1 if the upload failed.

        '''

        print("Timing tape is not supported yet.")
        return False

        length = len(tape)
        tape_addr_width = 9
        entry_length = 9 + 3 + 1

        # Check out of bounds
        if awg_nr < 0 or awg_nr > 2:
            raise ValueError
        if tape_nr < 0 or tape_nr > 6:
            raise ValueError
        if length < 1 or length > 512:
            raise ValueError("The conditional tape only supports a length from 1 to 512.")

        cmd = defHeaders.AwgCondionalTape
        data_bytes = bytearray()

        data_bytes.extend(c.encode_byte(awg_nr, 4,      # add AWG number
                                        expected_number_of_bytes = 1))
        data_bytes.extend(c.encode_byte(tape_nr, 4,     # add the tape number
                                        expected_number_of_bytes = 1))
        data_bytes.extend(c.encode_byte(length-1, 7,    # add the tape length
                          expected_number_of_bytes=np.ceil(tape_addr_width/7.0)))
        data_bytes.extend(c.encode_array(               # add the tape entries
                          convert_arrary_to_signed(tape, entry_length),
                          data_bits_per_byte = 7,
                          expected_number_of_bytes = np.ceil(entry_length/7.0)))

        message = c.create_message(cmd, data_bytes)
        (stat, mesg) = self.serial_write(message)
        return (stat, mesg)

    def set_segmented_tape(self, awg_nr, tape):
        '''
        NOTE: The ControlBox does not support timing tape till now.

        set the conditional tape content for an awg

        @param awg : the awg of the dac, (0,1,2).
        @param tape : the array of entries, with a maximum number of entries 29184.
            Every entry is an integer has the following structure:
                |WaitingTime (9bits) | PUlse number (3 bits) | EndofSegment marker (1bit)|
            WaitingTime: The waiting time before the end of last pulse or trigger, in ns.
            Pulse number: 0~7, indicating which pulse to be output
            EndofSegment marker: 1 if the entry is the last entry of a segment, otherwise 0.
        @return stat : 0 if the upload succeeded and 1 if the upload failed.

        '''

        print("Timing tape is not supported yet.")
        return False


        length = len(tape)
        tape_addr_width = 15
        entry_length = 9 + 3 + 1

        # Check out of bounds
        if awg_nr < 0 or awg_nr > 2:
            raise ValueError("Awg number error!")
        if length < 1 or length > 29184:
            raise ValueError("The segemented tape only supports a length from 1 to 29184.")

        cmd = defHeaders.AwgSegmentedTape
        data_bytes = bytearray()
        data_bytes.extend(c.encode_byte(awg_nr, 4,
                                        expected_number_of_bytes = 1))
        data_bytes.extend(c.encode_byte(length-1, 7,
                          expected_number_of_bytes=np.ceil(tape_addr_width / 7.0)))
        data_bytes.extend(c.encode_array(
                          convert_arrary_to_signed(tape, entry_length),
                          data_bits_per_byte = 7,
                          expected_number_of_bytes = np.ceil(entry_length/7.0)))

        message = self.create_message(cmd, data_bytes)
        (stat, mesg) = self.serial_write(message)
        return (stat, mesg)

    def create_entry(self, interval, pulse_num, end_of_marker):
        '''
        @param interval : The waiting time before the end of last pulse or trigger in ns,
                          ranging from 0ns to 2560ns with minimum step of 5ns.
        @param pulse_num : 0~7, indicating which pulse to be output
        @param end_of_marker : 1 if the entry is the last entry of a segment, otherwise 0.
        '''

        if interval < 0 or interval > 2560:
            raise ValueError
        if pulse_num < 0 or pulse_num > 7:
            raise ValueError
        if end_of_marker < 0 or end_of_marker > 1:
            raise ValueError

        entry_bits = BitArray(Bits(uint=interval, length=9))
        entry_bits.append(BitArray(Bits(uint=pulse_num, length=3)))
        entry_bits.append(BitArray(Bits(uint=end_of_marker, length=1)))
        # print "The entry generated is: ",
        # print entry_bits.uint

        return entry_bits.uint

    def convert_arrary_to_signed(self, unsigned_array, bit_width):
        '''
        Inteprete the input unsinged number array into a signed number array
        based on the given bitwidth.

        @param unsigned_array: the unsigned number array.
        @param bit_width: Bit width of the output signed number.
        '''

        signed_array = []
        for sample in unsigned_array:
            signed_array.append(self.convert_to_signed(sample, bit_width))

        return signed_array

    def convert_to_signed(self, unsigned_number, bit_width):
        '''
        Inteprete the input unsinged number into a signed number given the bitwidth.

        @param unsigned_number: the unsigned number.
        @param bit_width: Bit width of the output signed number.
        '''

        if (not isinstance(unsigned_number, (int))) or (unsigned_number < 0):
            raise ValueError("The number %d should be a positive integer." \
                              % unsigned_number)

        if unsigned_number < 0 or unsigned_number >= 2**bit_width:
            raise ValueError("Given number %d is too large in terms of the \
                              given bit_width %d." % unsigned_number, bit_width)

        if unsigned_number >= 2**(bit_width-1):
            signed_number = unsigned_number - 2**bit_width;
        else:
            signed_number = unsigned_number;

        return signed_number
