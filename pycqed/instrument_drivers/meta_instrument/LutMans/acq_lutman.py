from .base_lutman import Base_LutMan
from pycqed.measurement.waveform_control_CC import waveform as wf
from qcodes.instrument.parameter import ManualParameter
from qcodes.utils import validators as vals
import numpy as np
import copy as copy
from collections import OrderedDict
from pycqed.measurement import detector_functions as det
from collections import defaultdict

class Base_Acq_LutMan(Base_LutMan):
	    def __init__(self, name, feedline_nrs=tuple:(0,1), **kw):
        self.add_parameter('weight_nr_samples', unit='samples', vals=vals.Ints(0,2**20), #taking arbitrarily large max val
                           parameter_class=ManualParameter,
                           initial_value=4096)
        self.add_parameter('averages', unit='samples', vals=vals.Ints(0,2**20), #taking arbitrarily large max val
                           parameter_class=ManualParameter,
                           initial_value=2*11)
        self.add_parameter('digitized', vals=vals.Bool(),
                           initial_value=False,
                           parameter_class=ManualParameter)
        self.add_parameter('integration_length'.format(res), unit='s',
                           vals=vals.Numbers(1e-9, 8000e-9),
                           parameter_class=ManualParameter,
                           initial_value=2000e-9)
	    self.add_parameter('weight_type',
        	               initial_value='DSB',
        	               vals=vals.Enum('SSB', 'DSB', 'optimal', 'optimal IQ'),
        	               docstring=ro_acq_docstr,
        	               parameter_class=ManualParameter)
        self._feedline_nrs = feedline_nrs
        if self._feedline_nrs==(0):
            self._resonator_codeword_bit_mapping=[[0,2,3,5,6]]
        elif self._feedline_nrs==(1): 
            self._resonator_codeword_bit_mapping=[[1,4]]
        elif self._feedline_nrs==(0,1):
            self._resonator_codeword_bit_mapping=[[0,2,3,5,6],[1,4]]
        else:
            raise NotImplementedError(
              'hardcoded for feedline 0 and/or 1 of Surface-7')
                logging.info(__name__ + ' : Initializing instrument')
        super().super().__init__(name, **kw)
        self.add_parameter('sampling_rate', unit='Hz',
                           vals=vals.Numbers(1, 1e10),
                           initial_value=1e9,
                           parameter_class=ManualParameter)
        for feedline_nr in self._feedline_nrs():
            self.add_parameter(
                'instr_feed_{}'.format(feedline_nr), parameter_class=InstrumentRefParameter, docstring=(
                    "Name of the acquisition instrument used"),
                vals=vals.Strings())

    def _add_waveform_parameters(self):
        """
        mixer corrections are done globally per feedline
        """
        #demodulation mixer predistortion no yet implemented, parameters are there
        # for feedline_nr in self._feedline_nrs():
        #     # self.add_parameter('mixer_apply_correction_matrix_feed_{}'.format(feedline_nr),
        #     #                    vals=vals.Bool(),
        #     #                    parameter_class=ManualParameter,
        #     #                    initial_value=False)
        #     # self.add_parameter('mixer_alpha_feed_{}'.format(feedline_nr), vals=vals.Numbers(),
        #     #                    parameter_class=ManualParameter,
        #     #                    initial_value=1.0)
        #     # self.add_parameter('mixer_phi_feed_{}'.format(feedline_nr), vals=vals.Numbers(), unit='deg',
        #     #                    parameter_class=ManualParameter,
        #     #                    initial_value=0.0)
        #     # self.add_parameter('mixer_offs_I_feed_{}'.format(feedline_nr), unit='V',
        #     #                    parameter_class=ManualParameter, initial_value=0)
        #     # self.add_parameter('mixer_offs_Q_feed_{}'.format(feedline_nr), unit='V',
        #     #                    parameter_class=ManualParameter, initial_value=0)  
        #Generating integration parameters for each qubit
        for ress in self._resonator_codeword_bit_mapping:
            for res in ress:
                self.add_parameter('R{}_modulation'.format(res),
                               vals=vals.Numbers(), unit='Hz',
                               parameter_class=ManualParameter,
                               initial_value=20.0e6)
                self.add_parameter('R{}_opt_weights_I',
                               vals=vals.Arrays(),
                               label='Optimized weights for I channel',
                               parameter_class=ManualParameter)
            	self.add_parameter('R{}_opt_weights_Q',
                               vals=vals.Arrays(),
                               label='Optimized weights for Q channel',
                               parameter_class=ManualParameter)
            	self.add_parameter('R{}_digitized_threshold', unit='V',
                               initial_value=0,
                               parameter_class=ManualParameter)
    
    def generate_standard_waveforms(self):
        """
        Generating base cosines and sines for demodulation
        """
        self._wave_dict = {}
        for ress in self._resonator_codeword_bit_mapping:
            for res in ress:
            	IF = self.get('R{}_modulation'.format(res))
            	# 1. Generate weigt envelopes
                ## DSB weights
                trace_length = self.weight_nr_samples()
    	        tbase = np.arange(0, trace_length/self.sampling_rate(), 
                            1/self.sampling_rate())
    	        cos = np.array(np.cos(2*np.pi*IF*tbase))
    	        sin = np.array(np.sin(2*np.pi*IF*tbase))
    	        self._wave_dict['R{}_cos'.format(res)] = cos
    	        self._wave_dict['R{}_sin'.format(res)] = sin

    def get_acq_instr(self, qubit_nr):
    	'''gets the right acquisition instrument depending on the qubit number''' 
    	i, res_nr = self._resonator_codeword_bit_mapping.index(qubit_nr)
        feed_line_nr = self.feedline_nrs[i]
        instr = self.parameters['instr_feed_{}'.format(feedline_nr)].get_instr()
        return instr
    
    def assign_weight_channels_hard_coded(self, qubit_nrs):
    	"""
    	static weight assginment for single or multi-qubit experiments. 
    	Assginment follows he CC hardcoded bitmapping.
    	This mapping is required when using CC for feedback as the DIO bits are hardcoded
    	to the weight channels.

    	input
    		qubit_nrs
    	returns
            feedline_nrs        : list of feedline numbers to assign the instrument
            ro_ch_idx           : channel indices for acquisition
            value_names         : convenient labels
    	"""
    	channels_list = []
        for qubit_nr in qubit_nrs:
	    	i, weight_I_channel = self._resonator_codeword_bit_mapping.index(qubit_nr)
    		feed_line_nr = self.feedline_nrs[i]

	    	if self.weight_type()=='optimal':
	            channels_list.append((feed_line_nr, weight_I_channel,
	                'feed{} w{} q{}'.format(feed_line_nr, weight_I_channel, qubit_nr)))
		    else:
		    	if len(qubit_nrs)==1:
		    		last_weight_channel = self._weight_channels_per_instr-1
			        if weight_I_channel==last_weight_channel:
			    		weight_Q_channel = 0
			    	else:
			    		weight_Q_channel = weight_I_channel + 1 #borrowing integration channel from the next resonator
			        channels_list.append((feed_line_nr, weight_I_channel,
			            'feed{} w{} q{} I'.format(feed_line_nr, weight_I_channel, qubit_nr)))
			        channels_list.append((feed_line_nr, weight_I_channel,
			            'feed{} w{} q{} Q'.format(feed_line_nr, weight_I_channel, qubit_nr)))
			    else:
			     	NotImplementedError('multi-qubit hardcoded assignemnt is only possible optimal weights')
        
        feed_line_nrs = [feed_line_nr for feed_line_nr, _, _ in channels_list]
        weight_channels = [ch for _, ch, _ in channels_list]
        value_names = [n for _, _, n in channels_list]
    	return feed_line_nrs, weight_channels, value_names

    def assign_weight_channels_dynamic(self, qubit_nrs):
        """
        dynamic weight assginment for single or multi-qubit experiments. 
        single or two quadrature is possible. Dynamic assignment is not compatible
        with feedback.

    	input
    		qubit_nrs
    	returns
            feedline_nrs        : list of feedline numbers to assign the instrument
            ro_ch_idx           : channel indices for acquisition
            value_names         : convenient labels
    	"""
    	if self.weight_type() == 'optimal':
            nr_of_weight_channels_per_qubit = 1
        else:
            nr_of_weight_channels_per_qubit = 2

        next_unused_weight_channels = defaultdict(int)
    	channels_list = []
        for qubit_nr in reversed(qubit_nrs):
        	# ensures that the LSQ (last one) get's assigned the lowest ch_idx
	    	i, weight_I_channel = self._resonator_codeword_bit_mapping.index(qubit_nr)
    		feed_line_nr = self.feedline_nrs[i]

            # allocate different acquisition channels
            # optimal weight channel assignment
          	if nr_of_acquisition_channels_per_qubit==1:
          		weight_channel = next_unused_weight_channels[feed_line_nr]
            	next_unused_weight_channels[feed_line_nr] += 1
	            channels_list.append((feed_line_nr, weight_channel,
	                'feed{} w{} q{}'.format(feed_line_nr, weight_channel, qubit_nr)))
            # I and Q channel assignment
            elif nr_of_acquisition_channels_per_qubit==2:
            	weight_channel = next_unused_weight_channels[feed_line_nr]
            	next_unused_weight_channels[feed_line_nr] += 1
	            channels_list.append((feed_line_nr, weight_channel,
	                'feed{} w{} q{} I'.format(feed_line_nr, weight_channel, qubit_nr)))
	            weight_channel = next_unused_weight_channels[feed_line_nr]
            	next_unused_weight_channels[feed_line_nr] += 1
	            channels_list.append((feed_line_nr, weight_channel,
	                'feed{} w{} q{} Q'.format(feed_line_nr, weight_channel, qubit_nr)))
        
        feed_line_nrs = [feed_line_nr for feed_line_nr, _, _ in channels_list]
        weight_channels = [ch for _, ch, _ in channels_list]
        value_names = [n for _, _, n in channels_list]
    	return feed_line_nrs, weight_channels, value_names

	def get_input_average_detector(self, CC, qubit_nrs,  **kw):
    	self.prepare_weights_and_detectors(CC=CC, qubit_nrs=qubit_nrs)
        return self._input_average_detector

    def get_int_avg_det(self, CC, qubit_nrs, **kw):
    	self.prepare_weights_and_detectors(CC=CC, qubit_nrs=qubit_nrs) 
    	return self._int_avg_det

    def get_int_avg_det_single(self, CC, qubit_nrs, **kw):
    	self.prepare_weights_and_detectors(CC=CC, qubit_nrs=qubit_nrs)  
        return self._int_avg_det_single

    def get_int_log_det(self, CC, qubit_nrs, **kw):
    	self.prepare_weights_and_detectors(CC=CC, qubit_nrs=qubit_nrs)  
        return self._int_log_det

class UHFQC_Acq_LutMan(Base_Acq_LutMan):
	def __init__(self, name, **kw):
        super().__init__(name, **kw)
        self.sampling_rate(1.8e9)
        self.weight_nr_samples(4096)
        self._weight_channels_per_instr = 9 

    def load_waveform_and_threshold_onto_instr_lookuptable(self, qubit_nr, weight_channels, regenerate_waveforms: bool=False):
    	if regenerate_waveforms:
    		generate_standard_waveforms()
        cos = self._wave_dict['R{}_cos'.format(qubit_nr)]
    	sin = self._wave_dict['R{}_sin'.format(qubit_nr)]

    	weight_I_channel = weight_channels[0]

    	#setting the threshold for digitized readout
    	if self.digitized():
	    	threshold = self.get('R{}_digitized_threshold'.format(qubit_nr))
	    	instr.set('quex_thres_{}_level'.format(weight_I_channel), threshold)

	   	# getting the right instrument
	    instr = self.get_acq_instr(qubit_nr)
    	if self.weight_type() == 'DSB':
    		weight_Q_channel = weight_channels[1]
    		instr.set('quex_wint_weights_{}_real'.format(weight_I_channel),
	                 np.array(cos))
	    	instr.set('quex_wint_weights_{}_real'.format(weight_Q_channel), 
	                 np.array(sin))
	    	instr.set('quex_rot_{}_real'.format(weight_I_channel), 2.0)
        	instr.set('quex_rot_{}_imag'.format(weight_I_channel), 0.0)
        	instr.set('quex_rot_{}_real'.format(weight_Q_channel), 2.0)
        	instr.set('quex_rot_{}_imag'.format(weight_Q_channel), 0.0)
        elif self.weight_type() == 'SSB':
        	weight_Q_channel = weight_channels[1]
        	instr.set('quex_wint_weights_{}_real'.format(weight_I_channel),
                 np.array(cos))
        	instr.set('quex_wint_weights_{}_imag'.format(weight_I_channel),
                 np.array(sin))
        	instr.set('quex_wint_weights_{}_real'.format(weight_Q_channel),
                 np.array(sin))
        	instr.set('quex_wint_weights_{}_imag'.format(weight_Q_channel),
                 np.array(cos))
			instr.set('quex_rot_{}_real'.format(weight_I_channel), 1.0)
			instr.set('quex_rot_{}_imag'.format(weight_I_channel), 1.0)
			instr.set('quex_rot_{}_real'.format(weight_Q_channel), 1.0)
			instr.set('quex_rot_{}_imag'.format(weight_Q_channel), -1.0)
		elif 'optimal' in self.weight_type():
			instr.set('quex_wint_weights_{}_real'.format(
                        weight_I_channel),
                        self.get('R{}_opt_weights_I'.format(qubit_nr)))
	        instr.set('quex_wint_weights_{}_imag'.format(
	            		weight_I_channel),
	            		self.get('R{}_opt_weights_Q'.format(qubit_nr)))
	        instr.set('quex_rot_{}_real'.format(
	            weight_I_channel), 1.0)
	        instr.set('quex_rot_{}_imag'.format(
	            weight_I_channel), -1.0) 
	        if self.weight_type() == 'optimal IQ':
	        	weight_Q_channel = weight_channels[1]
	        	instr.set('quex_wint_weights_{}_real'.format(
                    weight_Q_channel),
                        self.get('R{}_opt_weights_I'.format(qubit_nr)))
	        	instr.set('quex_wint_weights_{}_imag'.format(
	            	weight_Q_channel),
	            		self.get('R{}_opt_weights_Q'.format(qubit_nr)))
	        	instr.set('quex_rot_{}_real'.format(
	            	weight_Q_channel), 1.0)
	        	instr.set('quex_rot_{}_imag'.format(
	            	weight_Q_channel), 1.0) 
	    

   	def prepare_weights_and_detectors(self, CC, qubit_nrs, hardcoded_assignment=True):
   		#assigning the weight channels
   		if hardcoded_assignment:
	   		feed_line_nrs, weight_channels, value_names = 
	   					assign_weight_channels_hard_coded(qubit_nrs=qubit_nrs)
	   	else:
	        feed_line_nrs, weight_channels, value_names = 
	        			assign_weight_channels_dynamic(qubit_nrs=qubit_nrs)

	    #loading the integration weights
        for i, qubit_nr in enumerate(qubit_nrs):
        	self.generate_standard_waveforms()
        	if self.weight_type()=='optimal':
        		qubit_weight_channels = list(weight_channels[i])
        	else: 
        		qubit_weight_channels = weight_channels[2*i:2*i+1]
		    self.load_waveform_onto_instr_lookuptable(qubit_nr=qubit_nr, 
        					weight_channels=qubit_weight_channels, 
        					regenerate_waveforms=False)
        
        #rearranging detector values per feedline to allow combined detector functions
        sort_args = numpy.argsort(feed_line_nrs)
        invert_sort_args=sort_args.argsort()#used later for unsorting
		sorted_feed_line_nrs = [feed_line_nrs[i] for i in sort_args]
		sorted_qubits = [qubits[i] for i in sort_args]
		sorted_channel_names = [channel_names[i] for i in sort_args]

		sorted_feed_line_set=list(set(sorted_feed_line_nrs))

		input_avg_dets = []
		input_avg_dets = []
		int_avg_det_singles = []
		int_log_dets = []
		single_qubit_statistics_logging_dets = []

		# prepare individual detectors
        for feedline_nr in sorted_feed_line_set:
	        #picking the right instruments based on the feedline
	        instr = self.parameters['instr_feed_{}'.format(feedline_nr)].get_instr()
        	#selecting the qubits and value names for this feedline
        	for i, sorted_feed_line_nr in enumerate(sorted_feed_line_nrs):
        		sorted_feedline_qubits = []
        		sorted_feedline_channel_names = []
        		if sorted_feed_line_nr==feedline_nr:
        			sorted_feedline_qubits.append(sorted_qubits[i])
        			sorted_feedline_channel_names.append(sorted_channel_names[i])

        	#choosing the right integration output paths
   			if self.weight_type() == 'optimal':
            	result_logging_mode = 'lin_trans'
            	if self.digitized():
                	result_logging_mode = 'digitized'
        	else:
            	result_logging_mode = 'raw'

	        input_avg_det = det.UHFQC_input_average_detector(
	            UHFQC=instr,
	            AWG=CC,
	            nr_averages=self.averages(),
	            nr_samples=self.weight_nr_samples(), 
	            **kw)
	        input_avg_dets.append(input_avg_det)

	        int_avg_det = det.UHFQC_integrated_average_detector(
	        	UHFQC=instr,
	        	AWG=CC,
	            channels=weight_channels,
	            result_logging_mode=result_logging_mode,
	            nr_averages=self.averages(),
	            integration_length=self.integration_length(), **kw)
	        int_avg_det.value_names = value_names
	        int_avg_dets.append(int_avg_det)

	        int_avg_det_single = det.UHFQC_integrated_average_detector(
	            UHFQC=instr, 
	            AWG=CC,
	            channels=weight_channels,
	            result_logging_mode=result_logging_mode,
	            nr_averages=self.averages(),
	            real_imag=True, single_int_avg=True,
	            integration_length=self.integration_length(), **kw)
	        int_avg_det_single.value_names = value_names
	        int_avg_det_singles.append(int_avg_det_single)

	        int_log_det = det.UHFQC_integration_logging_det(
	            UHFQC=instr, 
	            AWG=CC,
	            channels=weight_channels,
	            result_logging_mode=result_logging_mode,
	            integration_length=self.integration_length(), **kw)
	        int_log_det.value_names = value_names
	        int_log_dets.append(int_log_det)

	        single_qubit_statistics_logging_det = det.UHFQC_single_qubit_statistics_logging_det(
	            UHFQC=instr,
	            AWG=CC, 
	            nr_shots=4*4095,
	            integration_length=self.integration_length(),
	            channel=weight_channels,
	            statemap={'0': '1', '1': '0'})
	        single_qubit_statistics_logging_dets.append(single_qubit_statistics_logging_det)

	    #combine detector lists to individual multi-feedline detector functions
	    self._input_avg_det = det.Multi_detector(detectors=input_avg_dets, 
	    	det_idx_suffix=False, name=input_avg_dets[0].name, sort_args=invert_sort_args)
	    self._int_avg_det = det.Multi_detector(detectors=int_avg_dets, 
	    	det_idx_suffix=False, name=int_avg_det[0].name, sort_args=invert_sort_args)
	    self._int_avg_det_single = det.Multi_detector(detectors=int_avg_det_singles, 
	    	det_idx_suffix=False, name=int_avg_det_single[0].name, sort_args=invert_sort_args)
	    self._int_log_det = det.Multi_detector(detectors=int_log_dets, 
	    	det_idx_suffix=False, name=int_log_det[0].name, sort_args=invert_sort_args)
	    self._single_qubit_statistics_logging_det = det.Multi_detector(detectors=single_qubit_statistics_logging_dets, 
	    	det_idx_suffix=False, name=single_qubit_statistics_logging_dets[0].name,, sort_args=invert_sort_args)

















    
