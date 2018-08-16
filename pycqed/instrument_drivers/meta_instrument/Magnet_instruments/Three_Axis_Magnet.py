
import time
import logging
import numpy as np
from qcodes.utils import validators as vals

from pycqed.measurement import detector_functions as det

from pycqed.measurement import sweep_functions as swf
from pycqed.analysis import analysis_toolbox as a_tools
import logging
from copy import deepcopy,copy

import qcodes as qc
from qcodes.instrument.base import Instrument
from qcodes.utils import validators as vals
from qcodes.instrument.parameter import ManualParameter
from pycqed.instrument_drivers.pq_parameters import InstrumentParameter




class Three_Axis_Magnet(Instrument):


    def __init__(self, name,
                 Coil_X,
                 field_step_x,
                 Coil_Y,
                 field_step_y,
                 Coil_Z,
                 field_step_z,
                 MC_inst,**kw):

        super().__init__(name, **kw)

        self.protection_state=False
        # Set instrumentss
        self.add_parameter('Coil_X', parameter_class=InstrumentParameter)
        self.coil_x = Coil_X
        self.add_parameter('Coil_Y', parameter_class=InstrumentParameter)
        self.coil_y = Coil_Y
        self.add_parameter('Coil_Z', parameter_class=InstrumentParameter)
        self.coil_z = Coil_Z
        self.MC = MC_inst

        self.coil_list = [self.coil_x,self.coil_y,self.coil_z]

        self.add_parameter('field',
                           get_cmd=self.get_field,
                           set_cmd=self.set_field,
                           label='Persistent Field',
                           unit='T',
                           vals=vals.Numbers(min_value=0.),
                           docstring='Persistent absolute magnetic field')
        self.add_parameter('angle',
                           get_cmd=self.get_angle,
                           set_cmd=self.set_angle,
                           label='In-plane angle',
                           unit='deg',
                           vals=vals.Numbers(min_value=-180.,max_value=180.),
                           docstring='Angle of the field wrt. v1')
        self.add_parameter('field_vector',
                           get_cmd=self.get_field_vector,
                           set_cmd=self._set_field_vector,
                           label='Field vector',
                           unit='T',
                           vals=vals.Arrays(),
                           docstring='Field vector in coil frame')
        self.add_parameter('field_steps',
                           parameter_class=ManualParameter,
                           initial_value=np.array([field_step_x,field_step_y,field_step_z]),
                           label='Field steps',
                           vals=vals.Arrays(),
                           docstring='Maximum stepsize when changing field')
        self.add_parameter('v1',
                           initial_value=np.array([0.,0.,1.]),
                           parameter_class=ManualParameter,
                           label='In-plane vector 1',
                           vals=vals.Arrays(),
                           docstring='Vector that defines 0 angle in the sample plane')
        self.add_parameter('v2',
                           initial_value=np.array([0.,1.,0.]),
                           parameter_class=ManualParameter,
                           label='In-plane vector 2',
                           vals=vals.Arrays(),
                           docstring='Second vector in the sample plane')
        self.add_parameter('switch_state',
                           get_cmd=self.get_switch_state,
                           set_cmd=self.set_switch_state,
                           label='Switch State',
                           vals=vals.Enum('SuperConducting','NormalConducting'),
                           docstring='Indicating whether at least one of the\
                                persistent current switches is superconducting')

        self.get_all()

    def get_all(self):
        timestamp = a_tools.latest_data(contains='Define_vector_v1',
                                       return_timestamp=True)[0]
        params_dict = {'v1':'Magnet.v1'}
        data = a_tools.get_data_from_timestamp_list([timestamp], params_dict,
                             filter_no_analysis=False)
        v1 = np.fromstring(data['v1'][0][1:-1], dtype=float, sep=' ')
        self.v1(v1)
        timestamp = a_tools.latest_data(contains='Define_vector_v2',
                               return_timestamp=True)[0]
        params_dict = {'v2':'Magnet.v2'}
        data = a_tools.get_data_from_timestamp_list([timestamp], params_dict,
                             filter_no_analysis=False)
        v2 = np.fromstring(data['v2'][0][1:-1], dtype=float, sep=' ')
        self.switch_state()
        self.v2(v2)
        self.field_vector()
        self.field()
        self.angle()
        self.field_steps()

    def set_current_pos_as_v1(self):
        current_field = self.field_vector()
        self.v1(np.array(current_field))
        self.fake_folder('Define_vector_v1')


    def set_current_pos_as_v2(self):
        current_field = self.field_vector()
        self.v2(np.array(current_field))
        self.fake_folder('Define_vector_v2')

    def get_unit_normal(self):
        '''
        Returns the normal vector of the plane, where:
            v1: is used as the x-axis in the new polar coord system
            v2: defines the plane of the in-plane angle
        '''
        if self.v1 is None:
            raise ValueError('v1 is not defined')
        elif self.v2 is None:
            raise ValueError('v2 is not defined')
        else:
            v_normal = np.cross(self.v1(), self.v2())
            self.unit_normal = v_normal/np.linalg.norm(v_normal)
            return self.unit_normal

    def rotation_matrix(self, theta):

        """
        Return matrix to rotate about axis defined by unit_vector.
        """
        angle = theta*np.pi/180.
        sina = np.sin(angle)
        cosa = np.cos(angle)
        unit_normal = self.get_unit_normal()
        ux, uy, uz = unit_normal
        # rotation matrix around unit vector
        R = np.diag([cosa, cosa, cosa])
        R += np.outer(unit_normal, unit_normal) * (1.0 - cosa)
        R += np.array([[0.0,-uz*sina,uy*sina],
                          [uz*sina, 0.0,-ux*sina],
                          [-uy*sina, ux*sina,  0.0]])
        return R

    def set_field(self,field):
        '''
        Changes the field magnitude
        '''
        current_angle = self.angle()
        target_coords = self.rotation_matrix(current_angle).dot(self.v1())/np.linalg.norm(self.v1())*field
        self._set_field_vector(target_coords)
        return 0.

    def get_field_vector(self):
        x_coord = self.coil_x.field()
        y_coord = self.coil_y.field()
        z_coord = self.coil_z.field()
        return np.array([x_coord,y_coord,z_coord])

    def get_field(self):
        coords = self.field_vector()
        return np.linalg.norm(coords)

    def set_angle(self,angle):
        '''
        Changes the in-plane angle of the field
        '''
        current_field = self.field()
        target_coords = self.rotation_matrix(angle).dot(self.v1())/np.linalg.norm(self.v1())*current_field
        self._set_field_vector(target_coords)
        return 0.

    def get_angle(self):
        coords = self.field_vector()
        phi = np.arccos(np.dot(coords,self.v1())/(np.linalg.norm(coords)*np.linalg.norm(self.v1())))
        return phi*180./np.pi

    def _set_field_vector(self,target_coords):
        '''
        Gets the current coordinates, defines step sizes for all directions,
        steps in a straight line to the target coordinates.
        '''
        for coil in self.coil_list:
            if coil.persistent is True:
                if coil.switch_state()=='SuperConducting':
                    raise ValueError('Coil +'+coil.name+' has a superconducting switch!')
        current_coords = self.field_vector()
        delta = target_coords - current_coords
        max_num_steps = int(max(np.ceil(np.abs(delta)/self.field_steps())))
        stepsizes = delta/max_num_steps
        for tt in range(max_num_steps):
            for ii, coil in enumerate(self.coil_list):
                current_coords[ii]+=stepsizes[ii]
                coil.field(current_coords[ii])
        for ii, coil in enumerate(self.coil_list):
            coil.field(target_coords[ii])

    def get_switch_state(self):
        for coil in self.coil_list:
            if coil.persistent is True:
                if coil.switch_state()=='SuperConducting':
                    return 'SuperConducting'
        return 'NormalConducting'

    def set_switch_state(self,desired_state):
        if desired_state=='NormalConducting':
            for coil in self.coil_list:
                if coil.persistent is True:
                    print('Turning coil '+coil.name+' normal...')
                    coil.bring_source_to_field()
                    coil.switch_state('NormalConducting')
            print('Wait 30 s until the switches are normalconducting')
            time.sleep(30)
            return 0.
        elif desired_state=='SuperConducting':
            for coil in self.coil_list:
                if coil.persistent is True:
                    print('Turning coil '+coil.name+' persistent...')
                    coil.switch_state('SuperConducting')
            print('Wait 120 s until the switches are superconucting')
            time.sleep(120)
            return 0.
        else:
            raise ValueError('Switch state is either "NormalConducting"\
                             or "SuperConducting"!')

    def compensate_for_drift(self,comnpensation_factor=1.05):
        for coil in self.coil_list:
            if coil.persistent is True:
                coil.compensate_for_drift(
                    compensation_factor=comnpensation_factor)
        return 0.

    def fake_folder(self,folder_name):
        '''
        Give folder_name in the form of a string.
        Create a fake folder in the nanowire experiments folder with the desired folder name.
        This is usefull to use when magnetic fields change or when you start a new measurement cycle.
        Such that you can distinguish the different measurement sets.
        '''
        if isinstance(folder_name,str):
            sweep_pts = np.linspace(0, 10, 3)
            self.MC.set_sweep_function(swf.None_Sweep())
            self.MC.set_sweep_points(sweep_pts)
            self.MC.set_detector_function(det.Dummy_Detector_Soft())
            dat = self.MC.run(folder_name)
        else:
            raise ValueError('Please enter a string as the folder name!')


