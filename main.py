import numpy as np

def pol_to_comps(mag, ang, unit='rad'):
    """
    polar form to components


    inputs:
        mag - magnitude of the vector
        ang - angle of the vector
        unit - rad for radians deg for degrees"""
    if unit=='deg':
        ang = np.deg2rad(ang)

    unit_x = np.array([1.0, 0.0]) # i vector

    rot = np.array([[np.cos(ang), np.sin(ang)],
                   [-np.sin(ang), np.cos(ang)]]) # basic rotation matrix
    return unit_x @ rot

def comps_to_pol(vec, out_unit='rad'):
    # type checking
    assert isinstance(vec, np.ndarray)
    assert len(vec) == 2

    if out_unit == 'deg':
        ang = np.rad2deg(np.atan2(vec[1], vec[0]))
    else:
        ang = np.atan2(vec[1], vec[0])

    mag = np.sqrt(np.dot(vec, vec))

    return mag, ang

class V_triangle:
    # should be eliminated because you do the transform first or second depending on whether it is a rotor or stator.
    # do it stage based instead idiot
    # also look back at your 351 notes for this, idiot
    # also look back at the class notes for this, idiot
    # stop being so stupid.
    def __innit__(self, c1, U, beta):
        """c is the vector velocity of the incoming flow. Type must be array of dimension 2
            u is the speed of the rotor"""
        
        # check types
        assert isinstance(c1, np.ndarray)
        assert isinstance(U, float)
        assert isinstance(beta, float)

        # assign input values
        self.c1 = c1 # absolute velocity
        self.U = U
        self.beta = beta

        # calculate values in frame a
        self.v1 = c1-U
        self.

 
    def plot(self):
        # for debugging/viewing
        assert False, "stub"
        return -1

class Stage:
    def __innit__(self, c, what_info_do_i_need):
        # rotor
        rotor = V_triangle(c=c)
        rotor.out_velocity()

        # stator
    def get_dhn(self):
        print(f'need to get the de hallard number')
    def get

class Compressor:
    # stack the stages together
    # input parameters:
    #   overall pressure ratio
    #   number of stages
    #   
    def __innit__(self):
        pass

def main():
    # calculate dimensions
    # 
pi_oc = 4.15 # compressure overall pressure ratio
m_dot = 20 # air mass flow, kg/sec
TIT = 1100 # turbine inlet temperature, kelvin

