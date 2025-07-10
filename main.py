import numpy as np
import matplotlib.pyplot as plt


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
    # process:
    #   take the inlet velocity
    #   adjust it to the frame of reference of the cascade (be it rotor or stator)
    #   calculate deflection
    #   flow leaves (all outflows are relative to the cascade)
    def __init__(self, v_inlet, v_blade, turn_angle):
        """c is the vector velocity of the incoming flow. Type must be array of dimension 2
            u is the speed of the rotor"""
        
        # velocity of the blade (positive means right to left motion of the blade)
        # I expect rotors to be positive and stators to be negative
        assert isinstance(v_blade, float)
        self.v_blade = v_blade # 52ish?

        # v_inlet = 0, 30
        # velocity of the air at the inlet, adjusted for blade velocity
        assert isinstance(v_inlet, np.ndarray), type(v_inlet)
        self.rel_v_inlet = np.array([v_inlet[0] + v_blade,  
                                     v_inlet[1]])

        # angle of the turning angle in radians
        assert isinstance(turn_angle, float)
        self.turn_angle = turn_angle

        # calculate relative exit velocity
        rotation_induced_angle = np.atan2(self.rel_v_inlet[0], self.rel_v_inlet[1]) # angle of attack relative to the blade
        outlet_angle = rotation_induced_angle + turn_angle
        self.v_outlet = np.array([v_inlet[1]*np.tan(np.pi/2 - outlet_angle),
                                  v_inlet[1]])
 
        print(f'initial angle: 0\neffective angle after transform: { 180/3.14*np.atan2(self.rel_v_inlet[0], self.rel_v_inlet[1])}\nplus turning angle: {180/3.14*outlet_angle}\noutlet velocity: {self.v_outlet}')

    def plot(self, title='velocity triangle', verbose=True):
        # for debugging/viewing

        print(f'absolute inlet velocity: {'amogus'}\n\
              relative inlet velocity: {self.rel_v_inlet}\n\
              outlet velocity: {self.v_outlet}')

        fig, ax = plt.subplots()

        ax.set_xscale('linear')
        ax.set_yscale('linear')

        # arrow because quiver is comically broken
        ax.arrow(0, 0,                     0,                      30,                     color='k', head_width=0.5, head_length=1.0)
        ax.arrow(0, 0,                     self.rel_v_inlet[0],    self.rel_v_inlet[1],    color='r', head_width=0.5, head_length=1.0)
        ax.arrow(0, self.rel_v_inlet[1],   self.v_blade,           0,                      color='g', head_width=0.5, head_length=1.0)
        ax.arrow(0, 0,                     self.v_outlet[0],       self.v_outlet[1],       color='b', head_width=0.5, head_length=1.0)
        fig.suptitle(f'{title}') # broken

        if verbose:
            plt.show()
        elif not verbose:
            plt.savefig(f'{title}.png')
            plt.clf()

class Stage:
    def __init__(self, c, rpm, r):
        # rotor
        v_blade = (rpm*6.28/3600)*(r) # v = omega*r

        rotor = V_triangle(c, v_blade, np.deg2rad(20))
        v2 = rotor.v_outlet()

        # stator
        stator = V_triangle(v2, -v_blade, -np.deg2rad(20))
        v3 = stator.v_outlet()

    def get_dhn(self):
        print(f'need to get the de hallard number')
    def get_stuff(self):
        print(f'cuh')

class Compressor:
    # stack the stages together
    # input parameters:
    #   overall pressure ratio
    #   number of stages
    #   
    def __init__(self):
        pass

def main():
    # calculate dimensions
    # do something

    # velocity triangles!
    v_blade = (20000*6.28/3600)*(1.5) # v = omega*r

    rotor = V_triangle(np.array([0.0, 30.0]), v_blade, np.deg2rad(20))
    v2 = rotor.v_outlet
    rotor.plot()

if __name__ == "__main__":
    main()
pi_oc = 4.15 # compressure overall pressure ratio
m_dot = 20 # air mass flow, kg/sec
TIT = 1100 # turbine inlet temperature, kelvin

