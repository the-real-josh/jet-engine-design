import numpy as np
import matplotlib.pyplot as plt

# constants and other parameters
gamma = 1.4
R = 287.05 # J/kg K
cp = (gamma*R/(gamma-1))
cv = (cp-R)

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

        assert isinstance(v_inlet, np.ndarray)
        self.abs_v_inlet = v_inlet

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
 
    def plot(self, title='velocity triangle', verbose=True):
        # for debugging/viewing

        print(f'absolute inlet velocity: {np.sqrt(np.dot(self.abs_v_inlet, self.abs_v_inlet))}\n\
        relative inlet velocity: {self.rel_v_inlet}\n\
        outlet velocity: {self.v_outlet}')

        fig, ax = plt.subplots()

        ax.set_xscale('linear')
        ax.set_yscale('linear')

        # arrow because quiver is comically broken
        ax.arrow(0,                         0,                     self.abs_v_inlet[0],    self.abs_v_inlet[1],    color='k', head_width=0.5, head_length=1.0)
        ax.arrow(0,                         0,                     self.rel_v_inlet[0],    self.rel_v_inlet[1],    color='r', head_width=0.5, head_length=1.0)
        ax.arrow(self.abs_v_inlet[0],       self.rel_v_inlet[1],   self.v_blade,           0,                      color='g', head_width=0.5, head_length=1.0)
        ax.arrow(0,                         0,                     self.v_outlet[0],       self.v_outlet[1],       color='b', head_width=0.5, head_length=1.0)
        fig.suptitle(f'{title}') 

        plt.legend([f'absolute inlet v',
                    f'relative inlet v',
                    f'blade velocity (1d)',
                    f'outlet velocity'])

        if verbose:
            plt.show()
        elif not verbose:
            plt.savefig(f'{title}.png')
            plt.clf()

class Stage:
    """inputs: 
    c - inlet velocity
    rpm - rpm of the stage
    r - radius of the section you are examining
    T_inlet - temperature of the air in the inlet
    
    public values:
        output velocity (self.v2)
        output 
        output velocity"""

    # NOTE: please add real imperfect gas behavior as it can be very easy in python    
    # BUG: need to account for inefficiency
    def __init__(self,
                 v_inlet: np.ndarray, # meters per second
                  rpm: float, # revolutions per minute
                    r: float, # meter
                      T_inlet:float, # kelvin
                        p_inlet:float, # pascals
                         rot_defl_ang: float, # radians
                           stat_defl_ang: float): # radians
        
        omega_blade = (rpm*6.28/3600.0)       # speed of rotation in rad/sec
        v_blade = omega_blade*(r)           # assumes constant r
        h_inlet = cp*T_inlet                # specific enthalpy for the gas 
        self.T_inlet = T_inlet

        norm = np.linalg.norm

        # rotor
        self.rotor = V_triangle(v_inlet, v_blade, rot_defl_ang)
        v1_5  = self.rotor.v_outlet

        # stator
        self.stator = V_triangle(v1_5, -v_blade, -stat_defl_ang)
        self.v2 = self.stator.v_outlet

        # euler's equation for turbomachinery
        self.w = (omega_blade*r*v_inlet - omega_blade*r*v1_5) # specific work (energy per mass flow); NOTE: assuming that this includes all enthalpy added (including velocity)

        # outlet enthalpy and temp
        h_outlet = self.w + h_inlet - 0.5*self.v2**2 # account for static, not stagnation
        self.T_outlet = h_outlet/cv

        # flow coefficient
        self.phi = v_inlet/v_blade

        # worst-case mach number
        a = np.sqrt(gamma*R*T_inlet) # maybe T_inlet is wrong, but it certainly will result in a worst-case M
        M = v1_5 / a # v1.5 is probably wrong

        # degree of reaction (called Λ, but python and I hate non-roman variable names)
        # tan(alpha2) - tan(alpha1) = tan(beta1) - tan(beta2)
        self.DRXN = 1 - (norm(v_inlet) / (2*v_blade))*(np.tan(rot_defl_ang) - np.tan(stat_defl_ang))

        # de haller number
        self.DHN = norm(self.rotor.v_outlet) / norm(self.rotor.rel_v_inlet)

        # outlet temperature
        self.T_outlet = h_outlet/cp

        # assumes polytropic efficiency = 0.90
        # no comprendo
        self.isentropic_p_outlet = (self.T_outlet/T_inlet)**(gamma/(gamma-1))
        self.poly_n = 0.90
        self.p_outlet = (self.T_outlet/T_inlet)**(self.poly_n/(self.poly_n-1))

    def plot_triangles(self):
        self.rotor.plot()
        self.stator.plot()
    def print_stats(self):
        status_text = [f'stats for my stage:',
        f'De Hallard number: {self.DHN}',
        f'Degree of reaction: {self.DRXN}',
        f'flow coefficient: {self.phi}',
        f'stage work: {self.w}']
        print('\n'.join(status_text))
    def plot_stats(self):
        # use this to plot a graph of all the numbers across the stage or something (idk)
        assert False, 'function stub'
        return -1
    def plot_mollier(self):
        # use this to plot a mollier diagram of the compresion process

        # need to get s as a function of temperature or enthalpy
        # need to do this for isentropic and polytropic processes
        # s_polytrope = lambda T: (self.T_outlet/T_inlet)**(gamma/(gamma-1))
        # s_isentrope = lambda s:

        # s_points = np.linspace(something)
        # h_points = s_isentrope(s_points) # or something

        # plt.plot(s_points, h_points)
        # plt.show()


class Compressor:
    
    # stack the stages together
    # input parameters:
    #   overall pressure ratio
    #   number of stages
    #   
    def __init__(self):
        stage_rpm = 15000.0
        stage1 = Stage(v_inlet=np.array([0.0, 30.0]),
                    rpm=stage_rpm,
                     r=1.5,
                       T_inlet=300.0,
                        p_inlet=101325.0,
                         rot_defl_ang=np.deg2rad(20),
                          stat_defl_ang=np.deg2rad(10))
        stage1.print_stats()
def main():
    c = Compressor()

    # terms to define compressor
    pi_oc = 4.15 # compressure overall pressure ratio
    m_dot = 20 # air mass flow, kg/sec
    TIT = 1100 # turbine inlet temperature, kelvin
    pass

if __name__ == "__main__":
    main()


