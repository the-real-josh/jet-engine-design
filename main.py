import numpy as np
import matplotlib.pyplot as plt

# gas constants for air
gamma = 1.4
R = 287.05 # J/kg K
cp = (gamma*R/(gamma-1))
cv = (cp-R)


# helper functions
norm = np.linalg.norm
def pol_to_comps(mag, ang, unit='rad'):
    """
    vector in polar form to 2D vector in components

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
    """2D vector in components into polar form"""

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
    """ process:
       1) take the inlet velocity
       2) adjust it to the frame of reference of the cascade (be it rotor or stator)
       3) calculate deflection
       4) flow leaves (all outflows are relative to the cascade) """
    def __init__(self, v_inlet: np.ndarray,
                  v_blade: float,
                    turn_angle: float):
        """ v_inlet is the vector velocity of the incoming flow. Type must be array of dimension 2
            v_blade is the scalar speed of the rotor
            turn_angle is the cascade's turning angle in radians"""
        
        # velocity of the blade (positive means right to left motion of the blade)
        # I expect rotors to be positive and stators to be negative

        # type checking
        assert isinstance(v_blade, float)
        self.v_blade = v_blade
        assert isinstance(v_inlet, np.ndarray)
        self.abs_v_inlet = v_inlet
        assert isinstance(v_inlet, np.ndarray), type(v_inlet)
        assert isinstance(turn_angle, float)

        # class input
        self.turn_angle = turn_angle

        # account for the speed of the blade AND the inlet velocity (vector sum)
        self.rel_v_inlet = np.array([v_inlet[0] + v_blade,  
                                     v_inlet[1]])

        # calculate relative exit velocity
        rotation_induced_angle = np.atan2(self.rel_v_inlet[0], self.rel_v_inlet[1]) # angle of attack relative to the blade
        outlet_angle = rotation_induced_angle + turn_angle
        self.v_outlet = np.array([v_inlet[1]*np.tan(np.pi/2 - outlet_angle),
                                  v_inlet[1]])
 
    def plot(self, title='velocity triangle', verbose=True):
        """ for debugging/viewing
            Note that velocity triangles are drawn upside down in here for convenience."""

        print(f'absolute inlet velocity: {np.sqrt(np.dot(self.abs_v_inlet, self.abs_v_inlet))}\n\
        relative inlet velocity: {self.rel_v_inlet}\n\
        outlet velocity: {self.v_outlet}')

        # plot velocity vectors
        # using plt.arrow because plt.quiver is comically broken
        fig, ax = plt.subplots()
        ax.arrow(0,                         0,                     self.abs_v_inlet[0],    self.abs_v_inlet[1],    color='k', head_width=0.5, head_length=1.0)
        ax.arrow(0,                         0,                     self.rel_v_inlet[0],    self.rel_v_inlet[1],    color='r', head_width=0.5, head_length=1.0)
        ax.arrow(self.abs_v_inlet[0],       self.rel_v_inlet[1],   self.v_blade,           0,                      color='g', head_width=0.5, head_length=1.0)
        ax.arrow(0,                         0,                     self.v_outlet[0],       self.v_outlet[1],       color='b', head_width=0.5, head_length=1.0)
        fig.suptitle(f'{title}') 

        # legend (in order)
        plt.legend([f'absolute inlet v',
                    f'relative inlet v',
                    f'blade velocity (1d)',
                    f'outlet velocity'])

        # save graphs or view graphs depending on parameter "verbose"
        if verbose:
            plt.show()
        elif not verbose:
            plt.savefig(f'{title}.png')
            plt.clf()


class Stage:
    """inputs: 
        v_inlet                 vector velocity of inlet, m/s
        rpm                     float compressor's rotational speed in revolutions per minute
        r                       float average radius of compressor stage in meters
        T_inlet                 float stage inlet temperature in k
        p_inlet                 float stage inlet pressure in pa
        rot_defl_ang            float rotor deflection angle in radians
        stat_defl_ang           float stator deflection angle in radian
    
    class values (public):
        T_inlet                 stage inlet temperature in k
        p_inlet                 stage inlet pressure in pa
        stator                  class instance of the stator
        rotor                   class instance of the rotor
        w                       specific work from the compressor
        T_outlet                emperature of stage outlet gas, Kelvin
        isentropic_p_outlet     pressure of the stage outlet in Pa, if the stage was isentropic
        poly_n                  polytropic efficiency (assumed to be 0.9)
        p_outlet                stage outlet pressure based on the polytropic efficiency, in Pa
        phi                     stage flow coefficient
        DRXN                    stage degree of reaction
        DHN                     stage de Haller number
        """

    # NOTE: please add real imperfect gas behavior as it can be very easy in python    
    def __init__(self,
                 v_inlet: np.ndarray, # meters per second
                  rpm: float, # revolutions per minute
                    r: float, # meter
                      T_inlet:float, # kelvin
                        p_inlet:float, # pascals
                         rot_defl_ang: float, # radians
                           stat_defl_ang: float): # radians
        
        omega_blade = (rpm*6.28/3600.0)         # speed of rotation in rad/sec
        v_blade = omega_blade*(r)               # assumes constant r
        h_inlet = cp*T_inlet                    # specific enthalpy for the gas 
        self.T_inlet = T_inlet                  # stage inlet temperature in k
        self.p_inlet = p_inlet                  # stage inlet pressure in pa

        # rotor of the stage
        self.rotor = V_triangle(v_inlet, v_blade, rot_defl_ang)
        v1_5  = self.rotor.v_outlet

        # stator of the stage
        self.stator = V_triangle(v1_5, -v_blade, -stat_defl_ang)
        self.v2 = self.stator.v_outlet

        # euler's equation for turbomachinery - since the compressor q
        # specific work (energy per mass flow); 
        # NOTE: assuming that this includes all enthalpy added (including velocity)
        # this will be negative, as the compressor REQUIRES energy to operate 
        self.w = (omega_blade*r*norm(v_inlet) - omega_blade*r*norm(v1_5)) 

        # outlet enthalpy and temp
        # negative self.w because the work the work that is put into the system (negative) ADDS to the energy in the system
        # why 0.5*(norm(v_inlet)**2 - norm(self.v2))**2? To account for static change, not stagnation
        h_outlet = -self.w + h_inlet + 0.5*(norm(v_inlet)**2 - norm(self.v2))**2 

        # get outlet temperature by perfect gas laws
        self.T_outlet = h_outlet/cp

        # get outlet pressures by polytropic gas laws
        # assumes polytropic efficiency = 0.90
        # I don't quite understand the reasoning behind polytropic efficiencies
        self.isentropic_p_outlet = (self.T_outlet/T_inlet)**(gamma/(gamma-1))
        self.poly_n = 0.90
        self.p_outlet = (self.T_outlet/T_inlet)**(self.poly_n/(self.poly_n-1))

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

    def plot_triangles(self):
        """Plots the two velocity triangles for the stage"""
        self.rotor.plot()
        self.stator.plot()

    def print_stats(self, verbose=True):
        """ prints out important data for the current stage
            returns a string containing all the stats"""
        status_text = [f'stats for my stage:',
        f'De Hallard number: {self.DHN}',
        f'Degree of reaction: {self.DRXN}',
        f'flow coefficient: {self.phi}',
        f'stage work: {self.w}',
        f'stage pressure ratio: {self.p_outlet/self.p_inlet}',
        f'stage temperature ratio: {self.T_outlet/self.T_inlet}']
        if verbose:
            print('\n'.join(status_text))
        return status_text

    def plot_mollier(self, verbose=True):
        """ Plots a mollier diagram of the stage's compression process

            verbose:
                If verbose, then the plots are shown and then discarded
                if not verbose, then plots will be written on the global plt object, allowing other methods to annotate after this method."""
        # use this to plot a mollier diagram of the compresion process
        s0 = lambda p,T: cp*np.log(T/273.16) - R*np.log(p/101325.0) + 6608.1 # returns initial "absolute" entropy value in kj/kg K for a perfect gas

        # need to get s as a function of temperature or enthalpy
        # need to do this for isentropic and polytropic processes
        p_polytropic_ratio = lambda T: (T/self.T_inlet)**(self.poly_n/(self.poly_n-1))
        p_isentropic_ratio = lambda T: (T/self.T_inlet)**(gamma/(gamma-1))

        s_isentrope = lambda T: cp*np.log(T/self.T_inlet) - R*np.log(p_isentropic_ratio(T)) + s0(self.p_inlet, self.T_inlet)
        s_polytrope = lambda T: cp*np.log(T/self.T_inlet) - R*np.log(p_polytropic_ratio(T)) + s0(self.p_inlet, self.T_inlet)

        # y axis points (usually temperature or enthalpy)
        T_points = np.linspace(self.T_inlet, self.T_outlet, num=50)

        # plot, label, show
        plt.plot(s_isentrope(T_points), T_points); plt.plot(s_polytrope(T_points), T_points)
        plt.ylabel(f'Temperature (K)'); plt.xlabel(f'Entropy (Kj/Kg K)'); plt.title(f'Mollier diagram for compression process')
        if verbose:
            plt.legend([f'isentropic', f'polytropic process with n={self.poly_n}'])
            plt.show()


class Compressor:
    # just one stage for now
    # work in progress

    # my eventual plan:
        # stack the stages together
        # input parameters:
        #       overall pressure ratio
        #       number of stages
        # a while loop (or something) to generate stages until the pressure ratio has been reached and flow has not been reversed (or any other funny business) 
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
        stage1.plot_mollier()
        stage1.plot_triangles()

def main():
    c = Compressor()

    # terms to define compressor (incorporate later), work in progress
    pi_oc = 4.15 # compressure overall pressure ratio
    m_dot = 20 # air mass flow, kg/sec
    TIT = 1100 # turbine inlet temperature, kelvin
    pass

# program entry point
if __name__ == "__main__":
    main()


