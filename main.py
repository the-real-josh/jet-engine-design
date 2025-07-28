import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math

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
        rotation_induced_angle = np.atan2(self.rel_v_inlet[0], self.rel_v_inlet[1]) # angle of attack relative to the bladehd600 sennheiser
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
        ax.arrow(0,                         0,                     self.abs_v_inlet[0],    self.abs_v_inlet[1],    color='k', head_width=5.0, head_length=5.0)
        ax.arrow(0,                         0,                     self.rel_v_inlet[0],    self.rel_v_inlet[1],    color='r', head_width=5.0, head_length=5.0)
        ax.arrow(self.abs_v_inlet[0],       self.rel_v_inlet[1],   self.v_blade,           0,                      color='g', head_width=5.0, head_length=5.0)
        ax.arrow(0,                         0,                     self.v_outlet[0],       self.v_outlet[1],       color='b', head_width=5.0, head_length=5.0)
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


class Stage_1D:
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
    __total_instance_counter = 0

    # NOTE: please add real imperfect gas behavior as it can be very easy in python    
    def __init__(self,
                 v_inlet: np.ndarray, # meters per second
                  rpm: float, # revolutions per minute
                    r: float, # meter
                      T_inlet:float, # kelvin
                        p_inlet:float, # pascals
                         rot_defl_ang: float, # radians
                           stat_defl_ang: float, # radians
                            s_inlet=-1.0):
        
        Stage_1D.__total_instance_counter += 1 # increase the counter of the protected value
        self.instance_number = Stage_1D.__total_instance_counter

        omega_blade = (rpm*6.28/60.0)         # speed of rotation in rad/sec
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
        # NOTE: this is only valid for constant velocity along the length of the blade (only radially narrow streamtubes) 
        self.w = (omega_blade*r*norm(v_inlet) - omega_blade*r*norm(v1_5)) 

        # outlet enthalpy and temp (static), assumes roughly constant axial velocity
        # negative self.w because the work the work that is put into the system (negative) ADDS to the energy in the system
        # why 0.5*(norm(v_inlet)**2 - norm(self.v2))**2? To account for static change, not stagnation
        h_outlet = -self.w + h_inlet # + 0.5*(norm(v_inlet)**2 - norm(self.v2))**2 

        # get outlet temperature by perfect gas laws
        self.T_outlet = h_outlet/cp

        # get outlet pressures by polytropic gas laws
        # assumes polytropic efficiency = 0.90
        # I don't quite understand the reasoning behind polytropic efficiencies
        self.isentropic_p_outlet = self.p_inlet*(self.T_outlet/T_inlet)**(gamma/(gamma-1))
        self.poly_n = 0.90
        self.p_outlet = self.p_inlet*(self.T_outlet/T_inlet)**((gamma/(gamma-1))*(self.poly_n))

        # flow coefficient
        self.phi = norm(v_inlet)/v_blade

        # worst-case mach number
        a = np.sqrt(gamma*R*T_inlet) # maybe T_inlet is wrong, but it certainly will result in a worst-case M
        M = v1_5 / a # v1.5 is probably wrong

        # degree of reaction (called Λ, but python and I hate non-roman variable names)
        # tan(alpha2) - tan(alpha1) = tan(beta1) - tan(beta2)
        self.DRXN = 1 - (norm(v_inlet) / (2*v_blade))*(np.tan(rot_defl_ang) - np.tan(stat_defl_ang))

        # de haller number
        self.DHN = norm(self.rotor.v_outlet) / norm(self.rotor.rel_v_inlet)

        # recommended outlet area (ratio)
        rho_1 = self.p_inlet  / (R*self.T_inlet)
        rho_2 = self.p_outlet / (R*self.T_outlet)
        self.A_ratio = (rho_1 * norm(v_inlet)) / (rho_2 * norm(self.v2))


        if (s_inlet < 0.0):
            self.s_inlet = cp*np.log(self.T_inlet/273.16) - R*np.log(self.p_inlet/101325.0) + 6608.1
        else:
            self.s_inlet = s_inlet
        self.s_outlet =  cp*np.log(self.T_outlet/self.T_inlet) - R*np.log(self.p_outlet/self.p_inlet) + self.s_inlet        

    def plot_triangles(self):
        """Plots the two velocity triangles for the stage"""
        self.rotor.plot(title=f'Velocity triangle: {self._inst_name()} rotor')
        self.stator.plot(title=f'Velocity triangle: {self._inst_name()} stator')

    def _inst_name(self):
        # for var_name, var_val in globals().items():
        #     if var_val is self:
        #         return str(var_name)
        return f'stage {self.instance_number}'


    def print_stats(self, verbose=True):
        """ prints out important data for the current stage
            returns a string containing all the stats"""
        status_text = [f'stats for {self._inst_name()}:',
        f'----------------------------',
        f'De Haller number: {self.DHN:.2f}',
        f'Degree of reaction: {self.DRXN:.2f}',
        f'flow coefficient: {self.phi:.2f}',
        f'stage work: {self.w:.1f} J',
        f'stage pressure: {self.p_inlet:.1f} Pa -> {self.p_outlet:.1f} Pa (Ratio of {self.p_outlet/self.p_inlet:.3f})',
        f'temperature (inlet->outlet): {self.T_inlet:.2f} K -> {self.T_outlet:.2f} K  (Ratio of {self.T_outlet/self.T_inlet:.3f})',
        f'entropy (inlet->outlet): {self.s_inlet:.2f}->{self.s_outlet:.2f} Generated {self.s_outlet - self.s_inlet:.2f} J/kg K of entropy',
        f'area ratio: (out/in): {self.A_ratio:.4f}']
        if verbose:
            print('\n'.join(status_text), end='\n\n')
        return status_text

    def plot_mollier(self, verbose=True):
        """ Plots a mollier diagram of the stage's compression process

            verbose:
                If verbose, then the plots are shown and then discarded
                if not verbose, then plots will be written on the global plt object, allowing other methods to annotate after this method."""
        # use this to plot a mollier diagram of the compresion process

        # need to get s as a function of temperature or enthalpy
        # need to do this for isentropic and polytropic processes
        p_polytropic_ratio = lambda T: (T/self.T_inlet)**((gamma/(gamma-1))*(0.9))
        p_isentropic_ratio = lambda T: (T/self.T_inlet)**(gamma/(gamma-1))

        s_isentrope = lambda T: cp*np.log(T/self.T_inlet) - R*np.log(p_isentropic_ratio(T)) + self.s_inlet
        s_polytrope = lambda T: cp*np.log(T/self.T_inlet) - R*np.log(p_polytropic_ratio(T)) + self.s_inlet

        # y axis points (usually temperature or enthalpy)
        T_points = np.linspace(self.T_inlet, self.T_outlet, num=50)

        # plot, label, show
        plt.plot(s_isentrope(T_points), T_points); plt.plot(s_polytrope(T_points), T_points)
        plt.ylabel(f'Temperature (K)'); plt.xlabel(f'Entropy (J/Kg K)'); plt.title(f'Mollier diagram for compression process')
        if verbose:
            plt.legend([f'isentropic', f'polytropic process with n={self.poly_n}'])
            plt.show()


class Annulus:
    # implement later??
    def __init__(self):
        # get radii (minimum and maximum

        # calculate streamtube inlet area

        # calculate 
        pass

class TrueCompressor:
    # incorporate the blading and 3d effects
    def __init__(self):
        # define inlet specifications
        inlet_area = np.pi*(0.2**2 - 0.1**2)
        inlet_hub_radius = 0.1
        inlet_tip_radius = 0.2

        compressor_rpm = 16000

        # lists to hold the radii of all the stages
        hub_radii = [inlet_hub_radius]
        tip_radii = [inlet_tip_radius]

        # define annuli (streamtubes)
        strealine_radii = np.array([0.15]) # change to an arbitrarily large np.linspace() later :)
        for radius in strealine_radii:
            pass


class PrelimCompressor:
    # a compressor designed with only 1D calculations in mind
    # useful to get the preliminary area calculations!!


    def __init__(self,
                 inlet_hub_r : float,
                    inlet_tip_r : float,
                        comp_rpm : float,
                            comp_T_inlet : float,
                                comp_p_inlet : float,
                                    comp_v_inlet : np.ndarray,
                                        rotor_defl_angles : np.ndarray,
                                            stator_defl_angles : np.ndarray):

        # radius ranges from 0.1 m to 0.2 m. keep the external radius constant and increase internal radius as indicated by Stage.A2
        
        # first stage dimensions
        inlet_area = np.pi*(inlet_tip_r**2 - inlet_hub_r**2)
        hub_radii = [inlet_hub_r]
        tip_radii = [inlet_tip_r]

        # do constant outer diameter
        hub_radius_from_area = lambda A: np.sqrt((-A+np.pi*tip_radii[0]**2)/np.pi) # assumes constant exterior
        
        assert len(rotor_defl_angles) == len(stator_defl_angles)

        stages = []

        stages.append(Stage_1D(v_inlet=comp_v_inlet,
                rpm=comp_rpm,
                    r=0.5*(hub_radii[-1]+tip_radii[-1]),
                    T_inlet=comp_T_inlet,
                    p_inlet=comp_p_inlet,
                        rot_defl_ang=rotor_defl_angles[0],
                        stat_defl_ang=stator_defl_angles[0]))

        # remove 1st entry (scuffed, ik)
        rotor_defl_angles = rotor_defl_angles[1:]
        stator_defl_angles = stator_defl_angles[1:]
        # area calcs for second stage
        hub_radii.append(hub_radius_from_area(inlet_area*stages[-1].A_ratio))
        tip_radii.append(tip_radii[0])

        # generate rest of the stages
        for rotor_defl_ang, stator_defl_ang in zip(rotor_defl_angles, stator_defl_angles):
            stages.append(
                Stage_1D(v_inlet=comp_v_inlet,
                    rpm=comp_rpm,
                        r=0.5*(hub_radii[-1]+tip_radii[-1]),
                            T_inlet=stages[-1].T_outlet,
                                p_inlet=stages[-1].p_outlet,
                                    rot_defl_ang=rotor_defl_ang,
                                        stat_defl_ang=stator_defl_ang,
                                            s_inlet=stages[-1].s_outlet)
                                            )

            # area calcs for second stage
            hub_radii.append(hub_radius_from_area(np.pi*(tip_radii[-1]**2 - hub_radii[-1]**2)*stages[-1].A_ratio))
            tip_radii.append(tip_radii[0])

            self.stages = stages

        def print_stats(self):
            # prinout stats and triangles
            for stage in stages:
                stage.print_stats()
                
                print(f'overall compressor stats:\n'
                f'total work: {sum([s.w for s in stages]):.2f} J \n'
                f'pressure ratio: {stages[-1].p_outlet / stages[0].p_inlet} ({stages[0].p_inlet/1000:.2f} kPa --> {stages[-1].p_outlet/1000:.2f})\n'
                )


        def print_triangles(self):
            for stage in stages:
                stage.plot_triangles()
        
        def print_mollier_triangles(self):
            # printout mollier diagrams
            for stage in stages:
                stage.plot_mollier(verbose=False)
            plt.show()

        def print_illustrations(self):
            # 1D stage illustrations
            fig, ax = plt.subplots()
            for i in range(len(hub_radii)):
                ax.add_patch(patches.Rectangle((i/20, hub_radii[i]), 1/40, (tip_radii[i]-hub_radii[i])))
                ax.add_patch(patches.Rectangle((i/20, -hub_radii[i]-tip_radii[i]), 1/40, tip_radii[i]))
            m = max(tip_radii + [(0.05*(len(hub_radii)+1))])
            ax.set_xlim((0, 2*m))
            ax.set_ylim((-m, m))
            plt.title(f'1D stages (in meters)')
            plt.show()
            plt.clf()





def main():
    c = PrelimCompressor(
        inlet_hub_r = 0.1,
        inlet_tip_r = 0.2,
        comp_rpm = 16000.0,
        comp_T_inlet = 288.0,
        comp_p_inlet = 101325.0,
        comp_v_inlet = np.array([0.0, 150.0]),
        rotor_defl_angles = np.deg2rad(np.array([12, 12])),
        stator_defl_angles = np.deg2rad(np.array([35, 12])))


# program entry point
if __name__ == "__main__":
    main()


