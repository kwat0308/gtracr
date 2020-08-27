from scipy import interpolate
from gtracr.lib.constants import EARTH_RADIUS, SPEED_OF_LIGHT, ELEMENTARY_CHARGE, KG_PER_GEVC2
from gtracr.lib.trajectorypoint import TrajectoryPoint
from gtracr.lib.magnetic_field import MagneticField, IGRF13
'''
Class that traces the trajectory of the particle

This is required to allow the pyIGRF model to be used within our code.
This also adds flexibility to our code incase we want to switch back to Python at some point

This is directly taken from the C++ version, but vectorized to some extent
with the aid of numpy arrays

'''
import os
import sys
import numpy as np
from datetime import datetime as dt

# CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# sys.path.append(CURRENT_DIR)


class pTrajectoryTracer:
    '''
    A class that traces the trajectory of the particle. 

    Members
    --------

    - charge: the particle's charge in electrons
    - mass: the particle's mass in GeV
    - escape_radius: the radius in which the particle effectively escaped (default at 10RE). Defined relative to Earth's center.
    - stepsize: the stepsize of the integrator (default:1e-5)
    - max_time: the maximum time for the trajectory tracing procedure (default=10s)
    - bfield_type: the type of magnetic field to use for trajectory tracing, either "d" for dipole or "i" for igrf (default="d")
    - max_step: the maximum number of steps, evaluated from max_time and stepsize
    - particle_escaped: boolean to check if particle has escaped or not
    - bfield: the magnetic field model that will be used
    - igrf_params : the path to the data directory and the current date used in the IGRF model

    Note:
    The "p" in front of pTrajectoryTracer
    indicates that this is the Python version. Such naming is required since we want to distinguish between
    the TrajectoryTracer class in Python vs C++, and the Python version is mainly used as a 
    tester for the C++ version.
    '''

    def __init__(self,
                 charge,
                 mass,
                 escape_radius=10. * EARTH_RADIUS,
                 stepsize=1e-5,
                 max_step=10000,
                 bfield_type="d",
                 igrf_params=None):
        self.charge = charge * ELEMENTARY_CHARGE  # convert to coulombs
        self.mass = mass * KG_PER_GEVC2  # convert to kg
        self.escape_radius = escape_radius
        self.stepsize = stepsize
        self.max_step = max_step
        self.particle_escaped = False  # check if particle escaped or not

        # initialize magnetic field
        if bfield_type.find("d") != -1:
            self.bfield = MagneticField()
        elif bfield_type.find("i") != -1:
            curr_year = igrf_params[1]  # the current date
            nmax = 13  # should be able to vary in future versions
            self.bfield = IGRF13(curr_year, nmax=nmax)
        else:
            raise Exception("Only modes 'dipole' and 'igrf' are allowed!")

        # the final coordinates
        self.final_time = 0.
        self.final_sixvector = np.zeros(6)

    def ode_lrz(self, t, vec):
        '''
        The system of ordinary differential equations that describe the motion of charged
        particles in Earth's magnetic field via the Lorentz force in spherical coordinates.
        This version differs from the C++ code in that we perform this in a vectorized fashion.

        Parameters
        -----------

        - t : the time
        - vec: the six-vector (r, theta, phi, pr, ptheta, pphi) at time t

        Returns
        --------
        - ode_lrz : the ordinary differential equation for the six vector based on the Lorentz force equation
        '''
        # unpack vector for readability
        (r, theta, phi, pr, ptheta, pphi) = vec

        # get lorentz factor
        # pmag = np.linalg.norm(vec[3:])  # momentum magnitude
        pmag = np.sqrt(pr**2. + ptheta**2. + pphi**2.)
        gamma = np.sqrt(1. + (pmag / (self.mass * SPEED_OF_LIGHT))**2.)
        rel_mass = self.mass * gamma

        # evaluate B-field
        # bf_r = self.bfield.Br(r, theta, phi)
        # bf_theta = self.bfield.Btheta(r, theta, phi)
        # bf_phi = self.bfield.Bphi(r, theta, phi)
        bf_r, bf_theta, bf_phi = self.bfield.values(r, theta, phi)

        # print(bf_r, bf_theta, bf_phi)

        # define the momentum odes first
        # invert charge for back tracking
        # _lrz indicates the lorentz force term and
        # _sphcmp indicates the auxiliary terms that
        # come from using spherical coordinates
        # dprdt
        dprdt_lrz = -1. * self.charge * (ptheta * bf_phi - bf_theta * pphi)
        dprdt_sphcmp = ((ptheta**2. + pphi**2.) / r)
        dprdt = dprdt_lrz + dprdt_sphcmp

        # dpthetadt
        dpthetadt_lrz = self.charge * (pr * bf_phi - bf_r * pphi)
        dpthetadt_sphcmp = ((pphi**2. * np.cos(theta)) /
                            (r * np.sin(theta))) - ((pr * ptheta) / r)
        dpthetadt = dpthetadt_lrz + dpthetadt_sphcmp

        # dpphidt
        dpphidt_lrz = -1. * self.charge * (pr * bf_theta - bf_r * ptheta)
        dpphidt_sphcmp = ((pr * pphi) / r) + ((ptheta * pphi * np.cos(theta)) /
                                              (r * np.sin(theta)))
        dpphidt = dpphidt_lrz - dpphidt_sphcmp

        # define the position odes while creating the vector for the ODE
        # drdt = pr, dthetadt = ptheta/r, dphidt = pphi/r*sin(theta)
        # divide each ODE by relativistic mass to account for using momentum
        # instead of velocity
        ode_lrz = np.array([
            pr, (ptheta / r),
            (pphi / (r * np.sin(theta))), dprdt, dpthetadt, dpphidt
        ]) / rel_mass

        return ode_lrz

    def evaluate(self, t0, vec0):
        '''
        Evaluate the trajectory by performing a 4th order Runge Kutta integration.

        Parameters
        ----------

        - t0: the initial time
        - vec0 : the initial six-vector (r0, theta0, phi0, pr0, pphi0, ptheta0)
        - get_data: whether to extract the trajectory arrays as a dictionary or not (default False)

        Returns
        --------

        - None
        '''
        # set initial conditions
        t = t0
        vec = vec0
        h = self.stepsize  # for more compact writing
        # start the loop
        for i in range(self.max_step):
            # evaluate k-coefficients
            k1_vec = h * self.ode_lrz(t, vec)
            k2_vec = h * self.ode_lrz(t + (0.5 * h), vec + (0.5 * k1_vec))
            k3_vec = h * self.ode_lrz(t + (0.5 * h), vec + (0.5 * k2_vec))
            k4_vec = h * self.ode_lrz(t + h, vec + k3_vec)

            # increment by weighted sum of k-coeffs
            vec += (1. / 6.) * (k1_vec + (2. * k2_vec) +
                                (2. * k3_vec) + k4_vec)
            t += h  # increment time

            # print(vec, t)

            # breaking conditions based on value of r
            r = vec[0]  # for readability
            # if particle has reached escape radius or not
            # then trajectory is allowed
            if r > EARTH_RADIUS + self.escape_radius:
                self.particle_escaped = True
                break
            # if particle has came back to earth
            # then trajectory is forbidden
            if r < EARTH_RADIUS:
                break
        
        # set the final time / vector
        self.final_time = t
        self.final_sixvector = vec

        # return None

    def evaluate_and_get_trajectory(self, t0, vec0):
        '''
        Evaluate the trajectory by performing a 4th order Runge Kutta integration and get the corresponding trajectory data, that is, the duration of the trajectory and the six-vector of the trajectory as a numpy array. 

        Parameters
        ----------

        - t0: the initial time
        - vec0 : the initial six-vector (r0, theta0, phi0, pr0, pphi0, ptheta0)
        - get_data: whether to extract the trajectory arrays as a dictionary or not (default False)

        Returns
        --------

        - trajectory_data (dict<str, numpy array>) : 
            The dictionary that contains the trajectory information as well as the final six-vector of the trajectory. 
            - keys are ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]
        '''
        # set initial conditions
        t = t0
        vec = vec0
        h = self.stepsize  # for more compact writing
        # print(t, vec)
        t_arr = []
        # vec_arr = []
        r_arr = []
        theta_arr = []
        phi_arr = []
        pr_arr = []
        ptheta_arr = []
        pphi_arr = []
        # start the loop
        for i in range(self.max_step):

            # print(vec)

            # append data
            t_arr.append(t)
            r_arr.append(vec[0])
            theta_arr.append(vec[1])
            phi_arr.append(vec[2])
            pr_arr.append(vec[3])
            ptheta_arr.append(vec[4])
            pphi_arr.append(vec[5])

            # evaluate k-coefficients
            k1_vec = h * self.ode_lrz(t, vec)
            k2_vec = h * self.ode_lrz(t + (0.5 * h), vec + (0.5 * k1_vec))
            k3_vec = h * self.ode_lrz(t + (0.5 * h), vec + (0.5 * k2_vec))
            k4_vec = h * self.ode_lrz(t + h, vec + k3_vec)

            # increment by weighted sum of k-coeffs
            vec += (1. / 6.) * (k1_vec + (2. * k2_vec) +
                                (2. * k3_vec) + k4_vec)
            t += h  # increment time

            # print(vec, t)

            # breaking conditions based on value of r
            r = vec[0]  # for readability
            # if particle has reached escape radius or not
            # then trajectory is allowed
            if r > EARTH_RADIUS + self.escape_radius:
                self.particle_escaped = True
                break
            # if particle has came back to earth
            # then trajectory is forbidden
            if r < EARTH_RADIUS:
                break

        # set the final time / vector
        self.final_time = t
        self.final_sixvector = vec

        # create the dictionary that contains the trajectory data
        trajectory_data = {
            "t": t_arr,
            "r": r_arr,
            "theta": theta_arr,
            "phi": phi_arr,
            "pr": pr_arr,
            "ptheta": ptheta_arr,
            "pphi": pphi_arr
        }

        # print(t_arr, vec_arr)

        return trajectory_data
