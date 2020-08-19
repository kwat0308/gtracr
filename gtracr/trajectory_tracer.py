'''
Class that traces the trajectory of the particle

This is required to allow the pyIGRF model to be used within our code.
This also adds flexibility to our code incase we want to switch back to Python at some point

This is directly taken from the C++ version, but vectorized to some extent
with the aid of numpy arrays

'''
import os, sys
import numpy as np
from datetime import datetime as dt

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(CURRENT_DIR)

from gtracr.trajectorypoint import TrajectoryPoint
from gtracr.constants import EARTH_RADIUS, SPEED_OF_LIGHT, ELEMENTARY_CHARGE, KG_PER_GEVC2
from gtracr.magnetic_field import MagneticField, IGRF13
from scipy import interpolate


class pTrajectoryTracer:
    '''
    A class that traces the trajectory of the particle. 
    
    Members
    --------

    - charge: the particle's charge in electrons
    - mass: the particle's mass in GeV
    - escape_radius: the radius in which the particle effectively escaped (default at 10RE). Defined relative to sea level (at Earth's surface).
    - stepsize: the stepsize of the integrator (default:1e-5)
    - max_time: the maximum time for the trajectory tracing procedure (default=10s)
    - bfield_type: the type of magnetic field to use for trajectory tracing, either "dipole" or "igrf" (default="dipole")
    - max_step: the maximum number of steps, evaluated from max_time and stepsize
    - particle_escaped: boolean to check if particle has escaped or not
    - bfield: the magnetic field model that will be used

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
                 max_time=10,
                 max_step=None,
                 bfield_type="dipole"):
        self.charge = charge * ELEMENTARY_CHARGE  # convert to coulombs
        self.mass = mass * KG_PER_GEVC2  # convert to kg
        self.escape_radius = escape_radius
        self.stepsize = stepsize
        self.max_time = max_time  # default 10s
        self.max_step = int(max_time / stepsize) \
            if max_step is None else max_step  # N = (tf - t0) / h
        self.particle_escaped = False  # check if particle escaped or not

        # initialize magnetic field
        if bfield_type.find("dipole") != -1:
            self.bfield = MagneticField()
        elif bfield_type.find("igrf") != -1:
            curr_year = dt.now().year
            nmax = 13  # should be able to vary in future versions
            self.bfield = IGRF13(curr_year, nmax=nmax)
        else:
            raise Exception("Only modes 'dipole' and 'igrf' are allowed!")

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

        #

        # evaluate B-field
        # bf_r = self.bfield.Br(r, theta, phi)
        # bf_theta = self.bfield.Btheta(r, theta, phi)
        # bf_phi = self.bfield.Bphi(r, theta, phi)
        bf_r, bf_theta, bf_phi = self.bfield.values(r, theta, phi)

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

    def evaluate(self, t0, vec0, get_data=False):
        '''
        Evaluate the trajectory by performing a 4th order Runge Kutta integration.
        
        Parameters
        ----------

        - t0: the initial time
        - vec0 : the initial six-vector (r0, theta0, phi0, pr0, pphi0, ptheta0)
        - get_data: whether to extract the trajectory arrays as a dictionary or not (default False)
        '''
        # set initial conditions
        t = t0
        vec = vec0

        # print(t, vec)
        # create empty arrays
        # if get_data is false it wont be appended so its fine to do this
        # without conditional case
        t_arr = np.zeros(self.max_step)
        vec_arr = np.zeros((self.max_step, 6))
        # start the loop
        for i in range(self.max_step):

            # print(vec)

            # append only if get_data is true
            if get_data:
                # print(get_data)
                # print(t, vec)
                # np.append(t_arr, t)
                # np.append(vec_arr, vec.reshape((1, 6)), axis=0)
                # np.append(t_arr, t)
                # np.append(vec_arr, vec.reshape((1, 6)))
                t_arr[i] = t
                vec_arr[i] = vec

            h = self.stepsize  # for more compact writing
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

        # define the last point as a trajectory point
        final_tp = TrajectoryPoint(*vec)

        # print(t_arr, vec_arr)

        # return the arrays regardless of the conditions
        # the upper interface will deal with the cases anyways
        return np.array(t_arr), np.array(vec_arr), final_tp