'''
Keeps track of particle trajectory with considerations to cutoffs and E-W effects etc
'''

import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

# from gtracr.utils import EARTH_RADIUS, g10, DEG_TO_RAD, get_sphcomp_momentum
from gtracr.utils import *
from gtracr.runge_kutta import runge_kutta
from gtracr.add_particles import particleDict

# vertical rigidity cutoff defined in Baldini's paper
def VertRigidityCutoff(r, theta):
    return (14.5 * np.sin(theta)**4.) / r**2.

# def get_initial_values(trajectory, energy):
#     # get initial position in spherical coordinates
#     (r0, theta0, phi0) = trajectory.sphericalCoord()
#     # set particle momenta and get their spherical coordinate equivalents
#     p0 = self.particle.setMomentum(energy)
#     (pr0, pth0, pph0) = get_sphcomp_momentum(p0, r0, theta0)
#     # initial time
#     t0 = 0.
#     # put all information into some initial value array
#     ival = np.array(t0, r0, theta0, phi0, pr0, pth0, pph0)
#     return ival

key_list = ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]


class ParticleTrajectory:
    '''
    Traces trajectory of a single particle, given the following information:
        - particleName: the name of the particle of interest, obtained from Particle class
        - initialAltitude : the starting altitude of the particle
        - finalAltitude : the altitude in which the particle trajectory ends
        - maxStep : the maximum number of steps to integrate for (default=100000)

    '''
    def __init__(self,
                 particleName,
                 startAltitude=565,
                 stopAltitude=2,
                 maxStep=100000,
                 escapeRadius=2):
        self.particle = particleDict[
            particleName]  # should be obtained from some dictionary
        self.startAltitude = startAltitude
        self.stopAltitude = stopAltitude
        self.escaperadius = escapeRadius
        self.maxStep = maxStep
        self.stepSize = np.abs(stopAltitude - startAltitude) / maxStep
        self.results = {key:np.zeros(maxStep) for key in key_list}
    
    

    # get the particle trajectory given:
    # - energy: the energy of the cosmic ray as measured from detector
    # - startLongitude: the magnetic longitude
    # - startLatitude: the magnetic latitude
    # - startZenith: the angle from the local zenith
    # - startAzimuth: the angle from the local North
    # default all to zero
    def getTrajectory(self,
                      energy,
                      startLongitude=0.,
                      startLatitude=0.,
                      startZenith=0.,
                      startAzimuth=0.):
        startTraj = TrajectoryPoint(startLongitude, startLatitude,
                                    self.startAltitude, startZenith,
                                    startAzimuth)

        # get initial position in spherical coordinates
        (r0, theta0, phi0) = startTraj.sphericalCoord()
        # set particle momenta and get their spherical coordinate equivalents
        self.particle.setMomentum(energy)
        self.particle.setRigidity(energy)
        (pr0, pth0, pph0) = get_sphcomp_momentum(self.particle.momentum, r0, theta0)
        # initial time
        t0 = 0.
        # put all information into some initial value array
        initial_value = np.array([t0, r0, theta0, phi0, pr0, pth0, pph0])

        print(initial_value)

        # initial_value = get_initial_values(startTraj, energy)
        # initialize array
        # t_arr = np.zeros(self.maxStep)
        # r_arr = np.zeros(self.maxStep)
        # th_arr = np.zeros(self.maxStep)
        # ph_arr = np.zeros(self.maxStep)
        # pr_arr = np.zeros(self.maxStep)
        # pth_arr = np.zeros(self.maxStep)
        # pph_arr = np.zeros(self.maxStep)

        # perform runge_kutta
        # this should be in some while loop
        i = 0
        while i < self.maxStep:
            # print(initial_value)
            # while i < self.maxStep: # (or something like this)
            (t, r, theta, phi, pr, ptheta,
             pphi) = runge_kutta(self.particle, self.stepSize, initial_value)

            # print(t, r, theta, phi, pr, ptheta,pphi)
            # t = result[0]
            # r = result[1]
            # theta = result[2]
            # phi = result[3]
            # pr = result[4]
            # ptheta = result[5]
            # pphi = result[6]
            # set this to particle's new position and momentum
            self.particle.momentum = get_momentum_magnitude(pr, ptheta, pphi, r, theta)
            
            # append
            # t_arr[i] = t
            # r_arr[i] = r
            # th_arr[i] = theta
            # ph_arr[i] = phi
            # pr_arr[i] = pr
            # pth_arr[i] = pth
            # pph_arr[i] = pph
            self.results["t"][i] = t
            self.results["r"][i] = r
            self.results["theta"][i] = theta
            self.results["phi"][i] = phi
            self.results["pr"][i] = pr
            self.results["ptheta"][i] = ptheta
            self.results["pphi"][i] = pphi


            # check for some conditions over here
            # if r <= self.escaperadius:
            #     break
            # if r <= 1. and i > 0:  # r*RE <= RE
            #     break
            
            # if self.particle.rigidity < VertRigidityCutoff(r, theta):
            #     break

            # update initial_value
            initial_value = np.array([t, r, theta, phi, pr, ptheta, pphi])

            i += 1

        # # add to dictionary
        # self.coord["t"] = t_arr
        # self.coord["r"] = r_arr
        # self.coord["theta"] = th_arr
        # self.coord["phi"] = ph_arr

        # return (t_arr, r_arr, th_arr, ph_arr)
        print(self.results)

        


class TrajectoryPoint:
    '''
    Constructs a single point of the trajectory given the following information:
        - longitude: the magnetic longitude (0 = prime meridian)
        - latitude: the magnetic latitude (0 = magnetic equator)
        - altitude: the distance from Earth's surface [km]
        - zenithAngle: the angle from the local zenith
        - azimuthAngle: the angle from the local North
    '''
    def __init__(self, longitude, latitude, altitude, zenithAngle,
                 azimuthAngle):
        self.longitude = longitude
        self.latitude = latitude
        self.altitude = altitude
        self.zenithAngle = zenithAngle
        self.azimuthAngle = azimuthAngle

    # converts into spherical coordinates
    # r [EARTH_RADIUS], theta [rad], phi [rad]
    def sphericalCoord(self):
        r = (self.altitude + EARTH_RADIUS) / EARTH_RADIUS
        theta = (90. - self.latitude) * DEG_TO_RAD
        phi = (180. + self.longitude) * DEG_TO_RAD
        return np.array([r, theta, phi])
