'''
Keeps track of particle trajectory with considerations to cutoffs and E-W effects etc
'''

import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.utils import EARTH_RADIUS, g10, DEG_TO_RAD
from gtracr.runge_kutta import runge_kutta
from gtracr.add_particles import particleDict


def VertRigidityCutoff(r, theta):
    return (14.5*np.sin(theta)**4.) / r**2.

class ParticleTrajectory:
    '''
    Traces trajectory of a single particle, given the following information:
        - particleName: the name of the particle of interest, obtained from Particle class
        - initialAltitude : the starting altitude of the particle
        - finalAltitude : the altitude in which the particle trajectory ends
        - maxStep : the maximum number of steps to integrate for (default=100000)

    '''
    def __init__(self, particleName, startAltitude, stopAltitude=100, maxStep=100000, escapeRadius=10):
        self.particle = particleDict[particleName]  # should be obtained from some dictionary
        self.startAltitude = startAltitude
        self.stopAltitude = stopAltitude
        self.escaperadius = escapeRadius
        self.maxStep = maxStep
        self.stepSize = np.abs(stopAltitude - startAltitude) / maxStep
        self.coord = dict.fromkeys(["t", "r", "theta", "phi"], np.zeros(maxStep))

    # get the particle trajectory given:
    # - energy: the energy of the cosmic ray as measured from detector
    # - startLongitude: the magnetic longitude
    # - startLatitude: the magnetic latitude
    # - startZenith: the angle from the local zenith
    # - startAzimuth: the angle from the local North
    # default all to zero
    def getTrajectory(self, energy, startLongitude=0., startLatitude=0., startZenith=0., startAzimuth=0.):
        startTraj = TrajectoryPoint(startLongitude, startLatitude, self.startAltitude, startZenith, startAzimuth)
        # get initial position in spherical coordinates
        (r0,theta0,phi0) = startTraj.sphericalCoord()
        # set particle momenta and get their spherical coordinate equivalents
        p0 = self.particle.momentum(energy)
        pr0 = p0
        pth0 = p0 / r0
        pph0 = p0 / (r0*np.sin(theta0))
        # initial time
        t0 = 0.
        # put all information into some initial value array
        ival = (t0, r0, theta0, phi0, pr0, pth0, pph0)

        # initialize array
        t_arr = np.zeros(self.maxStep)
        r_arr = np.zeros(self.maxStep)
        th_arr = np.zeros(self.maxStep)
        ph_arr = np.zeros(self.maxStep)
        pr_arr = np.zeros(self.maxStep)
        pth_arr = np.zeros(self.maxStep)
        pph_arr = np.zeros(self.maxStep)

        
        # perform runge_kutta
        # this should be in some while loop
        i = 0
        while i < self.maxStep:
        # while i < self.maxStep: # (or something like this)
            result = runge_kutta(self.particle, self.stepSize, ival)
            t = result[0]
            r = result[1]
            theta = result[2]
            phi = result[3]
            pr = result[4]
            pth = result[5]
            pph = result[6]
            # set this to particle's new position and momentum 
            # self.particle.momentum = np.sqrt(pr**2. + (r*pth)**2. + (r*np.sin(theta)*pph)**2.)

            # append
            t_arr[i] = t
            r_arr[i] = r
            th_arr[i] = theta
            ph_arr[i] = phi
            pr_arr[i] = pr
            pth_arr[i] = pth
            pph_arr[i] = pph

            # check for some conditions over here
            if r > self.escaperadius:
                break
            elif r <= 0.:
                break

            # update ival
            ival = (t, r, theta, phi, pr, pth, pph)

            i += 1

        # # add to dictionary
        # self.coord["t"] = t_arr
        # self.coord["r"] = r_arr
        # self.coord["theta"] = th_arr
        # self.coord["phi"] = ph_arr

        return (t_arr, r_arr, th_arr, ph_arr)
            

class TrajectoryPoint:
    '''
    Constructs a single point of the trajectory given the following information:
        - longitude: the magnetic longitude (0 = prime meridian)
        - latitude: the magnetic latitude (0 = magnetic equator)
        - altitude: the distance from Earth's surface [km]
        - zenithAngle: the angle from the local zenith
        - azimuthAngle: the angle from the local North
    '''
    def __init__(self, longitude, latitude, altitude, zenithAngle, azimuthAngle):
        self.longitude = longitude
        self.latitude = latitude
        self.altitude = altitude
        self.zenithAngle = zenithAngle
        self.azimuthAngle = azimuthAngle

    # converts into spherical coordinates
    # r [EARTH_RADIUS], theta [rad], phi [rad]
    def sphericalCoord(self):
        r = (self.altitude + EARTH_RADIUS) / EARTH_RADIUS
        theta = (90. - self.latitude)*DEG_TO_RAD
        phi = (180. + self.longitude)*DEG_TO_RAD
        return (r,theta,phi)
