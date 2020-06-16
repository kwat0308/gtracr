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
from gtracr.add_particle import particleDict


# vertical rigidity cutoff defined in Baldini's paper
def VertRigidityCutoff(r, theta):
    return (14.5 * np.sin(theta)**4.) / r**2.


key_list = ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]


class ParticleTrajectory:
    '''
    Traces trajectory of a single particle, given the following information:
    - particleName: the name of the particle of interest, obtained from Particle class
    - energy : the energy of particle / cosmic ray
    - startLatitude : initial latitude of location in decimal notation
    - startLongitude : initial longitude of location in decimal notation
    - startAltitude : the starting altitude of the particle [km]
    - stopAltitude : the altitude in which the particle trajectory ends [km]
    - maxStep : the maximum number of steps to integrate for (default=10000)

    '''
    def __init__(self,
                 particleName,
                 energy,
                 startLatitude=0.,
                 startLongitude=0.,
                 startAltitude=1.,
                 stopAltitude=500.,
                 maxStep=10000):
        self.particle = particleDict[
            particleName]  # should be obtained from some dictionary
        self.energy = energy
        self.startLatitude = startLatitude
        self.startLongitude = startLongitude
        self.startAltitude = startAltitude
        self.stopAltitude = stopAltitude
        self.maxStep = maxStep
        self.stepSize = np.abs(stopAltitude - startAltitude) / maxStep
        self.results = {key: np.zeros(maxStep) for key in key_list}

    # def __init__(self,
    #              particleName,
    #              startAltitude=1,
    #              stopAltitude=500,
    #              maxStep=10000):
    #     self.particle = particleDict[
    #         particleName]  # should be obtained from some dictionary
    #     self.startAltitude = startAltitude
    #     self.stopAltitude = stopAltitude
    #     self.maxStep = maxStep
    #     self.stepSize = np.abs(stopAltitude - startAltitude) / maxStep
    #     self.results = {key:np.zeros(maxStep) for key in key_list}

    # get the trajectory from a provided zenith and azimuthal angle

    def getTrajectory(self, zenith=0., azimuth=0.):
        # create trajectory point for initial point
        startTraj = TrajectoryPoint(self.startLatitude, self.startLongitude,
                                    self.startAltitude, zenith, azimuth)

        # get initial value
        initial_value = self.get_initial_values(startTraj, self.energy)
        print(initial_value)

        i = 0
        while i < self.maxStep:
            (t, r, theta, phi, vr, vtheta,
             vphi) = runge_kutta(self.particle, self.stepSize, initial_value)

             # check for some conditions over here
            # if r <= EARTH_RADIUS:  # particle already on / within surface of earth
            #     break
            
            # convert velocities to momenta
            vmag = vmag_spherical(vr, vtheta, vphi, r, theta)
            # print(vmag)
            # pr = gamma(vmag)*self.particle.mass*vr
            # ptheta = gamma(vmag)*self.particle.mass*vtheta
            # pphi = gamma(vmag)*self.particle.mass*vphi
            # print(t, r, theta, phi, pr, ptheta,pphi)
            # set this to particle's new position and momentum
            self.particle.momentum = gamma(vmag) * self.particle.mass * vmag
            # self.particle.momentum = vmag_spherical(pr, ptheta, pphi, r, theta)
            self.particle.set_rigidity_from_momentum()

            (pr, ptheta,
             pphi) = vp_components_spherical(self.particle.momentum, r, theta)

            # append
            self.results["t"][i] = t
            self.results["r"][i] = r
            self.results["theta"][i] = theta
            self.results["phi"][i] = phi
            self.results["pr"][i] = pr
            self.results["ptheta"][i] = ptheta
            self.results["pphi"][i] = pphi

             # check for some conditions over here
            if r <= EARTH_RADIUS:  # particle already on / within surface of earth
                break

            # update initial_value
            # initial_value = np.array([t, r, theta, phi, pr, ptheta, pphi])
            initial_value = np.array([t, r, theta, phi, vr, vtheta, vphi])

            i += 1

        # get the end trajectory from the last spherical coordinate values
        endTraj = TrajectoryPoint()
        endTraj.set_from_sphericalCoord(r, theta, phi)

        print(self.results)
        # print(len(self.results["r"]))

        return (startTraj, endTraj)

    # # get the particle trajectory given:
    # # - energy: the energy of the cosmic ray as measured from detector
    # # - startLongitude: the magnetic longitude
    # # - startLatitude: the magnetic latitude
    # # - startZenith: the angle from the local zenith
    # # - startAzimuth: the angle from the local North
    # # default all to zero
    # # return tuple of start and end trajectories
    # def getTrajectory(self,
    #                   energy,
    #                   startLongitude=0.,
    #                   startLatitude=0.,
    #                   startZenith=0.,
    #                   startAzimuth=0.):
    #     startTraj = TrajectoryPoint(startLongitude, startLatitude,
    #                                 self.startAltitude, startZenith,
    #                                 startAzimuth)

    #     # set rigidity
    #     self.particle.set_rigidity_from_energy(energy)
    #     # get initial position in spherical coordinates
    #     initial_value = self.get_initial_values(startTraj, energy)
    #     # perform runge_kutta
    #     i = 0
    #     while i < self.maxStep:
    #         (t, r, theta, phi, pr, ptheta,
    #          pphi) = runge_kutta(self.particle, self.stepSize, initial_value)

    #         # print(t, r, theta, phi, pr, ptheta,pphi)
    #         # set this to particle's new position and momentum
    #         self.particle.momentum = get_momentum_magnitude(pr, ptheta, pphi, r, theta)
    #         self.particle.set_rigidity_from_momentum()
    #         # append
    #         self.results["t"][i] = t
    #         self.results["r"][i] = r
    #         self.results["theta"][i] = theta
    #         self.results["phi"][i] = phi
    #         self.results["pr"][i] = pr
    #         self.results["ptheta"][i] = ptheta
    #         self.results["pphi"][i] = pphi

    #         # check for some conditions over here
    #         # if r >= self.escaperadius:  # if particle has escaped
    #         #     break
    #         if r <= EARTH_RADIUS: # particle already on / within surface of earth
    #             break

    #         # # print(self.particle.rigidity)
    #         # if self.particle.rigidity < VertRigidityCutoff(r, theta):
    #         #     break

    #         # update initial_value
    #         initial_value = np.array([t, r, theta, phi, pr, ptheta, pphi])

    #         i += 1

    #     # get the end trajectory from the last spherical coordinate values
    #     endTraj = TrajectoryPoint()
    #     endTraj.set_from_sphericalCoord(r, theta, phi)

    #     print(self.results)

    #     return (startTraj, endTraj)

    # get initial position and velocities for preparation
    # of integration
    def get_initial_values(self, trajectory, energy):
        # get initial position in spherical coordinates
        (r0, theta0, phi0) = trajectory.sphericalCoord()
        # set particle momenta and get their spherical coordinate equivalents
        self.particle.set_momentum(energy)
        v0 = self.particle.set_velocity()
        # (pr0, pth0, pph0) = get_sphcomp_momentum(p0, r0, theta0)
        (vr0, vth0, vph0) = vp_components_spherical(v0, r0, theta0)
        # initial time
        t0 = 0.
        # put all information into some initial value array
        # ival = np.array([t0, r0, theta0, phi0, pr0, pth0, pph0])
        ival = np.array([t0, r0, theta0, phi0, vr0, vth0, vph0])
        return ival


class TrajectoryPoint:
    '''
    Constructs a single point of the trajectory given the following information:
        - latitude: the magnetic latitude (0 = magnetic equator)
        - longitude: the magnetic longitude (0 = prime meridian)
        - altitude: the distance from Earth's surface [km]
        - zenithAngle: the angle from the local zenith
        - azimuthAngle: the angle from the local North
    '''
    def __init__(self,
                latitude=0.,
                 longitude=0.,
                 altitude=1.,
                 zenithAngle=0.,
                 azimuthAngle=0.):
        self.latitude = latitude
        self.longitude = longitude
        
        self.altitude = altitude
        self.zenithAngle = zenithAngle
        self.azimuthAngle = azimuthAngle

    # converts into spherical coordinates
    # r [EARTH_RADIUS], theta [rad], phi [rad]
    def sphericalCoord(self):
        # r = (self.altitude + EARTH_RADIUS) / EARTH_RADIUS

        # r = self.altitude + EARTH_RADIUS
        # longitude and latitude to (r0, theta0, phi0)
        r0 = EARTH_RADIUS
        theta0 = (90. - self.latitude) * DEG_TO_RAD
        phi0 = (180. + self.longitude) * DEG_TO_RAD

        # print(r0, theta0, phi0)

        # 3-vector for altitude, zenith angle and azimuthal angle
        # just for my own reference
        l = self.altitude
        xi = self.zenithAngle * DEG_TO_RAD
        alpha = self.azimuthAngle * DEG_TO_RAD

        # some convenient expression
        d = self.altitude * np.tan(xi) * np.cos(
            alpha)

        # print(d, l, xi, alpha)

        phi = phi0 - np.arctan2(d, EARTH_RADIUS * np.tan(theta0))
        # rcos(phi), obtained from cosine law
        rcosphi = np.sqrt((EARTH_RADIUS * np.cos(phi0))**2. +
                        ((self.altitude * np.cos(alpha)) /
                         np.cos(xi))**2. +
                        2 * EARTH_RADIUS * self.altitude *
                        np.cos(alpha) * np.cos(phi0))
        r = rcosphi / np.cos(phi)
        theta = theta0 - (d / rcosphi)

        # print(r, theta, phi, rcosphi)
        return np.array([r, theta, phi])

    def set_from_sphericalCoord(self, r, theta, phi):
        self.latitude = 90. - (theta * RAD_TO_DEG)
        self.longitude = (phi * RAD_TO_DEG) - 180.
        self.altitude = r - EARTH_RADIUS


    def __str__(self):
        return "Latitude: {0}, Longitude: {1}, Altitude: {3}".format(self.latitude, self.longitude, self.altitude)

    # def __eq__(self, other):
    #     return self.latitude == other.longitude and \
    #         self.longitude == other.longitude and \
    #             self.altitude == other.altitude

'''

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
        self.particle.set_rigidity_from_energy(energy)
        # (pr0, pth0, pph0) = get_sphcomp_momentum(self.particle.momentum, r0, theta0)
        # # initial time
        # t0 = 0.
        # # put all information into some initial value array
        # initial_value = np.array([t0, r0, theta0, phi0, pr0, pth0, pph0])

        # print(initial_value)

        initial_value = self.get_initial_values(startTraj, energy)
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
            self.particle.set_rigidity_from_momentum()
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
            if r <= self.escaperadius:
                break
            if r <= EARTH_RADIUS and i > 0:  # r*RE <= RE
                break
            
            # print(self.particle.rigidity)
            if self.particle.rigidity < VertRigidityCutoff(r, theta):
                break

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

'''