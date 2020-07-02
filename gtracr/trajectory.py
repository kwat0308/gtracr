'''
Keeps track of particle trajectory with considerations to cutoffs and E-W effects etc
'''

import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

# from gtracr.utils import EARTH_RADIUS, g10, DEG_TO_RAD, get_sphcomp_momentum
from gtracr.constants import EARTH_RADIUS, DEG_TO_RAD, RAD_TO_DEG
from gtracr.utils import CarCoord_to_SphCoord, CarVel_to_SphVel
# from gtracr.runge_kutta import runge_kutta
# from gtracr.runge_kutta import RungeKutta
from RungeKutta import RungeKutta
# import RungeKutta
from gtracr.add_particle import particle_dict

key_list = ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]


class Trajectory:
    '''
    Class that controls the trajectory of a particle at some given energy / rigidity
    Members:
    - particleName: the label of the particle
    - latitude: the geographic latitude, with 0 defined at the equator in degrees
    - longitude: the geographic longitude, with 0 defined at the Prime Meridian in degrees
    - altitude: the height from sea level (0=sea level) in km
    - zenithAngle: the angle from the local zenith, with 0 being at the local zenith
    - azimuthAngle: the angle with 0 being in the direction of the East hemisphere from the Prime Meridian in the local tangent plane
    - energy: the particle energy
    - rigidity: the particle rigidity (momentum / charge)
    - escapeAltitude: the altitude in which the particle has "escaped" Earth (default 1000km)
    - maxBuffer: maximum length of array
    '''
    def __init__(self,
                 particleName,
                 latitude,
                 longitude,
                 altitude,
                 zenithAngle,
                 azimuthAngle,
                    energy=None,
                 rigidity=None,
                 escapeAltitude=1000.,
                 maxBuffer=10000):
        self.particle = particle_dict[particleName]
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.zenithAngle = zenithAngle
        self.azimuthAngle = azimuthAngle
        self.escapeAltitude = escapeAltitude

        # define rigidity and energy only if they are provided, evaluate for the other member
        # also set momentum in each case

        if rigidity is None:
            self.particle.get_rigidity_from_energy(energy)
            self.particle.set_momentum_from_energy(energy)
            self.rigidity = self.particle.rigidity
            self.energy = energy
            # self.particle.set_momentum_from_energy(energy)
        elif energy is None:
            self.particle.set_momentum_from_rigidity(rigidity)
            self.rigidity = rigidity
            self.particle.rigidity = rigidity
            self.energy = self.particle.get_energy_from_rigidity(rigidity)
            # self.particle.set_momentum_from_rigidity(rigidity)
        # elif rigidity is None and energy is None:
        else:
            raise Exception(
                "Provide either energy or rigidity as input, not both!")

        # print(self.rigidity)
        # if rigidity is None and energy is None:
        #     raise Exception("Provide either energy or rigidity as input!")

        # self.energy = self.particle.get_energy_from_rigidity(rigidity) if energy is None else energy
        # self.rigidity = self.particle.get_rigidity_from_energy(energy) if rigidity is None else rigidity

        self.particle.set_velocity()  # set velocity here
        # print(self.particle.momentum, self.particle.velocity)

        self.particleEscaped = False  # check if trajectory is allowed or not

        # initialize required arrays here
        self.maxBuffer = maxBuffer
        # self.time_array = np.zeros(maxBuffer)
        # self.TJP_array = np.zeros(maxBuffer)
        self.time_array = [None] * maxBuffer
        self.TJP_array = [None] * maxBuffer  # list to append TJP objects
        # self.time_array = []
        # self.TJP_array = []

    # get the initial trajectory points based on the latitude, longitude, altitude, zenith, and azimuth
    # returns tuple of 2 trajectory points (the initial one and the first one relating to that of the zenith and azimuth one)
    def getInitTJP(self):

        origin_TJP = TrajectoryPoint(self.latitude, self.longitude,
                                     self.altitude)

        print(origin_TJP)
        # print(origin_TJP.getSphericalCoord())

        # transformation process for coordinate
        originCoord = origin_TJP.getCartesianCoord()
        LTPCoord = self.getLTPCoord(mag=1e-10)
        print(originCoord, LTPCoord)
        print(self.tf_matrix())
        (x, y, z) = self.LTP_to_ECEF(originCoord, LTPCoord)
        (r, theta, phi) = CarCoord_to_SphCoord(x, y, z)

        print(x, y, z)
        print(r, theta, phi)

        # transformation for velocity
        originVel = np.array([0., 0., 0.])
        LTPVel = self.getLTPCoord(mag=self.particle.velocity)
        (vx, vy, vz) = self.LTP_to_ECEF(originVel, LTPVel)

        (vr, vtheta, vphi) = CarVel_to_SphVel(vx, vy, vz, r, theta, phi)

        # create new trajectory point and set the new coordinate and velocity
        init_TJP = TrajectoryPoint(vr=vr, vtheta=vtheta, vphi=vphi)
        init_TJP.setSphericalCoord(r, theta, phi)
        # init_TJP.setCarVelocity(vr, vtheta, vphi)

        print(init_TJP)

        return (origin_TJP, init_TJP)

    # evaluates the trajectory using Runge-Kutta methods
    def getTrajectory(self, maxStep=10000, stepSize=0.1):

        # check if maxStep > maxBuffer, if so then update this
        # there is a better way to do this, im sure
        if maxStep > self.maxBuffer:
            # self.time_array = np.zeros(maxStep)
            # self.TJP_array = np.zeros(maxStep)
            self.time_array.extend([None] * ((maxStep) - self.maxBuffer))
            self.TJP_array.extend([None] * ((maxStep) - self.maxBuffer))
        elif maxStep < self.maxBuffer:
            self.time_array = self.time_array[:maxStep]
            self.TJP_array = self.TJP_array[:maxStep]

        # print(self.time_array, self.TJP_array)

        # self.time_array = [None] * maxStep
        # self.TJP_array = [None] * maxStep

        # initialize array
        # time_array = np.zeros(maxStep)
        # TJP_array = np.zeros(maxStep)

        # get the initial trajectory points
        (origin_TJP, init_TJP) = self.getInitTJP()

        # append both time array and TJP
        self.time_array[0:2] = [0., stepSize]
        self.TJP_array[0:2] = [origin_TJP, init_TJP]
        # self.time_array[0] = stepSize
        # self.TJP_array[0] = init_TJP

        # print(self.TJP_array[0], self.TJP_array[1])

        # start iteration process
        # RKI = RKIntegrator(self.particle.mass, self.particle.charge)
        RKI = RungeKutta(self.particle.charge, self.particle.mass, stepSize)
        i = 2
        curr_TJP = init_TJP
        t = stepSize
        (r, theta, phi, vr, vtheta, vphi) = curr_TJP.spherical()
        # valtup required for integration process due to conversions between theta values and latitude
        # valtup = self.valtup(t, curr_TJP)
        
        # valtup = (t, r, theta, phi, vr, vtheta, vphi)
        # ivals = (t, init_TJP.spherical())
        # ivals = init_TJP.spherical().insert(0, t)
        while i < maxStep:
            # print(t, curr_TJP)
            # (t, curr_TJP, valtup) = self.evalTrajectory(t, curr_TJP, stepSize, valtup)
            [t, r, theta, phi, vr, vtheta,
            vphi] = RKI.evaluate([t, r, theta, phi, vr, vtheta, vphi])

            print(t, r, theta, phi, vr, vtheta, vphi, '\n')
            # print(phi)
            # print(theta)

            # (t, r, theta, phi, vr, vtheta, vphi) = RungeKutta.evaluate(self.particle.mass, self.particle.charge, stepSize, t, r, theta, phi, vr, vtheta, vphi)

            new_TJP = TrajectoryPoint(vr=vr, vtheta=vtheta, vphi=vphi)
            new_TJP.vr = vr
            new_TJP.vtheta = vtheta
            new_TJP.vphi = vphi
            new_TJP.setSphericalCoord(r, theta, phi)
            # curr_TJP.vr = vr
            # curr_TJP.vtheta = vtheta
            # curr_TJP.vphi = vphi
            # curr_TJP.setSphericalCoord(r, theta, phi)
            # valtup = (t, r, theta, phi, vr, vtheta,
            # vphi)
            curr_TJP = new_TJP

            self.time_array[i] = t
            self.TJP_array[i] = new_TJP

            # conditions
            if curr_TJP.altitude > self.escapeAltitude:
                self.particleEscaped = True
                self.time_array = self.time_array[:i]
                self.TJP_array = self.TJP_array[:i]
                break

            if curr_TJP.altitude < 0.:
                self.time_array = self.time_array[:i]
                self.TJP_array = self.TJP_array[:i]
                break
            
            if (i-2) % (maxStep // 10) == 0 and (i-2) != 0:
                print("{0} iterations completed".format(i-2))

            # ivals = init_TJP.spherical().insert(0, t)          
            # valtup = self.valtup(t, curr_TJP)

            # some looping checker
            i += 1

        # trim zero values at the end
        # self.trim_arrays()
        # np.trim_zeros(self.time_array, trim="b")
        # np.trim_zeros(self.TJP_array, trim="b")

        # print(self.time_array, self.TJP_array)
        # print(self.time_array)

        print("All done!\n")

    # evaluate the trajectory at some time, position, and velocity using TrajectoryPoints
    # # returns a new time and new TrajPoint
    # def evalTrajectory(self, t0, TJP, stepSize, valtup):

    #     # (r0, theta0, phi0) = TJP.getSphericalCoord()

    #     # print(theta0)

    #     # valtup = (t0, r0, theta0, phi0, TJP.vr, TJP.vtheta, TJP.vphi)

    #     # print(valtup)

    #     (t, r, theta, phi, vr, vtheta,
    #      vphi) = runge_kutta(self.particle.mass, self.particle.charge, stepSize, *valtup)

    #     # print(t, r, theta, phi, vr, vtheta, vphi, '\n')
    #     # print(phi)
    #     # print(theta)

    #     # new_TJP = TrajectoryPoint(vr=vr, vtheta=vtheta, vphi=vphi)
    #     TJP.vr = vr
    #     TJP.vtheta = vtheta
    #     TJP.vphi = vphi
    #     TJP.setSphericalCoord(r, theta, phi)
    #     valtup = (t, r, theta, phi, vr, vtheta,
    #      vphi)

    #     # print(TJP, '\n')

    #     return np.array([t, TJP, valtup])

    # trim the unnecessary values at the end of the array
    # def trim_arrays(self):
    #     # np.trim_zeros(self.time_array, trim="b")
    #     # np.trim_zeros(self.TJP_array, trim="b")
    #     for i, val in enumerate(self.time_array):
    #         if val == None:
    #             self.time_array = self.time_array[:i]
    #             break
    #     for i, val in enumerate(self.TJP_array):
    #         if val == None:
    #             self.TJP_array = self.TJP_array[:i]
    #             break

    # get the cartesian coordinates from the array of trajectory points for plotting purposes
    def getPlotter(self):

        n = len(self.TJP_array)
        data_dict = {
            "t": self.time_array,
            "x": np.zeros(n),
            "y": np.zeros(n),
            "z": np.zeros(n),
            "vx": np.zeros(n),
            "vy": np.zeros(n),
            "vz": np.zeros(n)
        }

        for i, TJP in enumerate(self.TJP_array):
            # print(TJP)
            (x, y, z) = TJP.getCartesianCoord()
            (vx, vy, vz) = TJP.getCartesianVelocity()

            data_dict["x"][i] = x 
            data_dict["y"][i] = y 
            data_dict["z"][i] = z 
            data_dict["vx"][i] = vx
            data_dict["vy"][i] = vy
            data_dict["vz"][i] = vz

        return data_dict

    # convert between local tangent plane coordinates to ECEF coordinates
    def LTP_to_ECEF(self, originCoord, LTPCoord):
        return originCoord + np.dot(self.tf_matrix(), LTPCoord)

    # get the local tangent plane coordinates (in Cartesian) from zenith and azimuthal angles
    def getLTPCoord(self, mag):
        xi = self.zenithAngle * DEG_TO_RAD
        alpha = self.azimuthAngle * DEG_TO_RAD

        print(xi, alpha)

        xt = mag * np.sin(xi) * np.cos(alpha)
        yt = mag * np.sin(xi) * np.sin(alpha)
        zt = mag * np.cos(xi)

        return np.array([xt, yt, zt])

    # the transformation matrix from the local tangent plane to ECEF (Cartesian) coordinates
    # source: http://archive.psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf
    def tf_matrix(self):
        lmbda = self.latitude * DEG_TO_RAD
        eta = self.longitude * DEG_TO_RAD

        # print(lmbda, eta)

        row1 = np.array([
            -np.sin(eta), -np.cos(eta) * np.sin(lmbda),
            np.cos(lmbda) * np.cos(eta)
        ])
        row2 = np.array([
            np.cos(eta), -np.sin(lmbda) * np.sin(eta),
            np.cos(lmbda) * np.sin(eta)
        ])
        row3 = np.array([0., np.cos(lmbda), np.sin(lmbda)])

        return np.array([row1, row2, row3])

    # tuple that contains time and the values at some trajectory point
    def valtup(self, t, TJP):
        (r, theta, phi, vr, vtheta, vphi) = TJP.spherical()
        return (t, r, theta, phi, vr, vtheta, vphi)


class TrajectoryPoint:
    '''
    Class that records a single point in the particle trajectory
    Members:
    - latitude: the geographic latitude, with 0 at the equator (+90 at the North Pole, -90 at South pole) in degrees
    - longitude: the geographic longitude, with 0 at the Prime Meridian (negative degrees towards Western hemisphere) in degrees
    - altitude: the altitude above sea level
    - vr: the radial velocity
    - vtheta: the velocity in the polar direction
    - vphi: the velocity in the azimuthal direction
    '''
    def __init__(self,
                 latitude=0.,
                 longitude=0.,
                 altitude=0.,
                 vr=0.,
                 vtheta=0.,
                 vphi=0.):
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.vr = vr
        self.vtheta = vtheta
        self.vphi = vphi

    # get spherical coordinate equivalents of the latitude, longitude, and altitude
    def getSphericalCoord(self):
        r = EARTH_RADIUS + self.altitude
        theta = (
            90. - self.latitude
        ) * DEG_TO_RAD  # theta defined in [0, pi], theta = 90 at equator
        phi = self.longitude * DEG_TO_RAD  # phi defined in [-pi, pi], phi = 0 at prime meridian

        return np.array([r, theta, phi])

    # set latitude, longitude, altitude from spherical coordiantes
    def setSphericalCoord(self, r, theta, phi):
        self.latitude = 90. - (theta * RAD_TO_DEG)
        # self.longitude = (phi * RAD_TO_DEG) - 180.
        self.longitude = phi * RAD_TO_DEG
        self.altitude = r - EARTH_RADIUS

    # get cartesian coordinate equivalents of the latitude, longitude, and altitude
    def getCartesianCoord(self):
        lmbda = self.latitude * DEG_TO_RAD
        eta = self.longitude * DEG_TO_RAD
        x = (EARTH_RADIUS + self.altitude) * np.cos(lmbda) * np.cos(eta)
        y = (EARTH_RADIUS + self.altitude) * np.cos(lmbda) * np.sin(eta)
        z = (EARTH_RADIUS + self.altitude) * np.sin(lmbda)

        return np.array([x, y, z])

    # set latitude, longitude, altitude from Cartesian coordiantes
    def setCartesianCoord(self, x, y, z):
        pass

    # set vr, vtheta, vphi from cartesian velocities
    def setCartesianVelocity(self, vx, vy, vz):
        (r, theta, phi) = self.getSphericalCoord()

        self.vr = vx * np.sin(theta) * np.cos(phi) + vy * np.sin(
            theta) * np.sin(phi) + vz * np.cos(theta)
        self.vtheta = (vx * np.cos(theta) * np.cos(phi) + vy * np.cos(theta) *
                       np.sin(phi) - vz * np.sin(theta)) / r
        self.vphi = (-vx * np.sin(phi) + vy * np.cos(phi)) / (r *
                                                              np.sin(theta))

        # return np.array([vr, vtheta, vphi])

    # return cartesian velocities from spherical ones
    def getCartesianVelocity(self):
        (r, theta, phi) = self.getSphericalCoord()
        vx = self.vr * np.sin(theta) * np.cos(phi) + r * self.vtheta * np.cos(
            theta) * np.cos(phi) - r * self.vphi * np.sin(theta) * np.sin(phi)
        vy = self.vr * np.sin(theta) * np.sin(phi) + r * self.vtheta * np.cos(
            theta) * np.sin(phi) + r * self.vphi * np.sin(theta) * np.cos(phi)
        vz = self.vr * np.cos(theta) - r * self.vtheta * np.sin(theta)

        return np.array([vx, vy, vz])

    def __str__(self):
        return "Latitude: {:.6f}, Longitude: {:.6f}, Altitude: {:.6f}, Velocity (vr, vtheta, vphi): ({:.6f}, {:.6f}, {:.6f})".format(
            self.latitude, self.longitude, self.altitude, self.vr, self.vtheta,
            self.vphi)

    # both position and velocity in spherical coordinates
    def spherical(self):
        return np.concatenate((self.getSphericalCoord(), np.array([self.vr, self.vtheta, self.vphi])), axis=0)

    # both position and velocity in Cartesian coordinates
    def cartesian(self):
        return np.concatenate((self.getCartesianCoord(), self.getCartesianVelocity()), axis=0)


# class ParticleTrajectory:
#     '''
#     Traces trajectory of a single particle, given the following information:
#     - particleName: the name of the particle of interest, obtained from Particle class
#     - energy : the energy of particle / cosmic ray
#     - startLatitude : initial latitude of location in decimal notation
#     - startLongitude : initial longitude of location in decimal notation
#     - startAltitude : the starting altitude of the particle [km]
#     - stopAltitude : the altitude in which the particle trajectory ends [km]
#     - maxStep : the maximum number of steps to integrate for (default=10000)

#     '''
#     def __init__(self,
#                  particleName,
#                  energy,
#                  startLatitude=0.,
#                  startLongitude=0.,
#                  startAltitude=1.,
#                  stopAltitude=1000.,
#                  maxStep = 10000,
#                  stepSize=0.1):
#         self.particle = particle_dict[
#             particleName]  # should be obtained from some dictionary
#         self.energy = energy
#         self.startLatitude = startLatitude
#         self.startLongitude = startLongitude
#         self.startAltitude = startAltitude
#         self.stopAltitude = stopAltitude
#         self.maxStep = maxStep
#         # self.stepSize = np.abs(stopAltitude - startAltitude) / maxStep
#         self.stepSize = stepSize
#         # self.maxStep = int(np.abs(stopAltitude - startAltitude) / stepSize)
#         self.results = {key: np.zeros(maxStep) for key in key_list}
#         # self.results = {key: [] for key in key_list}

#     # def __init__(self,
#     #              particleName,
#     #              startAltitude=1,
#     #              stopAltitude=500,
#     #              maxStep=10000):
#     #     self.particle = particleDict[
#     #         particleName]  # should be obtained from some dictionary
#     #     self.startAltitude = startAltitude
#     #     self.stopAltitude = stopAltitude
#     #     self.maxStep = maxStep
#     #     self.stepSize = np.abs(stopAltitude - startAltitude) / maxStep
#     #     self.results = {key:np.zeros(maxStep) for key in key_list}

#     # get the trajectory from a provided zenith and azimuthal angle

#     def getTrajectory(self, zenith=0., azimuth=0.):
#         # create trajectory point for initial point
#         startTraj = TrajectoryPoint(self.startLatitude, self.startLongitude,
#                                     self.startAltitude)
#         startTraj.set_from_horizontal(zenith, azimuth)

#         # get initial value
#         (t, r, theta, phi, vr, vtheta,
#              vphi) = self.get_initial_values(startTraj, self.energy)
#         # info(1, initial_value)
#         # print(initial_value)
#         # print(startTraj)

#         # results = {key: np.zeros(self.maxStep) for key in key_list}

#         i = 0
#         # r=EARTH_RADIUS
#         # while i < self.maxStep:
#         while r < EARTH_RADIUS + self.stopAltitude:
#             (t, r, theta, phi, vr, vtheta,
#              vphi) = runge_kutta(self.particle, self.stepSize, (t, r, theta, phi, vr, vtheta,
#              vphi))
#             # info(3, r)
#             # print(r)
#             # check for some conditions over here
#             # if r <= EARTH_RADIUS:  # particle already on / within surface of earth
#             #     break

#             # convert velocities to momenta
#             vmag = vmag_spherical(vr, vtheta, vphi, r, theta)
#             # print(vmag)
#             # pr = gamma(vmag)*self.particle.mass*vr
#             # ptheta = gamma(vmag)*self.particle.mass*vtheta
#             # pphi = gamma(vmag)*self.particle.mass*vphi
#             # print(t, r, theta, phi, pr, ptheta,pphi)
#             # set this to particle's new position and momentum
#             self.particle.momentum = gamma(vmag) * self.particle.mass * vmag
#             # self.particle.momentum = vmag_spherical(pr, ptheta, pphi, r, theta)
#             self.particle.set_rigidity_from_momentum()
#             # print(self.particle.rigidity)

#             (pr, ptheta,
#              pphi) = vp_components_spherical(self.particle.momentum, r, theta)

#             # append
#             self.results["t"][i] = t
#             self.results["r"][i] = r
#             self.results["theta"][i] = theta
#             self.results["phi"][i] = phi
#             self.results["pr"][i] = pr
#             self.results["ptheta"][i] = ptheta
#             self.results["pphi"][i] = pphi

#             # self.results["t"].append(t)
#             # self.results["r"].append(r)
#             # self.results["theta"].append(theta)
#             # self.results["phi"].append(phi)
#             # self.results["pr"].append(pr)
#             # self.results["ptheta"].append(ptheta)
#             # self.results["pphi"].append(pphi)

#             # check for some conditions over here
#             # if r <= EARTH_RADIUS:  # particle already on / within surface of earth
#             #     break

#             # update initial_value
#             # initial_value = np.array([t, r, theta, phi, pr, ptheta, pphi])
#             # initial_value = np.array([t, r, theta, phi, vr, vtheta, vphi])

#             i += 1

#         # get the end trajectory from the last spherical coordinate values
#         # print(r, theta, phi)
#         endTraj = TrajectoryPoint()
#         endTraj.set_from_sphericalCoord(r, theta, phi)
#         # endTraj.set_from_sphericalCoord(r, theta / np.pi, phi / (2.*np.pi))
#         # trim the zeros from the arrays if there is a break
#         self.trim_zeros()
#         # info(2, results)
#         # print(self.results)
#         # print(len(self.results["r"]))

#         # return (startTraj, endTraj)
#         return self.results

#     # # get the particle trajectory given:
#     # # - energy: the energy of the cosmic ray as measured from detector
#     # # - startLongitude: the magnetic longitude
#     # # - startLatitude: the magnetic latitude
#     # # - startZenith: the angle from the local zenith
#     # # - startAzimuth: the angle from the local North
#     # # default all to zero
#     # # return tuple of start and end trajectories
#     # def getTrajectory(self,
#     #                   energy,
#     #                   startLongitude=0.,
#     #                   startLatitude=0.,
#     #                   startZenith=0.,
#     #                   startAzimuth=0.):
#     #     startTraj = TrajectoryPoint(startLongitude, startLatitude,
#     #                                 self.startAltitude, startZenith,
#     #                                 startAzimuth)

#     #     # set rigidity
#     #     self.particle.set_rigidity_from_energy(energy)
#     #     # get initial position in spherical coordinates
#     #     initial_value = self.get_initial_values(startTraj, energy)
#     #     # perform runge_kutta
#     #     i = 0
#     #     while i < self.maxStep:
#     #         (t, r, theta, phi, pr, ptheta,
#     #          pphi) = runge_kutta(self.particle, self.stepSize, initial_value)

#     #         # print(t, r, theta, phi, pr, ptheta,pphi)
#     #         # set this to particle's new position and momentum
#     #         self.particle.momentum = get_momentum_magnitude(pr, ptheta, pphi, r, theta)
#     #         self.particle.set_rigidity_from_momentum()
#     #         # append
#     #         self.results["t"][i] = t
#     #         self.results["r"][i] = r
#     #         self.results["theta"][i] = theta
#     #         self.results["phi"][i] = phi
#     #         self.results["pr"][i] = pr
#     #         self.results["ptheta"][i] = ptheta
#     #         self.results["pphi"][i] = pphi

#     #         # check for some conditions over here
#     #         # if r >= self.escaperadius:  # if particle has escaped
#     #         #     break
#     #         if r <= EARTH_RADIUS: # particle already on / within surface of earth
#     #             break

#     #         # # print(self.particle.rigidity)
#     #         # if self.particle.rigidity < VertRigidityCutoff(r, theta):
#     #         #     break

#     #         # update initial_value
#     #         initial_value = np.array([t, r, theta, phi, pr, ptheta, pphi])

#     #         i += 1

#     #     # get the end trajectory from the last spherical coordinate values
#     #     endTraj = TrajectoryPoint()
#     #     endTraj.set_from_sphericalCoord(r, theta, phi)

#     #     print(self.results)

#     #     return (startTraj, endTraj)

#     # get initial position and velocities for preparation
#     # of integration
#     def get_initial_values(self, trajectory, energy):
#         # get initial position in spherical coordinates
#         (r0, theta0, phi0) = trajectory.sphericalCoord()
#         # info(4, (r0, theta0, phi0))
#         print(r0, theta0, phi0)
#         # set particle momenta and get their spherical coordinate equivalents
#         self.particle.set_momentum(energy)
#         v0 = self.particle.set_velocity()
#         # (pr0, pth0, pph0) = get_sphcomp_momentum(p0, r0, theta0)
#         (vr0, vth0, vph0) = vp_components_spherical(v0, r0, theta0)
#         # initial time
#         t0 = 0.
#         # put all information into some initial value array
#         # ival = np.array([t0, r0, theta0, phi0, pr0, pth0, pph0])
#         ival = np.array([t0, r0, theta0, phi0, vr0, vth0, vph0])
#         return ival

#     # trim zeros at end of list
#     # an indicator for geomagnetic cutoff (array length comparison)
#     def trim_zeros(self):
#         for key, arr in list(self.results.items()):
#             self.results.update({key:np.trim_zeros(arr, 'b')})
#         # for key, arr in list(self.results.items()):
#         #     # np.trim_zeros(arr, trim="b")
#         #     nz_ind = np.nonzero(arr)
#         #     arr = arr[nz_ind]

# class TrajectoryPoint:
#     '''
#     Constructs a single point of the trajectory given the following information:
#         - latitude: the magnetic latitude (0 = magnetic equator)
#         - longitude: the magnetic longitude (0 = prime meridian)
#         - altitude: the distance from Earth's surface [km]
#         # - zenithAngle: the angle from the local zenith
#         # - azimuthAngle: the angle from the local North
#     '''
#     def __init__(self, latitude=0., longitude=0., altitude=0.):
#         self.latitude = latitude
#         self.longitude = longitude
#         self.altitude = altitude
#         # self.zenithAngle = zenithAngle
#         # self.azimuthAngle = azimuthAngle

#     # set a new latitude and longitude from the provided
#     # zenith and azimuthal angles
#     def set_from_horizontal(self, zenithAngle, azimuthAngle):
#         # 3-vector for altitude, zenith angle and azimuthal angle
#         # just for my own reference
#         # l = self.altitude
#         xi = zenithAngle * DEG_TO_RAD
#         alpha = azimuthAngle * DEG_TO_RAD
#         # print(self.latitude, self.longitude)
#         print(xi, alpha)
#         print(np.tan(xi), np.cos(alpha))
#         print(self.altitude * np.tan(xi) * np.cos(alpha))
#         # if altitude is zero make latitude and longitude same as before
#         if np.abs(self.altitude) < 1e-7:
#             pass
#         else:   # otherwise compute based on altitude, zenith, and azimuth
#             self.latitude += np.abs(self.altitude) * np.tan(xi) * np.sin(alpha)
#             self.longitude += np.abs(self.altitude) * np.tan(xi) * np.cos(alpha)
#         print(self.latitude, self.longitude)

#         # some convenient expression
#         # d = self.altitude * np.tan(xi) * np.cos(
#         #     alpha)

#         # # print(d, l, xi, alpha)

#         # phi = phi0 - np.arctan2(d, EARTH_RADIUS * np.tan(theta0))
#         # # rcos(phi), obtained from cosine law
#         # rcosphi = np.sqrt((EARTH_RADIUS * np.cos(phi0))**2. +
#         #                 ((self.altitude * np.cos(alpha)) /
#         #                  np.cos(xi))**2. +
#         #                 2 * EARTH_RADIUS * self.altitude *
#         #                 np.cos(alpha) * np.cos(phi0))
#         # r = rcosphi / np.cos(phi)
#         # theta = theta0 - (d / rcosphi)

#     # set the latitude, longitude, and altitude from
#     # spherical coordinates
#     def set_from_sphericalCoord(self, r, theta, phi):
#         print(phi*RAD_TO_DEG)
#         self.latitude = (theta * RAD_TO_DEG) - 90.
#         # self.longitude = (phi * RAD_TO_DEG) - 180.
#         self.longitude = phi*RAD_TO_DEG
#         self.altitude = r - EARTH_RADIUS

#     # converts into spherical coordinates
#     # r [km], theta [rad], phi [rad]
#     def sphericalCoord(self):
#         # r = (self.altitude + EARTH_RADIUS) / EARTH_RADIUS

#         r = self.altitude + EARTH_RADIUS
#         # longitude and latitude to (r0, theta0, phi0)
#         # r = EARTH_RADIUS + sel
#         theta = ((90. + self.latitude) * DEG_TO_RAD)
#         # phi = ((180. + self.longitude) * DEG_TO_RAD)
#         phi = ( self.longitude * DEG_TO_RAD)

#         # print(r0, theta0, phi0)

#         # # 3-vector for altitude, zenith angle and azimuthal angle
#         # # just for my own reference
#         # l = self.altitude
#         # xi = self.zenithAngle * DEG_TO_RAD
#         # alpha = self.azimuthAngle * DEG_TO_RAD

#         # # some convenient expression
#         # d = self.altitude * np.tan(xi) * np.cos(
#         #     alpha)

#         # # print(d, l, xi, alpha)

#         # phi = phi0 - np.arctan2(d, EARTH_RADIUS * np.tan(theta0))
#         # # rcos(phi), obtained from cosine law
#         # rcosphi = np.sqrt((EARTH_RADIUS * np.cos(phi0))**2. +
#         #                 ((self.altitude * np.cos(alpha)) /
#         #                  np.cos(xi))**2. +
#         #                 2 * EARTH_RADIUS * self.altitude *
#         #                 np.cos(alpha) * np.cos(phi0))
#         # r = rcosphi / np.cos(phi)
#         # theta = theta0 - (d / rcosphi)

#         # print(r, theta, phi, rcosphi)
#         return np.array([r, theta, phi])

#     def __str__(self):
#         return "Latitude: {0}, Longitude: {1}, Altitude: {2}".format(
#             self.latitude, self.longitude, self.altitude)

#     def __eq__(self, other):
#         return self.latitude == other.longitude and \
#             self.longitude == other.longitude and \
#                 self.altitude == other.altitude

# '''

# # get the particle trajectory given:
#     # - energy: the energy of the cosmic ray as measured from detector
#     # - startLongitude: the magnetic longitude
#     # - startLatitude: the magnetic latitude
#     # - startZenith: the angle from the local zenith
#     # - startAzimuth: the angle from the local North
#     # default all to zero
#     def getTrajectory(self,
#                       energy,
#                       startLongitude=0.,
#                       startLatitude=0.,
#                       startZenith=0.,
#                       startAzimuth=0.):
#         startTraj = TrajectoryPoint(startLongitude, startLatitude,
#                                     self.startAltitude, startZenith,
#                                     startAzimuth)

#         # get initial position in spherical coordinates
#         self.particle.set_rigidity_from_energy(energy)
#         # (pr0, pth0, pph0) = get_sphcomp_momentum(self.particle.momentum, r0, theta0)
#         # # initial time
#         # t0 = 0.
#         # # put all information into some initial value array
#         # initial_value = np.array([t0, r0, theta0, phi0, pr0, pth0, pph0])

#         # print(initial_value)

#         initial_value = self.get_initial_values(startTraj, energy)
#         # initialize array
#         # t_arr = np.zeros(self.maxStep)
#         # r_arr = np.zeros(self.maxStep)
#         # th_arr = np.zeros(self.maxStep)
#         # ph_arr = np.zeros(self.maxStep)
#         # pr_arr = np.zeros(self.maxStep)
#         # pth_arr = np.zeros(self.maxStep)
#         # pph_arr = np.zeros(self.maxStep)

#         # perform runge_kutta
#         # this should be in some while loop
#         i = 0
#         while i < self.maxStep:
#             # print(initial_value)
#             # while i < self.maxStep: # (or something like this)
#             (t, r, theta, phi, pr, ptheta,
#              pphi) = runge_kutta(self.particle, self.stepSize, initial_value)

#             # print(t, r, theta, phi, pr, ptheta,pphi)
#             # t = result[0]
#             # r = result[1]
#             # theta = result[2]
#             # phi = result[3]
#             # pr = result[4]
#             # ptheta = result[5]
#             # pphi = result[6]
#             # set this to particle's new position and momentum
#             self.particle.momentum = get_momentum_magnitude(pr, ptheta, pphi, r, theta)
#             self.particle.set_rigidity_from_momentum()
#             # append
#             # t_arr[i] = t
#             # r_arr[i] = r
#             # th_arr[i] = theta
#             # ph_arr[i] = phi
#             # pr_arr[i] = pr
#             # pth_arr[i] = pth
#             # pph_arr[i] = pph
#             self.results["t"][i] = t
#             self.results["r"][i] = r
#             self.results["theta"][i] = theta
#             self.results["phi"][i] = phi
#             self.results["pr"][i] = pr
#             self.results["ptheta"][i] = ptheta
#             self.results["pphi"][i] = pphi

#             # check for some conditions over here
#             if r <= self.escaperadius:
#                 break
#             if r <= EARTH_RADIUS and i > 0:  # r*RE <= RE
#                 break

#             # print(self.particle.rigidity)
#             if self.particle.rigidity < VertRigidityCutoff(r, theta):
#                 break

#             # update initial_value
#             initial_value = np.array([t, r, theta, phi, pr, ptheta, pphi])

#             i += 1

#         # # add to dictionary
#         # self.coord["t"] = t_arr
#         # self.coord["r"] = r_arr
#         # self.coord["theta"] = th_arr
#         # self.coord["phi"] = ph_arr

#         # return (t_arr, r_arr, th_arr, ph_arr)
#         print(self.results)

# vertical rigidity cutoff defined in Baldini's paper
# def VertRigidityCutoff(r, theta):
#     return (14.5 * np.sin(theta)**4.) / r**2.

# '''