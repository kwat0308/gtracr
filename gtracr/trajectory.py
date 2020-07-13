'''
Keeps track of particle trajectory with considerations to cutoffs and E-W effects etc
'''

import os, sys
import numpy as np

# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.constants import EARTH_RADIUS, DEG_TO_RAD, RAD_TO_DEG
from gtracr.trajectory_point import TrajectoryPoint
from RungeKutta import RungeKutta
from gtracr.add_particle import particle_dict

KEY_LIST = ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]


class Trajectory:
    '''
    Class that controls the trajectory of a particle at some given energy / rigidity
    Members:
    - particle_name: the label of the particle
    - latitude: the geographic latitude, with 0 defined at the equator in degrees
    - longitude: the geographic longitude, with 0 defined at the Prime Meridian in degrees
    - altitude: the height from sea level (0=sea level) in km
    - zenith_angle: the angle from the local zenith, with 0 being at the local zenith
    - azimuth_angle: the angle with 0 being in the direction of the East hemisphere from the Prime Meridian in the local tangent plane
    - energy: the particle energy
    - rigidity: the particle rigidity (momentum / charge)
    - escape_altitude: the altitude in which the particle has "escaped" Earth (default 1000km)
    - max_buffer: maximum length of array
    '''
    def __init__(self,
                 particle_name,
                 latitude,
                 longitude,
                 altitude,
                 zenith_angle,
                 azimuth_angle,
                 energy=None,
                 rigidity=None,
                 escape_altitude=10. * EARTH_RADIUS,
                 max_buffer=10000):
        self.particle = particle_dict[particle_name]
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.zenith_angle = zenith_angle
        self.azimuth_angle = azimuth_angle
        self.escape_altitude = escape_altitude

        # define rigidity and energy only if they are provided, evaluate for the other member
        # also set momentum in each case
        self.particle.print()

        if rigidity is None:
            # self.particle.get_rigidity_from_energy(energy)
            # self.particle.set_momentum_from_energy(energy)
            self.particle.set_from_energy(energy)
            self.rigidity = self.particle.rigidity
            self.energy = energy
            # self.particle.set_momentum_from_energy(energy)
        elif energy is None:
            # self.particle.set_momentum_from_rigidity(rigidity)
            # self.rigidity = rigidity
            # self.particle.rigidity = rigidity
            self.particle.set_from_rigidity(rigidity)
            self.rigidity = rigidity
            self.energy = self.particle.get_energy_rigidity()
        # elif rigidity is None and energy is None:
        else:
            raise Exception(
                "Provide either energy or rigidity as input, not both!")

        # print(self.particle.momentum, self.particle.velocity, self.partic)
        self.particle.print()
        self.particle_escaped = False  # check if trajectory is allowed or not

        # initialize required arrays here
        self.max_buffer = max_buffer
        # self.time_array = np.zeros(max_buffer)
        # self.tp_array = np.zeros(max_buffer)
        self.time_array = [None] * max_buffer
        self.tp_array = [None] * max_buffer  # list to append TJP objects

    # get the initial trajectory points based on the latitude, longitude, altitude, zenith, and azimuth
    # returns tuple of 2 trajectory points (the initial one and the first one relating to that of the zenith and azimuth one)
    def detector_to_particle(self):

        detector_tp = TrajectoryPoint()
        detector_tp.set_geodesic_coord(self.latitude, self.longitude)

        print(detector_tp)
        # print(detector_tp.getSphericalCoord())

        # transformation process for coordinate
        detector_coord = detector_tp.cartesian_coord()
        particle_coord = self.get_particle_coord(mag=self.altitude)
        print(detector_coord, particle_coord)
        # print(self.tf_matrix())
        (part_x, part_y, part_z) = self.transform(detector_coord,
                                                  particle_coord)

        print(part_x, part_y, part_z)

        # transformation for velocity
        detector_vel = np.zeros(3)
        particle_vel = self.get_particle_coord(mag=self.particle.velocity)
        print(detector_vel, particle_vel)
        (part_vx, part_vy, part_vz) = self.transform(detector_vel,
                                                     particle_vel)

        print(part_vx, part_vy, part_vz)

        # create new trajectory point and set the new coordinate and velocity
        particle_tp = TrajectoryPoint()
        particle_tp.set_cartesian_coord(part_x, part_y, part_z)
        particle_tp.set_cartesian_velocity(part_vx, part_vy, part_vz)
        # particle_tp.set_spherical_coord(part_r, part_theta, part_phi)

        print(particle_tp)

        return (detector_tp, particle_tp)

    # evaluates the trajectory using Runge-Kutta methods
    def get_trajectory(self, max_step=10000, step_size=0.01):

        # check if max_step > max_buffer, if so then update this
        # there is a better way to do this, im sure
        if max_step > self.max_buffer:
            # self.time_array = np.zeros(max_step)
            # self.tp_array = np.zeros(max_step)
            self.time_array.extend([None] * ((max_step) - self.max_buffer))
            self.tp_array.extend([None] * ((max_step) - self.max_buffer))
        elif max_step < self.max_buffer:
            self.time_array = self.time_array[:max_step]
            self.tp_array = self.tp_array[:max_step]

        # get the initial trajectory points
        (detector_tp, particle_tp) = self.detector_to_particle()

        # append both time array and TJP
        self.time_array[0:2] = [0., step_size]
        self.tp_array[0:2] = [detector_tp, particle_tp]
        # self.time_array[0] = step_size
        # self.tp_array[0] = particle_tp

        # print(self.tp_array[0], self.tp_array[1])

        # start iteration process
        # rk_integrator = rk_integratorntegrator(self.particle.mass, self.particle.charge)
        rk_integrator = RungeKutta(self.particle.charge, self.particle.mass,
                                   step_size)
        i = 2
        part_t = step_size
        # (part_r, part_theta, part_phi, part_vr, part_vtheta,
        #  part_vphi) = particle_tp.spherical()
        (part_r, part_theta, part_phi, part_vr, part_vtheta,
         part_vphi) = tuple(vars(particle_tp).values())
        while i < max_step:
            [
                part_t, part_r, part_theta, part_phi, part_vr, part_vtheta,
                part_vphi
            ] = rk_integrator.evaluate([
                part_t, part_r, part_theta, part_phi, part_vr, part_vtheta,
                part_vphi
            ])

            print(part_t, part_r, part_theta, part_phi, part_vr, part_vtheta,
                  part_vphi, '\n')
            # print(phi)
            # print(theta)

            next_tp = TrajectoryPoint(r=part_r,
                                      theta=part_theta,
                                      phi=part_phi,
                                      vr=part_vr,
                                      vtheta=part_vtheta,
                                      vphi=part_vphi)

            # print(next_tp)
            # next_tp.vr = vr
            # next_tp.vtheta = vtheta
            # next_tp.vphi = vphi
            # next_tp.set_spherical_coord(part_r, part_theta, part_phi)

            self.time_array[i] = part_t
            self.tp_array[i] = next_tp

            # conditions
            # if next_tp.altitude > self.escape_altitude:
            #     self.particle_escaped = True
            #     self.time_array = self.time_array[:i]
            #     self.tp_array = self.tp_array[:i]
            #     break

            # if next_tp.altitude < 0.:
            #     self.time_array = self.time_array[:i]
            #     self.tp_array = self.tp_array[:i]
            #     break

            if next_tp.r > EARTH_RADIUS + self.escape_altitude:
                self.particle_escaped = True
                self.time_array = self.time_array[:i]
                self.tp_array = self.tp_array[:i]
                break

            if next_tp.r < EARTH_RADIUS:
                self.time_array = self.time_array[:i]
                self.tp_array = self.tp_array[:i]
                break

            if (i - 2) % (max_step // 10) == 0 and (i - 2) != 0:
                print("{0} iterations completed".format(i - 2))

            # some looping checker
            i += 1

        print("All done!\n")

    # get the cartesian coordinates from the array of trajectory points for plotting purposes
    def get_plotting_variables(self):

        n = len(self.tp_array)
        data_dict = {
            "t": self.time_array,
            "x": np.zeros(n),
            "y": np.zeros(n),
            "z": np.zeros(n),
            "vx": np.zeros(n),
            "vy": np.zeros(n),
            "vz": np.zeros(n)
        }

        for i, tp in enumerate(self.tp_array):
            # print(tp)
            (x, y, z) = tp.cartesian_coord()
            (vx, vy, vz) = tp.cartesian_velocity()

            data_dict["x"][i] = x
            data_dict["y"][i] = y
            data_dict["z"][i] = z
            data_dict["vx"][i] = vx
            data_dict["vy"][i] = vy
            data_dict["vz"][i] = vz

        return data_dict

    # convert between local tangent plane coordinates to ECEF coordinates
    def transform(self, detector_coord, particle_coord):
        return detector_coord + np.dot(self.tf_matrix(), particle_coord)

    # get the local tangent plane coordinates (in Cartesian) from zenith and azimuthal angles
    def get_particle_coord(self, mag):
        xi = self.zenith_angle * DEG_TO_RAD
        alpha = self.azimuth_angle * DEG_TO_RAD

        # print(xi, alpha)

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


# class ParticleTrajectory:
#     '''
#     Traces trajectory of a single particle, given the following information:
#     - particle_name: the name of the particle of interest, obtained from Particle class
#     - energy : the energy of particle / cosmic ray
#     - startLatitude : initial latitude of location in decimal notation
#     - startLongitude : initial longitude of location in decimal notation
#     - startAltitude : the starting altitude of the particle [km]
#     - stopAltitude : the altitude in which the particle trajectory ends [km]
#     - max_step : the maximum number of steps to integrate for (default=10000)

#     '''
#     def __init__(self,
#                  particle_name,
#                  energy,
#                  startLatitude=0.,
#                  startLongitude=0.,
#                  startAltitude=1.,
#                  stopAltitude=1000.,
#                  max_step = 10000,
#                  step_size=0.1):
#         self.particle = particle_dict[
#             particle_name]  # should be obtained from some dictionary
#         self.energy = energy
#         self.startLatitude = startLatitude
#         self.startLongitude = startLongitude
#         self.startAltitude = startAltitude
#         self.stopAltitude = stopAltitude
#         self.max_step = max_step
#         # self.step_size = np.abs(stopAltitude - startAltitude) / max_step
#         self.step_size = step_size
#         # self.max_step = int(np.abs(stopAltitude - startAltitude) / step_size)
#         self.results = {key: np.zeros(max_step) for key in KEY_LIST}
#         # self.results = {key: [] for key in KEY_LIST}

#     # def __init__(self,
#     #              particle_name,
#     #              startAltitude=1,
#     #              stopAltitude=500,
#     #              max_step=10000):
#     #     self.particle = particleDict[
#     #         particle_name]  # should be obtained from some dictionary
#     #     self.startAltitude = startAltitude
#     #     self.stopAltitude = stopAltitude
#     #     self.max_step = max_step
#     #     self.step_size = np.abs(stopAltitude - startAltitude) / max_step
#     #     self.results = {key:np.zeros(max_step) for key in KEY_LIST}

#     # get the trajectory from a provided zenith and azimuthal angle

#     def get_trajectory(self, zenith=0., azimuth=0.):
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

#         # results = {key: np.zeros(self.max_step) for key in KEY_LIST}

#         i = 0
#         # r=EARTH_RADIUS
#         # while i < self.max_step:
#         while r < EARTH_RADIUS + self.stopAltitude:
#             (t, r, theta, phi, vr, vtheta,
#              vphi) = runge_kutta(self.particle, self.step_size, (t, r, theta, phi, vr, vtheta,
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
#     # def get_trajectory(self,
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
#     #     while i < self.max_step:
#     #         (t, r, theta, phi, pr, ptheta,
#     #          pphi) = runge_kutta(self.particle, self.step_size, initial_value)

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
#         # - zenith_angle: the angle from the local zenith
#         # - azimuth_angle: the angle from the local North
#     '''
#     def __init__(self, latitude=0., longitude=0., altitude=0.):
#         self.latitude = latitude
#         self.longitude = longitude
#         self.altitude = altitude
#         # self.zenith_angle = zenith_angle
#         # self.azimuth_angle = azimuth_angle

#     # set a new latitude and longitude from the provided
#     # zenith and azimuthal angles
#     def set_from_horizontal(self, zenith_angle, azimuth_angle):
#         # 3-vector for altitude, zenith angle and azimuthal angle
#         # just for my own reference
#         # l = self.altitude
#         xi = zenith_angle * DEG_TO_RAD
#         alpha = azimuth_angle * DEG_TO_RAD
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
#         # xi = self.zenith_angle * DEG_TO_RAD
#         # alpha = self.azimuth_angle * DEG_TO_RAD

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
#     def get_trajectory(self,
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
#         # t_arr = np.zeros(self.max_step)
#         # r_arr = np.zeros(self.max_step)
#         # th_arr = np.zeros(self.max_step)
#         # ph_arr = np.zeros(self.max_step)
#         # pr_arr = np.zeros(self.max_step)
#         # pth_arr = np.zeros(self.max_step)
#         # pph_arr = np.zeros(self.max_step)

#         # perform runge_kutta
#         # this should be in some while loop
#         i = 0
#         while i < self.max_step:
#             # print(initial_value)
#             # while i < self.max_step: # (or something like this)
#             (t, r, theta, phi, pr, ptheta,
#              pphi) = runge_kutta(self.particle, self.step_size, initial_value)

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