import os, sys
import numpy as np

# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.constants import EARTH_RADIUS, DEG_TO_RAD, RAD_TO_DEG

# KEY_LIST = ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]


class TrajectoryPoint:
    '''
    Class that records a single point in the particle trajectory
    Members:
    - r: the radial component [m]
    - theta: the polar component, with 0 defined at the North Pole
    - phi: the azimuthal component, with 0 defined at the Prime Meridian (domain: [-pi, pi])
    - vr: the radial velocity
    - vtheta: the velocity in the polar direction
    - vphi: the velocity in the azimuthal direction
    '''
    def __init__(self,
                 r=EARTH_RADIUS,
                 theta=0.,
                 phi=0.,
                 vr=0.,
                 vtheta=0.,
                 vphi=0.):
        self.r = r
        self.theta = theta
        self.phi = phi
        self.vr = vr
        self.vtheta = vtheta
        self.vphi = vphi

    # get geodesic coordinate equivalents of spherical ones
    def geodesic_coordinate(self):
        latitude = 90. - (self.theta * RAD_TO_DEG)
        # longitude = (phi * RAD_TO_DEG) - 180.
        longitude = self.phi * RAD_TO_DEG
        return np.array([latitude, longitude])

    # set latitude, longitude, altitude from geodesic coordiantes
    def set_geodesic_coord(self, latitude, longitude):
        self.theta = (
            90. - latitude
        ) * DEG_TO_RAD  # theta defined in [0, pi], theta = 90 at equator
        self.phi = longitude * DEG_TO_RAD  # phi defined in [-pi, pi], phi = 0 at prime meridian

    # get cartesian coordinate equivalents of spherical ones
    def cartesian_coord(self):
        x = self.r * np.sin(self.theta) * np.cos(self.phi)
        y = self.r * np.sin(self.theta) * np.sin(self.phi)
        z = self.r * np.cos(self.theta)

        return np.array([x, y, z])

    # set latitude, longitude, altitude from Cartesian coordiantes
    def set_cartesian_coord(self, x, y, z):
        self.r = np.sqrt(x**2. + y**2. + z**2.)
        self.theta = np.arccos(z / self.r)
        self.phi = np.arctan2(y, x)

    # set vr, vtheta, vphi from cartesian velocities
    def set_cartesian_velocity(self, vx, vy, vz):
        self.vr = vx * np.sin(self.theta) * np.cos(self.phi) + vy * np.sin(
            self.theta) * np.sin(self.phi) + vz * np.cos(self.theta)
        self.vtheta = (vx * np.cos(self.theta) * np.cos(self.phi) +
                       vy * np.cos(self.theta) * np.sin(self.phi) -
                       vz * np.sin(self.theta)) / self.r
        self.vphi = (-vx * np.sin(self.phi) +
                     vy * np.cos(self.phi)) / (self.r * np.sin(self.theta))

        # return np.array([vr, vtheta, vphi])

    # return cartesian velocities from spherical ones
    def cartesian_velocity(self):
        vx = self.vr * np.sin(self.theta) * np.cos(
            self.phi) + self.r * self.vtheta * np.cos(self.theta) * np.cos(
                self.phi) - self.r * self.vphi * np.sin(self.theta) * np.sin(
                    self.phi)
        vy = self.vr * np.sin(self.theta) * np.sin(
            self.phi) + self.r * self.vtheta * np.cos(self.theta) * np.sin(
                self.phi) + self.r * self.vphi * np.sin(self.theta) * np.cos(
                    self.phi)
        vz = self.vr * np.cos(self.theta) - self.r * self.vtheta * np.sin(
            self.theta)

        return np.array([vx, vy, vz])

    def __str__(self):
        return "Coordinate (r, theta, phi): ({:.12f}, {:.12f}, {:.12f}), \
            \nVelocity (vr, vtheta, vphi): ({:.12f}, {:.12f}, {:.12f})\n".format(
            *tuple(vars(self).values()))

    # both position and velocity in Cartesian coordinates
    def cartesian(self):
        return np.concatenate(
            (self.cartesian_coord(), self.cartesian_velocity()), axis=0)


# class TrajectoryPoint:
#     '''
#     Class that records a single point in the particle trajectory
#     Members:
#     - latitude: the geographic latitude, with 0 at the equator (+90 at the North Pole, -90 at South pole) in degrees
#     - longitude: the geographic longitude, with 0 at the Prime Meridian (negative degrees towards Western hemisphere) in degrees
#     - altitude: the altitude above sea level
#     - vr: the radial velocity
#     - vtheta: the velocity in the polar direction
#     - vphi: the velocity in the azimuthal direction
#     '''
#     def __init__(self,
#                  latitude=0.,
#                  longitude=0.,
#                  altitude=0.,
#                  vr=0.,
#                  vtheta=0.,
#                  vphi=0.):
#         self.latitude = latitude
#         self.longitude = longitude
#         self.altitude = altitude
#         self.vr = vr
#         self.vtheta = vtheta
#         self.vphi = vphi

#     # get spherical coordinate equivalents of the latitude, longitude, and altitude
#     def get_spherical_coord(self):
#         r = EARTH_RADIUS + self.altitude
#         theta = (
#             90. - self.latitude
#         ) * DEG_TO_RAD  # theta defined in [0, pi], theta = 90 at equator
#         phi = self.longitude * DEG_TO_RAD  # phi defined in [-pi, pi], phi = 0 at prime meridian

#         return np.array([r, theta, phi])

#     # set latitude, longitude, altitude from spherical coordiantes
#     def set_spherical_coord(self, r, theta, phi):
#         self.latitude = 90. - (theta * RAD_TO_DEG)
#         # self.longitude = (phi * RAD_TO_DEG) - 180.
#         self.longitude = phi * RAD_TO_DEG
#         self.altitude = r - EARTH_RADIUS

#     # get cartesian coordinate equivalents of the latitude, longitude, and altitude
#     def get_cartesian_coord(self):
#         lmbda = self.latitude * DEG_TO_RAD
#         eta = self.longitude * DEG_TO_RAD
#         x = (EARTH_RADIUS + self.altitude) * np.cos(lmbda) * np.cos(eta)
#         y = (EARTH_RADIUS + self.altitude) * np.cos(lmbda) * np.sin(eta)
#         z = (EARTH_RADIUS + self.altitude) * np.sin(lmbda)

#         return np.array([x, y, z])

#     # set latitude, longitude, altitude from Cartesian coordiantes
#     def set_cartesian_coord(self, x, y, z):
#         pass

#     # set vr, vtheta, vphi from cartesian velocities
#     def set_cartesian_velocity(self, vx, vy, vz):
#         (r, theta, phi) = self.get_spherical_coord()

#         self.vr = vx * np.sin(theta) * np.cos(phi) + vy * np.sin(
#             theta) * np.sin(phi) + vz * np.cos(theta)
#         self.vtheta = (vx * np.cos(theta) * np.cos(phi) + vy * np.cos(theta) *
#                        np.sin(phi) - vz * np.sin(theta)) / r
#         self.vphi = (-vx * np.sin(phi) + vy * np.cos(phi)) / (r *
#                                                               np.sin(theta))

#         # return np.array([vr, vtheta, vphi])

#     # return cartesian velocities from spherical ones
#     def get_cartesian_velocity(self):
#         (r, theta, phi) = self.get_spherical_coord()
#         vx = self.vr * np.sin(theta) * np.cos(phi) + r * self.vtheta * np.cos(
#             theta) * np.cos(phi) - r * self.vphi * np.sin(theta) * np.sin(phi)
#         vy = self.vr * np.sin(theta) * np.sin(phi) + r * self.vtheta * np.cos(
#             theta) * np.sin(phi) + r * self.vphi * np.sin(theta) * np.cos(phi)
#         vz = self.vr * np.cos(theta) - r * self.vtheta * np.sin(theta)

#         return np.array([vx, vy, vz])

#     def __str__(self):
#         return "Latitude: {:.6f}, Longitude: {:.6f}, Altitude: {:.6f}, Velocity (vr, vtheta, vphi): ({:.6f}, {:.6f}, {:.6f})".format(
#             self.latitude, self.longitude, self.altitude, self.vr, self.vtheta,
#             self.vphi)

#     # both position and velocity in spherical coordinates
#     def spherical(self):
#         return np.concatenate((self.get_spherical_coord(),
#                                np.array([self.vr, self.vtheta, self.vphi])),
#                               axis=0)

#     # both position and velocity in Cartesian coordinates
#     def cartesian(self):
#         return np.concatenate(
#             (self.get_cartesian_coord(), self.get_cartesian_velocity()),
#             axis=0)