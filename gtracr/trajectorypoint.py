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
    - pr: the radial momentum
    - ptheta: the momentum in the polar direction
    - pphi: the momentum in the azimuthal direction
    '''
    def __init__(self,
                 r=EARTH_RADIUS,
                 theta=0.,
                 phi=0.,
                 pr=0.,
                 ptheta=0.,
                 pphi=0.):
        self.r = r
        self.theta = theta
        self.phi = phi
        self.pr = pr
        self.ptheta = ptheta
        self.pphi = pphi

    # get geodesic coordinate equivalents of spherical ones
    def geodesic_coordinate(self):
        latitude = 90. - (self.theta * RAD_TO_DEG)
        # longitude = (phi * RAD_TO_DEG) - 180.
        longitude = self.phi * RAD_TO_DEG
        altitude = self.r - EARTH_RADIUS
        return np.array([latitude, longitude, altitude])

    # set latitude, longitude, altitude from geodesic coordiantes
    def set_geodesic_coord(self, latitude, longitude, altitude):
        self.r = EARTH_RADIUS + altitude
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

    # set pr, ptheta, pphi from cartesian momentum
    def set_cartesian_momentum(self, px, py, pz):
        self.pr = px * np.sin(self.theta) * np.cos(self.phi) + py * np.sin(
            self.theta) * np.sin(self.phi) + pz * np.cos(self.theta)
        self.ptheta = (px * np.cos(self.theta) * np.cos(self.phi) +
                       py * np.cos(self.theta) * np.sin(self.phi) -
                       pz * np.sin(self.theta))
        self.pphi = (-px * np.sin(self.phi) + py * np.cos(self.phi))

        # return np.array([vr, vtheta, vphi])

    # return cartesian momenta from spherical ones
    def cartesian_momentum(self):
        # px = self.pr * np.sin(self.theta) * np.cos(
        #     self.phi) + self.r * self.ptheta * np.cos(self.theta) * np.cos(
        #         self.phi) - self.r * self.pphi * np.sin(self.theta) * np.sin(
        #             self.phi)
        # py = self.pr * np.sin(self.theta) * np.sin(
        #     self.phi) + self.r * self.ptheta * np.cos(self.theta) * np.sin(
        #         self.phi) + self.r * self.pphi * np.sin(self.theta) * np.cos(
        #             self.phi)
        # pz = self.pr * np.cos(self.theta) - self.r * self.ptheta * np.sin(
        #     self.theta)

        tfmat_sph2car = np.array([[
            [
                np.sin(self.theta) * np.cos(self.phi),
                self.r * np.cos(self.theta) * np.cos(self.phi),
                -self.r * np.sin(self.theta) * np.sin(self.phi)
            ],
            [
                np.sin(self.theta) * np.sin(self.phi),
                self.r * np.cos(self.theta) * np.sin(self.phi),
                self.r * np.sin(self.theta) * np.cos(self.phi),
            ],
            [np.cos(self.theta), -self.r * np.sin(self.theta), 0.],
        ]])

        cart_momentum_arr = np.dot(tfmat_sph2car,
                                   np.array([self.pr, self.ptheta, self.pphi]))

        return cart_momentum_arr

    def __str__(self):
        return "Coordinate (r, theta, phi): ({:.5e}, {:.5e}, {:.5e}), \
            \nMomentum (pr, ptheta, pphi): ({:.5e}, {:.5e}, {:.5e})\n".format(
            *tuple(vars(self).values()))

    # both position and momentum in Cartesian coordinates
    def cartesian(self):
        return np.concatenate(
            (self.cartesian_coord(), self.cartesian_momentum()), axis=0)
