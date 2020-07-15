'''
Keeps track of particle trajectory with considerations to cutoffs and E-W effects etc
'''

import os, sys
import numpy as np

# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.constants import *
from gtracr.trajectory_point import TrajectoryPoint
from RungeKutta import RungeKutta
from gtracr.add_particle import particle_dict

KEY_LIST = ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]


class Trajectory:
    '''
    Class that controls the trajectory of a particle at some given energy / rigidity
    Members:
    - plabel: the label of the particle defined in particle_dict in add_particle.py
    - latitude: the geographic latitude of the detector, with 0 defined at the equator in degrees
    - longitude: the geographic longitude of the detector, with 0 defined at the Prime Meridian in degrees
    - detector_altitude: the height of the detector from sea level (0=sea level) in km
    - zenith_angle: the angle of the cosmic ray trajectory from the local zenith, with 0 being at the local zenith
    - azimuth_angle: the angle of the cosmic ray trajectory with 0 being in the direction of the geographic North in the local tangent plane
    - particle_altitude: the altitude in which the cosmic ray hits Earth's atmosphere and creates showers
    - energy: the cosmic ray energy
    - rigidity: the cosmic ray rigidity 
    - escape_altitude: the altitude in which the particle has "escaped" Earth (default 10 * RE)
    - max_buffer: maximum length of array 
    '''
    def __init__(self,
                 plabel,
                 latitude,
                 longitude,
                 detector_altitude,
                 zenith_angle,
                 azimuth_angle,
                 particle_altitude,
                 energy=None,
                 rigidity=None,
                 escape_altitude=10. * EARTH_RADIUS,
                 max_buffer=10000):
        self.particle = particle_dict[plabel]
        self.latitude = latitude
        self.longitude = longitude
        self.detector_altitude = detector_altitude * (1e3)  # convert to meters
        self.zenith_angle = zenith_angle
        self.azimuth_angle = azimuth_angle
        self.particle_altitude = particle_altitude * (1e3)  # convert to meters
        self.escape_altitude = escape_altitude

        # define rigidity and energy only if they are provided, evaluate for the other member
        # also set momentum in each case
        # self.particle.print()

        if rigidity is None:
            self.particle.set_from_energy(energy)
            self.rigidity = self.particle.rigidity
            self.energy = energy
        elif energy is None:
            self.particle.set_from_rigidity(rigidity)
            self.rigidity = rigidity
            self.energy = self.particle.get_energy_rigidity()
        # elif rigidity is None and energy is None:
        else:
            raise Exception(
                "Provide either energy or rigidity as input, not both!")

        # self.particle.print()
        self.particle_escaped = False  # check if trajectory is allowed or not

        # initialize required arrays here
        self.max_buffer = max_buffer
        # self.time_array = np.zeros(max_buffer)
        # self.tp_array = np.zeros(max_buffer)
        self.time_array = [None] * max_buffer
        self.tp_array = [None] * max_buffer  # list to append TJP objects

    # get the initial trajectory points based on the latitude, longitude, altitude, zenith, and azimuth
    # returns tuple of 2 trajectory points (the initial one and the first one relating to that of the zenith and azimuth one)
    def detector_to_geocentric(self, detector_tp):

        # transformation process for coordinate
        detector_coord = detector_tp.cartesian_coord()
        particle_coord = self.get_particle_coord(
            altitude=self.particle_altitude, magnitude=1e-10)
        print(detector_coord, particle_coord)
        # print(self.tf_matrix())
        (part_x, part_y, part_z) = self.transform(detector_coord,
                                                  particle_coord)

        # print(part_x, part_y, part_z)

        # transformation for velocity
        detector_momentum = np.zeros(3)
        particle_momentum = self.particle.momentum * KG_M_S_PER_GEVC  # convert from natural units to SI units
        particle_momentum = self.get_particle_coord(
            altitude=0., magnitude=particle_momentum)
        print(detector_momentum, particle_momentum)
        (part_px, part_py, part_pz) = self.transform(detector_momentum,
                                                     particle_momentum)

        # print(part_px, part_py, part_pz)

        # create new trajectory point and set the new coordinate and velocity
        particle_tp = TrajectoryPoint()
        particle_tp.set_cartesian_coord(part_x, part_y, part_z)
        particle_tp.set_cartesian_momentum(part_px, part_py, part_pz)
        # particle_tp.set_spherical_coord(part_r, part_theta, part_phi)

        # print(particle_tp)

        return particle_tp

    # evaluates the trajectory using Runge-Kutta methods
    def get_trajectory(self, max_step=10000, step_size=1e-5):

        # check if max_step > max_buffer, if so then update this
        # there is a better way to do this, im sure
        if max_step > self.max_buffer:
            self.time_array.extend([None] * ((max_step) - self.max_buffer))
            self.tp_array.extend([None] * ((max_step) - self.max_buffer))
        elif max_step < self.max_buffer:
            self.time_array = self.time_array[:max_step]
            self.tp_array = self.tp_array[:max_step]

        # get the initial trajectory points
        detector_tp = TrajectoryPoint()
        detector_tp.set_geodesic_coord(self.latitude, self.longitude,
                                       self.detector_altitude)

        # print(detector_tp)

        particle_tp = self.detector_to_geocentric(detector_tp)

        # append both time array and TJP
        # self.time_array[0:2] = [0., step_size]
        # self.tp_array[0:2] = [detector_tp, particle_tp]
        self.time_array[0] = 0.
        self.tp_array[0] = particle_tp

        # start iteration process
        # runge kutta integrator
        rk_integrator = RungeKutta(self.particle.charge, self.particle.mass,
                                   step_size)

        # print(rk_integrator.charge, rk_integrator.mass)
        i = 1
        part_t = step_size
        (part_r, part_theta, part_phi, part_pr, part_ptheta,
         part_pphi) = tuple(vars(particle_tp).values())

        print(part_r, part_theta, part_phi, part_pr, part_ptheta, part_pphi)
        while i < max_step:
            [
                part_t, part_r, part_theta, part_phi, part_pr, part_ptheta,
                part_pphi
            ] = rk_integrator.evaluate([
                part_t, part_r, part_theta, part_phi, part_pr, part_ptheta,
                part_pphi
            ])

            # print(part_t, part_r, part_theta, part_phi, part_pr, part_ptheta,
            #       part_pphi, '\n')

            next_tp = TrajectoryPoint(r=part_r,
                                      theta=part_theta,
                                      phi=part_phi,
                                      pr=part_pr,
                                      ptheta=part_ptheta,
                                      pphi=part_pphi)

            self.time_array[i] = part_t
            self.tp_array[i] = next_tp

            # conditions

            if next_tp.r > EARTH_RADIUS + self.escape_altitude:
                print("Allowed Trajectory!")
                self.particle_escaped = True
                self.time_array = self.time_array[:i]
                self.tp_array = self.tp_array[:i]
                break

            if next_tp.r < EARTH_RADIUS:
                print("Forbidden Trajectory!")
                self.time_array = self.time_array[:i]
                self.tp_array = self.tp_array[:i]
                break

            # some looping checker
            if (i - 2) % (max_step // 10) == 0 and (i - 2) != 0:
                print("{0} iterations completed".format(i - 2))

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
            "px": np.zeros(n),
            "py": np.zeros(n),
            "pz": np.zeros(n)
        }

        for i, tp in enumerate(self.tp_array):
            # print(tp)
            (x, y, z, px, py, pz) = tp.cartesian()

            data_dict["x"][i] = x
            data_dict["y"][i] = y
            data_dict["z"][i] = z
            data_dict["px"][i] = px
            data_dict["py"][i] = py
            data_dict["pz"][i] = pz

        return data_dict

    # convert between detector coordinates to geocentric coordinates
    def transform(self, detector_coord, particle_coord):
        return detector_coord + np.dot(self.transform_matrix(), particle_coord)

    # get the detector coordinates (in Cartesian) from zenith and azimuthal angles
    def get_particle_coord(self, altitude, magnitude):
        xi = self.zenith_angle * DEG_TO_RAD
        alpha = self.azimuth_angle * DEG_TO_RAD

        # xt and yt are flipped from usual conversions from spherical coordinates
        # to allow azimuth = 0 to point to the geographic north pole
        # (if we use normal spherical coordinate conversion, azimuth = 0
        #  means pointing west in detector coordinates)

        xt = magnitude * np.sin(xi) * np.sin(alpha)
        yt = magnitude * np.sin(xi) * np.cos(alpha)
        zt = magnitude * np.cos(xi) + altitude

        return np.array([xt, yt, zt])

    # the transformation matrix from the detector to geocentric (Cartesian) coordinates
    # source: http://archive.psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf
    def transform_matrix(self):
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