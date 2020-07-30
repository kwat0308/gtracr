'''
Keeps track of particle trajectory with considerations to cutoffs and E-W effects etc
'''

import os, sys
import numpy as np
# from _trajectorytracer import TrajectoryTracer
from _gtracr import TrajectoryTracer
# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.constants import *
from gtracr.trajectorypoint import TrajectoryPoint
from gtracr.add_particle import particle_dict


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
    - particle_altitude: the altitude in which the cosmic ray hits Earth's atmosphere and creates showers (default 100km)
    - energy: the cosmic ray energy
    - rigidity: the cosmic ray rigidity 
    - escape_altitude: the altitude in which the particle has "escaped" Earth (default 10 * RE)
    '''
    def __init__(self,
                 plabel,
                 latitude,
                 longitude,
                 detector_altitude,
                 zenith_angle,
                 azimuth_angle,
                 particle_altitude=100.,
                 energy=None,
                 rigidity=None,
                 escape_altitude=10. * EARTH_RADIUS):
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
        self.particle_escaped = 0  # check if trajectory is allowed or not

    # evaluates the trajectory using Runge-Kutta methods
    def get_trajectory(self, max_step=10000, step_size=1e-5, get_data=False):

        # get the 6-vector for the detector location
        detector_tp = TrajectoryPoint()
        detector_tp.set_geodesic_coord(self.latitude, self.longitude,
                                       self.detector_altitude)

        # print(detector_tp)
        # get the 6-vector for the particle, initially defined in
        # detector frame, and transform it to geocentric
        # coordinates
        particle_tp = self.detector_to_geocentric(detector_tp)

        # start iteration process

        # initialize the trajectory tracer
        traj_tracer = TrajectoryTracer(self.particle.charge,
                                       self.particle.mass,
                                       self.escape_altitude, step_size,
                                       max_step)

        # get the initial values
        part_t = 0.
        (part_r, part_theta, part_phi, part_pr, part_ptheta,
         part_pphi) = tuple(vars(particle_tp).values())

        initial_values = [
            part_t, part_r, part_theta, part_phi, part_pr, part_ptheta,
            part_pphi
        ]

        if get_data:
            # evaluate the trajectory tracer
            # get data dictionary of the trajectory
            trajectory_datadict = traj_tracer.evaluate_and_get_trajectories(
                initial_values)

            # print(trajectory_datadict)

            # get the final point of the trajectory
            # and make it into a trajectory point
            # not sure if we would use this, but we might...
            particle_final_sixvector = tuple(
                trajectory_datadict.pop("final_values"))

            particle_finaltp = TrajectoryPoint(*particle_final_sixvector)

            # convert all data to numpy arrays for computations etc
            # this should be done within C++ in future versions
            for key, arr in list(trajectory_datadict.items()):
                trajectory_datadict[key] = np.array(arr)

            # lastly get the boolean of if the particle has escaped or not
            # in binary format
            # this helps with the geomagnetic cutoff procedure
            # alternatively this can be inside the geomagnetic things
            self.particle_escaped = int(traj_tracer.particle_escaped)

            return trajectory_datadict

        else:
            # simply evaluate without returning the dictionary
            traj_tracer.evaluate(initial_values)
            # lastly get the boolean of if the particle has escaped or not
            # in binary format
            # this helps with the geomagnetic cutoff procedure
            # alternatively this can be inside the geomagnetic things
            self.particle_escaped = int(traj_tracer.particle_escaped)

            return None

        # print("All done!\n")

        # return trajectory_datadict

    # get the initial trajectory points based on the latitude, longitude, altitude, zenith, and azimuth
    # returns tuple of 2 trajectory points (the initial one and the first one relating to that of the zenith and azimuth one)
    def detector_to_geocentric(self, detector_tp):

        # transformation process for coordinate
        detector_coord = detector_tp.cartesian_coord()

        # change particle initial location if zenith angle is > 90
        # so that we only consider upward moving particles
        if self.zenith_angle > 90.:
            # here we count both altitude and magnitude as a whole
            # for ease of computation
            # particle_altitude = 0.
            # particle_magnitude = -(2. * EARTH_RADIUS + self.particle_altitude
            #                        ) * np.cos(self.zenith_angle * DEG_TO_RAD)
            # particle_magnitude = (2. * EARTH_RADIUS + self.particle_altitude)

            particle_coord = self.get_particle_coord(
                altitude=0.,
                magnitude=-(2. * EARTH_RADIUS + self.particle_altitude) *
                np.cos(self.zenith_angle * DEG_TO_RAD))

        elif self.zenith_angle <= 90.:
            particle_coord = self.get_particle_coord(
                altitude=self.particle_altitude, magnitude=1e-10)

        # print(detector_coord, particle_coord)
        # print(self.tf_matrix())
        (part_x, part_y, part_z) = self.transform(detector_coord,
                                                  particle_coord)

        # print(part_x, part_y, part_z)

        # transformation for momentum
        # need to convert from natural units to SI units
        detector_momentum = np.zeros(3)
        # particle_momentum = self.particle.momentum * KG_M_S_PER_GEVC
        particle_momentum = self.get_particle_coord(
            altitude=0., magnitude=self.particle.momentum * KG_M_S_PER_GEVC)
        # print(detector_momentum, particle_momentum)
        (part_px, part_py, part_pz) = self.transform(detector_momentum,
                                                     particle_momentum)

        # print(part_px, part_py, part_pz)

        # create new trajectory point and set the new coordinate and momentum
        # convert from natural units to SI units
        particle_tp = TrajectoryPoint()
        particle_tp.set_cartesian_coord(part_x, part_y, part_z)
        particle_tp.set_cartesian_momentum(part_px, part_py, part_pz)
        # particle_tp.set_spherical_coord(part_r, part_theta, part_phi)

        # print(particle_tp)

        return particle_tp

    # convert between detector coordinates to geocentric coordinates
    def transform(self, detector_coord, particle_coord):
        return detector_coord + np.dot(self.transform_matrix(), particle_coord)

    # get the detector coordinates (in Cartesian) from zenith and azimuthal angles
    def get_particle_coord(self, altitude, magnitude):
        xi = self.zenith_angle * DEG_TO_RAD
        alpha = self.azimuth_angle * DEG_TO_RAD
        # alpha = (self.azimuth_angle + 90.) * DEG_TO_RAD

        # xt and yt are flipped from usual conversions from spherical coordinates
        # to allow azimuth = 0 to point to the geographic north pole
        # (if we use normal spherical coordinate conversion, azimuth = 0
        #  means pointing west in detector coordinates)

        # xt = magnitude * np.sin(xi) * np.sin(alpha)
        # yt = magnitude * np.sin(xi) * np.cos(alpha)
        # zt = magnitude * np.cos(xi) + altitude

        xt = -magnitude * np.sin(xi) * np.sin(alpha)
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