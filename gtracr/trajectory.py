'''
Keeps track of particle trajectory with considerations to cutoffs and E-W effects etc
'''

# from add_particle import particle_dict
from gtracr.utils import get_particledict, get_locationdict
from gtracr.lib.trajectorypoint import TrajectoryPoint
from gtracr.lib.constants import EARTH_RADIUS, DEG_PER_RAD, RAD_PER_DEG, KG_M_S_PER_GEVC
from gtracr.lib.trajectory_tracer import pTrajectoryTracer

import os
import sys
import numpy as np
import pickle
from datetime import datetime as dt
from gtracr.lib._libgtracr import TrajectoryTracer, uTrajectoryTracer

# from trajectory_tracer import TrajectoryTracer
# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(CURRENT_DIR, "data")


class Trajectory:
    '''
    Class that controls the trajectory of a particle at some given energy / rigidity

    Members
    -------

    - plabel: the label of the particle defined in particle_dict
    - latitude: the geographic latitude of the detector, with 0 defined at the equator in degrees
    - longitude: the geographic longitude of the detector, with 0 defined at the Prime Meridian in degrees
    - detector_altitude: the height of the detector from sea level (0=sea level) in km
    - zenith_angle: the angle of the cosmic ray trajectory from the local zenith, with 0 being at the local zenith
    - azimuth_angle: the angle of the cosmic ray trajectory with 0 being in the direction of the geographic North in the local tangent plane
    - particle_altitude: the altitude in which the cosmic ray hits Earth's atmosphere and creates showers (default 100km)
    - energy: the cosmic ray energy
    - rigidity: the cosmic ray rigidity 
    - escape_altitude: the altitude in which the particle has "escaped" Earth (default 10 * RE)
    - bfield_type: the type of bfield to evaluate the trajectory with (either 'dipole' or 'igrf', default: dipole)
    - date: the date in which the field is evaluated in.
    '''
    def __init__(self,
                 plabel,
                 zenith_angle,
                 azimuth_angle,
                 particle_altitude=100.,
                 latitude=0.,
                 longitude=0.,
                 detector_altitude=0.,
                 location_name=None,
                 energy=None,
                 rigidity=None,
                 escape_altitude=10. * EARTH_RADIUS,
                 bfield_type="dipole",
                 date=dt.now().year):
        self.zenith_angle = zenith_angle
        self.azimuth_angle = azimuth_angle
        self.particle_altitude = particle_altitude * (1e3)  # convert to meters
        self.escape_altitude = escape_altitude

        # define particle from particle_dict
        particle_dict = get_particledict()
        self.particle = particle_dict[plabel]

        # only import location dictionary and use those values if location_name is not None
        if location_name is not None:
            location_dict = get_locationdict()
            loc = location_dict[location_name]

            latitude = loc.latitude
            longitude = loc.longitude
            detector_altitude = loc.altitude

        self.latitude = latitude
        self.longitude = longitude
        self.detector_altitude = detector_altitude * (1e3)  # convert to meters
        # define rigidity and energy only if they are provided, evaluate for the other member
        # also set momentum in each case
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

        self.particle_escaped = False  # check if trajectory is allowed or not
        # type of bfield to use
        # take only first character for compatibility with char in c++
        self.bfield_type = bfield_type[0]

        # find the path to the data and set current date for igrf bfield
        datapath = os.path.abspath(os.path.join(CURRENT_DIR, "data"))
        # print(datapath)
        #         curr_year = dt.now(
        #         ).year  # should be in decimal years with mm/dd implemented in future
        # print(curr_year)
        self.igrf_params = (datapath, date)

        # final time and six-vector, used for testing purposes
        self.final_time = 0.
        self.final_sixvector = np.zeros(6)

        # get the 6-vector for the detector location
        detector_tp = TrajectoryPoint()
        detector_tp.set_geodesic_coord(self.latitude, self.longitude,
                                       self.detector_altitude)

        # print(detector_tp)
        # get the 6-vector for the particle, initially defined in
        # detector frame, and transform it to geocentric
        # coordinates
        self.particle_tp = self.detector_to_geocentric(detector_tp)

    def get_trajectory(self,
                       dt=1e-5,
                       max_time=1,
                       max_step=None,
                       get_data=False,
                       use_python=False,
                       use_unvectorized=False):
        '''
        Evaluate the trajectory of the particle within Earth's magnetic field
        and determines whether particle has escaped or not.
        Optionally also returns the information of the trajectory (the duration
        and the six-vector in spherical coordinates) if `get_data == True`.

        Parameters
        ----------

        dt : float
            the time step between each iteration of the integration (default: 1e-5)
        max_time : float
            the maximum duration in which the integration would occur in seconds (default: 10)
        max_step : int, optional
            maximum number of steps to integrate for (default None). If `max_step` is not `None`,
            then `max_step` will override the evaluation of maximum number of steps based on `max_time`.
        get_data : bool, optional
            decides whether we want to extract the information (time and six vector)
            for the whole trajectory for e.g. debugging purposes (default: False)
        use_python : bool, optional
            decides whether to use the python implementation for the TrajectoryTracer class instead of
            that implemented in C++. This is mainly enabled for debugging purposes (default: False)
        use_unvectorized : bool, optional
            decides whether to evaluate the Runge Kutta integration in the C++ version in its 
            unvectorized or vectorized form. This is mainly enabled for debugging purposes (default: False)

        Returns
        ---------

        - trajdata_dict : dict
            a dictionary that contains the information of the whole trajectory in
            spherical coordinates.
            Keys are ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]
            - only returned when `get_data` is True
        '''

        # evaluate max_step only when max_time is given, else use the user-given
        # max step
        max_step = int(np.ceil(max_time /
                               dt)) if max_step is None else max_step

        # raise issues if both use python and use unvectorized form is True
        if use_python and use_unvectorized:
            raise Exception("Unvectorized Python version does not exist!")

        # start iteration process

        # initialize trajectory tracer
        if use_python:
            # the python trajectory tracer version
            traj_tracer = pTrajectoryTracer(self.particle.charge,
                                            self.particle.mass,
                                            self.escape_altitude, dt, max_step,
                                            self.bfield_type, self.igrf_params)
        elif use_unvectorized:
            # the unvectorized trajectory tracer version
            traj_tracer = uTrajectoryTracer(self.particle.charge,
                                            self.particle.mass,
                                            self.escape_altitude, dt, max_step,
                                            self.bfield_type, self.igrf_params)
        else:
            # the vectorized trajectory tracer version
            traj_tracer = TrajectoryTracer(self.particle.charge,
                                           self.particle.mass,
                                           self.escape_altitude, dt, max_step,
                                           self.bfield_type, self.igrf_params)

        # set initial values
        particle_t0 = 0.
        particle_vec0 = self.particle_tp.asarray()

        if get_data:
            # evaluate the trajectory tracer
            # get data dictionary of the trajectory
            trajectory_datadict = traj_tracer.evaluate_and_get_trajectory(
                particle_t0, particle_vec0)

            # convert all data to numpy arrays for computations etc
            # this should be done within C++ in future versions
            for key, arr in list(trajectory_datadict.items()):
                trajectory_datadict[key] = np.array(arr)

            # lastly get the boolean of if the particle has escaped or not
            # in binary format
            # this helps with the geomagnetic cutoff procedure
            # alternatively this can be inside the geomagnetic things
            self.particle_escaped = traj_tracer.particle_escaped

            # set the final time and six-vector from the evaluator
            self.final_time = traj_tracer.final_time
            self.final_sixvector = np.array(traj_tracer.final_sixvector)

            return trajectory_datadict

        else:
            # simply evaluate without returning the dictionary
            traj_tracer.evaluate(particle_t0, particle_vec0)
            # lastly get the boolean of if the particle has escaped or not
            # in binary format
            # this helps with the geomagnetic cutoff procedure
            # alternatively this can be inside the geomagnetic things
            self.particle_escaped = traj_tracer.particle_escaped

            # set the final time and six-vector from the evaluator
            self.final_time = traj_tracer.final_time
            self.final_sixvector = np.array(traj_tracer.final_sixvector)

            return None

    # get the initial trajectory points based on the latitude, longitude, altitude, zenith, and azimuth
    # returns tuple of 2 trajectory points (the initial one and the first one relating to that of the zenith and azimuth one)
    def detector_to_geocentric(self, detector_tp):
        '''
        Convert the coordinates defined in the detector frame (the coordinate system 
        defined in the local tangent plane to Earth's surface at some specified
        latitude and longitude) to geocentric (Earth-centered, Earth-fixed) coordinates.

        Parameters
        -----------
        detector_tp : TrajectoryPoint
            The six vector of the detector location, defined as a TrajectoryPoint
            This is defined based on the geodesic coordinates (latitude, longitude, altitude).

        Returns
        -------
        particle_tp : TrajectoryPoint
            The six vector of the particle defined as a TrajectoryPoint, evaluated
            based on the location of the detector, the direction in which the particle
            comes from, and the altitude in which a shower occurs.

        '''

        # transformation process for coordinate
        detector_coord = detector_tp.cartesian_coord()

        # change particle initial location if zenith angle is > 90
        # so that we only consider upward moving particles
        if self.zenith_angle > 90.:
            # here we count both altitude and magnitude as a whole
            # for ease of computation

            particle_coord = self.get_particle_coord(
                altitude=0.,
                magnitude=-(2. * EARTH_RADIUS + self.particle_altitude) *
                np.cos(self.zenith_angle * DEG_PER_RAD))

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
        xi = self.zenith_angle * DEG_PER_RAD
        alpha = self.azimuth_angle * DEG_PER_RAD
        # alpha = (self.azimuth_angle + 90.) * DEG_TO_RAD

        # xt and yt are flipped from usual conversions from spherical coordinates
        # to allow azimuth = 0 to point to the geographic north pole
        # (if we use normal spherical coordinate conversion, azimuth = 0
        #  means pointing west in detector coordinates)

        # the below coordinate convention for detector coordinates are used
        # to follow the convention used in Honda's 2002 paper, in which
        # an azimuth angle of 0 correlates to direction of the
        # geographic South Pole, and + azimuth will indicate particles
        # coming from the west
        xt = magnitude * np.sin(xi) * np.sin(alpha)
        yt = -magnitude * np.sin(xi) * np.cos(alpha)
        zt = magnitude * np.cos(xi) + altitude

        return np.array([xt, yt, zt])

    # the transformation matrix from the detector to geocentric (Cartesian) coordinates
    # source: http://archive.psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf
    def transform_matrix(self):
        lmbda = self.latitude * DEG_PER_RAD
        eta = self.longitude * DEG_PER_RAD

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
