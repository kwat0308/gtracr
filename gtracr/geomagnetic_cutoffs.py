import sys
import os
import numpy as np
from scipy.interpolate import griddata
from tqdm import tqdm
from p_tqdm import p_map
from datetime import date

from gtracr.trajectory import Trajectory
from gtracr.utils import location_dict

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)


class GMRC():
    '''
    Evaluates the geomagnetic cutoff rigidities associated to a specific location on the globe for each zenith and azimuthal angle (a zenith angle > 90 degrees are for upward-moving particles, that is, for cosmic rays coming from the other side of Earth).

    The cutoff rigidities are evaluated using a Monte-Carlo sampling scheme, combined with a 2-dimensional linear interpolation using `scipy.interpolate`. 

    The resulting cutoffs can be plotted as 2-dimensional heatmap.

    Parameters
    -----------

    - location : str
        The location in which the geomagnetic cutoff rigidities are evaluated (default = "Kamioka"). 
        The names must be one of the locations contained in `location_dict`, which is configured in `gtracr.utils`.
        If set to None, location may be set manually with latitude, longitude
    - latitude: float
        the geographic latitude of the detector, with 0 defined at the equator in degrees. Overridden
        if the location name is specified.
    - longitude: float
        the geographic longitude of the detector, with 0 defined at the Prime Meridian in degrees. Overridden
        if the location name is specified.
    - particle_altitude : float
        The altitude in which the cosmic ray interacts with the atmosphere in km (default = 100).
    - detector_altitude: float
        the height of the detector from sea level in km (default = 0km). Overridden if location is specified
    - iter_num : int
        The number of iterations to perform for the Monte-Carlo sampling routine (default = 10000) 
    - bfield_type : str
        The type of magnetic field model to use for the evaluation of the cutoff rigidities (default = "igrf"). Set to "dipole" to use the dipole approximation of the geomagnetic field instead.
    - particle_type : str
        The type of particle of the cosmic ray (default  ="p+").
    - date : str
        The specific date in which the geomagnetic rigidity cutoffs are evaluated. Defaults to the current date.
    - min_rigidity : float
        The minimum rigidity to which we evaluate the cutoff rigidities for (default = 5 GV).
    - max_rigidity : float
        The maximum rigidity to which we evaluate the cutoff rigidities for (default = 55 GV).
    - delta_rigidity : float
        The spacing between each rigidity (default = 5 GV). Sets the coarseness of the rigidity sample space.
    - dt : float
        The stepsize of each trajectory evaluation (default = 1e-5)
    - max_time : float
        The maximal time of each trajectory evaluation (default = 1.).
    '''
    def __init__(self,
                 location="Kamioka",
                 latitude=0.,
                 longitude=0.,
                 particle_altitude=100.,
                 detector_altitude=0.,
                 iter_num=10000,
                 bfield_type="igrf",
                 particle_type="p+",
                 date=str(date.today()),
                 min_rigidity=5.,
                 max_rigidity=55.,
                 delta_rigidity=1.,
                 dt=1e-5,
                 max_time=1,
                 method='serial'):
        # set class attributes
        self.location = location

        # only import location dictionary and use those values if location is not None
        if location is not None:
            if location in location_dict:
                loc = location_dict[location]

                latitude = loc.latitude
                longitude = loc.longitude
                detector_altitude = loc.altitude

        self.lat = latitude
        self.lon = longitude
        self.dalt = detector_altitude
        self.palt = particle_altitude
        self.iter_num = iter_num
        self.bfield_type = bfield_type
        self.plabel = particle_type
        self.date = date
        '''
        Rigidity configurations
        '''
        self.rmin = min_rigidity
        self.rmax = max_rigidity
        self.rdelta = delta_rigidity
        self.dt = dt
        self.max_time = max_time
        self.method = method

        # generate list of rigidities
        self.rigidity_list = np.arange(self.rmin, self.rmax, self.rdelta)

        # initialize container for rigidity cutoffs
        # along with the zenith and azimuthal arrays
        self.data_dict = {
            "azimuth": np.zeros(self.iter_num),
            "zenith": np.zeros(self.iter_num),
            "rcutoff": np.zeros(self.iter_num)
        }

    def evaluate_angle(self, azimuth, zenith):
        # iterate through each rigidity, and break the loop
        # when particle is able to escape earth
        for rigidity in self.rigidity_list:

            traj = Trajectory(plabel=self.plabel,
                              location_name=None,
                              latitude = self.lat,
                              longitude = self.lon,
                              zenith_angle=zenith,
                              azimuth_angle=azimuth,
                              particle_altitude=self.palt,
                              detector_altitude=self.dalt,
                              rigidity=rigidity,
                              bfield_type=self.bfield_type,
                              date=self.date)

            traj.get_trajectory(dt=self.dt, max_time=self.max_time)

            # break loop and return current rigidity if particle has escaped
            if traj.particle_escaped == True:
                return rigidity

        # didn't escape through all the rigidities in the list, return None
        return None


    def evaluate_serial(self):
        '''
        Evaluate the rigidity cutoff value at some provided location
        on Earth for a given cosmic ray particle.
        '''

        # perform Monte Carlo integration to get cutoff rigidity
        for i in tqdm(range(self.iter_num)):
            # get a random zenith and azimuth angle
            # zenith angles range from 0 to 180
            # azimuth angles range from 0 to 360
            [azimuth, zenith] = np.random.rand(2)
            azimuth *= 360.
            zenith *= 180.

            rigidity = self.evaluate_angle(azimuth, zenith)

            if rigidity:
                self.data_dict["azimuth"][i] = azimuth
                self.data_dict["zenith"][i] = zenith
                self.data_dict["rcutoff"][i] = rigidity

    def evaluate_parallel(self):
        # generate lists of random zenith and azimuth angles
        azimuth = np.random.rand(self.iter_num) * 360.0
        zenith = np.random.rand(self.iter_num) * 180.0

        rigidity = p_map(self.evaluate_angle, azimuth, zenith)

        # insert the non-None's into the data_dict
        for i in range(self.iter_num):
            if rigidity[i]:
                self.data_dict['azimuth'][i] = azimuth[i]
                self.data_dict['zenith'][i] = zenith[i]
                self.data_dict['rcutoff'][i] = rigidity[i]


    def evaluate(self):
        if self.method == 'serial':
            self.evaluate_serial()
        elif self.method == 'parallel':
            self.evaluate_parallel()


    def interpolate_results(self,
                            method="linear",
                            ngrid_azimuth=70,
                            ngrid_zenith=70):
        '''
        Interpolate the rigidity cutoffs using `scipy.interpolate.griddata`

        Parameters
        ----------
        - method : str
            The type of linear interpolation used for `griddata` (default = "linear"). Choices are between "nearest", "linear", and "cubic".
        - ngrid_azimuth, ngrid_zenith : int
            The number of grids for the azimuth and zenith angles used for the interpolation (default = 70).

        Returns
        --------

        Returns a tuple of the following objects:

        - azimuth_grid : np.array(float), size ngrid_azimuth
            The linearly spaced values of the azimuthal angle
        - zenith_grid : np.array(float), size ngrid_zenith
            The linearly spaced values of the zenith angle
        - rcutoff_grid : np.array(float), size ngrid_azimuth x ngrid_zenith
            The interpolated geomagnetic cutoff rigidities.
        '''

        azimuth_grid = np.linspace(np.min(self.data_dict["azimuth"]),
                                   np.max(self.data_dict["azimuth"]),
                                   ngrid_azimuth)
        zenith_grid = np.linspace(np.max(self.data_dict["zenith"]),
                                  np.min(self.data_dict["zenith"]),
                                  ngrid_zenith)

        rcutoff_grid = griddata(points=(self.data_dict["azimuth"],
                                        self.data_dict["zenith"]),
                                values=self.data_dict["rcutoff"],
                                xi=(azimuth_grid[None, :], zenith_grid[:,
                                                                       None]),
                                method=method)

        return (azimuth_grid, zenith_grid, rcutoff_grid)


