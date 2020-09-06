import gtracr.lib.igrf_utils as iuf
from gtracr.lib.constants import EARTH_RADIUS, G10
'''
Library that controls the equation for the Earth's magnetic field.
'''

import os
import sys
import os.path as p
import numpy as np
# import csv
from scipy.interpolate import interp1d

# sys.path.append(os.getcwd())
# sys.path.append(p.join(os.getcwd(), "gtracr"))

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
# sys.path.append(PARENT_DIR)


class MagneticField:
    '''
    The base class that contains the properties of 
    Earth's geomagnetic field defined in spherical coordinates.

    The base class uses the ideal dipole approximation
    in which the first-order approximation of the series expansion
    of the magnetic potential of Earth's geomagnetic field
    is used. 

    Notes
    ------
    - The Earth's magnetic field is assumed to consist of 
    no external currents, i.e. curl(B) = 0

    '''
    def __init__(self):
        pass

    def values(self, r, theta, phi):
        '''
        The values of the magnetic field in spherical 
        coordinates at a given location (r, theta, phi)

        Parameters
        ----------
        - r : the radial component in meters
        - theta: the polar component in radians
        - phi : the phi component in radians

        Returns
        -------
        - Br, Btheta, Bphi: The spherical components
                            of the magnetic field
        '''

        Br = 2. * (EARTH_RADIUS / r)**3. * G10 * np.cos(theta)
        Btheta = (EARTH_RADIUS / r)**3. * G10 * np.sin(theta)
        Bphi = 0.

        return np.array([Br, Btheta, Bphi])

    # def Br(self, r, theta, phi):
    #     return 2. * (EARTH_RADIUS / r)**3. * G10 * np.cos(theta)

    # def Btheta(self, r, theta, phi):
    #     return (EARTH_RADIUS / r)**3. * G10 * np.sin(theta)

    # def Bphi(self, r, theta, phi):
    #     return 0.


class IGRF13(MagneticField):
    '''
    Derived class from MagneticField class that 
    describes Earth's geomagnetic field using the 
    International Geomagnetic Reference Field (IGRF)
    of the 13th Edition (2020).

    The scripts that produces the IGRF magnetic field
    are obtained from the Python scripts that are 
    officially constructed by the V-MOD group in the 
    International Union of Geodesy and Geophysics.

    More information about the IGRF-13 model can be obtained 
    here: https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html 

    Members
    -------
    - curr_year: the current year in years
    - igrf_coeffs : the Gauss coefficients for the year 
                    specified in curr_year.
    - nmax : the number of truncation

    Initialization Parameters
    -------------------------
    - curr_year: the current year in years. This is required as the 
            magnetic field values vary in a linear fashion 
            per year.

    - nmax (optional): the number to truncate the series expansion of
                        the magnetic field. The default value is obtained
                        from the .shc file.

    '''
    def __init__(self, curr_year, nmax=None):
        # override MagneticField __init__ function
        # not necessary, but for decorative sake
        # super().__init__(self)

        self.curr_year = curr_year
        fpath = os.path.join(PARENT_DIR, "data", "IGRF13.shc")

        leap_year = False

        # check if current year is a leap year or not
        if self.curr_year % 4 == 0:
            leap_year = True

        # import the coefficients from .shc file
        igrf = iuf.load_shcfile(fpath, leap_year=leap_year)

        # interpolate to account for any year
        # linear interpolation since time variation in
        # coefficients is a linear one
        interp_coeffs = interp1d(
            igrf.time,
            igrf.coeffs,
            kind="linear",
        )

        # get the number of truncation based on the
        # igrf .shc file
        self.nmax = igrf.parameters["nmax"] if nmax is None else nmax

        self.igrf_coeffs = interp_coeffs(self.curr_year).T

    def values(self, r, theta, phi):
        # obtaining the values can easily be done by using
        # the function synth_values from igrf_utils.py
        Br, Btheta, Bphi = iuf.synth_values(self.igrf_coeffs,
                                            r,
                                            theta,
                                            phi,
                                            nmax=self.nmax)

        return np.array([Br, Btheta, Bphi])

    # def get_igrfcoeffs(self, fpath):
    #     '''
    #     Get the IGRF Gauss coeffiecients for all years
    #     by performing an interpolation

    #     Parameters
    #     ----------
    #     - fpath: the file path in which the .shc file is in

    #     Returns
    #     -------
    #     - coeffs_curr : the Gauss coefficients for the current year
    #     - nmax : the number of truncation for the series expansion
    #             of the IGRF model. Obtained from the .shc file.
    #     '''
    #     # check if current year is a leap year or not
    #     if self.curr_year % 4 == 0:
    #         leap_year = True

    #     # import the coefficients from .shc file
    #     igrf = iuf.load_shcfile(fpath, leap_year=leap_year)

    #     # interpolate to account for any year
    #     # linear interpolation since time variation in
    #     # coefficients is a linear one
    #     interp_coeffs = interp1d(
    #         igrf.time,
    #         igrf.coeffs,
    #         kind="linear",
    #     )

    #     # get the number of truncation based on the
    #     # igrf .shc file
    #     nmax = igrf.parameters["nmax"]

    #     coeffs_curr = interp_coeffs(self.curr_year).T

    #     return coeffs_curr, nmax


# if __name__ == "__main__":
#     bfield = MagneticField()
#     bfield.import_coeffs()

#     print(bfield.gcoeffs, bfield.hcoeffs)
