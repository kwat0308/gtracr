'''
Library that controls the equation for the Earth's magnetic field.
'''

import os, sys
import os.path as p
import numpy as np
# import csv

sys.path.append(os.getcwd())
sys.path.append(p.join(os.getcwd(), "gtracr"))

from gtracr.constants import EARTH_RADIUS, G10


class MagneticField:
    def __init__(self):
        pass

    # def __init__(self):
    #     self.gcoeffs = np.zeros((13,13))  # N is truncation number
    #     self.hcoeffs = np.zeros((13,13))

    # # import the IGRF gauss coefficients
    # def import_coeffs(self):
    #     # np.genfromtxt(p.join(os.getcwd(), "data", "IGRF13coeffs.csv"), dtype=float, skip_header=4, )
    #     fname = p.join(os.getcwd(), "data", "IGRF13coeffs.csv")
    #     with open(fname, newline='') as csvfile:
    #         reader = csv.reader(csvfile)
    #         # get the coefficients here
    #         # skip 4 rows / headers
    #         for row in reader:
    #             if 'g' == row[0]:
    #                 n = int(row[1])
    #                 m = int(row[2])
    #                 self.gcoeffs[n][m] = row[4]
    #                 print(row[4])
    #             elif 'h' == row[0]:
    #                 n = int(row[1])
    #                 m = int(row[2])
    #                 self.hcoeffs[n][m] = row[4]
    #                 print(row[4])
    #             else:
    #                 pass

    # print(self.gcoeffs, self.hcoeffs)

    # below are the IDEAL DIPOLE APPROX.
    def Br(self, r, theta, phi):
        return 2. * (EARTH_RADIUS / r)**3. * G10 * np.cos(theta)

    def Btheta(self, r, theta, phi):
        return (EARTH_RADIUS / r)**3. * G10 * np.sin(theta)

    def Bphi(self, r, theta, phi):
        return 0.


# # magnetic field (ideal dipole)
# def B_r(r, theta, phi):
#     return 2. * (EARTH_RADIUS / r)**3. * g10 * np.cos(theta)
#     # return 2.*(1./r)**3.*g10*np.cos(theta)

# #
# def B_theta(r, theta, phi):
#     return (EARTH_RADIUS / r)**3. * g10 * np.sin(theta)
#     # return (1./r)**3.*g10*np.sin(theta)

# def B_phi(r, theta, phi):
#     return 0.

# if __name__ == "__main__":
#     bfield = MagneticField()
#     bfield.import_coeffs()

#     print(bfield.gcoeffs, bfield.hcoeffs)