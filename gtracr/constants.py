'''
Library that contain relavant constants for this package
'''

import scipy.constants as sc
import numpy as np

SPEED_OF_LIGHT = sc.c

ELEMENTARY_CHARGE = 1.602e-19

EARTH_RADIUS = 6371.2 * (1e3)  # earth radius (meters)
G10 = -29404.8 * (1e-9)  # B-field parameter from IGRF 2020 (in Teslas)

DEG_PER_RAD = np.pi / 180.
RAD_PER_DEG = 180. / np.pi

KG_PER_GEVC2 = 1.78e-27  # conversion factor between kg and GeV/c^2
KG_M_S_PER_GEVC = 5.36e-19  # conversion factor between kg m/s (momentum SI units) and GeV/c
