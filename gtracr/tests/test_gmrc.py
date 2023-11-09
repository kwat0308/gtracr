import time
import sys
import numpy as np
from p_tqdm import p_map
from gtracr.plotting import plot_gmrc_heatmap

# using the local GMRC class instead of from the package
sys.path.append("..")
from geomagnetic_cutoffs import GMRC

# a main guard is required to be multiprocessing safe in Windows
# see eg. https://stackoverflow.com/questions/20360686/compulsory-usage-of-if-name-main-in-windows-while-using-multiprocessi
if __name__ == "__main__":

    # default GMRC at Kamioka
    gmrc_serial = GMRC(iter_num= 5000)

    start = time.time()
    gmrc_serial.evaluate()
    end = time.time()

    print("GMRC serial evaluated in", end - start, "secs")

    interpd_gmrc_data = gmrc_serial.interpolate_results()

    plot_gmrc_heatmap(interpd_gmrc_data,
                      gmrc_serial.rigidity_list,
                      locname=gmrc_serial.location,
                      plabel=gmrc_serial.plabel,
                      bfield_type = '',
                      show_plot = True,
                      plotdir_path = '')

    gmrc_parallel = GMRC(iter_num = 5000, method='parallel')

    start = time.time()
    # azimuth = np.random.rand(gmrc_parallel.iter_num) * 360.0
    # zenith = np.random.rand(gmrc_parallel.iter_num) * 180.0
    # rigidity = p_map(gmrc_parallel.evaluate_angle, azimuth, zenith)
    gmrc_parallel.evaluate()
    end = time.time()

    print("GMRC parallel evaluated in", end - start, "secs")

    interpd_gmrc_data = gmrc_parallel.interpolate_results()

    plot_gmrc_heatmap(interpd_gmrc_data,
                    gmrc_parallel.rigidity_list,
                    locname=gmrc_parallel.location,
                    plabel=gmrc_parallel.plabel,
                    bfield_type = '',
                    show_plot = True,
                    plotdir_path = '')