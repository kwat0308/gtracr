import time
from gtracr.geomagnetic_cutoffs import GMRC

# default GMRC at Kamioka
gmrc_serial = GMRC()

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