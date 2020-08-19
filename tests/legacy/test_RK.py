import os, sys
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

sys.path.append(os.path.join(os.getcwd(), "gtracr"))
sys.path.append(os.path.join(os.getcwd(), "..", "gtracr"))

# import runge kutta
from runge_kutta import runge_kutta

# global variables (yikes!)
R_E = 6.371e6  # earth radius
g10 = -29404.8 * (1e-9) # B-field parameter from IGRF 2020


# magnetic field (ideal dipole)

def B_r(r, theta):
    return 2.*(R_E/r)**3.*g10*np.cos(theta)

def B_theta(r, theta):
    return (R_E/r)**3.*g10*np.sin(theta)


# here we test for the magnetic field
# radial momentum DE
def dprdt(t,r,th,ph,vr,vth,vph):
    return -(B_theta(r,th)*vph) / 0.937

# theta component momentum DE
def dpthdt(t,r,th,ph,vr,vth,vph):
    return (B_r(r,th)*vph) / 0.937

# phi comp mom. DE
def dpphdt(t,r,th,ph,vr,vth,vph):
    return (B_theta(r,th)*vr - B_r(r,th)*vth) / 0.937

# radial coord DE
def drdt(t,r,th,ph,vr,vth,vph):
    return vr

# theta component coord DE
def dthdt(t,r,th,ph,vr,vth,vph):
    return vth

# phi comp mom. DE
def dphdt(t,r,th,ph,vr,vth,vph):
    return vph



if __name__ == "__main__":
    # set initial values
    ival = (0.001, 1e7, 0., 0., 0.1, 0.1, 0.1)

    # perform RK
    result_tup = runge_kutta(drdt, dthdt, dphdt, dprdt, dpthdt, dpphdt, h=0.01, n=1000, ival=ival)

    t = result_tup[0]
    r = result_tup[1]
    theta = result_tup[2]
    phi =  result_tup[3]

    # convert back to cartesian
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    # plot
    ax = plt.axes(projection="3d")
    ax.plot3D(x, y, z, color="k")
    plt.show()

    fig = go.Figure(data=[go.Scatter3d(x=x,y=y,z=z,mode="markers",
                    marker=dict(size=2.0),
                    line=dict(color='darkblue', width=2))])

    fig.show()















# # define some diff equation
# def dxdt(t,x,y,z):
#     return -np.sin(t)

# def dydt(t,x,y,z):
#     return np.cos(t)

# def dzdt(t,x,y,z):
#     return t

# # initial values
# ival = (0.001, 0.1, 0.1, 0.1)
# # stepsize
# h = 5e-2
# # max iterations
# n = 1000

# t, x, y, z = runge_kutta(dxdt, dydt, dzdt, h, n, ival)

# # plot
# ax = plt.axes(projection="3d")
# ax.plot3D(x,y,z,color="k")
# plt.show()

# fig = go.Figure(data=[go.Scatter3d(x=x,y=y,z=z,mode="markers",
#                 marker=dict(size=2.0),
#                 line=dict(color='darkblue', width=2))])

# fig.show()