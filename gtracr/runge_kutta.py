'''
Performs Runge-Kutta integration to the 6 coupled differential equations for
the Lorenz Force equation.

Edit (June 9th, 2020): I feel like there should be a much nicer implementation...
'''

import os, sys
import numpy as np

sys.path.append(os.getcwd())
sys.path.append(os.path.join(os.getcwd(), "gtracr"))

from gtracr.utils import EARTH_RADIUS, g10, SPEED_OF_LIGHT, gamma, vmag_spherical
from gtracr.magnetic_field import B_r, B_theta, B_phi

# magnetic field (ideal dipole)

# def B_r(r, theta):
#     return 2.*(EARTH_RADIUS/r)**3.*g10*np.cos(theta)
#     # return 2.*(1./r)**3.*g10*np.cos(theta)
# #
# def B_theta(r, theta):
#     return (EARTH_RADIUS/r)**3.*g10*np.sin(theta)
#     # return (1./r)**3.*g10*np.sin(theta)

# def B_phi(r, theta):
#     return 0.

# # here we test for the magnetic field (ideal dipole)
# # radial momentum DE
# def dprdt(t,r,th,ph,vr,vth,vph,particle):
#     return -particle.charge*(B_theta(r,th)*vph) / particle.mass

# # theta component momentum DE
# def dpthdt(t,r,th,ph,vr,vth,vph,particle):
#     return particle.charge*(B_r(r,th)*vph) / particle.mass

# # phi comp mom. DE
# def dpphdt(t,r,th,ph,vr,vth,vph, particle):
#     return particle.charge*(B_theta(r,th)*vr - B_r(r,th)*vth) / particle.mass


# here we test for the magnetic field (ideal dipole)
# radial velocity DE
def dvrdt(t, r, theta, phi, vr, vtheta, vphi, particle):
    term1 = particle.charge * (
        vtheta * B_phi(r, theta) - B_theta(r, theta) * vphi) / (
            particle.mass * gamma(vmag_spherical(vr, vtheta, vphi, r, theta)))
    term2 = vtheta**2. / r
    term3 = vphi**2. / r
    return term1 + term2 + term3


# theta component velocity DE
def dvthetadt(t, r, theta, phi, vr, vtheta, vphi, particle):
    term1 = particle.charge * (vphi * B_r(r, theta) - B_phi(r, theta) * vr) / (
        particle.mass * gamma(vmag_spherical(vr, vtheta, vphi, r, theta)))
    term2 = (vr * vtheta) / r
    term3 = vphi**2. / (r * np.tan(theta))
    return term1 - term2 + term3


# phi comp mom. DE
def dvphidt(t, r, theta, phi, vr, vtheta, vphi, particle):
    term1 = particle.charge * (
        B_theta(r, theta) * vr - B_r(r, theta) * vtheta) / (
            particle.mass * gamma(vmag_spherical(vr, vtheta, vphi, r, theta)))
    term2 = (vr * vphi) / r
    term3 = (vtheta * vphi) / (r * np.tan(theta))
    return term1 - term2 - term3


def wsum(n1, n2, n3, n4):
    return (1. / 6.) * n1 + (1. / 3.) * n2 + (1. / 3.) * n3 + (1. / 6.) * n4


# evaluate 4th-order Runge Kutta with 6 coupled differential equations
# this code can be reduced greatly by removing certain arguments, however by
# making the code more general those arguments have to stay there.
# Edit 2: this only performs one iteration of RK
# the position DEs are replaced with the spherical definition for velocity
def runge_kutta(particle, h, ival):

    #set initial conditions to array
    t = ival[0]
    r = ival[1]
    th = ival[2]
    ph = ival[3]
    vr = ival[4]
    vth = ival[5]
    vph = ival[6]

    # get the RK terms
    # k,l,m are for drdt, dthdt, dphdt
    # p,q,s are for momenta
    k1 = h * vr
    l1 = h * vth / r
    m1 = h * vph / (r * np.sin(th))
    a1 = h * dvrdt(t, r, th, ph, vr, vth, vph, particle)
    b1 = h * dvthetadt(t, r, th, ph, vr, vth, vph, particle)
    c1 = h * dvphidt(t, r, th, ph, vr, vth, vph, particle)

    k2 = h * vr + 0.5 * a1
    l2 = h * vth + 0.5 * b1 / (r + 0.5 * k1)
    m2 = h * vph + 0.5 * c1 / ((r + 0.5 * k1) * (np.sin(th + 0.5 * l1)))
    a2 = h * dvrdt(t + 0.5 * h, r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                   vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1, particle)
    b2 = h * dvthetadt(t + 0.5 * h, r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                       vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1, particle)
    c2 = h * dvphidt(t + 0.5 * h, r + 0.5 * k1, th + 0.5 * l1, ph + 0.5 * m1,
                     vr + 0.5 * a1, vth + 0.5 * b1, vph + 0.5 * c1, particle)

    k3 = h * vr + 0.5 * a2
    l3 = h * vth + 0.5 * b2 / (r + 0.5 * k2)
    m3 = h * vph + 0.5 * c2 / ((r + 0.5 * k2) * (np.sin(th + 0.5 * l2)))
    a3 = h * dvrdt(t + 0.5 * h, r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                   vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2, particle)
    b3 = h * dvthetadt(t + 0.5 * h, r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                       vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2, particle)
    c3 = h * dvphidt(t + 0.5 * h, r + 0.5 * k2, th + 0.5 * l2, ph + 0.5 * m2,
                     vr + 0.5 * a2, vth + 0.5 * b2, vph + 0.5 * c2, particle)

    k4 = h * vr + a3
    l4 = h * vth + b3 / (r + k3)
    m4 = h * vph + c3 / ((r + k3) * (np.sin(th + l3)))
    a4 = h * dvrdt(t + h, r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                   vph + c3, particle)
    b4 = h * dvthetadt(t + h, r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                        vph + c3, particle)
    c4 = h * dvphidt(t + h, r + k3, th + l3, ph + m3, vr + a3, vth + b3,
                     vph + c3, particle)

    # get the weighted sum of each component
    k = wsum(k1, k2, k3, k4)
    l = wsum(l1, l2, l3, l4)
    m = wsum(m1, m2, m3, m4)
    a = wsum(a1, a2, a3, a4)
    b = wsum(b1, b2, b3, b4)
    c = wsum(c1, c2, c3, c4)
    # iterate by weighted sum of stepsize
    r = r + k
    th = th + l
    ph = ph + m
    vr = vr + a
    vth = vth + b
    vph = vph + c
    t = t + h

    return np.array([t, r, th, ph, vr, vth, vph])


# # evaluate 4th-order Runge Kutta with 6 coupled differential equations
# # this code can be reduced greatly by removing certain arguments, however by
# # making the code more general those arguments have to stay there.
# # for this version, the drdt / dthdt / dphdt arguments are replaced with
# # pr, pth, pph in the RK steps
# def runge_kutta(particle, h, n, ival):
#     i = 1

#     #create arrays that will be our final solution
#     t = np.zeros(n)
#     r = np.zeros(n)
#     th = np.zeros(n)
#     ph = np.zeros(n)
#     pr = np.zeros(n)
#     pth = np.zeros(n)
#     pph = np.zeros(n)

#     #set initial conditions to array
#     t[0] = ival[0]
#     r[0] = ival[1]
#     th[0] = ival[2]
#     ph[0] = ival[3]
#     pr[0] = ival[4]
#     pth[0] = ival[5]
#     pph[0] = ival[6]

#     while i < n:
#         j = i - 1  # this is the term before "n+1th" term (ie nth term)

#         # get the RK terms
#         # k,l,m are for drdt, dthdt, dphdt
#         # p,q,s are for momenta
#         k1 = h * pr[j]
#         l1 = h * pth[j]
#         m1 = h * pph[j]
#         a1 = h * dprdt(t[j], r[j], th[j], ph[j], pr[j], pth[j], pph[j], particle)
#         b1 = h * dpthdt(t[j], r[j], th[j], ph[j], pr[j], pth[j], pph[j], particle)
#         c1 = h * dpphdt(t[j], r[j], th[j], ph[j], pr[j], pth[j], pph[j], particle)

#         k2 = h * pr[j] + 0.5 * a1
#         l2 = h * pth[j] + 0.5 * b1
#         m2 = h * pph[j] + 0.5 * c1
#         a2 = h * dprdt(t[j] + 0.5 * h, r[j] + 0.5 * k1, th[j] + 0.5 * l1,
#                        ph[j] + 0.5 * m1, pr[j] + 0.5 * a1, pth[j] + 0.5 * b1,
#                        pph[j] + 0.5 * c1, particle)
#         b2 = h * dpthdt(t[j] + 0.5 * h, r[j] + 0.5 * k1, th[j] + 0.5 * l1,
#                        ph[j] + 0.5 * m1, pr[j] + 0.5 * a1, pth[j] + 0.5 * b1,
#                        pph[j] + 0.5 * c1, particle)
#         c2 = h * dpphdt(t[j] + 0.5 * h, r[j] + 0.5 * k1, th[j] + 0.5 * l1,
#                        ph[j] + 0.5 * m1, pr[j] + 0.5 * a1, pth[j] + 0.5 * b1,
#                        pph[j] + 0.5 * c1, particle)

#         k3 = h * pr[j] + 0.5 * a2
#         l3 = h * pth[j] + 0.5 * b2
#         m3 = h * pph[j] + 0.5 * c2
#         a3 = h * dprdt(t[j] + 0.5 * h, r[j] + 0.5 * k2, th[j] + 0.5 * l2,
#                        ph[j] + 0.5 * m2, pr[j] + 0.5 * a2, pth[j] + 0.5 * b2,
#                        pph[j] + 0.5 * c2, particle)
#         b3 = h * dpthdt(t[j] + 0.5 * h, r[j] + 0.5 * k2, th[j] + 0.5 * l2,
#                        ph[j] + 0.5 * m2, pr[j] + 0.5 * a2, pth[j] + 0.5 * b2,
#                        pph[j] + 0.5 * c2, particle)
#         c3 = h * dpphdt(t[j] + 0.5 * h, r[j] + 0.5 * k2, th[j] + 0.5 * l2,
#                        ph[j] + 0.5 * m2, pr[j] + 0.5 * a2, pth[j] + 0.5 * b2,
#                        pph[j] + 0.5 * c2, particle)

#         k4 = h * pr[j] + a3
#         l4 = h * pth[j] + b3
#         m4 = h * pph[j] + c3
#         a4 = h * dprdt(t[j] + h, r[j] + k3, th[j] + l3, ph[j] + m3, pr[j] + a3,
#                        pth[j] + b3, pph[j] + c3, particle)
#         b4 = h * dpthdt(t[j] + h, r[j] + k3, th[j] + l3, ph[j] + m3, pr[j] + a3,
#                        pth[j] + b3, pph[j] + c3, particle)
#         c4 = h * dpphdt(t[j] + h, r[j] + k3, th[j] + l3, ph[j] + m3, pr[j] + a3,
#                        pth[j] + b3, pph[j] + c3, particle)

#         # get the weighted sum of each component
#         k = wsum(k1, k2, k3, k4)
#         l = wsum(l1, l2, l3, l4)
#         m = wsum(m1, m2, m3, m4)
#         a = wsum(a1, a2, a3, a4)
#         b = wsum(b1, b2, b3, b4)
#         c = wsum(c1, c2, c3, c4)
#         # iterate by weighted sum of stepsize
#         r[i] = r[j] + k
#         th[i] = th[j] + l
#         ph[i] = ph[j] + m
#         pr[i] = pr[j] + a
#         pth[i] = pth[j] + b
#         pph[i] = pph[j] + c
#         t[i] = t[j] + h
#         i += 1  # increment

#     return np.array(t, r, th, ph, pr, pth, pph)

# # evaluate 4th-order Runge Kutta with 6 coupled differential equations
# # this code can be reduced greatly by removing certain arguments, however by
# # making the code more general those arguments have to stay there.
# # for this version, the drdt / dthdt / dphdt arguments are replaced with
# # pr, pth, pph in the RK steps
# def runge_kutta(dprdt, dpthdt, dpphdt, particle, h, n, ival):
#     i = 1

#     #create arrays that will be our final solution
#     t = np.zeros(n)
#     r = np.zeros(n)
#     th = np.zeros(n)
#     ph = np.zeros(n)
#     pr = np.zeros(n)
#     pth = np.zeros(n)
#     pph = np.zeros(n)

#     #set initial conditions to array
#     t[0] = ival[0]
#     r[0] = ival[1]
#     th[0] = ival[2]
#     ph[0] = ival[3]
#     pr[0] = ival[4]
#     pth[0] = ival[5]
#     pph[0] = ival[6]

#     while i < n:
#         j = i - 1  # this is the term before "n+1th" term (ie nth term)

#         # get the RK terms
#         # k,l,m are for drdt, dthdt, dphdt
#         # p,q,s are for momenta
#         k1 = h * pr[j]
#         l1 = h * pth[j]
#         m1 = h * pph[j]
#         a1 = h * dprdt(t[j], r[j], th[j], ph[j], pr[j], pth[j], pph[j], particle)
#         b1 = h * dpthdt(t[j], r[j], th[j], ph[j], pr[j], pth[j], pph[j], particle)
#         c1 = h * dpphdt(t[j], r[j], th[j], ph[j], pr[j], pth[j], pph[j], particle)

#         k2 = h * pr[j] + 0.5 * a1
#         l2 = h * pth[j] + 0.5 * b1
#         m2 = h * pph[j] + 0.5 * c1
#         a2 = h * dprdt(t[j] + 0.5 * h, r[j] + 0.5 * k1, th[j] + 0.5 * l1,
#                        ph[j] + 0.5 * m1, pr[j] + 0.5 * a1, pth[j] + 0.5 * b1,
#                        pph[j] + 0.5 * c1, particle)
#         b2 = h * dpthdt(t[j] + 0.5 * h, r[j] + 0.5 * k1, th[j] + 0.5 * l1,
#                        ph[j] + 0.5 * m1, pr[j] + 0.5 * a1, pth[j] + 0.5 * b1,
#                        pph[j] + 0.5 * c1, particle)
#         c2 = h * dpphdt(t[j] + 0.5 * h, r[j] + 0.5 * k1, th[j] + 0.5 * l1,
#                        ph[j] + 0.5 * m1, pr[j] + 0.5 * a1, pth[j] + 0.5 * b1,
#                        pph[j] + 0.5 * c1, particle)

#         k3 = h * pr[j] + 0.5 * a2
#         l3 = h * pth[j] + 0.5 * b2
#         m3 = h * pph[j] + 0.5 * c2
#         a3 = h * dprdt(t[j] + 0.5 * h, r[j] + 0.5 * k2, th[j] + 0.5 * l2,
#                        ph[j] + 0.5 * m2, pr[j] + 0.5 * a2, pth[j] + 0.5 * b2,
#                        pph[j] + 0.5 * c2, particle)
#         b3 = h * dpthdt(t[j] + 0.5 * h, r[j] + 0.5 * k2, th[j] + 0.5 * l2,
#                        ph[j] + 0.5 * m2, pr[j] + 0.5 * a2, pth[j] + 0.5 * b2,
#                        pph[j] + 0.5 * c2, particle)
#         c3 = h * dpphdt(t[j] + 0.5 * h, r[j] + 0.5 * k2, th[j] + 0.5 * l2,
#                        ph[j] + 0.5 * m2, pr[j] + 0.5 * a2, pth[j] + 0.5 * b2,
#                        pph[j] + 0.5 * c2, particle)

#         k4 = h * pr[j] + a3
#         l4 = h * pth[j] + b3
#         m4 = h * pph[j] + c3
#         a4 = h * dprdt(t[j] + h, r[j] + k3, th[j] + l3, ph[j] + m3, pr[j] + a3,
#                        pth[j] + b3, pph[j] + c3, particle)
#         b4 = h * dpthdt(t[j] + h, r[j] + k3, th[j] + l3, ph[j] + m3, pr[j] + a3,
#                        pth[j] + b3, pph[j] + c3, particle)
#         c4 = h * dpphdt(t[j] + h, r[j] + k3, th[j] + l3, ph[j] + m3, pr[j] + a3,
#                        pth[j] + b3, pph[j] + c3, particle)

#         # get the weighted sum of each component
#         k = wsum(k1, k2, k3, k4)
#         l = wsum(l1, l2, l3, l4)
#         m = wsum(m1, m2, m3, m4)
#         a = wsum(a1, a2, a3, a4)
#         b = wsum(b1, b2, b3, b4)
#         c = wsum(c1, c2, c3, c4)
#         # iterate by weighted sum of stepsize
#         r[i] = r[j] + k
#         th[i] = th[j] + l
#         ph[i] = ph[j] + m
#         pr[i] = pr[j] + a
#         pth[i] = pth[j] + b
#         pph[i] = pph[j] + c
#         t[i] = t[j] + h
#         i += 1  # increment

#     return np.array(t, r, th, ph, pr, pth, pph)

# # evaluate 4th-order Runge Kutta with 6 coupled differential equations
# # this code can be reduced greatly by removing certain arguments, however by
# # making the code more general those arguments have to stay there.
# def runge_kutta(dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt, h, n, ival):
#     i = 1

#     #create arrays that will be our final solution
#     t = np.zeros(n)
#     x = np.zeros(n)
#     y = np.zeros(n)
#     z = np.zeros(n)
#     vx = np.zeros(n)
#     vy = np.zeros(n)
#     vz = np.zeros(n)

#     #set initial conditions to array
#     t[0] = ival[0]
#     x[0] = ival[1]
#     y[0] = ival[2]
#     z[0] = ival[3]
#     vx[0] = ival[4]
#     vy[0] = ival[5]
#     vz[0] = ival[6]

#     while i < n:
#         j = i - 1  # this is the term before "n+1th" term (ie nth term)

#         # get the RK terms k,l,m
#         k1 = h * dxdt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
#         l1 = h * dydt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
#         m1 = h * dzdt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
#         p1 = h * dvxdt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
#         q1 = h * dvydt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
#         s1 = h * dvzdt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])

#         k2 = h * dxdt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
#                       z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
#                       vz[j] + 0.5 * s1)
#         l2 = h * dydt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
#                       z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
#                       vz[j] + 0.5 * s1)
#         m2 = h * dzdt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
#                       z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
#                       vz[j] + 0.5 * s1)
#         p2 = h * dvxdt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
#                        z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
#                        vz[j] + 0.5 * s1)
#         q2 = h * dvydt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
#                        z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
#                        vz[j] + 0.5 * s1)
#         s2 = h * dvzdt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
#                        z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
#                        vz[j] + 0.5 * s1)

#         k3 = h * dxdt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
#                       z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
#                       vz[j] + 0.5 * s2)
#         l3 = h * dydt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
#                       z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
#                       vz[j] + 0.5 * s2)
#         m3 = h * dzdt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
#                       z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
#                       vz[j] + 0.5 * s2)
#         p3 = h * dvxdt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
#                        z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
#                        vz[j] + 0.5 * s2)
#         q3 = h * dvydt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
#                        z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
#                        vz[j] + 0.5 * s2)
#         s3 = h * dvzdt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
#                        z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
#                        vz[j] + 0.5 * s2)

#         k4 = h * dxdt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
#                       vy[j] + q3, vz[j] + s3)
#         l4 = h * dydt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
#                       vy[j] + q3, vz[j] + s3)
#         m4 = h * dzdt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
#                       vy[j] + q3, vz[j] + s3)
#         p4 = h * dvxdt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
#                        vy[j] + q3, vz[j] + s3)
#         q4 = h * dvydt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
#                        vy[j] + q3, vz[j] + s3)
#         s4 = h * dvzdt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
#                        vy[j] + q3, vz[j] + s3)

#         # get the weighted sum of each component
#         k = wsum(k1, k2, k3, k4)
#         l = wsum(l1, l2, l3, l4)
#         m = wsum(m1, m2, m3, m4)
#         p = wsum(p1, p2, p3, p4)
#         q = wsum(q1, q2, q3, q4)
#         s = wsum(s1, s2, s3, s4)
#         # iterate by weighted sum of stepsize
#         x[i] = x[j] + k
#         y[i] = y[j] + l
#         z[i] = z[j] + m
#         vx[i] = vx[j] + p
#         vy[i] = vy[j] + q
#         vz[i] = vz[j] + s
#         t[i] = t[j] + h
#         i += 1  # increment

#     return (t, x, y, z, vx, vy, vz)

# performs 4th order Runge-Kutta until value x
# returns
# inputs:
#   - diff_eq: the DE that we are trying to solve for (assume time-independent)
#   - max_iter: the maximum iteration to perform RK for (not yet)
#   - h : the step size
#   - ftup : tuple of final values (tf, rf, thetaf, phif)
#   - itup : tuple of initial value (t0, r0, theta0, phi0)
# def runge_kutta(diff_eq, h, ftup, itup):
#     # get max iterations from step
#     max_iter = np.floor((ftup[0] - itup[0])/h)

#     # set initial values
#     t = itup[0]
#     r = itup[1]
#     theta = itup[2]
#     phi = itup[3]

#     for i in range(max_iter):
#         k1 = h*dxdt(taui, xi, yi, zi, params)
#         l1 = h*dydt(taui, xi, yi, zi, params)
#         m1 = h*dzdt(taui, xi, yi, zi, params)

#         k2 = h*dxdt(taui + 0.5*h, xi + 0.5*k1, yi + 0.5*l1, zi + 0.5*m1, params)
#         l2 = h*dydt(taui + 0.5*h, xi + 0.5*k1, yi + 0.5*l1, zi + 0.5*m1, params)
#         m2 = h*dzdt(taui + 0.5*h, xi + 0.5*k1, yi + 0.5*l1, zi + 0.5*m1, params)

#         k3 = h*dxdt(taui + 0.5*h, xi + 0.5*k2, yi + 0.5*l2, zi + 0.5*m2, params)
#         l3 = h*dydt(taui + 0.5*h, xi + 0.5*k2, yi + 0.5*l2, zi + 0.5*m2, params)
#         m3 = h*dzdt(taui + 0.5*h, xi + 0.5*k2, yi + 0.5*l2, zi + 0.5*m2, params)

#         k4 = h*dxdt(taui + h, xi + k3, yi + l3, zi + m3, params)
#         l4 = h*dydt(taui + h, xi + k3, yi + l3, zi + m3, params)
#         m4 = h*dzdt(taui + h, xi + k3, yi + l3, zi + m3, params)

#         # get the weighted sum of each component
#         k = wsum(k1, k2, k3, k4)
#         l = wsum(l1, l2, l3, l4)
#         m = wsum(m1, m2, m3, m4)

#         return k, l, m
