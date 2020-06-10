'''
Performs Runge-Kutta integration to the 6 coupled differential equations for
the Lorenz Force equation.

Edit (June 9th, 2020): I feel like there should be a much nicer implementation...
'''

import numpy as np


def wsum(n1, n2, n3, n4):
    return (1. / 6.) * n1 + (1. / 3.) * n2 + (1. / 3.) * n3 + (1. / 6.) * n4


# evaluate 4th-order Runge Kutta with 6 coupled differential equations
# this code can be reduced greatly by removing certain arguments, however by
# making the code more general those arguments have to stay there.
def runge_kutta(dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt, h, n, ival):
    i = 1

    #create arrays that will be our final solution
    t = np.zeros(n)
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    vx = np.zeros(n)
    vy = np.zeros(n)
    vz = np.zeros(n)

    #set initial conditions to array
    t[0] = ival[0]
    x[0] = ival[1]
    y[0] = ival[2]
    z[0] = ival[3]
    vx[0] = ival[4]
    vy[0] = ival[5]
    vz[0] = ival[6]

    while i < n:
        j = i - 1  # this is the term before "n+1th" term (ie nth term)

        # get the RK terms k,l,m
        k1 = h * dxdt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
        l1 = h * dydt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
        m1 = h * dzdt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
        p1 = h * dvxdt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
        q1 = h * dvydt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])
        s1 = h * dvzdt(t[j], x[j], y[j], z[j], vx[j], vy[j], vz[j])

        k2 = h * dxdt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
                      z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
                      vz[j] + 0.5 * s1)
        l2 = h * dydt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
                      z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
                      vz[j] + 0.5 * s1)
        m2 = h * dzdt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
                      z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
                      vz[j] + 0.5 * s1)
        p2 = h * dvxdt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
                       z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
                       vz[j] + 0.5 * s1)
        q2 = h * dvydt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
                       z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
                       vz[j] + 0.5 * s1)
        s2 = h * dvzdt(t[j] + 0.5 * h, x[j] + 0.5 * k1, y[j] + 0.5 * l1,
                       z[j] + 0.5 * m1, vx[j] + 0.5 * p1, vy[j] + 0.5 * q1,
                       vz[j] + 0.5 * s1)

        k3 = h * dxdt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
                      z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
                      vz[j] + 0.5 * s2)
        l3 = h * dydt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
                      z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
                      vz[j] + 0.5 * s2)
        m3 = h * dzdt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
                      z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
                      vz[j] + 0.5 * s2)
        p3 = h * dvxdt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
                       z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
                       vz[j] + 0.5 * s2)
        q3 = h * dvydt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
                       z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
                       vz[j] + 0.5 * s2)
        s3 = h * dvzdt(t[j] + 0.5 * h, x[j] + 0.5 * k2, y[j] + 0.5 * l2,
                       z[j] + 0.5 * m2, vx[j] + 0.5 * p2, vy[j] + 0.5 * q2,
                       vz[j] + 0.5 * s2)

        k4 = h * dxdt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
                      vy[j] + q3, vz[j] + s3)
        l4 = h * dydt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
                      vy[j] + q3, vz[j] + s3)
        m4 = h * dzdt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
                      vy[j] + q3, vz[j] + s3)
        p4 = h * dvxdt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
                       vy[j] + q3, vz[j] + s3)
        q4 = h * dvydt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
                       vy[j] + q3, vz[j] + s3)
        s4 = h * dvzdt(t[j] + h, x[j] + k3, y[j] + l3, z[j] + m3, vx[j] + p3,
                       vy[j] + q3, vz[j] + s3)

        # get the weighted sum of each component
        k = wsum(k1, k2, k3, k4)
        l = wsum(l1, l2, l3, l4)
        m = wsum(m1, m2, m3, m4)
        p = wsum(p1, p2, p3, p4)
        q = wsum(q1, q2, q3, q4)
        s = wsum(s1, s2, s3, s4)
        # iterate by weighted sum of stepsize
        x[i] = x[j] + k
        y[i] = y[j] + l
        z[i] = z[j] + m
        vx[i] = vx[j] + p
        vy[i] = vy[j] + q
        vz[i] = vz[j] + s
        t[i] = t[j] + h
        i += 1  # increment

    return (t, x, y, z, vx, vy, vz)


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
