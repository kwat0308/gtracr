import os, sys
import numpy as np
import matplotlib.pyplot as plt

# define global variables
SPEED_OF_LIGHT = 3e8
ELEMENTARY_CHARGE = 1.602e-19
EARTH_RADIUS = 6.371 * (1e6)
RAD_PER_DEG = np.pi / 180.
KG_PER_GEVC2 = 1.78e-27
KGMS_PER_GEVC = 5.36e-19

KEY_LIST = ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]

# define the spherical components of Earth magnetic field
G10 = 2.91e-5  # value of B-field at center of Earth


def Br(r, theta, phi):
    return -2. * (EARTH_RADIUS / r)**3. * G10 * np.cos(theta)


def Btheta(r, theta, phi):
    return -(EARTH_RADIUS / r)**3. * G10 * np.sin(theta)


def Bphi(r, theta, phi):
    return 0.


# ODEs in spherical coordinates,
# based on lorentz force and the form of
# acceleration in spherical coordinates
def drdt(pr):
    return pr


def dthetadt(r, ptheta):
    return ptheta / r


def dphidt(r, theta, pphi):
    # dphi_dt = 0.0 if abs(theta) < 1e-12 else pphi / (r * np.sin(theta))
    return pphi / (r * np.sin(theta))


def dprdt(charge, r, theta, phi, pr, ptheta, pphi):
    lorentz_term = charge * (ptheta * Bphi(r, theta, phi) -
                             pphi * Btheta(r, theta, phi))
    auxiliary_terms = (ptheta**2. / r) + (pphi**2. / r)
    return lorentz_term + auxiliary_terms


def dpthetadt(charge, r, theta, phi, pr, ptheta, pphi):
    lorentz_term = -charge * (pr * Bphi(r, theta, phi) -
                              pphi * Br(r, theta, phi))
    auxiliary_terms = (pr**2. / (r * np.tan(theta))) - ((pr * ptheta) / r)
    return lorentz_term + auxiliary_terms


def dpphidt(charge, r, theta, phi, pr, ptheta, pphi):
    lorentz_term = charge * (pr * Btheta(r, theta, phi) -
                             ptheta * Br(r, theta, phi))
    auxiliary_terms = ((pr * pphi) / r) + ((ptheta * pphi) /
                                           (r * np.tan(theta)))
    return lorentz_term - auxiliary_terms


# convert between coordinates defined in detector frame with geocentric coordinates
def detector_to_geocentric(lat, lng, detector_alt, zenith, azimuth,
                           particle_alt, momentum):

    # first convert angles (latitude, longitude, zenith, azimuth)
    # to radians
    lat *= RAD_PER_DEG
    lng *= RAD_PER_DEG
    zenith *= RAD_PER_DEG
    azimuth *= RAD_PER_DEG

    # now convert the coordinates / momentum
    # this can realistically be done in one step by utilizing matrix multiplication
    # with transformation matrices between cartesian, geodesic, and spherical

    # convert detector latitude, longitude, altitude into cartesian
    detector_x = (detector_alt + EARTH_RADIUS) * np.cos(lat) * np.cos(lng)
    detector_y = (detector_alt + EARTH_RADIUS) * np.cos(lat) * np.sin(lng)
    detector_z = (detector_alt + EARTH_RADIUS) * np.sin(lat)

    detector_geoc_coord = np.array([detector_x, detector_y, detector_z])

    # convert zenith, azimuth, and particle altitude to cartesian in detector frame
    particle_detector_x = (particle_alt +
                           EARTH_RADIUS) * np.sin(zenith) * np.cos(azimuth)
    particle_detector_y = (particle_alt +
                           EARTH_RADIUS) * np.sin(zenith) * np.sin(azimuth)
    particle_detector_z = (particle_alt + EARTH_RADIUS) * np.cos(zenith)

    particle_detector_coord = np.array(
        [particle_detector_x, particle_detector_y, particle_detector_z])

    # the transformation matrix between detector frame and geodesic coordinates
    row1 = np.array(
        [-np.sin(lng), -np.cos(lng) * np.sin(lat),
         np.cos(lat) * np.cos(lng)])
    row2 = np.array(
        [np.cos(lng), -np.sin(lat) * np.sin(lng),
         np.cos(lat) * np.sin(lng)])
    row3 = np.array([0., np.cos(lat), np.sin(lat)])

    transform_matrix = np.array([row1, row2, row3])

    # now get the transformed coordinate vector in geocentric frame
    (particle_geoc_x,
     particle_geoc_y, particle_geoc_z) = detector_geoc_coord + np.dot(
         transform_matrix, particle_detector_coord)

    particle_r = np.sqrt(particle_geoc_x**2. + particle_geoc_y**2. +
                         particle_geoc_z**2.)
    particle_theta = np.arccos(particle_geoc_z / particle_r)
    particle_phi = np.arctan2(particle_geoc_y, particle_geoc_x)

    # now for momentum
    # convert the direction of the momentum of the particle (defined by zenith and azimuth angles)
    # with the magnitude of the momentum (obtained by energy / rigidity)
    particle_detector_px = momentum * np.sin(zenith) * np.cos(azimuth)
    particle_detector_py = momentum * np.sin(zenith) * np.sin(azimuth)
    particle_detector_pz = momentum * np.cos(zenith)

    particle_detector_momentum = np.array(
        [particle_detector_px, particle_detector_py, particle_detector_pz])

    # now transform the momentum vector into geocentric frame
    (particle_geoc_px, particle_geoc_py,
     particle_geoc_pz) = np.dot(transform_matrix, particle_detector_momentum)

    # here we should have conversion between cartesian and spherical again
    # but for momenta its more complicated
    # double check if this is actually correct (like the /r , /rsin(theta) parts)
    #  -> this is fixed, it should be removed to keep units consistent
    particle_pr = particle_geoc_px * np.sin(particle_theta) * np.cos(
        particle_phi) + particle_geoc_py * np.sin(particle_theta) * np.sin(
            particle_phi) + particle_geoc_pz * np.cos(particle_theta)
    particle_ptheta = (
        particle_geoc_px * np.cos(particle_theta) * np.cos(particle_phi) +
        particle_geoc_py * np.cos(particle_theta) * np.sin(particle_phi) -
        particle_geoc_pz * np.sin(particle_theta))
    particle_pphi = (-particle_geoc_px * np.sin(particle_phi) +
                     particle_geoc_py * np.cos(particle_phi))

    particle_sixvector = np.array([
        particle_r, particle_theta, particle_phi, particle_pr, particle_ptheta,
        particle_pphi
    ])

    return particle_sixvector


# perform the Runge kutta sequence
def perform_runge_kutta(mass, charge, t0, initial_sixvec, dt, max_iter):
    # define array to append to
    # t_arr = np.zeros(max_iter)
    # r_arr = np.zeros(max_iter)
    # theta_arr = np.zeros(max_iter)
    # phi_arr = np.zeros(max_iter)
    # pr_arr = np.zeros(max_iter)
    # ptheta_arr = np.zeros(max_iter)
    # pphi_arr = np.zeros(max_iter)

    traj_dict = {key: [] for key in KEY_LIST}

    # define the initial values
    t = t0
    r = initial_sixvec[0]
    theta = initial_sixvec[1]
    phi = initial_sixvec[2]
    pr = initial_sixvec[3]
    ptheta = initial_sixvec[4]
    pphi = initial_sixvec[5]

    # get the relativistic mass (mass * lorentz factor)
    # to do this we get magnitude of momentum
    # since no acceleration should occur, lorentz factor will
    # remain constant, so it can be defined here
    p = np.sqrt(pr**2. + ptheta**2. + pphi**2.)
    print(p)
    # print(mass * SPEED_OF_LIGHT)
    gamma = np.sqrt(1. + (p / (mass * SPEED_OF_LIGHT))**2.)

    print("Lorentz factor: {:.10e}\n".format(gamma))
    rel_mass = mass * gamma

    # perform the integration process
    for i in range(max_iter):
        # first append the values to the array
        traj_dict["t"].append(t)
        traj_dict["r"].append(r)
        traj_dict["theta"].append(theta)
        traj_dict["phi"].append(phi)
        traj_dict["pr"].append(pr)
        traj_dict["ptheta"].append(ptheta)
        traj_dict["pphi"].append(pphi)

        # perform runge kutta steps

        # first runge kutta variable
        r_k1 = dt * drdt(pr)
        theta_k1 = dt * dthetadt(r, ptheta)
        phi_k1 = dt * dphidt(r, theta, pphi)
        pr_k1 = dt * dprdt(charge, r, theta, phi, pr, ptheta, pphi)
        ptheta_k1 = dt * dpthetadt(charge, r, theta, phi, pr, ptheta, pphi)
        pphi_k1 = dt * dpphidt(charge, r, theta, phi, pr, ptheta, pphi)

        # second runge kutta variable
        r_k2 = dt * drdt(pr + 0.5 * pr_k1)
        theta_k2 = dt * dthetadt(r + 0.5 * r_k1, ptheta + 0.5 * ptheta_k1)
        phi_k2 = dt * dphidt(r + 0.5 * r_k1, theta + 0.5 * theta_k1,
                             pphi + 0.5 * pphi_k1)
        pr_k2 = dt * dprdt(charge, r + 0.5 * r_k1, theta + 0.5 * theta_k1,
                           phi + 0.5 * phi_k1, pr + 0.5 * pr_k1,
                           ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1)
        ptheta_k2 = dt * dpthetadt(
            charge, r + 0.5 * r_k1, theta + 0.5 * theta_k1, phi + 0.5 * phi_k1,
            pr + 0.5 * pr_k1, ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1)
        pphi_k2 = dt * dpphidt(charge, r + 0.5 * r_k1, theta + 0.5 * theta_k1,
                               phi + 0.5 * phi_k1, pr + 0.5 * pr_k1,
                               ptheta + 0.5 * ptheta_k1, pphi + 0.5 * pphi_k1)

        # third runge kutta variable
        r_k3 = dt * drdt(pr + 0.5 * pr_k2)
        theta_k3 = dt * dthetadt(r + 0.5 * r_k2, ptheta + 0.5 * ptheta_k2)
        phi_k3 = dt * dphidt(r + 0.5 * r_k2, theta + 0.5 * theta_k2,
                             pphi + 0.5 * pphi_k2)
        pr_k3 = dt * dprdt(charge, r + 0.5 * r_k2, theta + 0.5 * theta_k2,
                           phi + 0.5 * phi_k2, pr + 0.5 * pr_k2,
                           ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2)
        ptheta_k3 = dt * dpthetadt(
            charge, r + 0.5 * r_k2, theta + 0.5 * theta_k2, phi + 0.5 * phi_k2,
            pr + 0.5 * pr_k2, ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2)
        pphi_k3 = dt * dpphidt(charge, r + 0.5 * r_k2, theta + 0.5 * theta_k2,
                               phi + 0.5 * phi_k2, pr + 0.5 * pr_k2,
                               ptheta + 0.5 * ptheta_k2, pphi + 0.5 * pphi_k2)

        # forth runge kutta variable
        r_k4 = dt * drdt(pr + pr_k3)
        theta_k4 = dt * dthetadt(r + r_k3, ptheta + ptheta_k3)
        phi_k4 = dt * dphidt(r + r_k3, theta + theta_k3, pphi + pphi_k3)
        pr_k4 = dt * dprdt(charge, r + r_k3, theta + theta_k3, phi + phi_k3,
                           pr + pr_k3, ptheta + ptheta_k3, pphi + pphi_k3)
        ptheta_k4 = dt * dpthetadt(charge, r + r_k3, theta + theta_k3,
                                   phi + phi_k3, pr + pr_k3,
                                   ptheta + ptheta_k3, pphi + pphi_k3)
        pphi_k4 = dt * dpphidt(charge, r + r_k3, theta + theta_k3,
                               phi + phi_k3, pr + pr_k3, ptheta + ptheta_k3,
                               pphi + pphi_k3)

        # finally get the weighted sum of them
        r_k = (1. / (6. * rel_mass)) * (r_k1 + 2. * r_k2 + 2. * r_k3 + r_k4)
        theta_k = (1. / (6. * rel_mass)) * (theta_k1 + 2. * theta_k2 +
                                            2. * theta_k3 + theta_k4)
        phi_k = (1. / (6. * rel_mass)) * (phi_k1 + 2. * phi_k2 + 2. * phi_k3 +
                                          phi_k4)
        pr_k = (1. /
                (6. * rel_mass)) * (pr_k1 + 2. * pr_k2 + 2. * pr_k3 + pr_k4)
        ptheta_k = (1. / (6. * rel_mass)) * (ptheta_k1 + 2. * ptheta_k2 +
                                             2. * ptheta_k3 + ptheta_k4)
        pphi_k = (1. / (6. * rel_mass)) * (pphi_k1 + 2. * pphi_k2 +
                                           2. * pphi_k3 + pphi_k4)

        # increment by weighted sums (and stepsize for time)
        r += r_k
        theta += theta_k
        phi += phi_k
        pr += pr_k
        ptheta += ptheta_k
        pphi += pphi_k
        t += dt

        # now check conditions for breaking the loop
        # if particle has effectively escaped
        if r > 10 * EARTH_RADIUS:
            print("Allowed trajectory!\n")
            print(r)
            break

        # if particle has reach back to earth
        if r < EARTH_RADIUS:
            print("Forbidden trajectory...\n")
            print(r)
            break

    print("All iterations complete!")

    return traj_dict


# convert the spherical coordinate versions of the trajectory obtained by integration
# into cartesian to plot them
def get_trajectory_cartesian(trajectory_dict):
    trajectory_t = trajectory_dict["t"]
    trajectory_r = trajectory_dict["r"]
    trajectory_theta = trajectory_dict["theta"]
    trajectory_phi = trajectory_dict["phi"]

    trajectory_x = trajectory_r * np.sin(trajectory_theta) * np.cos(
        trajectory_phi)
    trajectory_y = trajectory_r * np.sin(trajectory_theta) * np.sin(
        trajectory_phi)
    trajectory_z = trajectory_r * np.cos(trajectory_theta)

    # divide by radius of earth for better units
    trajectory_x /= EARTH_RADIUS
    trajectory_y /= EARTH_RADIUS
    trajectory_z /= EARTH_RADIUS

    return (trajectory_t, trajectory_x, trajectory_y, trajectory_z)


# plot the trajectories
def plot_trajectories(t_arr, x_arr, y_arr, z_arr):
    fig_3d = plt.figure()
    ax_3d = fig_3d.add_subplot(111, projection="3d")
    cm_3d = ax_3d.scatter(x_arr, y_arr, z_arr, c=t_arr, marker='o')
    # ax_3d.scatter(x2, y2, z2, c="r", marker='o')
    cbar_3d = fig_3d.colorbar(cm_3d, ax=ax_3d)
    # plt.colorbar()
    # plt.savefig("test.png")
    ax_3d.set_xlim([-6, 6])
    ax_3d.set_ylim([-6, 6])
    ax_3d.set_zlim([-6, 6])
    ax_3d.set_xlabel(r"x [$R_E$]")
    ax_3d.set_ylabel(r"y [$R_E$]")
    ax_3d.set_zlabel(r"z [$R_E$]")
    cbar_3d.ax.set_ylabel("Time [s]")
    plt.show()

    fig_xy, ax_xy = plt.subplots()
    cm_xy = ax_xy.scatter(x_arr, y_arr, c=t_arr)
    cbar_xy = fig_xy.colorbar(cm_xy, ax=ax_xy)
    ax_xy.set_xlim([-6, 6])
    ax_xy.set_ylim([-6, 6])
    ax_xy.set_xlabel(r"x [$R_E$]")
    ax_xy.set_ylabel(r"y [$R_E$]")
    cbar_xy.ax.set_ylabel("Time [s]")
    plt.show()

    fig_xz, ax_xz = plt.subplots()
    cm_xz = ax_xz.scatter(x_arr, z_arr, c=t_arr)
    cbar_xz = fig_xz.colorbar(cm_xz, ax=ax_xz)
    ax_xz.set_xlim([-6, 6])
    ax_xz.set_ylim([-6, 6])
    ax_xz.set_xlabel(r"x [$R_E$]")
    ax_xz.set_ylabel(r"z [$R_E$]")
    cbar_xz.ax.set_ylabel("Time [s]")
    plt.show()

    fig_yz, ax_yz = plt.subplots()
    cm_yz = ax_yz.scatter(y_arr, x_arr, c=t_arr)
    cbar_yz = fig_yz.colorbar(cm_yz, ax=ax_yz)
    ax_yz.set_xlim([-6, 6])
    ax_yz.set_ylim([-6, 6])
    ax_yz.set_xlabel(r"y [$R_E$]")
    ax_yz.set_ylabel(r"z [$R_E$]")
    cbar_yz.ax.set_ylabel("Time [s]")
    plt.show()


if __name__ == "__main__":
    # define particle properties
    mass = 0.938 * KG_PER_GEVC2  # particle (rest) mass in kg
    charge = 1 * ELEMENTARY_CHARGE  # particle charge in coulombs

    print(mass, charge)

    # define the detector location
    detector_lat = 0.
    detector_lng = 0.
    detector_alt = 0.

    # define the zenith, azimuthal angle, and altitude in which particle comes from
    # defined within the local tangent plane of the detector
    particle_zenith = 0.
    particle_azimuth = 0.
    particle_altitude = 1.  # in km

    # define the energy / rigidity of the particle
    # energy = 10.
    rigidity = 10.

    # define the momentum (in magnitude, in SI units)
    momentum = rigidity * np.abs(1) * KGMS_PER_GEVC

    # momentum = 10. * KGMS_PER_GEVC

    print(momentum)

    # get the initial values
    # this is obtained from the angles and location in geodesic coordinates
    # the returned value it the six-vector of the initial stage in the trajectory
    particle_initialvec = detector_to_geocentric(detector_lat, detector_lng,
                                                 detector_alt, particle_zenith,
                                                 particle_azimuth,
                                                 particle_altitude, momentum)

    # particle_initialvec = [EARTH_RADIUS * 2., np.pi / 2., 0., momentum, 0., 0.]

    print(particle_initialvec)

    # set the initial time
    t0 = 0.

    # perform the runge_kutta sequence
    dt = 1e-8  # stepsize
    max_iter = 10000  # maximum number of iterations
    trajectory_dict = perform_runge_kutta(mass, charge, t0,
                                          particle_initialvec, dt, max_iter)

    # convert spherical coordinate trajectory arrays back to cartesian for plotting
    (trajectory_t, trajectory_x, trajectory_y,
     trajectory_z) = get_trajectory_cartesian(trajectory_dict)

    print(trajectory_x, trajectory_y, trajectory_z)

    # plot the 3d plot and the projections onto each plane as well
    plot_trajectories(trajectory_t, trajectory_x, trajectory_y, trajectory_z)

# import plotly.graph_objects as go

# sys.path.append(os.getcwd())
# sys.path.append(os.path.join(os.getcwd(), "gtracr"))
# sys.path.append(os.path.join(os.getcwd(), "..", "gtracr"))

# # from gtracr.trajectory import ParticleTrajectory
# # from gtracr.utils import spherical_to_cartesian
# from gtracr.constants import EARTH_RADIUS, DEG_TO_RAD, RAD_TO_DEG
# from gtracr.utils import *
# from gtracr.runge_kutta import runge_kutta, euler
# from gtracr.add_particle import particle_dict

# def func(latitude, longitude, start_altitude, stop_altitude, zenith, azimuth, particle, energy):

#     # transform the variables into ECEF coordinates

#     r0 = EARTH_RADIUS + start_altitude
#     # longitude and latitude to (r0, theta0, phi0)
#     theta0 = ((90. -  latitude) * DEG_TO_RAD)  # theta defined in [0, pi], theta = 0 at equator
#     # phi = ((180. + longitude) * DEG_TO_RAD)
#     phi0 = ((longitude) * DEG_TO_RAD)  # phi defined in [-pi, pi], phi = 0 at prime meridian

#     x0 = r0*np.sin(theta0)*np.cos(phi0)
#     y0 = r0*np.sin(theta0)*np.sin(phi0)
#     z0 = r0*np.cos(theta0)

#     # print(r0, theta0, phi0)

#     # RK variables
#     maxStep = 10000
#     stepSize = 0.1
#     key_list = ["t", "r", "theta", "phi", "pr", "ptheta", "pphi"]
#     results = {key: np.zeros(maxStep) for key in key_list}

#     # now specify some zenith and azimuthal angle
#     # zenith = 160.
#     # azimuth = 270.

#     # Rotation matrix from tangent plane coordinate system to geocentric coordinate system
#     # source: http://archive.psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf
#     row1 = np.array([-np.sin(longitude*DEG_TO_RAD), -np.cos(longitude*DEG_TO_RAD)*np.sin(latitude*DEG_TO_RAD), np.cos(latitude*DEG_TO_RAD)*np.cos(longitude*DEG_TO_RAD)])
#     row2 = np.array([np.cos(longitude*DEG_TO_RAD), -np.sin(latitude*DEG_TO_RAD)*np.sin(longitude*DEG_TO_RAD), np.cos(latitude*DEG_TO_RAD)*np.sin(longitude*DEG_TO_RAD)])
#     row3 = np.array([0., np.cos(latitude*DEG_TO_RAD), np.sin(latitude*DEG_TO_RAD)])

#     # row1 = np.array([np.sin(latitude+90)*np.cos(), -np.cos(longitude)*np.sin(latitude), np.cos(latitude)*np.cos(longitude)])
#     # row2 = np.array([np.cos(longitude), -np.sin(latitude)*np.sin(longitude), np.cos(latitude)*np.sin(longitude)])
#     # row3 = np.array([0., np.cos(latitude), np.sin(latitude)])

#     tf_matrix = np.array([row1, row2, row3])

#     # get the cartesian coordinates in the local tangent plane coordinate system from zenith and azimuth
#     l = stepSize
#     xi = zenith * DEG_TO_RAD
#     alpha = azimuth * DEG_TO_RAD

#     xt = l*np.sin(xi)*np.cos(alpha)
#     yt = l*np.sin(xi)*np.sin(alpha)
#     zt = l*np.cos(xi)

#     # transform zenith and azimuthal vector coordinates to ECEF
#     # print(np.shape(np.dot(tf_matrix, np.array([xt, yt, zt]))), np.dot(tf_matrix, np.array([xt, yt, zt])))
#     ecef_coords = np.array([x0, y0, z0]) + np.dot(tf_matrix, np.array([xt, yt, zt]))

#     x = ecef_coords[0]
#     y = ecef_coords[1]
#     z = ecef_coords[2]
#     # print(EARTH_RADIUS)
#     print(ecef_coords)

#     r = np.sqrt(x**2. + y**2. + z**2.)
#     # theta = np.arctan((np.sqrt(x**2. + y**2.)) / z)
#     theta = np.arccos(z / r)
#     phi = np.arctan(y / x)

#     # then r1, the next part of the integration step, would be the transformed coordinate
#     # d = np.sin(xi) * np.cos(alpha)

#     # print(d, l, xi, alpha)

#     # phi = phi0 - np.arctan2(d, r0 * np.tan(theta0))
#     # # rcos(phi), obtained from cosine law
#     # rcosphi = np.sqrt((r0 * np.cos(phi0))**2. +
#     #                   (np.cos(alpha))**2. +
#     #                   2. * r0 * np.cos(alpha) * np.cos(phi0)*np.cos(xi))
#     # r = rcosphi / np.cos(phi)
#     # theta = theta0 - (d / rcosphi)

#     # r1 = np.sqrt(r0**2. + 1. + 2.*r0*np.cos(xi))
#     # theta1 = theta0 - np.arcsin(np.sin(xi) / r1)
#     # phi1 = phi0

#     # rtest = np.sqrt(r0**2. + 1.)
#     # thetatest = theta0 - np.arccos(r0 / r1)
#     # print(r0, theta0)
#     # print(r1, theta1)
#     # print(rtest, thetatest)
#     # l = 1.
#     # # if xi > (-np.pi / 2.) and xi < (np.pi / 2.):
#     # #     sign = 1.
#     # # elif xi > (np.pi / 2.) and xi < ((3.*np.pi) / 2.):
#     # #     sign = -1.
#     # # r = sign*np.sqrt(r0**2. + l**2. + 2.*r0*np.cos(xi))

#     # r = np.sqrt(r0**2. + l**2. + 2.*r0*np.cos(xi))

#     # theta = theta0 - np.arcsin((l*np.sin(xi)) / r)

#     # phi = phi0

#     # cos_ag = ((r*np.sin(theta))**2. - (r0*np.sin(theta0))**2. - (l*np.sin(xi)**2.)**2.) / (-2.*r0*np.sin(theta0)*l*np.sin(theta)**2.)
#     # while np.isnan(np.arccos(cos_ag)):
#     #     if cos_ag > 1:
#     #         cos_ag -= 1
#     #     elif cos_ag < -1:
#     #         cos_ag += 1
#     #     print(cos_ag)
#     # gmma = np.arccos(cos_ag) - alpha

#     # phi = phi0 + np.arcsin((l*np.sin(xi)**2.*np.sin(alpha+gmma)) / (r*np.sin(theta)))
#     # if np.sin(theta0) < 1e-10 or np.sin(theta) < 1e-10:
#     #     phi = 0.
#     # else:
#     #     # side1 = np.hypot(r0*np.sin(theta0), r*np.sin(theta))
#     #     # side2 = np.hypot(l*np.sin(xi)*np.sin(theta0), r*np.sin(theta))
#     #     # side3 = np.hypot(l*np.sin(xi)*np.sin(theta0), r0*np.sin(theta0))
#     #     # cos_delphi = (-(side1)**2. + (side2)**2. + (side3)**2.) / (2.*side2*side3)
#     #     cos_delphi = ((-l*np.sin(xi)*np.sin(theta0))**2. + (r0*np.sin(theta0))**2. + (r*np.sin(theta))**2.) / (2.*r0*np.sin(theta0)*l*np.sin(theta))
#     #     print(cos_delphi)
#     #     delphi = np.arccos(cos_delphi)
#     #     while np.isnan(np.arccos(cos_delphi)):
#     #         if cos_delphi > 1:
#     #             cos_delphi -= 1
#     #         elif cos_delphi < -1:
#     #             cos_delphi += 1
#     #     print(cos_delphi)
#     #     delphi = np.arccos(cos_delphi)
#     #     phi = phi0 + delphi

#     print(r, theta, phi)
#     # print(phi)

#     # print(r, theta, phi, rcosphi)

#     # set particle momenta and get their spherical coordinate equivalents
#     particle.set_momentum_from_energy(energy)
#     # # print(particle.momentum)

#     # particle momentum is known in local tangent plane coordinates
#     # so we should transform this into ECEF coordinates

#     # components in local tangent plane coordinates
#     # (pr0, ptheta0, pphi0) = vp_components_spherical(particle.momentum, r, theta)

#     # particle.set_velocity()

#     # (vr0, vtheta0, vphi0) = vp_components_spherical(particle.velocity, r, theta)

#     # vx0 = vr0*np.sin(theta0)*np.cos(phi0) + r0*vtheta0*np.cos(theta0)*np.cos(phi0) - r0*vphi0*np.sin(theta0)*np.sin(phi0)
#     # vy0 = vr0*np.sin(theta0)*np.sin(phi0) + r0*vtheta0*np.cos(theta0)*np.sin(phi0) +r0*vphi0*np.sin(theta0)*np.cos(phi0)
#     # vz0 = vr0*np.cos(theta0) - r0*vtheta0*np.sin(theta0)

#     # # transform zenith and azimuthal vector coordinates to ECEF
#     # ecef_coords = np.array([x0, y0, z0]) + np.dot(tf_matrix, np.array([vx0, vy0, vz0]))

#     # vx = ecef_coords[0]
#     # vy = ecef_coords[1]
#     # vz = ecef_coords[2]

#     # # r = np.sqrt(vx**2. + vy**2. + vz**2.)
#     # # theta = np.arctan((np.sqrt(vx**2. + vy**2.)) / vz)
#     # # phi = np.arctan(vy / vx)

#     # vr = vx*np.sin(theta)*np.cos(phi) + vy*np.sin(theta)*np.sin(phi) + vz*np.cos(theta)
#     # vtheta = (vx*np.cos(theta)*np.cos(phi) + vy*np.cos(theta)*np.sin(phi) - vz*np.sin(theta)) / r
#     # vphi = (-vx*np.sin(phi) + vy*np.cos(phi)) / (r*np.sin(theta))

#     # print(particle.velocity)
#     # (pr0, ptheta0, pphi0) = vp_components_spherical(particle.momentum, r0,
#     #                                                 theta0)
#     (pr0, ptheta0, pphi0) = (0., 0., 0.)
#     (pr, ptheta, pphi) = vp_components_spherical(particle.momentum, r, theta)
#     # # (vr, vtheta, vphi) = vp_components_spherical(particle.velocity, r, theta)
#     # l = 1.
#     # pr = np.sqrt(pr0**2. + l**2. + 2.*pr0*np.cos(xi))
#     # ptheta = ptheta0 - np.arcsin((l*np.sin(xi)) / pr)
#     # cos_ag = ((pr*np.sin(theta))**2. - (pr0*np.sin(theta0))**2. - (l*np.sin(xi)**2.)**2.) / (-2.*pr0*np.sin(theta0)*l*np.sin(theta)**2.)
#     # # print(cos_ag-2)
#     # ag = np.arccos(cos_ag)
#     # # print(ag)
#     # if np.isneginf(ag):
#     #     ag = 0.
#     # elif np.isneginf(ag):
#     #     ag = np.pi
#     # elif np.isnan(ag):
#     #     while np.isnan(np.arccos(cos_ag)):
#     #         if cos_ag > 1:
#     #             cos_ag -= 1
#     #         elif cos_ag < -1:
#     #             cos_ag += 1
#     #         # print(cos_ag)
#     #     ag = np.arccos(cos_ag)
#     # gmma = ag - alpha

#     # pphi = pphi0 + np.arcsin((l*np.sin(xi)**2.*np.sin(alpha+gmma)) / (pr*np.sin(theta)))

#     # print(pr, ptheta, pphi)

#     # particle.set_velocity()
#     # (vr0, vtheta0, vphi0) = vp_components_spherical(particle.velocity, r0,
#     #                                                 theta0)

#     # (vr, vtheta, vphi) = vp_components_spherical(particle.velocity, r,
#     #                                                 theta)

#     particle.set_velocity()

#     print(particle.velocity)

#     # vr0 = particle.velocity
#     # vtheta0 = particle.velocity / l
#     # vphi0 = particle.velocity / (l*np.sin(xi)) if np.abs(xi) < 1e-10 else 0.

#     # # vx0 = vr0*np.sin(theta0)*np.cos(phi0)
#     # # vy0 = vr0*np.sin(theta0)*np.sin(phi0)
#     # # vz0 = vr0*np.cos(theta0)

#     # vx0 = vr0*np.sin(theta0)*np.cos(phi0) + r0*vtheta0*np.cos(theta0)*np.cos(phi0) - r0*vphi0*np.sin(theta0)*np.sin(phi0)
#     # vy0 = vr0*np.sin(theta0)*np.sin(phi0) + r0*vtheta0*np.cos(theta0)*np.sin(phi0) +r0*vphi0*np.sin(theta0)*np.cos(phi0)
#     # vz0 = vr0*np.cos(theta0) - r0*vtheta0*np.sin(theta0)

#     vxt = particle.velocity*np.sin(xi)*np.cos(alpha)
#     vyt = particle.velocity*np.sin(xi)*np.sin(alpha)
#     vzt = particle.velocity*np.cos(xi)

#     print(np.sqrt(vxt**2. + vyt**2. + vzt**2.))

#     # l = np.sqrt(vx0**2. + vy0**2. + vz0**2.)
#     # xi = zenith * DEG_TO_RAD
#     # alpha = azimuth * DEG_TO_RAD

#     # vxt = l*np.sin(xi)*np.cos(alpha)
#     # vyt = l*np.sin(xi)*np.sin(alpha)
#     # vzt = l*np.cos(xi)

#     # transform zenith and azimuthal vector coordinates to ECEF
#     ecef_coords = np.dot(tf_matrix, np.array([vxt, vyt, vzt]))

#     print(tf_matrix)

#     # the thing above needs fixing, there should be some consideration between the local tangent plane velocities
#     # with them being velocities themselves, some time derivative stuff should be considered
#     # this could also resolve the translation issue

#     print(ecef_coords)

#     vx = ecef_coords[0]
#     vy = ecef_coords[1]
#     vz = ecef_coords[2]

#     # r = np.sqrt(vx**2. + vy**2. + vz**2.)
#     # theta = np.arctan((np.sqrt(vx**2. + vy**2.)) / vz)
#     # phi = np.arctan(vy / vx)

#     vr = vx*np.sin(theta)*np.cos(phi) + vy*np.sin(theta)*np.sin(phi) + vz*np.cos(theta)
#     vtheta = (vx*np.cos(theta)*np.cos(phi) + vy*np.cos(theta)*np.sin(phi) - vz*np.sin(theta)) / r
#     vphi = (-vx*np.sin(phi) + vy*np.cos(phi)) / (r*np.sin(theta))

#     print(vr, vtheta, vphi)

#     (pr, ptheta, pphi) = gamma(vmag_spherical(vr, vtheta, vphi, r, theta))*particle.mass*np.array([vr, vtheta, vphi])

#     particle.set_momentum_from_velocity()

#     # (pr, ptheta, pphi) = vp_components_spherical(particle.momentum, r, theta)

#     # (vr, vtheta, vphi) = vp_components_spherical(particle.velocity, r, theta)
#     # (vr, vtheta, vphi) = vp_components_spherical(particle.velocity, r, theta)
#     # l = 1.
#     # vr = np.sqrt(vr0**2. + l**2. + 2.*vr0*np.cos(xi))
#     # # vr = sign*np.sqrt(vr0**2. + l**2. + 2.*vr0*np.cos(xi))
#     # vtheta = vtheta0 - np.arcsin((l*np.sin(xi)) / vr)

#     # vphi = vphi0

#     # cos_ag = ((vr*np.sin(theta))**2. - (vr0*np.sin(theta0))**2. - (l*np.sin(xi)**2.)**2.) / (2.*vr0*np.sin(theta0)*l*np.sin(theta)**2.)
#     # # print(cos_ag-2)
#     # ag = np.arccos(cos_ag)
#     # print(ag)
#     # if np.isneginf(ag):
#     #     ag = 0.
#     # elif np.isneginf(ag):
#     #     ag = np.pi
#     # elif np.isnan(ag):
#     #     while np.isnan(np.arccos(cos_ag)):
#     #         if cos_ag > 1:
#     #             cos_ag -= 1
#     #         elif cos_ag < -1:
#     #             cos_ag += 1
#     #         # print(cos_ag)
#     #     ag = np.arccos(cos_ag)

#     # # print(ag)
#     # # while np.isnan(np.arccos(cos_ag)):
#     # #     if cos_ag > 1:
#     # #         cos_ag -= 1
#     # #     elif cos_ag < -1:
#     # #         cos_ag += 1
#     # #     print(cos_ag)
#     # gmma = ag - alpha

#     # vphi = vphi0 + np.arcsin(-(l*np.sin(xi)**2.*np.sin(alpha+gmma)) / (vr*np.sin(theta)))

#     # if np.sin(theta0) < 1e-10 or np.sin(theta) < 1e-10:
#     #     vphi = 0.
#     # else:
#     #     # side1 = np.hypot(vr0*np.sin(theta0), vr*np.sin(theta))
#     #     # side2 = np.hypot(l*np.sin(xi)*np.sin(theta0), vr*np.sin(theta))
#     #     # side3 = np.hypot(l*np.sin(xi)*np.sin(theta0), vr0*np.sin(theta0))
#     #     # cos_delphi = (-(side1)**2. + (side2)**2. + (side3)**2.) / (2.*side2*side3)
#     #     cos_delphi = (-(l*np.sin(xi)*np.sin(theta0))**2. + (vr0*np.sin(theta0))**2. + (vr*np.sin(theta))**2.) / (2.*vr0*np.sin(theta0)*l*np.sin(theta))
#     #     print(cos_delphi)
#     #     delphi = np.arccos(cos_delphi)
#     #     while np.isnan(np.arccos(cos_delphi)):
#     #         if cos_delphi > 1:
#     #             cos_delphi -= 1
#     #         elif cos_delphi < -1:
#     #             cos_delphi += 1
#     #     delphi = np.arccos(cos_delphi)
#     #     vphi = vphi0 + delphi

#     print(vr, vtheta, vphi)
#     # (pr0, pth0, pph0) = get_sphcomp_momentum(p0, r0, theta0)
#     # (vr, vtheta, vphi) = vp_components_spherical(particle.velocity, r, theta)

#     # initial time
#     t0 = 0.
#     t = t0 + stepSize
#     # print(pr0, ptheta0, pphi0)
#     # print(pr, ptheta, pphi)
#     # print(vr, vtheta, vphi)

#     # append the initial value (the origin in local reference frame) and the first value after zenith / azimuth transformation
#     results["t"][0:2] = [t0, t]
#     results["r"][0:2] = [r0, r]
#     results["theta"][0:2] = [theta0, theta]
#     results["phi"][0:2] = [phi0, phi]
#     results["pr"][0:2] = [pr0, pr]
#     results["ptheta"][0:2] = [ptheta0, ptheta]
#     results["pphi"][0:2] = [pphi0, pphi]

#     # perform the integration
#     initial_values = (t, r, theta, phi, vr, vtheta, vphi)

#     # initial_values = (t, r0, theta0, phi0, vr0, vtheta0,
#     #         vphi0)

#     # print(initial_values)

#     i = 2
#     while i < maxStep - 2:

#         # # get the RK variables
#         (t, r, theta, phi, vr, vtheta,
#          vphi) = runge_kutta(particle, stepSize, initial_values)

#         # get the RK variables
#         # (t, r, theta, phi, vr, vtheta,
#         #  vphi) = euler(particle, stepSize, initial_values)

#         # convert velocities to momenta
#         vmag = vmag_spherical(vr, vtheta, vphi, r, theta)
#         particle.momentum = gamma(vmag) * particle.mass * vmag
#         # particle.momentum = vmag_spherical(pr, ptheta, pphi, r, theta)
#         particle.set_rigidity_from_momentum()
#         # print(particle.rigidity)

#         # (pr, ptheta, pphi) = vp_components_spherical(particle.momentum, r, theta)
#         (pr, ptheta, pphi) = gamma(vmag)*particle.mass*np.array([vr, vtheta, vphi])

#         # print(pr, ptheta, pphi)

#         # append
#         results["t"][i] = t
#         results["r"][i] = r
#         results["theta"][i] = theta
#         results["phi"][i] = phi
#         results["pr"][i] = pr
#         results["ptheta"][i] = ptheta
#         results["pphi"][i] = pphi

#         # if particle has escaped
#         if r > EARTH_RADIUS + stop_altitude:
#             break

#         # if particle touched Earth again
#         if r < EARTH_RADIUS:
#             break

#         # initial_values = (t, r1, theta1, phi1, vr1, vtheta1,
#         #     vphi1)

#         initial_values = (t, r, theta, phi, vr, vtheta, vphi)

#         i += 1

#     # trim all zeros
#     for key, arr in list(results.items()):
#         results.update({key: np.trim_zeros(arr, 'b')})

#     # plot
#     t = results["t"]
#     (x, y, z) = spherical_to_cartesian(results["r"] / EARTH_RADIUS,
#                                        results["theta"], results["phi"])
#     return results

# if __name__ == "__main__":

#     # some initial variables
#     latitude = 0.
#     longitude = -60.
#     start_altitude = 20.
#     stop_altitude = 565.

#     particle = particle_dict["p-"]
#     energy = 12.

#     result = func(latitude, longitude, start_altitude, stop_altitude, 0., 0., particle, energy)
#     # result2 = func(0., 0., 20., 565., 70., 90., particle, energy)

#     # plot
#     t = result["t"]
#     (x, y, z) = spherical_to_cartesian(result["r"] / EARTH_RADIUS,
#                                        result["theta"], result["phi"])

#     # t2 = result2["t"]
#     # (x2, y2, z2) = spherical_to_cartesian(result2["r"] / EARTH_RADIUS,
#     #                                    result2["theta"], result2["phi"])

#     plt.scatter(x, y, c="k")
#     # plt.scatter(x2, y2, c="r")
#     # plt.xlim([-.5, .5])
#     # plt.ylim([-.5, .5])
#     plt.show()

#     plt.scatter(x, z, c="k")
#     # plt.scatter(x2, z2, c="r")
#     # plt.xlim([-.5, .5])
#     # plt.ylim([-.5, .5])
#     plt.show()

#     plt.scatter(y, z, c="k")
#     # plt.scatter(y2, z2, c="r")
#     # plt.xlim([-.5, 1.5])
#     # plt.ylim([-1.5, 1.5])
#     plt.show()

#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection="3d")
#     cm = ax.scatter(x, y, z, c="k", marker='o')
#     # cm = ax.scatter(x2, y2, z2, c="r", marker='o')
#     fig.colorbar(cm, ax=ax)
#     # plt.colorbar()
#     # plt.savefig("test.png")
#     plt.show()

#     result1 = func(0., 0., 20., 565., 90., 0., particle_dict["p+"], 30.)
#     result2 = func(0., 0., 20., 565., 90., 180., particle_dict["p+"], 30.)

#     # plot
#     t1 = result1["t"]
#     (x1, y1, z1) = spherical_to_cartesian(result1["r"] / EARTH_RADIUS,
#                                        result1["theta"], result1["phi"])

#     t2 = result2["t"]
#     (x2, y2, z2) = spherical_to_cartesian(result2["r"] / EARTH_RADIUS,
#                                        result2["theta"], result2["phi"])

#     plt.scatter(x1, y1, c="k")
#     plt.scatter(x2, y2, c="r")
#     # plt.xlim([-1.5, 1.5])
#     # plt.ylim([-1.5, 1.5])
#     plt.show()

#     plt.scatter(x1, z1, c="k")
#     plt.scatter(x2, z2, c="r")
#     # plt.xlim([-1.5, 1.5])
#     # plt.ylim([-1.5, 1.5])
#     plt.show()

#     plt.scatter(y1, z1, c="k")
#     plt.scatter(y2, z2, c="r")
#     # plt.xlim([-1.5, 1.5])
#     # plt.ylim([-1.5, 1.5])
#     plt.show()

#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection="3d")
#     cm = ax.scatter(x1, y1, z1, c="k", marker='o')
#     cm = ax.scatter(x2, y2, z2, c="r", marker='o')
#     fig.colorbar(cm, ax=ax)
#     # plt.colorbar()
#     # plt.savefig("test.png")
#     plt.show()

# '''
#         # print(pr, ptheta, pphi)

#         # particle.set_velocity()

#         # # set particle momenta and get their spherical coordinate equivalents
#         # # particle.set_momentum(energy)
#         # # particle.set_velocity()
#         # (pr1, ptheta1,pphi1) = vp_components_spherical(particle.momentum, r1, theta1)
#         # # (pr, ptheta,pphi) = vp_components_spherical(particle.momentum, r, theta)
#         # # (pr0, pth0, pph0) = get_sphcomp_momentum(p0, r0, theta0)
#         # (vr1, vtheta1, vphi1) = vp_components_spherical(particle.velocity, r1, theta1)
#         # # initial time
#         # # t0 = 0.
#         # t = t + stepSize

#         # print(vr1, vtheta1, vphi1)

#         # # convert velocities to momenta
#         # vmag = vmag_spherical(vr, vtheta, vphi, r, theta)
#         # particle.momentum = gamma(vmag) * particle.mass * vmag
#         # # particle.momentum = vmag_spherical(pr, ptheta, pphi, r, theta)
#         # particle.set_rigidity_from_momentum()
#         # # print(particle.rigidity)

#         # (pr, ptheta,
#         #     pphi) = vp_components_spherical(particle.momentum, r, theta)

#         # print(pr, ptheta, pphi)

#         # (t, r, theta, phi, vr, vtheta,
#         #     vphi) = euler(particle, stepSize, initial_values)

#         # print((t, r, theta, phi, vr, vtheta, vphi))

#         # # xi = np.arccos((r-EARTH_RADIUS))
#         # xi = (np.pi / 2.) - np.sin(theta0 - theta)

#         # r1 = np.sqrt(r**2. + 10.**2. + 2.*r*np.cos(xi))
#         # theta1 = theta - np.arcsin((10.*np.sin(xi)) / r)
#         # phi1 = phi

#         # print(r1, theta1, phi1)

#         # print(r, theta, phi, rcosphi)
# '''
