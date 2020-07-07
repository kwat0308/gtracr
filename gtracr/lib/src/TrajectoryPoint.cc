/*
    Class that records a single point in the particle trajectory
    Members:
    - latitude: the geographic latitude, with 0 at the equator (+90 at the North Pole, -90 at South pole) in degrees
    - longitude: the geographic longitude, with 0 at the Prime Meridian (negative degrees towards Western hemisphere) in degrees
    - altitude: the altitude above sea level
    - v_r: the radial velocity
    - v_theta: the velocity in the polar direction
    - v_phi: the velocity in the azimuthal direction
*/

// #include <array>
#include <vector>
#include <iostream>
#include <math.h>
#include "constants.h"
#include "Matrix.h"
#include "Location.h"
#include "TrajectoryPoint.h"

// Constructors
// default constructor
TrajectoryPoint::TrajectoryPoint()
    :lat{0.}, lng{0.}, alt{0.}, v_r{0.}, v_theta{0.}, v_phi{0.}
{

}
// variable location and velocity
TrajectoryPoint::TrajectoryPoint(const double &latitude, const double &longitude, const double &altitude, const double &rvel, const double &thvel, const double &phvel)
    : lat{latitude}, lng{longitude}, alt{altitude}, v_r{rvel}, v_theta{thvel}, v_phi{phvel}
{
}

// initialize with Location object
TrajectoryPoint::TrajectoryPoint(const Location &loc, const double &rvel, const double &thvel, const double &phvel)
    : v_r{rvel}, v_theta{thvel}, v_phi{phvel}
{
    lat = loc.latitude();
    lng = loc.longitude();
    alt = loc.altitude();
}

// copy constructor
TrajectoryPoint::TrajectoryPoint(const TrajectoryPoint &TJP)
    : lat{TJP.lat}, lng{TJP.lng}, alt{TJP.alt}, v_r{TJP.v_r}, v_theta{TJP.v_theta}, v_phi{TJP.v_phi}
{
}

// copy assignment operator
TrajectoryPoint &TrajectoryPoint::operator=(const TrajectoryPoint &TJP)
{
    lat = TJP.lat;
    lng = TJP.lng;
    alt = TJP.alt;
    v_r = TJP.v_r;
    v_theta = TJP.v_theta;
    v_phi = TJP.v_phi;
    return *this;
}

// get the magnitude of the coordinate vector in geocentric coordinates
double TrajectoryPoint::magnitude()
{
    const double &mag = alt - constants::RE;
    return mag;
}

// get the magnitude of the velocity in geocentric coordinates
double TrajectoryPoint::velocity()
{
    std::vector<double> sphvec = getSphericalCoord();
    double r = sphvec[0];
    double theta = sphvec[1];
    double phi = sphvec[2];
    // sphvec[0] = r, sphvec[1] = theta, sphvec[2] = phi
    double vel = sqrt(v_r * v_r + r * r * v_theta * v_theta + r * r * sin(phi) * sin(phi) * v_phi * v_phi);
    return vel;
}

// convert geodesic coordinates to spherical coordinates
std::vector<double> TrajectoryPoint::getSphericalCoord()
{
    double r = constants::RE + alt;
    double theta = (90. - lat) * constants::DEG_TO_RAD;
    double phi = lng * constants::DEG_TO_RAD;
    std::vector<double> sphcoordarr = {r, theta, phi};
    // std::vector<double> sphcoordvec = std::vector<double>(0, 3, sphcoordarr);
    return sphcoordarr;
}

// convert spherical coordinates to geodesic coordinates
void TrajectoryPoint::setSphericalCoord(const double &r, const double &theta, const double &phi)
{
    lat = 90. - (theta * constants::RAD_TO_DEG);
    lng = phi * constants::RAD_TO_DEG;
    alt = r - constants::RE;
}

// convert geodesic to cartesian
std::vector<double> TrajectoryPoint::getCartesianCoord()
{
    double lmbda = lat * constants::DEG_TO_RAD;
    double eta = lng * constants::DEG_TO_RAD;
    double x = (constants::RE + alt) * cos(lmbda) * cos(eta);
    double y = (constants::RE + alt) * cos(lmbda) * sin(eta);
    double z = (constants::RE + alt) * sin(lmbda);
    std::vector<double> carcoordarr = {x, y, z};
    // std::cout << carcoordarr[0] << carcoordarr[1] << std::endl;
    // std::vector<double> carcoordvec = std::vector<double>(0, 3, carcoordarr);
    // carcoordvec.print();
    return carcoordarr;
}

// convert cartesian to geodesic
// void TrajectoryPoint::setCartesianCoord(const double &x, const double &y, const double &z)
// {

// }

// convert spherical velocity to cartesian velocity
std::vector<double> TrajectoryPoint::getCartesianVelocity()
{
    std::vector<double> sphvec = getSphericalCoord();
    double r = sphvec[0];
    double theta = sphvec[1];
    double phi = sphvec[2];
    // sphvec[0] = r, sphvec[1] = theta, sphvec[2] = phi
    double vx = v_r * sin(theta) * cos(phi) + r * v_theta * cos(theta) * cos(phi) - r * v_phi * sin(theta) * sin(phi);
    double vy = v_r * sin(theta) * sin(phi) + r * v_theta * cos(theta) * sin(phi) + r * v_phi * sin(theta) * cos(phi);
    double vz = v_r * cos(theta) - r * v_theta * sin(theta);
    std::vector<double> carvelarr = {vx, vy, vz};
    // std::vector<double> carvelvec = std::vector<double>(0, 3, carvelarr);
    return carvelarr;
}

// convert cartesian velocity to spherical velocity
void TrajectoryPoint::setCartesianVelocity(const double &vx, const double &vy, const double &vz)
{
    std::vector<double> sphvec = getSphericalCoord();
    double r = sphvec[0];
    double theta = sphvec[1];
    double phi = sphvec[2];
    // sphvec[0] = r, sphvec[1] = theta, sphvec[2] = phi
    v_r = vx * sin(theta) * cos(phi) + vy * sin(theta) * sin(phi) + vz * cos(theta);
    v_theta = (vx * cos(theta) * cos(phi) + vy * cos(theta) * sin(phi) - vz * sin(theta)) / r;
    v_phi = (-vx * sin(phi) + vy * cos(phi)) / (r * sin(theta));
}

// get the spherical components (coordinate and velocity)
std::vector<double> TrajectoryPoint::spherical()
{
    std::vector<double> sphcoordvec = getSphericalCoord();
    // sphvec[0] = r, sphvec[1] = theta, sphvec[2] = phi
    std::vector<double> sphvecarr = {sphcoordvec[0], sphcoordvec[1], sphcoordvec[2], v_r, v_theta, v_phi};
    // std::vector<double> sphvec = std::vector<double>(0, 6, sphvecarr);
    return sphvecarr;
}

// get the cartesian components (coordinate and velocity)
std::vector<double> TrajectoryPoint::cartesian()
{
    std::vector<double> carcoordvec = getCartesianCoord();
    std::vector<double> carvelvec = getCartesianVelocity();
    std::vector<double> carvecarr = {carcoordvec[0], carcoordvec[1], carcoordvec[2], carvelvec[0], carvelvec[1], carvelvec[2]};
    // std::vector<double> carvec = std::vector<double>(0, 6, carvecarr);
    return carvecarr;
}

// equivalent to __str__ in Python
// print the class object
void TrajectoryPoint::print()
{
    std::cout << "Latitude: " << lat << ", "
              << "Longitude: " << lng << ", "
              << "Altitude: " << alt << ", "
              << "Velocity (vr, vtheta, vphi): "
              << "( " << v_r << ", "
              << v_theta << ", "
              << v_phi << " )"
              << std::endl;
}
