// Class that controls trajectory of a particle at some given energy and rigidity

/*
    Class that controls the trajectory of a particle at some given energy / rigidity
    Members:
    - particleName: the label of the particle
    - latitude: the geographic latitude, with 0 defined at the equator in degrees
    - longitude: the geographic longitude, with 0 defined at the Prime Meridian in degrees
    - altitude: the height from sea level (0=sea level) in km
    - zenithAngle: the angle from the local zenith, with 0 being at the local zenith
    - azimuthAngle: the angle with 0 being in the direction of the East hemisphere from the Prime Meridian in the local tangent plane
    - energy: the particle energy
    - rigidity: the particle rigidity (momentum / charge)
    - escapeAltitude: the altitude in which the particle has "escaped" Earth (default 1000km)
    - maxBuffer: maximum length of array
*/
#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <string>
#include <map>
#include <math.h>
// #include <pybind11/pybind11.h>
// #include <pybind11/numpy.h>

#include "constants.h"
#include "Matrix.h"
// #include "Vector.h"
#include "Particle.h"
#include "Location.h"
#include "particlemap.h"
#include "locationmap.h"
#include "RungeKutta.h"
#include "TrajectoryPoint.h"
#include "Trajectory.h"

// Constructors
// Default Constructor
Trajectory::Trajectory()
    : part{Particle()}, lat{0.}, lng{0.}, alt{0.}, zang{0.}, azang{0.}, en{10.}, rg{0.}, escalt{1000.}, maxbuf{100000}, particleEscaped{false}
{
    // void init_array();
    part.set_from_energy(en);
    // part.print();
}

// Using location name
Trajectory::Trajectory(const std::string &plabel, const std::string &locname, const double &zenithAngle=0.0, const double &azimuthAngle=0.0, const double &energy = 0., const double &rigidity = 0., const double &escapeAltitude = 1000., const int maxBuffer = 100000)
    : zang{zenithAngle}, azang{azimuthAngle}, escalt{escapeAltitude}, maxbuf{maxBuffer}
{
    // part.print();
    // partmap[plabel].print();
    // Particle &p = partmap[plabel];
    // p.print();
    // part = &p;
    part = partmap[plabel];
    // part.print();
    Location loc = locmap[locname];
    lat = loc.latitude();
    lng = loc.longitude();
    alt = loc.altitude();
    set_kinematics(energy, rigidity);
}

// Using manual longitude / latitude / altitude
Trajectory::Trajectory(const std::string &plabel, const double &latitude=0., const double &longitude=0., const double &altitude=0., const double &zenithAngle=0., const double &azimuthAngle=0., const double &energy = 0., const double &rigidity = 0., const double &escapeAltitude = 1000., const int maxBuffer = 100000)
    : zang{zenithAngle}, azang{azimuthAngle}, escalt{escapeAltitude}, maxbuf{maxBuffer}
{
    // part.print();
    // partmap[plabel].print();
    // Particle &p = partmap[plabel];
    // p.print();
    // part = &p;
    part = partmap[plabel];
    // part.print();
    lat = latitude;
    lng = longitude;
    alt = altitude;
    // std::cout << lat << lng << alt;
    set_kinematics(energy, rigidity);
}

// set initial kinematics based on energy / rigidity
void Trajectory::set_kinematics(const double& energy, const double& rigidity)
{
    part.print();
    // use tolerance values to compare val == 0
    if (abs(energy) < 1e-10) {
        // rigidity is given as input
        part.set_from_rigidity(rigidity);
        // part.print();
        en = part.get_energy_rigidity();
        // part.print();
    }
    else if (abs(rigidity) < 1e-10) {
        // energy is given as input
        part.set_from_energy(energy);
        rg = part.rigidity();
        // part.print();
    }
    else {
        throw std::runtime_error("Provide energy or rigidity as inputs!");
    }

    // part.print();
}
// initialize arrays
// void Trajectory::init_array()
// {
//     // for (int i = 0; i < maxbuf; ++i){
//     //     timearr[i] = 0.0;
//     //     TJParr[i] = TrajectoryPoint();
//     //     xarr[i] = 0.0;
//     //     yarr[i] = 0.0;
//     //     zarr[i] = 0.0;
//     // }
// }

// // function for initializing object
// void init(const char *plabel, const double &latitude, const double &longitude, const double &altitude, const double &zenithAngle, const double &azimuthAngle)
// {

// }

// Destructor
Trajectory::~Trajectory()
{
    // delete[] timearr;
    // delete[] TJParr;
    // delete[] xarr;
    // delete[] yarr;
    // delete[] zarr;
    // delete part;
}

// copy constructor
Trajectory::Trajectory(const Trajectory &traj)
    :lat{traj.lat}, lng{traj.lng}, alt{traj.alt}, zang{traj.zang}, azang{traj.azang}, en{traj.en}, rg{traj.rg}, escalt{traj.escalt}, maxbuf{traj.maxbuf}, particleEscaped{traj.particleEscaped}
{
    // delete part;
    // part = new Particle();
    part = traj.part;
    timearr = traj.timearr;
    TJParr = traj.TJParr;
}

// copy assignment operator
Trajectory &Trajectory::operator=(const Trajectory &traj)
{
    lat = traj.lat;
    lng = traj.lng;
    alt = traj.alt;
    zang = traj.zang;
    azang = traj.azang;
    en = traj.en;
    rg = traj.rg;
    escalt = traj.escalt;
    maxbuf = traj.maxbuf;
    particleEscaped = traj.particleEscaped;
    timearr = traj.timearr;
    TJParr = traj.TJParr;
    // delete part;
    // part = new Particle();
    part = traj.part;
    return *this;
}

// // Move constructor
// Trajectory::Trajectory(Trajectory &&traj)
// {
// }

// // Move assignment operator
// Trajectory &Trajectory::operator=(Trajectory &&traj)
// {
// }

// get the origin TJP and the initial TJP after considering zenith and azimuth direction
std::pair<TrajectoryPoint, TrajectoryPoint> Trajectory::getInitTJP()
{
    TrajectoryPoint origin_TJP = TrajectoryPoint(lat, lng, alt, 0., 0., 0.);
    // origin_TJP.print();
    std::vector<double> origin_coord = origin_TJP.getCartesianCoord();
    // origin_coord.print();
    std::vector<double> LTP_coord = get_LTPvec(1e-8);
    // LTP_coord.print();
    std::vector<double> init_cartcoord = LTP_to_ECEF(origin_coord, LTP_coord);

    // init_cartcoord.print();

    std::vector<double> origin_vel (3, 0.);
    // origin_vel.print();
    std::vector<double> LTP_vel = get_LTPvec(part.velocity());
    // LTP_vel.print();
    std::vector<double> init_cartvel = LTP_to_ECEF(origin_vel, LTP_vel);
    // init_cartvel.print();

    TrajectoryPoint init_TJP = convert_spherical(init_cartcoord, init_cartvel);
    // init_TJP.print();

    std::pair<TrajectoryPoint, TrajectoryPoint> init_tup = {origin_TJP, init_TJP};

    // init_tup.first.print();
    // init_tup.second.print();

    return init_tup;
}

void Trajectory::getTrajectory(const int maxStep=10000, const double& stepSize=0.01)
{  
    // resize vector if maxstep > or < maxbuf
    // at this stage doesnt matter, vector is dynamic array
    // if (maxStep > maxbuf) {
    //     timearr.resize(maxStep);
    // }
    // else if (maxStep < maxbuf) {
    //     timearr.resize()
    // }

    // part.print();

    std::pair<TrajectoryPoint, TrajectoryPoint> init_tup = getInitTJP();
    
    // append initial values to time, TJP array
    timearr.push_back(0.0);
    timearr.push_back(stepSize);
    TJParr.push_back(init_tup.first);
    TJParr.push_back(init_tup.second);

    // std::cout << timearr[0] << timearr[1] << std::endl;
    // TJParr[0].print();
    // TJParr[1].print();

    // start iteration process
    RungeKutta RKI = RungeKutta(part.charge(), part.mass(), stepSize);
    // TrajectoryPoint curr_TJP = init_tup.second;
    // curr_TJP.print();
    double t = stepSize;
    std::vector<double> vec = init_tup.second.spherical();
    auto it = vec.begin();
    vec.insert(it, t);
    // for (int i=0; i<7; ++i) {
    //     std::cout << vec[i] << '\t';
    // }
    std::cout << std::endl;

    for (int i=0; i<maxStep; ++i) {
        // evaluate
        vec = RKI.evaluate(vec);
        // for (double val:vec) {
        //     std::cout << val << '\t';
        // }
        // std::cout << std::endl;
        // create new TJP to append to
        TrajectoryPoint new_TJP = TrajectoryPoint(0., 0., 0., vec[4], vec[5], vec[6]);
        new_TJP.setSphericalCoord(vec[1], vec[2], vec[3]);
        // new_TJP.print();
        // curr_TJP = new_TJP;
        // curr_TJP.print();
       // append to arrays
        // timearr[i] = vec[0];
        // TJParr[i] = new_TJP;
        timearr.push_back(vec[0]);
        TJParr.push_back(new_TJP);

        // break conditions
        if (new_TJP.altitude() > escalt) {
            particleEscaped = true;
            break;
        }

        if (new_TJP.altitude() < 0.) {
            break;
        }
    }

    std::cout << "All iterations completed. " << std::endl;

}

std::map<std::string, std::vector<double> > Trajectory::getPlottingVariables()
{
    int n = timearr.size();
    // initialize vectors
    std::vector<double> xarr(n, 0.);
    std::vector<double> yarr(n, 0.);
    std::vector<double> zarr(n, 0.);
    std::vector<double> vxarr(n, 0.);
    std::vector<double> vyarr(n, 0.);
    std::vector<double> vzarr(n, 0.);

    // append variables
    for (int i=0; i<n; ++i) {
        std::vector<double> vec = TJParr[i].cartesian();

        xarr[i] = vec[0];
        yarr[i] = vec[1];
        zarr[i] = vec[2];
        vxarr[i] = vec[3];
        vyarr[i] = vec[4];
        vzarr[i] = vec[5];
    }

    // create map and insert vectors into there
    std::map<std::string, std::vector<double> > m = {
        {"t", timearr}, {"x", xarr}, {"y", yarr}, {"z", zarr}, {"vx", vxarr}, {"vy", vyarr}, {"vz", vzarr}
    };

    return m;
}

// convert between locatl tangent plane coordinates to geocentric coordinates (Earth-Centered, Earth-Fixed)
std::vector<double> Trajectory::LTP_to_ECEF(std::vector<double> origin, std::vector<double> LTP) 
{

    std::vector<double> prodvec = tf_matrix().dot(LTP);
    std::vector<double> newvec;
    for (int i=0; i<prodvec.size(); ++i) {
        newvec.push_back(prodvec[i] + origin[i]);
        std::cout << newvec[i] << '\t';
    }
    // std::vector<double> newvec = origin + tf_matrix().dot(LTP);
    // tf_matrix().print();
    // newvec.print();
    return newvec;
}

// get cartesian coordinates in local tangent plane
std::vector<double> Trajectory::get_LTPvec(const double& mag)
{
    double xi = zang * constants::DEG_TO_RAD;
    double alpha = azang * constants::DEG_TO_RAD;

    double xt = mag * sin(xi) * cos(alpha);
    double yt = mag * sin(xi) * sin(alpha);
    double zt = mag * cos(xi);

    std::cout << xt << yt << zt;

    std::vector<double> newvec = std::vector<double> {xt, yt, zt};
    // newvec.print();
    // newvec.set_value(0, 0, xt);
    // std::cout << newvec(0,0) << std::endl;
    // newvec.print();
    // newvec.set_value(0, 1, yt);
    // std::cout << newvec(0,1) << std::endl;
    // newvec.print();
    // newvec.set_value(0, 2, zt);
    // std::cout << newvec(0,2) << std::endl;
    // newvec.print();

    return newvec;
}

// transformation matrix between LTP and ECEF
Matrix Trajectory::tf_matrix()
{
    double lmbda = lat * constants::DEG_TO_RAD;
    double eta = lng * constants::DEG_TO_RAD;
    // initialize pointer to fill matrix values with
    // double *dp = new double[9];
    // // fill in the values
    // dp[0] = -sin(eta);
    // dp[1] = -cos(eta)*sin(lmbda);
    // dp[2] = cos(lmbda)*cos(eta);
    // dp[3] = cos(eta);
    // dp[4] = -sin(lmbda)*sin(eta);
    // dp[5] = cos(lmbda)*sin(eta);
    // dp[6] = 0.;
    // dp[7] = cos(lmbda);
    // dp[8] = sin(lmbda);

    // Matrix m = Matrix(3, 3, dp);
    // delete[] dp;

    std::vector<double> newvec = {
        -sin(eta),
        -cos(eta)*sin(lmbda),
        cos(lmbda)*cos(eta),
        cos(eta),
        -sin(lmbda)*sin(eta),
        cos(lmbda)*sin(eta),
        0.,
        cos(lmbda),
        sin(lmbda)
    };
    // fill in the values
    

    Matrix m = Matrix(3, 3, newvec);
    // m.print();
    return m;
}

// convert cartesian coordinates to spherical ones, and then into TJP
TrajectoryPoint Trajectory::convert_spherical(const std::vector<double>& coord, const std::vector<double>& vel) 
{
    double r = sqrt(coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2]);
    double theta = acos(coord[2] / r);
    double phi = atan2(coord[1], coord[0]);

    double vx = vel[0];
    double vy = vel[1];
    double vz = vel[2];

    double vr = vx * sin(theta) * cos(phi) + vy * sin(theta) * sin(phi) + vz * cos(theta);
    double vtheta = (vx * cos(theta) * cos(phi) +
                vy * cos(theta) * sin(phi) - vz * sin(theta)) / r;
    double vphi = (-vx * sin(phi) + vy * cos(phi)) / (r * sin(theta));

    TrajectoryPoint newTJP = TrajectoryPoint(0., 0., 0., vr, vtheta, vphi);

    newTJP.setSphericalCoord(r, theta, phi);

    return newTJP;
}

// print contents of trajectory
void Trajectory::print() 
{   
    std::cout << "Particle: " << part.name() << " "
                << "Particle Energy: " << en << " "
                << "Particle Rigidity: " << rg << std::endl
                << "Latitude: " << lat << " "
                << "Longitude: " << lng << " "
                << "Altitude: " << alt << " "
                << "Angle from zenith: " << zang << " "
                << "Azimuthal angle: " << azang << std::endl;
                
}



