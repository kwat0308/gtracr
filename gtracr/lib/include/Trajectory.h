// Header file for Trajectory class
#ifndef __TRAJECTORY_H_
#define __TRAJECTORY_H_

class Trajectory {
    private:
        // main members
        Particle part;   // particle
        double lat; // latitude
        double lng; //longitude
        double alt; //altitude
        double zang;    // zenith angle
        double azang;   // azimuth angle
        double en;    // energy
        double rg;    // rigidity
        double escalt;   // escape altitude
        int maxbuf;    // max buffer
        bool particleEscaped;  // if particle escaped or not
        // double *timearr;  // pointer to time array
        // TrajectoryPoint *TJParr;  // pointer to TrajectoryPoint array
        // // members to store plotting data
        // double *xarr;
        // double *yarr;
        // double *zarr;
        std::vector<double> timearr;
        std::vector<TrajectoryPoint> TJParr;

    public:
        // Constructors
        Trajectory();
        // using location name
        Trajectory(const std::string &, const std::string &, const double&, const double&, const double &, const double&, const double &, const int); 
        // using altitude / latitude / longitude
        Trajectory(const std::string &, const double&, const double&, const double&, const double&, const double&, const double &, const double&, const double &, const int); 
        // constructor utility functions
        void init_array();  // initialize array
        // constructor function, initializes object
        // void init(const char*, const char*, const double&, const double&, const double &);
        // set kinematics based on rigidity and energy
        void set_kinematics(const double&, const double&);
        // destructor
        ~Trajectory();
        // copy constructor / assignment operator
        Trajectory(const Trajectory&);
        Trajectory &operator=(const Trajectory&);
        // move constructor / assignment operator
        // Trajectory(Trajectory &&);
        // Trajectory &operator=(Trajectory &&);
        // main functions
        std::pair<TrajectoryPoint, TrajectoryPoint> getInitTJP(); // get origin and first TJP from unit vector
        void getTrajectory(const int, const double&); // get the trajectory
        std::map<std::string, std::vector<double> > getPlottingVariables();
        // utility functions
        // convert from local tangent plane to geocentric coordinates
        std::vector<double> LTP_to_ECEF(std::vector<double>, std::vector<double>);
        // get cartesian coordinates in local tangent plane
        std::vector<double> get_LTPvec(const double&);
        // the transformation matrix
        Matrix tf_matrix();
        // convert cartesian to spherical coordinates
        TrajectoryPoint convert_spherical(const std::vector<double>&, const std::vector<double>&);
        // print contents
        void print();
};

#endif //__TRAJECTORY_H_