// Header file for Trajectory Point class
#ifndef __TRAJECTORYPOINT_H_
#define __TRAJECTORYPOINT_H_

class TrajectoryPoint {
    private:
        double lat;
        double lng;
        double alt;
        double v_r;
        double v_theta;
        double v_phi;
    public:
        // constructor
        TrajectoryPoint();
        TrajectoryPoint(const double&, const double&, const double&, const double&, const double&, const double&);
        TrajectoryPoint(const Location&, const double&, const double&, const double&);
        // copy constructor
        TrajectoryPoint(const TrajectoryPoint&);
        TrajectoryPoint &operator=(const TrajectoryPoint&);
        // getters
        const double &latitude() {return lat;}
        const double &longitude() {return lng;}
        const double &altitude() {return alt;}
        const double &vr() {return v_r;}
        const double &vtheta() {return v_theta;}
        const double &vphi() {return v_phi;}
        // setters
        void set_latitude(const double &_lat) {lat = _lat;}
        void set_longitude(const double &_lng) {lng = _lng;}
        void set_altitude(const double &_alt) {alt = _alt;}
        void set_vr(const double &_v_r) {v_r = _v_r;}
        void set_vtheta(const double &_v_theta) {v_theta = _v_theta;}
        void set_vphi(const double &_v_phi) {v_phi = _v_phi;}
        // other getters and setters
        double magnitude();
        double velocity();
        // coordinate conversions
        std::vector<double> getSphericalCoord();
        void setSphericalCoord(const double&, const double&, const double&);
        std::vector<double> getCartesianCoord();
        // void setCartesianCoord(const double&, const double&, const double&);
        std::vector<double> getCartesianVelocity();
        void setCartesianVelocity(const double&, const double&, const double&);
        std::vector<double> spherical();
        std::vector<double> cartesian();
        // utility functions
        void print();
};

#endif //__TRAJECTORYPOINT_H_