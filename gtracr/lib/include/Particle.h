// Header file for Particle class
#ifndef __PARTICLE_H_
#define __PARTICLE_H_

class Particle
{
private:
    std::string nm;  // Name
    int pid;         // PDG ID
    double m;        // mass
    int ch;          // charge
    std::string lbl; // label

    double p; // momentum
    double v; // velocity
    double R; // rigidity

public:
    // constructor
    Particle();
    Particle(const std::string &, const int, const double &, const int, const std::string &);                                 // no energy
    Particle(const std::string &, const int, const double &, const int, const std::string &, const double &, const double &); // with energy / rigidity
    // destructor
    ~Particle();
    // copy constructor / operator
    Particle(const Particle &);
    Particle &operator=(const Particle &);
    // getters
    const std::string &name() const { return nm; }
    const double &mass() const { return m; }
    const int charge() const { return ch; }
    const int pdgid() const { return pid; }
    const std::string &label() const { return lbl; }
    const double &momentum() const { return p; }
    const double &velocity() const { return v; }
    const double &rigidity() const { return R; }
    // setters
    void set_name(const std::string &_nm) { nm = _nm; }
    void set_mass(const double &_m) { m = _m; }
    void set_charge(const double &_ch) { ch = _ch; }
    void set_pdgid(const int _pid) { pid = _pid; }
    void set_label(const std::string &_lbl) { lbl = _lbl; }
    void set_momentum(const double &_p) { p = _p; }
    void set_velocity(const double &_v) { v = _v; }
    void set_rigidity(const double &_R) { R = _R; }
    // setters (only for p, v, R for now)
    void set_from_energy(const double &);
    void set_from_momentum(const double &);
    void set_from_rigidity(const double &);
    void set_from_velocity(const double &);
    // utility function
    const double &gamma();               // using object velocity
    const double &gamma(const double &); // user input velocity
    void print();
    // other member functions
    double get_energy_rigidity();
};

#endif //__PARTICLE_H_