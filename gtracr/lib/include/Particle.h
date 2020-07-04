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
    Particle(const std::string &, const int, const double &, const int, const std::string &);                 // no energy
    Particle(const std::string &, const int, const double &, const int, const std::string &, const double &); // energy
    // destructor
    ~Particle();
    // copy constructor / operator
    Particle(const Particle &);
    Particle &operator=(const Particle &);
    // getters
    const std::string &name() { return nm; }
    const double &mass() { return m; }
    const int charge() { return ch; }
    const int pdgid() { return pid; }
    const std::string &label() { return lbl; }
    const double &momentum() { return p; }
    const double &velocity() { return v; }
    const double &rigidity() { return R; }
    // setters (only for p, v, R for now)
    void set_from_energy(const double&);
    void set_from_momentum(const double&);
    void set_from_rigidity(const double&);
    void set_from_velocity(const double&);
    // utility function
    const double &gamma();  // using object velocity
    const double &gamma(const double&); // user input velocity
    // other member functions 
    const double &get_energy_rigidity();
};

#endif //__PARTICLE_H_