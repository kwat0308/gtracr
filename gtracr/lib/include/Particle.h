// Header file for Particle class
#ifndef __PARTICLE_H_
#define __PARTICLE_H_

class Particle
{
private:
    const char *nm;  // Name
    int pid;         // PDG ID
    double m;        // mass
    int ch;          // charge
    const char *lbl; // label

    double p; // momentum
    double v; // velocity
    double R; // rigidity

public:
    // constructor
    Particle();
    Particle(const char *, const int, const double &, const int, const char *);                 // no energy
    Particle(const char *, const int, const double &, const int, const char *, const double &); // energy
    // destructor
    ~Particle();
    // copy constructor / operator
    Particle(const Particle &);
    Particle &operator=(const Particle &);
    // getters
    const char *name() { return nm; }
    const double &mass() { return m; }
    const int charge() { return ch; }
    const int pdgid() { return pid; }
    const char *label() { return lbl; }
    const double &momentum() { return p; }
    const double &velocity() { return v; }
    const double &rigidity() { return R; }
    // setters (only for p, v, R for now)
    // void set_momentum_energy(const double&);
    // void set_momentum_rigidity(const double&);
    // void set_momentum_velocity(const double&);
    // void set_velocity_energy(const double&);
    // void set_velocity_momentum(const double&);
    // void set_velocity_rigidity(const double&);
    // void set_rigidity_energy(const double&);
    // void set_rigidity_momentum(const double&);
    // void set_rigidity_velocity(const double&);
    // void set_from_energy(const double&);
    // void set_from_momentum(const double&);
    // void set_from_rigidity(const double&);
    // void set_from_velocity(const double&);
    // other member functions
    const double &get_energy_rigidity(const double &);
};

#endif //__PARTICLE_H_