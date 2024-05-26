#include <vector>
#include <string>
#include "const.hpp"
#include "particle_push.hpp"
#include "field_solver.hpp"
#include "current_calculater.hpp"


class PIC1D
{
private:
    std::vector<Particle> particlesIon;
    std::vector<Particle> particlesElectron;
    std::vector<std::vector<double>> E;
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> current;

    ParticlePush particlePush;
    FieldSolver fieldSolver;
    CurrentCalculater currentCalculater;

public:
    PIC1D() :
        particlesIon(totalNumIon), 
        particlesElectron(totalNumElectron), 
        E(3, std::vector<double>(nx, 0.0)), 
        B(3, std::vector<double>(nx, 0.0)), 
        current(3, std::vector<double>(nx, 0.0))
        {}
    
    virtual void initialize();
    
    void oneStep();

    void saveFields(
        std::string directoryname, 
        std::string filenameWithoutStep, 
        int step
    );

    void saveParticle(
        std::string directoryname, 
        std::string filenameWithoutStep, 
        int step
    );

    void getParticles();

    void getFields();

private:

};


