#include <vector>
#include "particle_struct.hpp"
#include "const.hpp"


class ParticlePush
{
private: 

public:

    void pushPosition(
        std::vector<Particle>& particlesIon, 
        std::vector<Particle>& particlesElectron, 
        double dt
    );

    void pushVelocityForWallBoundary(
        std::vector<Particle>& particlesIon, 
        std::vector<Particle>& particlesElectron, 
        const std::vector<std::vector<std::vector<double>>>& B, 
        const std::vector<std::vector<std::vector<double>>>& E, 
        double dt
    );

private:

    void pushPositionOfOneSpecies(
        std::vector<Particle>& particlesSpecies, 
        int nStart, int nEnd, 
        double dt
    );

    void pushVelocityOfOneSpeciesForWallBoundary(
        std::vector<Particle>& particlesSpecies, 
        const std::vector<std::vector<std::vector<double>>>& B, 
        const std::vector<std::vector<std::vector<double>>>& E, 
        double q, double m, 
        int nStart, int nEnd, double dt
    );

    ParticleField getParticleFieldsForWallBoundary(
        const std::vector<std::vector<std::vector<double>>>& B, 
        const std::vector<std::vector<std::vector<double>>>& E, 
        const Particle& particle
    );
};


