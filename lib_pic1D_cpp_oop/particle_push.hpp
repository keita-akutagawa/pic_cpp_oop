#include <vector>
#include "particle_struct.hpp"
#include "const.hpp"


class ParticlePush
{
private: 

public:

    void pushVelocity(
        const std::vector<std::vector<double>>& B, 
        const std::vector<std::vector<double>>& E, 
        std::vector<Particle>& particlesIon, 
        std::vector<Particle>& particlesElectron
    );
    void pushPosition(
        std::vector<Particle>& particlesIon, 
        std::vector<Particle>& particlesElectron
    );

private:

    void pushVelocityOfOneSpecies(
        const std::vector<std::vector<double>>& B, 
        const std::vector<std::vector<double>>& E, 
        std::vector<Particle>& particlesSpecies, 
        double q, double m, double totalNumSpecies
    );

    ParticleField getParticleFields(
        const std::vector<std::vector<double>>& B, 
        const std::vector<std::vector<double>>& E, 
        const Particle& particle
    );

    void pushPositionOfOneSpecies(
        std::vector<Particle>& particlesSpecies, 
        double totalNumSpecies
    );

};


