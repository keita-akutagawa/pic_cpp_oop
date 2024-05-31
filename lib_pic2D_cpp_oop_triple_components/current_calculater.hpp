#include <vector>
#include "const.hpp"
#include "particle_struct.hpp"


class CurrentCalculater
{
private: 

public: 
    void resetCurrent(
        std::vector<std::vector<std::vector<double>>>& current
    );

    void calculateCurrentForWallBoundary(
        std::vector<std::vector<std::vector<double>>>& current, 
        const std::vector<Particle>& particlesIon, 
        const std::vector<Particle>& particlesEleectron
    );

private:

    void calculateCurrentOfOneSpeciesForWallBoundary(
        std::vector<std::vector<std::vector<double>>>& current, 
        const std::vector<Particle>& particlesSpecies, 
        double q, int nStart, int nEnd
    );
};

