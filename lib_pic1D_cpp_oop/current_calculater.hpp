#include <vector>
#include "const.hpp"
#include "particle_struct.hpp"


class CurrentCalculater
{
private: 

public: 
    void resetCurrent(
        std::vector<std::vector<double>>& current
    );

    void calculateCurrent(
        const std::vector<Particle>& particlesIon, 
        const std::vector<Particle>& particlesEleectron, 
        std::vector<std::vector<double>>& current
    );

private:

    void calculateCurrentOfOneSpecies(
        const std::vector<Particle>& particlesSpecies, 
        std::vector<std::vector<double>>& current, 
        double q, double totalNumSpecies
    );
};

