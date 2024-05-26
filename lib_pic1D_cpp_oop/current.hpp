#include <vector>
#include "const.hpp"
#include "particle_struct.hpp"


class Current
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

};

