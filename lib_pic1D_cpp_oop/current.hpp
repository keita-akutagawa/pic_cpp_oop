#include <vector>
#include "const.hpp"
#include "particle_struct.hpp"


class Current
{
private: 
    std::vector<std::vector<double>> current;

public: 
    void resetCurrent();

    void calculateCurrent(
        const std::vector<Particle>& particlesIon, 
        const std::vector<Particle>& particlesEleectron
    );

    std::vector<std::vector<double>> getCurrent();

private:

};

