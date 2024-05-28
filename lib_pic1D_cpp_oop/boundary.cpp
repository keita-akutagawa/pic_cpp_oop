#include "boundary.hpp"


void Boundary::periodicBoundaryParticleX(
    std::vector<Particle>& particlesIon,
    std::vector<Particle>& particlesElectron
)
{
    for (int i = 0; i < totalNumIon; i++) {
        if (particlesIon[i].x < xmin) {
            particlesIon[i].x += xmax;
        }
        if (particlesIon[i].x > xmax) {
            particlesIon[i].x -= xmax;
        }
    }
}


void Boundary::periodicBoundaryBX()
{
    return;
}


void Boundary::periodicBoundaryEX()
{
    return;
}


void Boundary::periodicBoundaryCurrentX()
{
    return;
}


