#include "boundary.hpp"


void Boundary::periodicBoundaryParticleX(
    std::vector<Particle>& particlesIon,
    std::vector<Particle>& particlesElectron
)
{
    for (int i = 0; i < totalNumIon; i++) {
        if (particlesIon[i].x < xmin) {
            particlesIon[i].x += xmax - xmin;
        }
        if (particlesIon[i].x > xmax) {
            particlesIon[i].x -= xmax - xmin;
        }
    }

    for (int i = 0; i < totalNumElectron; i++) {
        if (particlesElectron[i].x < xmin) {
            particlesElectron[i].x += xmax - xmin;
        }
        if (particlesElectron[i].x > xmax) {
            particlesElectron[i].x -= xmax - xmin;
        }
    }
}


void Boundary::periodicBoundaryParticleY(
    std::vector<Particle>& particlesIon,
    std::vector<Particle>& particlesElectron
)
{
    for (int i = 0; i < totalNumIon; i++) {
        if (particlesIon[i].y < ymin) {
            particlesIon[i].y += ymax - ymin;
        }
        if (particlesIon[i].y > ymax) {
            particlesIon[i].y -= ymax - ymin;
        }
    }

    for (int i = 0; i < totalNumElectron; i++) {
        if (particlesElectron[i].y < ymin) {
            particlesElectron[i].y += ymax - ymin;
        }
        if (particlesElectron[i].y > ymax) {
            particlesElectron[i].y -= ymax - ymin;
        }
    }
}


void Boundary::periodicBoundaryBX()
{
    return;
}


void Boundary::periodicBoundaryBY()
{
    return;
}


void Boundary::periodicBoundaryEX()
{
    return;
}


void Boundary::periodicBoundaryEY()
{
    return;
}


void Boundary::periodicBoundaryCurrentX()
{
    return;
}


void Boundary::periodicBoundaryCurrentY()
{
    return;
}


