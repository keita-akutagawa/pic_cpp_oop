#include <cmath>
#include "current.hpp"


void Current::resetCurrent()
{
    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < nx; i++) {
            current[comp][i] = 0.0;
        }
    }
}


void Current::calculateCurrent(
    const std::vector<Particle>& particlesIon, 
    const std::vector<Particle>& particlesEleectron
)
{
    double cx1, cx2, xIndex1, xIndex2;
    double xOverDx;
    double qOverGamma, qVxOverGamma, qVyOverGamma, qVzOverGamma;

    for (int i = 0; i < totalNumIon; i++) {
        xOverDx = particlesIon[i].x / dx;

        xIndex1 = std::floor(xOverDx);
        xIndex2 = xIndex1 + 1;
        xIndex2 = (xIndex2 == nx) ? 0 : xIndex2;

        cx1 = xOverDx - xIndex1;
        cx2 = 1.0 - cx1;

        qOverGamma = qIon / particlesIon[i].gamma;
        qVxOverGamma = qOverGamma * particlesIon[i].vx;
        qVyOverGamma = qOverGamma * particlesIon[i].vy;
        qVzOverGamma = qOverGamma * particlesIon[i].vz;

        current[0][xIndex1] += qVxOverGamma * cx2;
        current[0][xIndex2] += qVxOverGamma * cx1;

        current[1][xIndex1] += qVyOverGamma * cx2;
        current[1][xIndex2] += qVyOverGamma * cx1;

        current[2][xIndex1] += qVzOverGamma * cx2;
        current[2][xIndex2] += qVzOverGamma * cx1;
    }
}


std::vector<std::vector<double>> Current::getCurrent()
{
    return current;
}

