#include <vector>
#include "const.hpp"
#include "particle_struct.hpp"


class Boundary
{
private:

public:

    void periodicBoundaryParticleX(
        std::vector<Particle>& particlesIon,
        std::vector<Particle>& particlesElectron
    );
    void conductingWallBoundaryParticleX(
        std::vector<Particle>& particlesIon,
        std::vector<Particle>& particlesElectron
    );
    void periodicBoundaryParticleY(
        std::vector<Particle>& particlesIon,
        std::vector<Particle>& particlesElectron
    );
    void conductingWallBoundaryParticleY(
        std::vector<Particle>& particlesIon,
        std::vector<Particle>& particlesElectron
    );

    void periodicBoundaryBX();
    void conductingWallBoundaryBX();
    void periodicBoundaryBY();
    void conductingWallBoundaryBY();

    void periodicBoundaryEX();
    void conductingWallBoundaryEX();
    void periodicBoundaryEY();
    void conductingWallBoundaryEY();

    void periodicBoundaryCurrentX();
    void conductingWallBoundaryCurrentX();
    void periodicBoundaryCurrentY();
    void conductingWallBoundaryCurrentY();

private:

};


