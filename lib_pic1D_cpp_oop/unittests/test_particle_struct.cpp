#include "../particle_push.hpp"
#include "../const.hpp"
#include <gtest/gtest.h>
#include <iostream>


TEST(ParticlePush, setParticle)
{
    const int numParticles = 10;
    std::vector<Particle> particles(numParticles);

    EXPECT_EQ(particles.size(), numParticles);

    for (int i = 0; i < numParticles; i++) {
        EXPECT_EQ(particles[i].x, 0.0);
    }
}


