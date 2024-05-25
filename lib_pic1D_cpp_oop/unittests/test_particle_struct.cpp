#include "../particle_push.hpp"
#include "../const.hpp"
#include <gtest/gtest.h>
#include <iostream>


TEST(ParticlePush, setParticle)
{
    const int totalNumParticles = 10;
    std::vector<Particle> particles(totalNumParticles);

    EXPECT_EQ(particles.size(), totalNumParticles);

    for (int i = 0; i < totalNumParticles; i++) {
        EXPECT_EQ(particles[i].x, 0.0);
    }
}


