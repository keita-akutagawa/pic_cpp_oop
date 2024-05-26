#include <vector>
#include "const.hpp"
#include "particle_push.hpp"


class PIC1D
{
private:
    std::vector<std::vector<double>> E;
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> current;

public:
    PIC1D() :
        E(3, std::vector<double>(nx, 0.0)), 
        B(3, std::vector<double>(nx, 0.0)), 
        current(3, std::vector<double>(nx, 0.0))
        {}
    
    void initializeParticles();

    void initializeFields();
    
    void oneStep();

    void save();

    void getParticles();

    void getFields();

private:

};


