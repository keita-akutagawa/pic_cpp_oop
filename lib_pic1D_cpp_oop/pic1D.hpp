#include <vector>
#include "const.hpp"
#include "particle_push.hpp"


class PIC1D
{
private:
    std::vector<std::vector<double>> E;
    std::vector<std::vector<double>> B;

public:
    PIC1D() :
        E(8, std::vector<double>(nx, 0.0)), 
        B(8, std::vector<double>(nx, 0.0))
        {}
    
    void initializeParticles();

    void initializeFields();
    
    void oneStep();

    void save();

    void getParticles();

    void getFields();

private:

}


