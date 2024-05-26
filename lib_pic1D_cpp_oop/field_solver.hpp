#include <vector>
#include "const.hpp"

class FieldSolver
{
private:
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> E;

public:
    void timeEvolutionB();

    void timeEvolutionE(
        const std::vector<std::vector<double>>& current
    );

    std::vector<std::vector<double>> getB();
    std::vector<std::vector<double>> getE();

private:

};


