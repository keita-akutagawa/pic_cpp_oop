#include <vector>
#include "const.hpp"

class FieldSolver
{
private:

public:
    void timeEvolutionB(
        const std::vector<std::vector<double>>& E, 
        std::vector<std::vector<double>>& B, 
        double dt
    );

    void timeEvolutionE(
        const std::vector<std::vector<double>>& B, 
        const std::vector<std::vector<double>>& current, 
        std::vector<std::vector<double>>& E, 
        double dt
    );

private:

};


