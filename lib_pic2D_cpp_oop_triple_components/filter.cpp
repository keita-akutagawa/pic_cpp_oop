#include "filter.hpp"


void Filter::resetRho()
{
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            rho[i][j] = 0.0;
        }
    }
}


void Filter::langdonMarderCorrection(
    std::vector<std::vector<double>>& F, 
    std::vector<std::vector<std::vector<double>>>& E
)
{
    for (int i = 1; i < nx; i++) {
        for (int j = 1; j < ny; j++) {
            F[i][j] = ((E[0][i][j] - E[0][i - 1][j])/dx + (E[1][i][j] - E[1][i][j - 1])/dy)
                    - rho[i][j] / epsilon0;
        }
    }
    for (int j = 0; j < ny; j++) {
        F[0][j] = ((E[0][0][j] - 0.0)/dx + (E[1][0][j] - 0.0)/dy)
                - rho[0][j] / epsilon0;
    }
    for (int i = 0; i < nx; i++) {
        F[i][0] = 0.0;
        F[i][ny - 1] = 0.0;
    }
}

