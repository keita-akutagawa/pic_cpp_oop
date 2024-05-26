#include <string>
#include <fstream>
#include "../lib_pic1D_cpp_oop/initialize_particle.hpp"
#include "../lib_pic1D_cpp_oop/pic1D.hpp"


const double c;
const double epsilon0;
const double mu0;

const int nx;
const double dx;
const double dt;

const int numberDensityIon;
const int numberDensityElectron;

const int totalNumIon;
const int totalNumElectron;
const int totalNumParticles;

const double mRatio;
const double mIon;
const double mElectron;
const double qRatio;
const double qIon;
const double qElectron;

const double tRatio;
const double tIon;
const double tElectron;

const double omegaPe;
const double omegaPi;
const double omegaCe;
const double omegaCi;

const double debyeLength;

const double bulkVxIon;
const double bulkVyIon;
const double bulkVzIon;
const double bulkVxElectron;
const double bulkVyElectron;
const double bulkVzElectron;
const double vThIon;
const double vThElectron;

const int totalStep;
const double totalTime;


void PIC1D::initialize()
{
    
}


int main()
{
    std::string directoryname = "results_CT";
    std::string filenameWithoutStep = "orszag_tang";
    std::ofstream logfile("log.txt");
    int recordStep = 100;

    PIC1D pIC1D;

    pIC1D.initialize();

    for (int step = 0; step < totalStep+1; step++) {
        if (step % recordStep) {
            pIC1D.saveFields(
                directoryname, filenameWithoutStep, step
            );
            pIC1D.saveParticle(
                directoryname, filenameWithoutStep, step
            );
        }

        pIC1D.oneStep();
    }

    return 0;
}



