#include <string>
#include <fstream>
#include <cmath>
#include "../../lib_pic1D_cpp_oop/pic1D.hpp"


const double c = 1.0;
const double epsilon0 = 1.0;
const double mu0 = 1.0;

const int nx = 512;
const double dx = 1.0;
extern const double xmin = 0.0; 
extern const double xmax = nx * dx;

const double dt = 0.5;

const int numberDensityIon = 10;
const int numberDensityElectron = 10;

const int totalNumIon = nx * numberDensityIon;
//追加
const int totalNumElectronBeam1 = nx * numberDensityElectron / 2;
const int totalNumElectronBeam2 = nx * numberDensityElectron / 2;
const int totalNumElectron = totalNumElectronBeam1 + totalNumElectronBeam2;
const int totalNumParticles = totalNumIon + totalNumElectron;

const double B0 = sqrt(static_cast<double>(numberDensityElectron)) / 10.0;

const double mRatio = 9.0;
const double mElectron = 1.0;
const double mIon = mRatio * mElectron;

const double tRatio = 1.0;
const double tElectron = 0.5 * mElectron * pow(0.01 * c, 2);
const double tIon = tRatio * tElectron;

const double qRatio = -1.0;
const double qElectron = -1.0 * sqrt(epsilon0 * tElectron / static_cast<double>(numberDensityElectron));
const double qIon = qRatio * qElectron;

const double omegaPe = sqrt(static_cast<double>(numberDensityElectron) * pow(qElectron, 2) / mElectron / epsilon0);
const double omegaPi = sqrt(static_cast<double>(numberDensityIon) * pow(qIon, 2) / mIon / epsilon0);
const double omegaCe = abs(qElectron * B0 / mElectron);
const double omegaCi = qIon * B0 / mIon;

const double debyeLength = sqrt(epsilon0 * tElectron / static_cast<double>(numberDensityElectron) / pow(qElectron, 2));

const double vThIon = sqrt(2.0 * tIon / mIon);
const double vThElectron = sqrt(2.0 * tElectron / mElectron);
const double bulkVxIon = 0.0;
const double bulkVyIon = 0.0;
const double bulkVzIon = 0.0;
const double bulkVxElectron = -10.0 * vThIon;
const double bulkVyElectron = 0.0;
const double bulkVzElectron = 0.0;
//追加
const double bulkVxElectronBeam = 10.0 * vThIon;
const double bulkVyElectronBeam = 0.0;
const double bulkVzElectronBeam = 0.0;

const int totalStep = 10000;
const double totalTime = 0.0;


void PIC1D::initialize()
{
    initializeParticle.uniformForPositionX(
        0, totalNumIon, 0, particlesIon
    );
    initializeParticle.uniformForPositionX(
        0, totalNumElectron, 100, particlesElectron
    );
    initializeParticle.maxwellDistributionForVelocity(
        bulkVxIon, bulkVyIon, bulkVzIon, vThIon, 
        0, totalNumIon, 200, particlesIon
    );
    initializeParticle.maxwellDistributionForVelocity(
        bulkVxElectron, bulkVyElectron, bulkVzElectron, vThElectron, 
        0, totalNumElectron / 2, 300, particlesElectron
    );
    initializeParticle.maxwellDistributionForVelocity(
        bulkVxElectronBeam, bulkVyElectronBeam, bulkVzElectronBeam, vThElectron, 
        totalNumElectron / 2, totalNumElectron, 400, particlesElectron
    );

    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < nx; i++) {
            B[comp][i] = 0.0;
            E[comp][i] = 0.0;
            current[comp][i] = 0.0;
        }
    }
}


int main()
{
    std::string directoryname = "results";
    std::string filenameWithoutStep = "two_stream_electron";
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



