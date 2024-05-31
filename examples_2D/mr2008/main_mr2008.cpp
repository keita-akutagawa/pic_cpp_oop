#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../../lib_pic2D_cpp_oop/pic2D_symX_conY.hpp"


const double c = 1.0;
const double epsilon0 = 1.0;
const double mu0 = 1.0;

const int numberDensityIon = 10;
const int numberDensityElectron = 10;

const double B0 = sqrt(static_cast<double>(numberDensityElectron)) / 1.0;

const double mRatio = 9.0;
const double mElectron = 1.0;
const double mIon = mRatio * mElectron;


const double tRatio = 1.0;
const double tElectron = (B0 * B0 / 2.0 / mu0) / (numberDensityIon + numberDensityElectron * tRatio);
const double tIon = tRatio * tElectron;

const double qRatio = -1.0;
const double qElectron = -1.0 * sqrt(epsilon0 * tElectron / static_cast<double>(numberDensityElectron));
const double qIon = qRatio * qElectron;

const double omegaPe = sqrt(static_cast<double>(numberDensityElectron) * pow(qElectron, 2) / mElectron / epsilon0);
const double omegaPi = sqrt(static_cast<double>(numberDensityIon) * pow(qIon, 2) / mIon / epsilon0);
const double omegaCe = abs(qElectron * B0 / mElectron);
const double omegaCi = qIon * B0 / mIon;

const double debyeLength = sqrt(epsilon0 * tElectron / static_cast<double>(numberDensityElectron) / pow(qElectron, 2));
//追加
const double ionInertialLength = c / omegaPi;

const int nx = int(100 * ionInertialLength);
const double dx = 1.0;
const double xmin = 0.0; 
const double xmax = nx * dx;

const int ny = int(30 * ionInertialLength);
const double dy = 1.0;
const double ymin = 0.0; 
const double ymax = ny * dy;

const double dt = 0.5;

const int totalNumIon = nx * ny * numberDensityIon;
const int totalNumElectron = nx * ny * numberDensityElectron;
const int totalNumParticles = totalNumIon + totalNumElectron;

const double vThIon = sqrt(tIon / mIon);
const double vThElectron = sqrt(tElectron / mElectron);
const double bulkVxIon = 0.0;
const double bulkVyIon = 0.0;
const double bulkVzIon = 0.0;
const double bulkVxElectron = 0.0;
const double bulkVyElectron = 0.0;
const double bulkVzElectron = 0.0;

const int totalStep = 100;
double totalTime = 0.0;


void PIC2DSymXConY::initialize()
{
    initializeParticle.uniformForPositionX(
        0, totalNumIon, 0, particlesIon
    );
    initializeParticle.uniformForPositionX(
        0, totalNumElectron, 100, particlesElectron
    );
    initializeParticle.uniformForPositionY(
        0, totalNumIon, 200, particlesIon
    );
    initializeParticle.uniformForPositionY(
        0, totalNumElectron, 300, particlesElectron
    );
    for (int i = 0; i < totalNumIon; i++) {
        particlesIon[i].z = 0.0;
    }
    for (int i = 0; i < totalNumElectron; i++) {
        particlesElectron[i].z = 0.0;
    }

    initializeParticle.maxwellDistributionForVelocity(
        bulkVxIon, bulkVyIon, bulkVzIon, vThIon, vThIon, vThIon, 
        0, totalNumIon, 400, particlesIon
    );
    initializeParticle.maxwellDistributionForVelocity(
        bulkVxElectron, bulkVyElectron, bulkVzElectron, vThElectron, vThElectron, vThElectron, 
        0, totalNumElectron, 500, particlesElectron
    );

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            B[0][i][j] = 0.0;
            B[1][i][j] = 0.0;
            B[2][i][j] = 0.0;
            E[0][i][j] = 0.0;
            E[1][i][j] = 0.0;
            E[2][i][j] = 0.0;
            current[0][i][j] = 0.0;
            current[1][i][j] = 0.0;
            current[2][i][j] = 0.0;
        }
    }
}


int main()
{
    std::string directoryname = "results";
    std::string filenameWithoutStep = "mr2008_forcefree";
    std::ofstream logfile("log.txt");
    int recordStep = 100;

    std::cout << "total number of partices is " << totalNumParticles << std::endl;
    std::cout << "box size is " << nx << " X " << ny << std::endl;
    std::cout << std::setprecision(4) 
              << "omega_pe * t = " << totalStep * dt * omegaPe << std::endl;

    PIC2DSymXConY pIC2DSymXConY;

    pIC2DSymXConY.initialize();

    for (int step = 0; step < totalStep+1; step++) {
        if (step % recordStep == 0) {
            std::cout << std::to_string(step) << " step done : total time is "
                      << std::setprecision(4) << totalTime
                      << std::endl;
            logfile << std::setprecision(6) << totalTime << std::endl;
            pIC2DSymXConY.saveFields(
                directoryname, filenameWithoutStep, step
            );
            pIC2DSymXConY.saveEnergy(
                directoryname, filenameWithoutStep, step
            );
            pIC2DSymXConY.saveMoments(
                directoryname, filenameWithoutStep, step
            );
            //pIC2DSymXConY.saveParticle(
            //    directoryname, filenameWithoutStep, step
            //);
        }

        pIC2DSymXConY.oneStep();

        totalTime += dt;
    }

    return 0;
}



