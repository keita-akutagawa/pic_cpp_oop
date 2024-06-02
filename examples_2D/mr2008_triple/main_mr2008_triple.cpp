#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../../lib_pic2D_cpp_oop_triple_components/pic2D_symX_conY.hpp"


std::string directoryname = "results_mr=200-9_n=2";
std::string filenameWithoutStep = "mr2008";
std::ofstream logfile("log.txt");

const int totalStep = 30000;
const int recordStep = 100;
double totalTime = 0.0;


const double c = 1.0;
const double epsilon0 = 1.0;
const double mu0 = 1.0;
const double dOfLangdonMarderCorrection = 0.01;

const int numberDensityIon = 10;
const int numberDensityElectron = 10;
const int numberDensityHeavyIon = 2;

const double B0 = sqrt(static_cast<double>(numberDensityElectron)) / 1.0;

const double mRatio = 9.0;
const double mElectron = 1.0;
const double mIon = mRatio * mElectron;
const double mHeavyIon = mIon * 200;


const double tRatio = 1.0;
const double tElectron = (B0 * B0 / 2.0 / mu0) / (numberDensityIon + numberDensityElectron * tRatio);
const double tIon = tRatio * tElectron;
const double tHeavyIon = tIon;

const double qRatio = -1.0;
const double qElectron = -1.0 * sqrt(epsilon0 * tElectron / static_cast<double>(numberDensityElectron));
const double qIon = qRatio * qElectron;
const double qHeavyIon = qIon;

const double omegaPe = sqrt(static_cast<double>(numberDensityElectron) * pow(qElectron, 2) / mElectron / epsilon0);
const double omegaPi = sqrt(static_cast<double>(numberDensityIon) * pow(qIon, 2) / mIon / epsilon0);
const double omegaCe = abs(qElectron * B0 / mElectron);
const double omegaCi = qIon * B0 / mIon;

const double debyeLength = sqrt(epsilon0 * tElectron / static_cast<double>(numberDensityElectron) / pow(qElectron, 2));
//追加
const double ionInertialLength = c / omegaPi;

const int nx = int(200 * ionInertialLength);
const double dx = 1.0;
const double xmin = 0.5 * dx; 
const double xmax = nx * dx - 1.0 * dx;

const int ny = int(50 * ionInertialLength);
const double dy = 1.0;
const double ymin = 0.5 * dy; 
const double ymax = ny * dy - 1.0 * dy;

const double dt = 0.5;

//追加
const double sheatThickness = 1.5 * ionInertialLength;
const double reconnectionTriggerRatio = 0.1;
const double xPointPosition = 20.0 * ionInertialLength;

//追加
const int harrisNumIon = int(nx * numberDensityIon * 2.0 * sheatThickness);
const int backgroundNumIon = int(0.2 * nx * ny * numberDensityIon);
const int backgroundNumHeavyIon = nx * ny * numberDensityHeavyIon;
const int totalNumIon = harrisNumIon + backgroundNumIon + backgroundNumHeavyIon;
const int harrisNumElectron = int(nx * numberDensityElectron * 2.0 * sheatThickness);
const int backgroundNumElectron = int(0.2 * nx * ny * numberDensityElectron)
                                + backgroundNumHeavyIon * std::abs(qHeavyIon / qElectron);
const int totalNumElectron = harrisNumElectron + backgroundNumElectron;
const int totalNumParticles = totalNumIon + totalNumElectron;

const double vThIon = sqrt(2.0 * tIon / mIon);
const double vThElectron = sqrt(2.0 * tElectron / mElectron);
const double bulkVxElectron = 0.0;
const double bulkVyElectron = 0.0;
const double bulkVzElectron = c * debyeLength / sheatThickness * sqrt(2 / (1.0 + 1.0/tRatio));
const double bulkVxIon = -bulkVxElectron / tRatio;
const double bulkVyIon = -bulkVyElectron / tRatio;
const double bulkVzIon = -bulkVzElectron / tRatio;

const double vThIonB = sqrt(2.0 * tIon / 10.0 / mIon);
const double vThElectronB = sqrt(2.0 * tElectron / 10.0 / mElectron);
const double vThHeavyIonB = sqrt(2.0 * tHeavyIon / 10.0 / mHeavyIon);
const double bulkVxElectronB = 0.0;
const double bulkVyElectronB = 0.0;
const double bulkVzElectronB = 0.0;
const double bulkVxIonB = 0.0;
const double bulkVyIonB = 0.0;
const double bulkVzIonB = 0.0;
const double bulkVxHeavyIonB = 0.0;
const double bulkVyHeavyIonB = 0.0;
const double bulkVzHeavyIonB = 0.0;


void PIC2DSymXConY::initialize()
{
    initializeParticle.uniformForPositionX(
        0, totalNumIon, 0, particlesIon
    );
    initializeParticle.uniformForPositionX(
        0, totalNumElectron, 100, particlesElectron
    );

    initializeParticle.harrisForPositionY(
        0, harrisNumIon, 200, particlesIon, sheatThickness
    );
    initializeParticle.uniformForPositionY(
        harrisNumIon, harrisNumIon + backgroundNumIon, 300, particlesIon
    );
    initializeParticle.uniformForPositionY(
        harrisNumIon + backgroundNumIon, totalNumIon, 350, particlesIon
    );
    initializeParticle.harrisForPositionY(
        0, harrisNumElectron, 400, particlesElectron, sheatThickness
    );
    initializeParticle.uniformForPositionY(
        harrisNumElectron, totalNumElectron, 500, particlesElectron
    );

    for (int i = 0; i < totalNumIon; i++) {
        particlesIon[i].z = 0.0;
    }
    for (int i = 0; i < totalNumElectron; i++) {
        particlesElectron[i].z = 0.0;
    }

    initializeParticle.maxwellDistributionForVelocity(
        bulkVxIon, bulkVyIon, bulkVzIon, vThIon, vThIon, vThIon, 
        0, harrisNumIon, 600, particlesIon
    );
    initializeParticle.maxwellDistributionForVelocity(
        bulkVxIonB, bulkVyIonB, bulkVzIonB, vThIonB, vThIonB, vThIonB, 
        harrisNumIon,harrisNumIon + backgroundNumIon, 700, particlesIon
    );
    initializeParticle.maxwellDistributionForVelocity(
        bulkVxHeavyIonB, bulkVyHeavyIonB, bulkVzHeavyIonB, vThHeavyIonB, vThHeavyIonB, vThHeavyIonB, 
        harrisNumIon + backgroundNumIon, totalNumIon, 750, particlesIon
    );
    initializeParticle.maxwellDistributionForVelocity(
        bulkVxElectron, bulkVyElectron, bulkVzElectron, vThElectron, vThElectron, vThElectron, 
        0, harrisNumElectron, 800, particlesElectron
    );
    initializeParticle.maxwellDistributionForVelocity(
        bulkVxElectronB, bulkVyElectronB, bulkVzElectronB, vThElectronB, vThElectronB, vThElectronB, 
        harrisNumElectron, totalNumElectron, 900, particlesElectron
    );

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double yCenter = 0.5 * (ymax - ymin);
            B[0][i][j] = B0 * tanh((j * dy - yCenter) / sheatThickness);
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

    // reconnection trigger
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double yCenter = 0.5 * (ymax - ymin);
            B[0][i][j] += -B0 * reconnectionTriggerRatio * (j * dy - yCenter) / sheatThickness
                        * exp(-(pow((i * dx - xPointPosition), 2) + pow((j * dy - yCenter), 2))
                        / pow(2.0 * sheatThickness, 2));
            B[1][i][j] += B0 * reconnectionTriggerRatio * (i * dx - xPointPosition) / sheatThickness
                        * exp(-(pow((i * dx - xPointPosition), 2) + pow((j * dy - yCenter), 2))
                        / pow(2.0 * sheatThickness, 2));  
            B[2][i][j] = 0.0;
        }
    }
}


//-------------------------------------------------------------


int main()
{
    std::cout << "total number of partices is " << totalNumParticles << std::endl;
    std::cout << "ion : " << harrisNumIon + backgroundNumIon << std::endl;
    std::cout << "heavy ion : " << backgroundNumHeavyIon << std::endl;
    std::cout << "electron : " << totalNumElectron << std::endl;
    std::cout << "box size is " << nx << " X " << ny << std::endl;
    std::cout << "sheat thickness is " 
              << std::setprecision(4) << sheatThickness << "[debye length]" << std::endl;
    std::cout << std::setprecision(4) 
              << "omega_ci * t = " << totalStep * dt * omegaCi << std::endl;

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



