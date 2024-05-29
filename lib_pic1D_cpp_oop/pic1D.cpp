#include <fstream>
#include <iomanip>
#include "pic1D.hpp"


void PIC1D::oneStep()
{
    fieldSolver.timeEvolutionB(B, E, dt/2.0);

    for (int i = 0; i < nx; i++) {
        tmpB[0][i] = B[0][i];
        tmpB[1][i] = 0.5 * (B[1][i] + B[1][(i-1+nx)%nx]);
        tmpB[2][i] = 0.5 * (B[2][i] + B[2][(i-1+nx)%nx]);
        tmpE[0][i] = 0.5 * (E[0][i] + E[0][(i-1+nx)%nx]);
        tmpE[1][i] = E[1][i];
        tmpE[2][i] = E[2][i];
    }

    particlePush.pushVelocity(
        particlesIon, particlesElectron, tmpB, tmpE, dt
    );

    particlePush.pushPosition(
        particlesIon, particlesElectron, dt/2.0
    );
    boundary.periodicBoundaryParticleX(
        particlesIon, particlesElectron
    );

    currentCalculater.resetCurrent(tmpCurrent);
    currentCalculater.calculateCurrent(
        tmpCurrent, particlesIon, particlesElectron
    );
    for (int i = 0; i < nx; i++) {
        current[0][i] = 0.5 * (tmpCurrent[0][i] + tmpCurrent[0][(i+1)%nx]);
        current[1][i] = tmpCurrent[1][i];
        current[2][i] = tmpCurrent[2][i];
    }

    fieldSolver.timeEvolutionB(B, E, dt/2.0);

    fieldSolver.timeEvolutionE(E, B, current, dt);

    particlePush.pushPosition(
        particlesIon, particlesElectron, dt/2.0
    );
    boundary.periodicBoundaryParticleX(
        particlesIon, particlesElectron
    );
}


void PIC1D::saveFields(
    std::string directoryname, 
    std::string filenameWithoutStep, 
    int step
)
{
    std::string filenameB, filenameE, filenameCurrent;

    filenameB = directoryname + "/"
             + filenameWithoutStep + "_B_" + std::to_string(step)
             + ".txt";
    filenameE = directoryname + "/"
             + filenameWithoutStep + "_E_" + std::to_string(step)
             + ".txt";
    filenameCurrent = directoryname + "/"
             + filenameWithoutStep + "_current_" + std::to_string(step)
             + ".txt";


    std::ofstream ofsB(filenameB);
    ofsB << std::setprecision(6);
    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < nx-1; i++) {
            ofsB << B[comp][i] << ",";
        }
        ofsB << B[comp][nx-1];
        ofsB << std::endl;
    }

    std::ofstream ofsE(filenameE);
    ofsE << std::setprecision(6);
    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < nx-1; i++) {
            ofsE << E[comp][i] << ",";
        }
        ofsE << E[comp][nx-1];
        ofsE << std::endl;
    }

    std::ofstream ofsCurrent(filenameCurrent);
    ofsCurrent << std::setprecision(6);
    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < nx-1; i++) {
            ofsCurrent << current[comp][i] << ",";
        }
        ofsCurrent << current[comp][nx-1];
        ofsCurrent << std::endl;
    }
}


void PIC1D::saveParticle(
    std::string directoryname, 
    std::string filenameWithoutStep, 
    int step
)
{
    std::string filenameXIon, filenameXElectron;
    std::string filenameVIon, filenameVElectron;

    filenameXIon = directoryname + "/"
             + filenameWithoutStep + "_x_ion_" + std::to_string(step)
             + ".txt";
    filenameXElectron = directoryname + "/"
             + filenameWithoutStep + "_x_electron_" + std::to_string(step)
             + ".txt";
    filenameVIon = directoryname + "/"
             + filenameWithoutStep + "_v_ion_" + std::to_string(step)
             + ".txt";
    filenameVElectron = directoryname + "/"
             + filenameWithoutStep + "_v_electron_" + std::to_string(step)
             + ".txt";


    std::ofstream ofsXIon(filenameXIon);
    ofsXIon << std::setprecision(6);
    for (int i = 0; i < totalNumIon; i++) {
        ofsXIon << particlesIon[i].x << "," 
                << particlesIon[i].y << "," 
                << particlesIon[i].z << std::endl ;
    }

    std::ofstream ofsXElectron(filenameXElectron);
    ofsXElectron << std::setprecision(6);
    for (int i = 0; i < totalNumElectron; i++) {
        ofsXElectron << particlesElectron[i].x << "," 
                     << particlesElectron[i].y << "," 
                     << particlesElectron[i].z << std::endl ;
    }

    std::ofstream ofsVIon(filenameVIon);
    ofsVIon << std::setprecision(6);
    for (int i = 0; i < totalNumIon; i++) {
        ofsVIon << particlesIon[i].vx << "," 
                << particlesIon[i].vy << "," 
                << particlesIon[i].vz << std::endl ;
    }

    std::ofstream ofsVElectron(filenameVElectron);
    ofsVElectron << std::setprecision(6);
    for (int i = 0; i < totalNumElectron; i++) {
        ofsVElectron << particlesElectron[i].vx << "," 
                     << particlesElectron[i].vy << "," 
                     << particlesElectron[i].vz << std::endl ;
    }
}

