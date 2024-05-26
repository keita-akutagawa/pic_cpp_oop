#include <fstream>
#include <iomanip>
#include "pic1D.hpp"


void PIC1D::oneStep()
{
    
}


void PIC1D::saveFields(
    std::string directoryname, 
    std::string filenameWithoutStep, 
    int step
)
{
    std::string filenameB, filenameE, filenameCurrent;

    filenameB = directoryname + "/"
             + filenameWithoutStep + "_" + std::to_string(step)
             + "_B.txt";
    filenameE = directoryname + "/"
             + filenameWithoutStep + "_" + std::to_string(step)
             + "_E.txt";
    filenameCurrent = directoryname + "/"
             + filenameWithoutStep + "_" + std::to_string(step)
             + "_current.txt";


    std::ofstream ofs(filenameB);
    ofs << std::fixed << std::setprecision(6);
    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < nx-1; i++) {
            ofs << B[comp][i] << ",";
        }
        ofs << B[comp][nx-1];
        ofs << std::endl;
    }

    std::ofstream ofs(filenameE);
    ofs << std::fixed << std::setprecision(6);
    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < nx-1; i++) {
            ofs << E[comp][i] << ",";
        }
        ofs << E[comp][nx-1];
        ofs << std::endl;
    }

    std::ofstream ofs(filenameCurrent);
    ofs << std::fixed << std::setprecision(6);
    for (int comp = 0; comp < 3; comp++) {
        for (int i = 0; i < nx-1; i++) {
            ofs << current[comp][i] << ",";
        }
        ofs << current[comp][nx-1];
        ofs << std::endl;
    }
}


void PIC1D::saveFields(
    std::string directoryname, 
    std::string filenameWithoutStep, 
    int step
)
{
    std::string filenameXIon, filenameXElectron;
    std::string filenameVIon, filenameVElectron;

    filenameXIon = directoryname + "/"
             + filenameWithoutStep + "_" + std::to_string(step)
             + "_x_ion.txt";
    filenameXElectron = directoryname + "/"
             + filenameWithoutStep + "_" + std::to_string(step)
             + "_x_electron.txt";
    filenameVIon = directoryname + "/"
             + filenameWithoutStep + "_" + std::to_string(step)
             + "_v_ion.txt";
    filenameVElectron = directoryname + "/"
             + filenameWithoutStep + "_" + std::to_string(step)
             + "_v_electron.txt";


    std::ofstream ofs(filenameXIon);
    ofs << std::fixed << std::setprecision(6);
    for (int i = 0; i < totalNumIon; i++) {
        ofs << particlesIon[i].x << "," 
            << particlesIon[i].y << "," 
            << particlesIon[i].z << std::endl ;
    }

    std::ofstream ofs(filenameXElectron);
    ofs << std::fixed << std::setprecision(6);
    for (int i = 0; i < totalNumElectron; i++) {
        ofs << particlesElectron[i].x << "," 
            << particlesElectron[i].y << "," 
            << particlesElectron[i].z << std::endl ;
    }

    std::ofstream ofs(filenameVIon);
    ofs << std::fixed << std::setprecision(6);
    for (int i = 0; i < totalNumIon; i++) {
        ofs << particlesIon[i].vx << "," 
            << particlesIon[i].vy << "," 
            << particlesIon[i].vz << std::endl ;
    }

    std::ofstream ofs(filenameVElectron);
    ofs << std::fixed << std::setprecision(6);
    for (int i = 0; i < totalNumElectron; i++) {
        ofs << particlesElectron[i].vx << "," 
            << particlesElectron[i].vy << "," 
            << particlesElectron[i].vz << std::endl ;
    }
}

