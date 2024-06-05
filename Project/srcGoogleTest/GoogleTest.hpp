#ifndef __TESTPOLYGONS_H
#define __TESTPOLYGONS_H

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include "Paraview.hpp"

using namespace Eigen;
using namespace std;
using namespace DFNLibrary;

// Test per la funzione letturaDatiFileFR
TEST(DFNTest, LetturaDatiFileFR){
    DFN Fract;
    unsigned int numFract = 0;
    string percorsoFileFR = "DFN/FR3_data.txt";
    bool result = letturaDatiFileFR(percorsoFileFR, Fract, numFract);
    ASSERT_TRUE(result);
}

// Test per la funzione calcoloBaricentriEDistBaricentroVertici
TEST(DFNTest, CalcoloBaricentriEDistBaricentroVertici){
    DFN Fract;
    unsigned int numFract = 0;
    unsigned int numIntersezioniFratture = 0;
    bool result = calcoloBaricentriEDistBaricentroVertici(Fract, numFract, numIntersezioniFratture);
    ASSERT_TRUE(result);
}

// Test per la funzione calcoloEqPianoEdEqRetteLati
TEST(DFNTest, calcoloEqPianoEdEqRetteLati){
    DFN Fract;
    bool result = calcoloEqPianoEdEqRetteLati(Fract);
    ASSERT_TRUE(result);
}

// Test per la funzione calcoloIntersezionePiani
TEST(DFNTest, CalcoloIntersezionePiani){
    DFN Fract;
    unsigned int numIntersezioniFratture = 0;
    bool result = calcoloIntersezionePiani(Fract, numIntersezioniFratture);
    ASSERT_TRUE(result);
}

// Test per la funzione calcoloIntersezioneRettaTracciaERettalati
TEST(DFNTest, CalcoloIntersezioneRettaTracciaERettalati){
    DFN Fract;
    unsigned int numeroTracceTotali = 0;
    bool result = calcoloIntersezioneRettaTracciaERettalati(Fract, numeroTracceTotali);
    ASSERT_TRUE(result);
}

// Test per la funzione stampaDatiSuiFileDiOutput
TEST(DFNTest, StampaDatiSuiFileDiOutput){
    DFN Fract;
    string percorsoFileOutputPuntiDiIntersezione = "DFN/puntiDiIntersezione.txt";
    string percorsoFileOutputLunghezzaTracce = "DFN/lunghezzaTracce.txt";
    unsigned int numeroTracceTotali = 0;
    bool result = stampaDatiSuiFileDiOutput(percorsoFileOutputPuntiDiIntersezione, percorsoFileOutputLunghezzaTracce, Fract, numeroTracceTotali);
    ASSERT_TRUE(result);
}

TEST(DFNTest, StampaDatiSulFileVTKDiParaview){
    DFN Fract;
    string fileVTK = "./intersezioniTracce.vtk";
    bool result = stampaDatiSulFileVTKDiParaview(fileVTK, Fract);
    ASSERT_TRUE(result);
}

#endif
