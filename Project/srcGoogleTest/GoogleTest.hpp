#ifndef __TESTPOLYGONS_H
#define __TESTPOLYGONS_H

#include <gtest/gtest.h>
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include "Paraview.hpp"
#include <map>

using std::map;
using namespace Eigen;
using namespace std;
using namespace DFNLibrary;

// Test per la funzione letturaDatiFileFR
TEST(TestFunzioneLetturaDatiFileFR, NonEsisteIlFile){

    Frattura Fratt;
    unsigned int numFract = 0;
    string percorsoFileFR = "DFN/non_esiste.txt";
    // Chiamata alla funzione
    bool result = letturaDatiFileFR(percorsoFileFR, Fratt, numFract);

    // Verifica del risultato
    ASSERT_FALSE(result);

    EXPECT_EQ(numFract, 0);
}

TEST(TestFunzioneLetturaDatiFileFR, FileApertoCorrettamente){

    Frattura Fratt;
    unsigned int numFract = 0;
    string percorsoFileFR = "DFN/FR3_data.txt";

    // Chiamata alla funzione
    bool result = letturaDatiFileFR(percorsoFileFR, Fratt, numFract);
    EXPECT_NE(numFract, 0);

    // Verifica del risultato
    ASSERT_TRUE(result);

    ASSERT_GT(Fratt.coordinateFratture.size(), 0); // verifico che la mappa non sia vuota
}

// TEST sull'equazione dei piani che contengono il poligono e l'equazione delle rette passanti per i lati
// verifichiamo che date le coordinate dei vertici di un poligono
// il codice resituisca correttamente l'eq del piano e le eq delle rette passanti per i lati
TEST(TestFunzioneCalcoloEqPianoEdEqRetteLati, TestEqPianiERette){

    Frattura Fratt;
    double tol = 1e+4 * numeric_limits<double>::epsilon();
    // Creazione di una matrice di coordinate per una frattura con 4 vertici
    MatrixXd matrCoordinateFratture1(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture1 << 0, 1, 1, 0,
                               0, 0, 1, 1,
                               0, 0, 0, 0;

    MatrixXd matrCoordinateFratture2(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture2 << 0.8, 0.8, 0.8, 0.8,
                                 0, 0, 1, 1,
                              -0.1, 0.299999, 0.299999, -0.1;

    Fratt.coordinateFratture.reserve(2);

    // Assegna le coordinate della frattura all'id 0
    Fratt.coordinateFratture.push_back(matrCoordinateFratture1);
    // Assegna le coordinate della frattura all'id 1
    Fratt.coordinateFratture.push_back(matrCoordinateFratture2);

    // Valori attesi
    Vector3d expectedabcPiano(0, 0, 1);
    double expectedtermineNotodPiano = 0;

    MatrixXd expectedCoeffDirettoriRettaLati(4, 3);
    expectedCoeffDirettoriRettaLati <<  1,  0, 0,
                                        0,  1, 0,
                                       -1,  0, 0,
                                        0, -1, 0;

    Fratt.idFrattureCheSiIntersecano.insert({0, make_pair(0, 1)});

    // Chiamata alla funzione
    bool result = calcoloEqPianoEdEqRetteLati(Fratt);

    // Verifica del risultato
    ASSERT_TRUE(result);

    // Verifica delle equazioni del piano
    ASSERT_EQ(Fratt.coeffabcPiano.size(), Fratt.coordinateFratture.size());
    ASSERT_EQ(Fratt.coeffdPiano.size(), Fratt.coordinateFratture.size());
    EXPECT_NEAR(Fratt.coeffabcPiano[0][0] , expectedabcPiano(0), tol);
    EXPECT_NEAR(Fratt.coeffabcPiano[0][1] , expectedabcPiano(1), tol);
    EXPECT_NEAR(Fratt.coeffabcPiano[0][2] , expectedabcPiano(2), tol);
    EXPECT_NEAR(Fratt.coeffdPiano[0], expectedtermineNotodPiano, tol);

    ASSERT_EQ(Fratt.coeffDirettoriRettaLati.size(), Fratt.coordinateFratture.size());
    EXPECT_EQ(Fratt.coeffDirettoriRettaLati[0].rows(), expectedCoeffDirettoriRettaLati.rows());
    EXPECT_EQ(Fratt.coeffDirettoriRettaLati[0].cols(), expectedCoeffDirettoriRettaLati.cols());

    for(unsigned int i = 0; i < expectedCoeffDirettoriRettaLati.rows(); i++){
        for(unsigned int j = 0; j < expectedCoeffDirettoriRettaLati.cols(); j++){
            EXPECT_NEAR(Fratt.coeffDirettoriRettaLati[0](i, j), expectedCoeffDirettoriRettaLati(i, j), tol);
        }
    }
}

// TEST sul calcolo del punto P, punto appartenente alla retta della traccia
// passiamo l'eq dei piani delle due fratture
TEST(TestFunzioneCalcoloIntersezionePiani, TestCoordinatePuntoP){

    Frattura Fratt;
    Traccia Trac;
    double tol = 1e+4 * numeric_limits<double>::epsilon();
    unsigned int numIntersezioniFratture = 1;
    Vector3d vett1 = Vector3d::Zero();
    vett1 << 0, 0, 1;
    Vector3d vett2 = Vector3d::Zero();
    vett2 << -0.4, 0, 0;

    Fratt.coeffabcPiano[0] = vett1;
    Fratt.coeffabcPiano[1] = vett2;
    Fratt.coeffdPiano[0] = 0;
    Fratt.coeffdPiano[1] = 0.32;

    Fratt.idFrattureCheSiIntersecano[0] = make_pair(0,1);
    Vector3d expectedcoordinatePuntoP(0.8, 0, 0);

    // Chiamata alla funzione
    bool result = calcoloIntersezionePiani(Fratt, Trac, numIntersezioniFratture);

    // Verifica del risultato
    ASSERT_TRUE(result);

    EXPECT_NEAR(Trac.coordinatePuntoP[make_pair(0,1)][0] , expectedcoordinatePuntoP(0), tol);
    EXPECT_NEAR(Trac.coordinatePuntoP[make_pair(0,1)][1] , expectedcoordinatePuntoP(1), tol);
    EXPECT_NEAR(Trac.coordinatePuntoP[make_pair(0,1)][2] , expectedcoordinatePuntoP(2), tol);
}

// TEST sull'eq della retta della traccia
// verifichiamo la correttezza dei coefficienti direttori della retta della traccia
TEST(TestFunzioneCalcoloIntersezionePiani, TestCoefficientiDirettoriRettaTraccia){

    Frattura Fratt;
    Traccia Trac;
    double tol = 1e+4 * numeric_limits<double>::epsilon();
    unsigned int numIntersezioniFratture = 1;
    Vector3d vett1 = Vector3d::Zero();
    vett1 << 0, 0, 1;
    Vector3d vett2 = Vector3d::Zero();
    vett2 << -0.4, 0, 0;

    Fratt.coeffabcPiano[0] = vett1;
    Fratt.coeffabcPiano[1] = vett2;
    Fratt.coeffdPiano[0] = 0;
    Fratt.coeffdPiano[1] = 0.32;

    Fratt.idFrattureCheSiIntersecano[0] = make_pair(0,1);
    Vector3d expectedcoefficientiDirettoriRettaTraccia(0, -0.4, 0);

    // Chiamata alla funzione
    bool result = calcoloIntersezionePiani(Fratt, Trac, numIntersezioniFratture);

    // Verifica del risultato
    ASSERT_TRUE(result);

    EXPECT_NEAR(Trac.coeffDirettoriRettaTraccia[make_pair(0,1)][0] , expectedcoefficientiDirettoriRettaTraccia(0), tol);
    EXPECT_NEAR(Trac.coeffDirettoriRettaTraccia[make_pair(0,1)][1] , expectedcoefficientiDirettoriRettaTraccia(1), tol);
    EXPECT_NEAR(Trac.coeffDirettoriRettaTraccia[make_pair(0,1)][2] , expectedcoefficientiDirettoriRettaTraccia(2), tol);
}

// TEST sul punto iniziale e finale del segmento della traccia
TEST(TestFunzioneCalcoloIntersezioneRettaTracciaERettalati, CoordinatePuntiTraccia){

    Frattura Fratt;
    Traccia Trac;
    double tol = 1e+4 * numeric_limits<double>::epsilon();
    unsigned int numeroTracceTotali = 1;

    MatrixXd matrCoordinateFratture1(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture1 << 0, 1, 1, 0,
                               0, 0, 1, 1,
                               0, 0, 0, 0;

    Fratt.coordinateFratture.reserve(2);

    // Assegna le coordinate della frattura all'id 0
    Fratt.coordinateFratture.push_back(matrCoordinateFratture1);

    MatrixXd matrCoordinateFratture2(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture2 << 0.8, 0.8, 0.8, 0.8,
                                0, 0, 1, 1,
                            -0.1, 0.299999999, 0.299999999, -0.1;

    // Assegna le coordinate della frattura all'id 1
    Fratt.coordinateFratture.push_back(matrCoordinateFratture2);

    Vector3d coordinateP(0.8, 0, 0);
    Vector3d coeffDirRettaTraccia(0, -0.4, 0);

    Trac.coordinatePuntoP[make_pair(0,1)] = coordinateP;
    Trac.coeffDirettoriRettaTraccia[make_pair(0,1)] = coeffDirRettaTraccia;

    MatrixXd matrLati1(4, 3); // 4 righe per x, y, z; 3 colonne per i 4 vertici
    matrLati1 << 1, 0, 0,
                 0, 1, 0,
                -1, 0, 0,
                 0, -1, 0;

    MatrixXd matrLati2(4, 3); // 4 righe per x, y, z; 3 colonne per i 4 vertici
    matrLati2 << 0, 0, 0.4,
                 0, 1, 0,
                 0, 0, -0.4,
                 0, -1, 0;

    Fratt.coeffDirettoriRettaLati[0] = matrLati1;
    Fratt.coeffDirettoriRettaLati[1] = matrLati2;

    MatrixXd expectedmatrTraccia(2, 3);
    expectedmatrTraccia << 0.8, 0, 0,
                           0.8, 1, 0;

    // Chiamata alla funzione
    bool result = calcoloIntersezioneRettaTracciaERettalati(Fratt, Trac, numeroTracceTotali);

    // Verifica del risultato
    ASSERT_TRUE(result);

    for(unsigned int i = 0; i < expectedmatrTraccia.rows(); i++){
        for(unsigned int j = 0; j < expectedmatrTraccia.cols(); j++){
            EXPECT_NEAR(Trac.coordinateIntersezioniTracce[make_pair(0,1)](i, j), expectedmatrTraccia(i, j), tol);
        }
    }
}

// TEST sulla lunghezza della traccia
TEST(TestFunzioneCalcoloIntersezioneRettaTracciaERettalati, lunghezzaTracce){

    Frattura Fratt;
    Traccia Trac;
    double tol = 1e+4 * numeric_limits<double>::epsilon();
    unsigned int numeroTracceTotali = 1;

    MatrixXd matrCoordinateFratture1(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture1 << 0, 1, 1, 0,
                               0, 0, 1, 1,
                               0, 0, 0, 0;

    Fratt.coordinateFratture.reserve(2);

    // Assegna le coordinate della frattura all'id 0
    Fratt.coordinateFratture.push_back(matrCoordinateFratture1);

    MatrixXd matrCoordinateFratture2(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture2 << 0.8, 0.8, 0.8, 0.8,
                                 0, 0, 1, 1,
                                -0.1, 0.299999999, 0.299999999, -0.1;

    // Assegna le coordinate della frattura all'id 1
    Fratt.coordinateFratture.push_back(matrCoordinateFratture2);

    Vector3d coordinateP(0.8, 0, 0);
    Vector3d coeffDirRettaTraccia(0, -0.4, 0);

    Trac.coordinatePuntoP[make_pair(0,1)] = coordinateP;
    Trac.coeffDirettoriRettaTraccia[make_pair(0,1)] = coeffDirRettaTraccia;

    MatrixXd matrLati1(4, 3); // 4 righe per x, y, z; 3 colonne per i 4 vertici
    matrLati1 << 1, 0, 0,
                 0, 1, 0,
                -1, 0, 0,
                 0, -1, 0;

    MatrixXd matrLati2(4, 3); // 4 righe per x, y, z; 3 colonne per i 4 vertici
    matrLati2 << 0, 0, 0.4,
                 0, 1, 0,
                 0, 0, -0.4,
                 0, -1, 0;

    Fratt.coeffDirettoriRettaLati[0] = matrLati1;
    Fratt.coeffDirettoriRettaLati[1] = matrLati2;

    double expectedLunghezzaTraccia = 1;

    // Chiamata alla funzione
    bool result = calcoloIntersezioneRettaTracciaERettalati(Fratt, Trac, numeroTracceTotali);

    // Verifica del risultato
    ASSERT_TRUE(result);

    EXPECT_NEAR(Trac.lunghezzaTracce[make_pair(0,1)], expectedLunghezzaTraccia, tol);
}

// TEST sulla lunghezza della traccia nulla
TEST(TestFunzioneCalcoloIntersezioneRettaTracciaERettalati, lunghezzaTracceNulla){

    Frattura Fratt;
    Traccia Trac;
    double tol = 1e+4 * numeric_limits<double>::epsilon();
    unsigned int numeroTracceTotali = 2;

    MatrixXd matrCoordinateFratture1(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture1 << 0, 1, 1, 0,
                               0, 0, 1, 1,
                               0, 0, 0, 0;

    Fratt.coordinateFratture.reserve(3);

    // Assegna le coordinate della frattura all'id 0
    Fratt.coordinateFratture.push_back(matrCoordinateFratture1);

    MatrixXd matrCoordinateFratture2(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture2 << 0.8, 0.8, 0.8, 0.8,
                                 0, 0, 1, 1,
                              -0.1, 0.299999999, 0.299999999, -0.1;

    // Assegna le coordinate della frattura all'id 1
    Fratt.coordinateFratture.push_back(matrCoordinateFratture2);

    MatrixXd matrCoordinateFratture3(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture3 << -0.23777799999999999, 0.31618370000000001, 0.31618370000000001, -0.23777799999999999,
                                0.50000000000000000, 0.50000000000000000, 0.50000000000000000, 0.50000000000000000,
                               -0.34444000000000002, -0.34444000000000002, 0.45283889999999999, 0.45283889999999999;

    // Assegna le coordinate della frattura all'id 2
    Fratt.coordinateFratture.push_back(matrCoordinateFratture3);

    Vector3d coordinateP1(0.8, 0, 0);
    Vector3d coeffDirRettaTraccia1(0, -0.4, 0);

    Vector3d coordinateP2(0, 0.5, 0);
    Vector3d coeffDirRettaTraccia2(0.4416619748, 0, 0);

    Trac.coordinatePuntoP[make_pair(0,1)] = coordinateP1;
    Trac.coeffDirettoriRettaTraccia[make_pair(0,1)] = coeffDirRettaTraccia1;

    Trac.coordinatePuntoP[make_pair(0,2)] = coordinateP2;
    Trac.coeffDirettoriRettaTraccia[make_pair(0,2)] = coeffDirRettaTraccia2;

    MatrixXd matrLati1(4, 3); // 4 righe per x, y, z; 3 colonne per i 4 vertici
    matrLati1 << 1, 0, 0,
                 0, 1, 0,
                -1, 0, 0,
                0, -1, 0;

    MatrixXd matrLati2(4, 3); // 4 righe per x, y, z; 3 colonne per i 4 vertici
    matrLati2 << 0, 0, 0.4,
                 0, 1, 0,
                 0, 0, -0.4,
                 0, -1, 0;

    MatrixXd matrLati3(4, 3); // 4 righe per x, y, z; 3 colonne per i 4 vertici
    matrLati3 << 0.55396170000000000, 0, 0,
                 0, 0, 0.79727890000000001,
                -0.55396170000000000, 0, 0,
                 0, 0, -0.79727890000000001;

    Fratt.coeffDirettoriRettaLati[0] = matrLati1;
    Fratt.coeffDirettoriRettaLati[1] = matrLati2;
    Fratt.coeffDirettoriRettaLati[2] = matrLati3;

    // Chiamata alla funzione
    bool result = calcoloIntersezioneRettaTracciaERettalati(Fratt, Trac, numeroTracceTotali);

    // Verifica del risultato
    ASSERT_TRUE(result);

    ASSERT_GT(Trac.lunghezzaTracce[make_pair(0,1)], tol);
    ASSERT_GT(Trac.lunghezzaTracce[make_pair(0,2)], tol);
}

// Test per la funzione stampaDatiSuiFileDiOutput
TEST(TestFunzioneStampaDatiSuiFileDiOutput, AperturaCorrettaFileDiOutput){

    Frattura Fratt;
    Traccia Trac;
    string percorsoFileOutputPuntiDiIntersezione = "DFN/puntiDiIntersezione.txt";
    string percorsoFileOutputLunghezzaTracce = "DFN/lunghezzaTracce.txt";

    string controlloFileOutputPuntiDiIntersezione = "DFN/punto";
    string controlloFileOutputLunghezzaTracce = "DFN/lunghezza";

    EXPECT_NE(percorsoFileOutputPuntiDiIntersezione, controlloFileOutputPuntiDiIntersezione);
    EXPECT_NE(percorsoFileOutputLunghezzaTracce, controlloFileOutputLunghezzaTracce);

    unsigned int numeroTracceTotali = 0;
    bool result = stampaDatiSuiFileDiOutput(percorsoFileOutputPuntiDiIntersezione, percorsoFileOutputLunghezzaTracce, Fratt, Trac, numeroTracceTotali);
    ASSERT_TRUE(result);
}

// Test per la funzione stampaDatiSulFileFrattureParaview
TEST(TestFunzioneStampaDatiSulFileFrattureParaview, AperturaCorrettaFileFrattureParaview){

    Frattura Fratt;
    string percorsoFileFrattureParaview = "DFN/fratture.vtk";
    string controlloFileFrattureParaview = "DFN/nonEsiste.vtk";

    EXPECT_NE(percorsoFileFrattureParaview, controlloFileFrattureParaview);

    bool result = stampaDatiSulFileFrattureParaview(percorsoFileFrattureParaview, Fratt);
    ASSERT_TRUE(result);
}

// Test per la funzione stampaDatiSulFileTracceParaview
TEST(TestFunzioneStampaDatiSulFileTracceParaview, AperturaCorrettaFileTracceParaview){

    Traccia Trac;
    string percorsoFileTracceParaview = "DFN/tracce.vtk";
    string controlloFileTracceParaview = "DFN/nonEsiste.vtk";

    EXPECT_NE(percorsoFileTracceParaview, controlloFileTracceParaview);

    bool result = stampaDatiSulFileTracceParaview(percorsoFileTracceParaview, Trac);
    ASSERT_TRUE(result);
}

#endif

