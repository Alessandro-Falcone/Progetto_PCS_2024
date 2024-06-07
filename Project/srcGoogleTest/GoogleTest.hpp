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

DFN Fract;

// Test per la funzione letturaDatiFileFR
TEST(TestFunzioneLetturaDatiFileFR, NonEsisteIlFile){
    unsigned int numFract = 0;
    string percorsoFileFR = "DFN/non_esiste.txt";
    bool result = letturaDatiFileFR(percorsoFileFR, Fract, numFract);
    EXPECT_EQ(numFract, 0);
    ASSERT_FALSE(result);
}

TEST(TestFunzioneLetturaDatiFileFR, FileApertoCorrettamente){
    unsigned int numFract = 0;
    string percorsoFileFR = "DFN/FR3_data.txt";;
    bool result = letturaDatiFileFR(percorsoFileFR, Fract, numFract);
    EXPECT_NE(numFract, 0);
    ASSERT_TRUE(result);
    ASSERT_GE(Fract.coordinateFratture.size(), 0); // verifico che la mappa non sia vuota
}

// Test per la funzione calcoloBaricentriEDistBaricentroVertici
TEST(TestFunzioneCalcoloBaricentriEDistBaricentroVertici, RestituisceTrue){
    unsigned int numFract = 0;
    unsigned int numIntersezioniFratture = 0;
    bool result = calcoloBaricentriEDistBaricentroVertici(Fract, numFract, numIntersezioniFratture);
    ASSERT_TRUE(result);
    ASSERT_GE(Fract.maxDistanzaBaricentri.size(), 0); // verifico che la mappa non sia vuota
    ASSERT_GE(Fract.idFrattureCheSiIntersecano.size(), 0); // verifico che la mappa non sia vuota
}

TEST(DistanzaMassimaBaricentroTest, DistanzaBaricentri){
    double tol = 1e+4 * numeric_limits<double>::epsilon();

    for(const auto& frattura : Fract.maxDistanzaBaricentri){
        double maxDistanza = Fract.maxDistanzaBaricentri[frattura.second];
        ASSERT_GT(maxDistanza, tol);
    }
}

// TEST SULLA DISTANZA MASSIMA TRA BARICENTRO E VERTICI e le coordinate baricentriche
TEST(DistanzaMassimaBaricentroTest, BasicTest){
    // DFN Fract;
    double tol = 1e+4 * numeric_limits<double>::epsilon();
    // Creazione di una matrice di coordinate per una frattura con 4 vertici
    MatrixXd matrCoordinateFratture(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture << 0, 1, 1, 0,
                              0, 0, 1, 1,
                              0, 0, 0, 0;

    // Assegna le coordinate della frattura all'id 0
    Fract.coordinateFratture[0] = matrCoordinateFratture;

    // Valori attesi
    double expectedMaxDistanza = 0.5; // Distanza tra il baricentro e l'angolo (0,1,0)
    // Numero di fratture e intersezioni
    unsigned int numFract = 1;
    unsigned int numIntersezioniFratture = 0;
    // Chiamata alla funzione
    bool result = calcoloBaricentriEDistBaricentroVertici(Fract, numFract, numIntersezioniFratture);

    // Verifica del risultato
    ASSERT_TRUE(result);
    EXPECT_NEAR(Fract.maxDistanzaBaricentri[0], expectedMaxDistanza, tol); // Verifica della massima distanza
}

// TEST(TestFunzioneCalcoloBaricentriEDistBaricentroVertici, MappaidFrattureCheSiIntersecanoNonVuota){
//     unsigned int numeroDiIntersezioni = 0;
//     for(unsigned int i = 0; i < Fract.idFrattureCheSiIntersecano.size(); i++){
//         numeroDiIntersezioni++;
//     }
//     ASSERT_GE(numeroDiIntersezioni, 0);
//     ASSERT_GE(Fract.idFrattureCheSiIntersecano.size(), 0);
// }

// Test per la funzione calcoloEqPianoEdEqRetteLati
TEST(TestFunzioneCalcoloEqPianoEdEqRetteLati, RestituisceTrue){
    bool result = calcoloEqPianoEdEqRetteLati(Fract);
    ASSERT_TRUE(result);
    ASSERT_GE(Fract.coeffabcPiano.size(), 0); // verifico che la mappa non sia vuota
    ASSERT_GE(Fract.coeffdPiano.size(), 0); // verifico che la mappa non sia vuota
    ASSERT_GE(Fract.coeffDirettoriRettaLati.size(), 0); // verifico che la mappa non sia vuota
}

// TEST sull'equazione dei piani che contengono il poligono e l'equazione delle rette passanti per i lati
TEST(TestFunzioneCalcoloEqPianoEdEqRetteLati, BasicTest) {
    // DFN Fract;
    double tol = 1e+4 * numeric_limits<double>::epsilon();
    // Creazione di una matrice di coordinate per una frattura con 4 vertici
    MatrixXd matrCoordinateFratture(3, 4); // 3 righe per x, y, z; 4 colonne per i 4 vertici
    matrCoordinateFratture << 0, 1, 1, 0,
                              0, 0, 1, 1,
                              0, 0, 0, 0;

    // Assegna le coordinate della frattura all'id 0
    Fract.coordinateFratture[0] = matrCoordinateFratture;

    // Valori attesi
    Vector3d expectedabcPiano(0, 0, 1);
    double expectedtermineNotodPiano = 0;

    MatrixXd expectedCoeffDirettoriRettaLati(4, 3);
    expectedCoeffDirettoriRettaLati <<  1,  0, 0,
                                        0,  1, 0,
                                       -1,  0, 0,
                                        0, -1, 0;

    // Chiamata alla funzione
    bool result = calcoloEqPianoEdEqRetteLati(Fract);

    // Verifica del risultato
    ASSERT_TRUE(result);

    // Verifica delle equazioni del piano
    ASSERT_EQ(Fract.coeffabcPiano.size(), Fract.coordinateFratture.size());
    ASSERT_EQ(Fract.coeffdPiano.size(), Fract.coordinateFratture.size());
    EXPECT_NEAR(Fract.coeffabcPiano[0][0] , expectedabcPiano(0), tol);
    EXPECT_NEAR(Fract.coeffabcPiano[0][1] , expectedabcPiano(1), tol);
    EXPECT_NEAR(Fract.coeffabcPiano[0][2] , expectedabcPiano(2), tol);
    EXPECT_NEAR(Fract.coeffdPiano[0], expectedtermineNotodPiano, tol);

    ASSERT_EQ(Fract.coeffDirettoriRettaLati.size(), Fract.coordinateFratture.size());
    EXPECT_EQ(Fract.coeffDirettoriRettaLati[0].rows(), expectedCoeffDirettoriRettaLati.rows());
    EXPECT_EQ(Fract.coeffDirettoriRettaLati[0].cols(), expectedCoeffDirettoriRettaLati.cols());

    for(unsigned int i = 0; i < expectedCoeffDirettoriRettaLati.rows(); i++){
        for(unsigned int j = 0; j < expectedCoeffDirettoriRettaLati.cols(); j++){
            EXPECT_NEAR(Fract.coeffDirettoriRettaLati[0](i, j), expectedCoeffDirettoriRettaLati(i, j), tol);
        }
    }
}

// Test per la funzione calcoloIntersezionePiani
TEST(TestFunzioneCalcoloIntersezionePiani, RestituisceTrue){
    unsigned int numIntersezioniFratture = 0;
    bool result = calcoloIntersezionePiani(Fract, numIntersezioniFratture);
    ASSERT_TRUE(result);
    ASSERT_GE(Fract.coordinatePuntoP.size(), 0); // verifico che la mappa non sia vuota
    ASSERT_GE(Fract.coeffDirettoriRettaTraccia.size(), 0); // verifico che la mappa non sia vuota
}

// Test per la funzione calcoloIntersezioneRettaTracciaERettalati
TEST(TestFunzioneCalcoloIntersezioneRettaTracciaERettalati, RestituisceTrue){
    unsigned int numeroTracceTotali = 0;
    bool result = calcoloIntersezioneRettaTracciaERettalati(Fract, numeroTracceTotali);
    ASSERT_TRUE(result);
    ASSERT_GE(Fract.coordinateIntersezioniTracce.size(), 0); // verifico che la mappa non sia vuota
    ASSERT_GE(Fract.traccePassantiONonPassanti1.size(), 0); // verifico che la mappa non sia vuota
    ASSERT_GE(Fract.traccePassantiONonPassanti2.size(), 0); // verifico che la mappa non sia vuota
    ASSERT_GE(Fract.lunghezzaTracce.size(), 0); // verifico che la mappa non sia vuota
}

TEST(TestFunzioneCalcoloIntersezioneRettaTracciaERettalati, lunghezzaTracce){

    double tol = 1e+4 * numeric_limits<double>::epsilon();

    for(const auto& traccia : Fract.lunghezzaTracce){
        unsigned int idFratt1 = traccia.first.first;
        unsigned int idFratt2 = traccia.first.second;
        double lunghezzaTraccia = Fract.lunghezzaTracce[make_pair(idFratt1, idFratt2)];
        ASSERT_GT(lunghezzaTraccia, tol);
    }
}

// Test per la funzione stampaDatiSuiFileDiOutput
TEST(TestFunzioneStampaDatiSuiFileDiOutput, AperturaCorrettaERestituisceTrue){
    string percorsoFileOutputPuntiDiIntersezione = "DFN/puntiDiIntersezione.txt";
    string percorsoFileOutputLunghezzaTracce = "DFN/lunghezzaTracce.txt";

    string controlloFileOutputPuntiDiIntersezione = "DFN/punto";
    string controlloFileOutputLunghezzaTracce = "DFN/lunghezza";

    EXPECT_NE(percorsoFileOutputPuntiDiIntersezione, controlloFileOutputPuntiDiIntersezione);
    EXPECT_NE(percorsoFileOutputLunghezzaTracce, controlloFileOutputLunghezzaTracce);

    unsigned int numeroTracceTotali = 0;
    bool result = stampaDatiSuiFileDiOutput(percorsoFileOutputPuntiDiIntersezione, percorsoFileOutputLunghezzaTracce, Fract, numeroTracceTotali);
    ASSERT_TRUE(result);
}

TEST(TestFunzioneStampaDatiSulFileVTKDiParaview, AperturaCorrettaERestituisceTrue){
    string percorsoFileVTK = "DFN/intersezioniTracce.vtk";

    string controlloFileVTK = "DFN/nonEsiste.vtk";

    EXPECT_NE(percorsoFileVTK, controlloFileVTK);

    bool result = stampaDatiSulFileVTKDiParaview(percorsoFileVTK, Fract);
    ASSERT_TRUE(result);
}

// TEST(CalcoloIntersezionePianiTest, BasicTest) {
//     DFN fract;

//     // Inizializzazione della struttura Fract con dati di esempio
//     fract.idFrattureCheSiIntersecano[0] = std::make_pair(1, 2); // Esempio di fratture che si intersecano

//     fract.coeffabcPiano[1] = Vector3d(1, 2, 3); // Coefficienti a, b, c del piano della frattura 1
//     fract.coeffdPiano[1] = 4.0; // Termine noto d del piano della frattura 1

//     fract.coeffabcPiano[2] = Vector3d(2, 1, 3); // Coefficienti a, b, c del piano della frattura 2
//     fract.coeffdPiano[2] = 5.0; // Termine noto d del piano della frattura 2

//     unsigned int numIntersezioniFratture = 1; // Numero di intersezioni

//     // Chiamata alla funzione
//     bool result = calcoloIntersezionePiani(fract, numIntersezioniFratture);

//     // Verifica del risultato
//     ASSERT_TRUE(result);

//     // Verifica se le intersezioni sono state calcolate correttamente
//     ASSERT_EQ(fract.coordinatePuntoP.size(), 1); // Dovrebbe esserci una sola intersezione
//     ASSERT_EQ(fract.coeffDirettoriRettaTraccia.size(), 1); // Dovrebbe esserci una sola intersezione

//     // Verifica delle coordinate del punto P
//     Vector3d expectedPuntoP(1.0, 1.0, -1.0); // Coordinate del punto P atteso
//     ASSERT_TRUE(fract.coordinatePuntoP.count(std::make_pair(1, 2)) > 0); // Verifica se l'intersezione è stata salvata
//     EXPECT_TRUE(fract.coordinatePuntoP[std::make_pair(1, 2)].isApprox(expectedPuntoP)); // Verifica se le coordinate del punto P sono corrette

//     // Verifica dei coefficienti direttori della retta di traccia
//     Vector3d expectedCoeffDirettoriRettaTraccia(3.0, -3.0, -5.0); // Coefficienti direttori attesi
//     ASSERT_TRUE(fract.coeffDirettoriRettaTraccia.count(std::make_pair(1, 2)) > 0); // Verifica se i coefficienti direttori sono stati salvati
//     EXPECT_TRUE(fract.coeffDirettoriRettaTraccia[std::make_pair(1, 2)].isApprox(expectedCoeffDirettoriRettaTraccia)); // Verifica se i coefficienti direttori sono corretti
// }





#endif

