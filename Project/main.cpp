#include <iostream>
#include <fstream>
#include <iomanip>
#include "srcUt/DFN.hpp"
#include "srcUt/Utils.hpp"
#include "Eigen/Eigen"

#include <vector>
#include <string>
#include "srcParaview/Paraview.hpp"

using namespace std;
using namespace Eigen;
using namespace DFNLibrary;

int main(){

    DFN Fract;

    string fileFR = "DFN/FR3_data.txt";
    string fileOutputPuntiDiIntersezione = "./puntiDiIntersezione.txt";
    string fileOutputLunghezzaTracce = "./lunghezzaTracce.txt";

    unsigned int numFract = 0; // numero di fratture
    unsigned int numIntersezioniFratture = 0; // si conta il numero di fratture che si intersecano
    unsigned int numeroTracceTotali = 0; // numero tracce totali

    if(!letturaDatiFileFR(fileFR, Fract, numFract)){
        cerr << "Errore: impossibile leggere la frattura" << endl;
        return 1;
    }

    if(!calcoloBaricentriEDistBaricentroVertici(Fract, numFract, numIntersezioniFratture)){
        cerr << "Errore: impossibile calcolare le coordinate dei baricentri o la distanza tra il baricentro e i vertici" << endl;
        return -1;
    }

    cout << "intersezioni: " << numIntersezioniFratture << endl;
    if(!calcoloEqPianoEdEqRetteLati(Fract)){
        cerr << "Errore: impossibile calcolare l'equazione del piano o l'equazione delle rette dei lati" << endl;
        return -1;
    }

    if(!calcoloIntersezionePiani(Fract, numIntersezioniFratture)){
        cerr << "Errore: impossibile calcolare l'equazione della retta della traccia intersezione tra i piani" << endl;
        return -1;
    }

    if(!calcoloIntersezioneRettaTracciaERettalati(Fract, numeroTracceTotali)){
        cerr << "Errore: impossibile calcolare l'intersezione della retta della traccia con i lati" << endl;
        return -1;
    }

    cout << endl;
    cout << "numero delle tracce: " << numeroTracceTotali << endl;
    if(!stampaDatiSuiFileDiOutput(fileOutputPuntiDiIntersezione, fileOutputLunghezzaTracce, Fract, numeroTracceTotali)){
        cerr << "Errore: impossibile stampare sui file di output" << endl;
        return -1;
    }

    // vector<intersezioni> tracce = readTracesFromFile("./puntiDiIntersezione.txt");
    string fileVTK = "./intersezioniTracce.vtk";
    if(!stampaDatiSulFileVTKDiParaview(fileVTK, Fract)){
        cerr << "Errore: impossibile stampare sul file da esportare su paraview" << endl;
        return -1;
    }else{
        cout << "il file VTK e' stato scritto con successo" << endl;
    }

    return 0;
}
