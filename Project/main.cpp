#include <iostream>
#include <fstream>
#include <iomanip>
#include "Eigen/Eigen"
#include "srcUt/DFN.hpp"
#include "srcUt/Utils.hpp"
#include "srcParaview/Paraview.hpp"
#include "srcParte2/PolygonalMesh.hpp"
#include "srcParte2/Utils2poligoni.hpp"

using namespace std;
using namespace Eigen;
using namespace DFNLibrary;
using namespace PolygonalLibrary;

int main(){

    // prima parte progetto
    DFN Fract;

    string percorsoFileFR = "DFN/FR3_data.txt";
    string percorsoFileOutputPuntiDiIntersezione = "./puntiDiIntersezione.txt";
    string percorsoFileOutputLunghezzaTracce = "./lunghezzaTracce.txt";
    string percorsoFileFrattureParaview = "./fratture.vtk";
    string percorsoFileTracceParaview = "./tracce.vtk";

    unsigned int numFract = 0; // numero di fratture
    unsigned int numIntersezioniFratture = 0; // si conta il numero di fratture che si intersecano
    unsigned int numeroTracceTotali = 0; // numero tracce totali

    if(!letturaDatiFileFR(percorsoFileFR, Fract, numFract)){
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

    // cout << endl;
    cout << "numero delle tracce: " << numeroTracceTotali << endl;
    if(!stampaDatiSuiFileDiOutput(percorsoFileOutputPuntiDiIntersezione, percorsoFileOutputLunghezzaTracce, Fract, numeroTracceTotali)){
        cerr << "Errore: impossibile stampare sui file di output" << endl;
        return 1;
    }

    if(!stampaDatiSulFileFrattureParaview(percorsoFileFrattureParaview, Fract)){
        cerr << "Errore: impossibile stampare sul file fratture da esportare su paraview" << endl;
        return 2;
    }

    if(!stampaDatiSulFileTracceParaview(percorsoFileTracceParaview, Fract)){
        cerr << "Errore: impossibile stampare sul file tracce da esportare su paraview" << endl;
        return 2;
    }

    // seconda parte progetto
    PolygonalMesh mesh;

    if(!letturaMesh(Fract, mesh)){
        cerr << "Errore: impossibile leggere la mesh" << endl;
        return 1;
    }else{
        cout << "tutto ok" << endl;
    }

    return 0;
}
