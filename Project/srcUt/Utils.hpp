#pragma once

#include <iostream>
#include "DFN.hpp"

using namespace std;

namespace DFNLibrary{

// funzione per la lettura dei dati dal file
bool letturaDatiFileFR(const string &percorsoFileFR, Frattura &Fratt, unsigned int &numFract);

// segue una serie di funzioni che servono per andare a escludere fratture lontane dall'intersezione, attraverso il calcolo di baricentri
bool calcoloBaricentriEDistBaricentroVertici(Frattura &Fratt, unsigned int &numFract, unsigned int &numIntersezioniFratture);

// questa funzione calcola l'equazione del piano dei vertici delle fratture e tutte le rette dei lati delle fratture
bool calcoloEqPianoEdEqRetteLati(Frattura &Fratt);

// questa funzione trova l'equazione della retta della traccia
bool calcoloIntersezionePiani(Frattura &Fratt, Traccia &Trac, unsigned int &numIntersezioniFratture);

// questa funzione le coordinate dei punti di intersezione della traccia
bool calcoloIntersezioneRettaTracciaERettalati(Frattura &Fratt, Traccia &Trac, unsigned int &numeroTracceTotali);

// stampa sui 2 file di output
bool stampaDatiSuiFileDiOutput(const string &percorsoFileOutputPuntiDiIntersezione, const string &percorsoFileOutputLunghezzaTracce, Frattura &Fratt, Traccia &Trac, unsigned int &numeroTracceTotali);

}
