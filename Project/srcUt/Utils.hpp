#pragma once

#include <iostream>
#include "DFN.hpp"

using namespace std;

namespace DFNLibrary{

// funzione per la lettura dei dati dal file
bool letturaDatiFileFR(const string &percorsoFileFR, Frattura &Fratt, unsigned int &numFract);

// segue una serie di funzioni che servono per andare a escludere fratture lontane dall'intersezione, attraverso il calcolo di baricentri
bool calcoloBaricentriEDistBaricentroVertici(Frattura &Fratt, unsigned int &numFract, unsigned int &numIntersezioniFratture);

// con questa funzione calcolo l'equazione del piano dei vertici delle fratture e tutte le rette dei lati delle fratture
bool calcoloEqPianoEdEqRetteLati(Frattura &Fratt);

bool calcoloIntersezionePiani(Frattura &Fratt, Traccia &Trac, unsigned int &numIntersezioniFratture);

bool calcoloIntersezioneRettaTracciaERettalati(Frattura &Fratt, Traccia &Trac, unsigned int &numeroTracceTotali);

bool stampaDatiSuiFileDiOutput(const string &percorsoFileOutputPuntiDiIntersezione, const string &percorsoFileOutputLunghezzaTracce, Frattura &Fratt, Traccia &Trac, unsigned int &numeroTracceTotali);

}
