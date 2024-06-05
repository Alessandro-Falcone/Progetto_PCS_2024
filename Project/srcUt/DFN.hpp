#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace DFNLibrary{

struct DFN{

    // definiamo una mappa che associa a una chiave (che coincide con l'id nel testo di lettura) la matrice delle coordinate dei vertici
    map<unsigned int, MatrixXd> coordinateFratture = {};

    // in base al numero di intersezioni, si tiene traccia di quali fratture si intersecano tramite il loro id
    map<unsigned int, pair<unsigned int, unsigned int>> idFrattureCheSiIntersecano = {};

    // tramite queste due mappe memorizziamo i coefficienti a, b, c e d, ax + by + cz + d = 0, delle equazioni dei piani che contengono le fratture
    map<unsigned int, Vector3d> coeffabcPiano = {};
    map<unsigned int, double> coeffdPiano = {};

    // in queste 2 mappe salviamo tutta l'informazione sulle rette, espresse in forma parametrica dei lati delle fratture,
    // nella prima salviamo i coefficienti direttori delle rette ovvero l, m, n, invece non ho bisogno di una seconda mappa
    // in quanto i punti iniziali delle rette ovvero x0, y0, z0 li ho già salvati nelle colonne della matrice delle coordinate fratture
    // x = x0 + l*t // y = y0 + m*t // z = z0 + n*t
    map<unsigned int, MatrixXd> coeffDirettoriRettaLati = {};

    // in queste 2 mappe salvo il punto iniziale P delle rette delle tracce e i coefficienti direttori delle rette tracce
    map<pair<unsigned int, unsigned int>, Vector3d> coordinatePuntoP = {};
    map<pair<unsigned int, unsigned int>, Vector3d> coeffDirettoriRettaTraccia = {};

    // in questa mappa salvo le coordinate dei 2 punti di intersezione che definiscono la traccia, identificata dalle 2 fratture
    map<pair<unsigned int, unsigned int>, MatrixXd> coordinateIntersezioniTracce = {};
    // in queste altre mappe che si possono cambiare in vettori??
    // salvo nelle prime 2 l'informazione che la traccia, identificata dalle 2 fratture, è passante o non passante
    // nella prima la prima colonna nella seconda la seconda colonna, sarà più chiaro nel corso del codice
    map<pair<unsigned int, unsigned int>, bool> traccePassantiONonPassanti1 = {};
    map<pair<unsigned int, unsigned int>, bool> traccePassantiONonPassanti2 = {};
    // in questa mappa invece salvo la lunghezza delle tracce, identificata dalle 2 fratture
    map<pair<unsigned int, unsigned int>, double> lunghezzaTracce = {};

    };
}

