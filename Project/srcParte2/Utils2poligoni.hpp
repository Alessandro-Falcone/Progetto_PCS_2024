#pragma once

#include <iostream>
#include "PolygonalMesh.hpp"
#include "DFN.hpp"

using namespace std;
using namespace DFNLibrary;

namespace PolygonalLibrary{

// Importa la mesh poligonale e verifica se la mesh è corretta
// mesh: una struct PolygonalMesh
// restituisce il risultato della lettura, vero se ha successo falso altrimenti
bool letturaMesh(Frattura &Fratt, Traccia &Trac, PolygonalMesh &mesh);

// Importa le proprietà della Cell0Ds dal file Cell0Ds.csv
// mesh: una struct PolygonalMesh
// restituisce il risultato della lettura, vero se ha successo falso altrimenti
bool letturaDatiCell0Ds(Frattura &Fratt, Traccia &Trac, PolygonalMesh &mesh);

// // Importa le proprietà della Cell1Ds dal file Cell1Ds.csv
// // mesh: una struct PolygonalMesh
// // restituisce il risultato della lettura, vero se ha successo falso altrimenti
// bool letturaDatiCell1Ds(DFN &Fract, PolygonalMesh &mesh);

// // Importa le proprietà della Cell2Ds dal file Cell2Ds.csv
// // mesh: una struct PolygonalMesh
// // restituisce il risultato della lettura, vero se ha successo falso altrimenti
// bool letturaDatiCell2Ds(DFN &Fract, PolygonalMesh &mesh);

}
