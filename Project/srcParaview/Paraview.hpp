#pragma once

#include "DFN.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;
namespace DFNLibrary{

bool stampaDatiSulFileVTKDiParaview(string &percorsoFileVTK, DFN &Fract){

    ofstream fileVTK;
    fileVTK.open(percorsoFileVTK);

    if(fileVTK.fail()){
        cerr << "errore: impossibile aprire il file " << percorsoFileVTK << endl;
        return false;
    }else{

        // Scrivo l'intestazione
        fileVTK << "# vtk DataFile Version 3.0" << endl;
        fileVTK << "VTK file for point and line visualization" << endl;
        fileVTK << "ASCII" << endl;
        fileVTK << "DATASET POLYDATA" << endl;

        // Scrivo i punti
        fileVTK << "POINTS " << 2 * Fract.coordinateIntersezioniTracce.size() << " double" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            fileVTK << fixed << scientific << setprecision(9) << coords(0, 0) << " " << coords(0, 1) << " " << coords(0, 2) << endl;
            fileVTK << fixed << scientific << setprecision(9) << coords(1, 0) << " " << coords(1, 1) << " " << coords(1, 2) << endl;
        }

        // Scrivo i segmenti
        fileVTK << "LINES " << Fract.coordinateIntersezioniTracce.size() << " " << 3 * Fract.coordinateIntersezioniTracce.size() << endl;
        for(unsigned int i = 0; i < Fract.coordinateIntersezioniTracce.size(); i++){
            fileVTK << "2 " << 2 * i << " " << 2 * i + 1 << endl;
        }

        // Aggiungo dati scalari ai punti
        fileVTK << "POINT_DATA " << 2 * Fract.coordinateIntersezioniTracce.size() << endl;

        fileVTK << "SCALARS TraceId int 1" << endl;
        fileVTK << "LOOKUP_TABLE default" << endl;
        for(unsigned int i = 0; i < Fract.coordinateIntersezioniTracce.size(); i++){
            fileVTK << i << endl;
            fileVTK << i << endl;
        }

        fileVTK << "SCALARS FractureId1 int 1" << endl;
        fileVTK << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            fileVTK << intersezioni.first.first << endl;
            fileVTK << intersezioni.first.first << endl;
        }

        fileVTK << "SCALARS FractureId2 int 1" << endl;
        fileVTK << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            fileVTK << intersezioni.first.second << endl;
            fileVTK << intersezioni.first.second << endl;
        }

        // Aggiungo vettori ai punti (direzioni dei segmenti)
        fileVTK << "VECTORS SegmentDirection double" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            double dx = coords(1, 0) - coords(0, 0);
            double dy = coords(1, 1) - coords(0, 1);
            double dz = coords(1, 2) - coords(0, 2);
            fileVTK << fixed << scientific << setprecision(9) << dx << " " << dy << " " << dz << endl;
            fileVTK << fixed << scientific << setprecision(9) << dx << " " << dy << " " << dz << endl;
        }

        // Aggiungo dati scalari ai segmenti
        fileVTK << "CELL_DATA " << Fract.coordinateIntersezioniTracce.size() << endl;

        fileVTK << "SCALARS SegmentLength double 1" << endl;
        fileVTK << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            Vector3d punto1(coords(0, 0), coords(0, 1), coords(0, 2));
            Vector3d punto2(coords(1, 0), coords(1, 1), coords(1, 2));
            double quadratoDistanza = (punto2 - punto1).squaredNorm();
            fileVTK << fixed << scientific << setprecision(9) << quadratoDistanza << endl;
        }

        fileVTK << "SCALARS SegmentId int 1" << endl;
        fileVTK << "LOOKUP_TABLE default" << endl;
        for(unsigned int i = 0; i < Fract.coordinateIntersezioniTracce.size(); i++){
            fileVTK << i << endl;
        }

        fileVTK.close();
        return true;
        }
    }
}

// void writeCSVFile(const string& filename, const vector<Trace>& traces){
//     ofstream csvFile(filename);
//     if(!csvFile.is_open()){
//         cerr << "Error: Unable to open file " << filename << endl;
//         return;
//     }

//     // Write header
//     // csvFile << "TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
//     csvFile << "TraceId;X1;Y1;Z1;X2;Y2;Z2" << endl;

//     // Write data
//     // for(const auto& trace : traces){
//     //     csvFile << trace.TraceId << ";" << trace.FractureId1 << ";" << trace.FractureId2 << ";"
//     //             << fixed << scientific << setprecision(9) << trace.X1 << ";" << trace.Y1 << ";" << trace.Z1 << ";"
//     //             << fixed << scientific << setprecision(9) << trace.X2 << ";" << trace.Y2 << ";" << trace.Z2 << endl;
//     // }
//     for(const auto& trace : traces){
//         csvFile << trace.TraceId << ";" << fixed << scientific << setprecision(9) << trace.X1 << ";" << trace.Y1 << ";" << trace.Z1 << ";"
//                 << fixed << scientific << setprecision(9) << trace.X2 << ";" << trace.Y2 << ";" << trace.Z2 << endl;
//     }

//     csvFile.close();
// }

