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

bool stampaDatiSulFileVTKDiParaview2(string &percorsoFileVTK2, DFN &Fract){

        ofstream fileVTK2;
        fileVTK2.open(percorsoFileVTK2);

        if(fileVTK2.fail()){
            cerr << "errore: impossibile aprire il file " << percorsoFileVTK2 << endl;
            return false;
        }else{
            // Scrivo l'intestazione
            fileVTK2 << "# vtk DataFile Version 3.0" << endl;
            fileVTK2 << "VTK file for point and line visualization" << endl;
            fileVTK2 << "ASCII" << endl;
            fileVTK2 << "DATASET POLYDATA" << endl;

            size_t totalPoints = 0;
            size_t totalLines = 0;
            size_t totalIndices = 0;

            // Calcolo il numero totale di punti e segmenti
            for(const auto& frattura : Fract.coordinateFratture){
                totalPoints += frattura.second.cols();
                totalLines += frattura.second.cols();
                totalIndices += 3 * frattura.second.cols();
            }

            // Scrivo i punti
            fileVTK2 << "POINTS " << totalPoints << " double" << endl;
            for(const auto& frattura : Fract.coordinateFratture){
                const MatrixXd& coordinateFrattura = frattura.second;
                for(unsigned int i = 0; i < coordinateFrattura.cols(); i++){
                    fileVTK2 << setprecision(4) << coordinateFrattura(0, i) << " "
                             << coordinateFrattura(1, i) << " " << coordinateFrattura(2, i) << endl;
                }
            }

            // Scrivo i segmenti
            fileVTK2 << "LINES " << totalLines << " " << totalIndices << endl;
            size_t pointIndex = 0;
            for(const auto& frattura : Fract.coordinateFratture){
                const MatrixXd& coordinateFrattura = frattura.second;
                for(unsigned int i = 0; i < coordinateFrattura.cols(); i++){
                    size_t startIdx = (pointIndex + i) % totalPoints; // Indice inizio del segmento
                    size_t endIdx = (pointIndex + (i + 1)) % totalPoints; // Indice fine del segmento
                    if(endIdx == (coordinateFrattura.cols() + pointIndex) || endIdx == 0){
                        unsigned int end = 0;
                        fileVTK2 << "2 " << startIdx << " " << end + pointIndex << endl;
                    }else{
                        fileVTK2 << "2 " << startIdx << " " << endIdx << endl;
                    }
                }
                pointIndex += coordinateFrattura.cols();
            }

            // Aggiungo dati scalari ai punti
            fileVTK2 << "POINT_DATA " << totalPoints << endl;
            fileVTK2 << "SCALARS FractureId int 1" << endl;
            fileVTK2 << "LOOKUP_TABLE default" << endl;
            for(const auto& frattura : Fract.coordinateFratture){
                for(unsigned int i = 0; i < frattura.second.cols(); i++){
                    fileVTK2 << frattura.first << endl;
                }
            }

            pointIndex = 0;
            // Aggiungo vettori ai punti (direzioni dei segmenti)
            fileVTK2 << "VECTORS SegmentDirection double" << endl;
            for(const auto& frattura : Fract.coordinateFratture){
                const MatrixXd& coordinateFrattura = frattura.second;
                for(unsigned int i = 0; i < coordinateFrattura.cols(); i++){
                    size_t startIdx = (i) % totalPoints; // Indice inizio del segmento
                    size_t endIdx = (i + 1) % totalPoints; // Indice fine del segmento
                    if(endIdx == coordinateFrattura.cols()){
                        unsigned int end = 0;
                        Vector3d segmenti = coordinateFrattura.col(end) - coordinateFrattura.col(startIdx);
                        fileVTK2 << setprecision(4) << segmenti(0) << " " << segmenti(1) << " " << segmenti(2) << endl;
                    }else{
                        Vector3d segmenti = coordinateFrattura.col(endIdx) - coordinateFrattura.col(startIdx);
                        fileVTK2 << setprecision(4) << segmenti(0) << " " << segmenti(1) << " " << segmenti(2) << endl;
                    }
                }
                pointIndex += coordinateFrattura.cols();
            }

            pointIndex = 0;
            // Aggiungo dati scalari ai segmenti
            fileVTK2 << "CELL_DATA " << totalLines << endl;
            fileVTK2 << "SCALARS SegmentLength double 1" << endl;
            fileVTK2 << "LOOKUP_TABLE default" << endl;
            for(const auto& frattura : Fract.coordinateFratture){
                const MatrixXd& coordinateFrattura = frattura.second;
                for(unsigned int i = 0; i < coordinateFrattura.cols(); i++){
                    size_t startIdx = (i) % totalPoints; // Indice inizio del segmento
                    size_t endIdx = (i + 1) % totalPoints; // Indice fine del segmento
                    if(endIdx == coordinateFrattura.cols()){
                        unsigned int end = 0;
                        Vector3d differenza = coordinateFrattura.col(end) - coordinateFrattura.col(startIdx);
                        double quadratoDistanza = differenza.squaredNorm();
                        fileVTK2 << setprecision(4) << quadratoDistanza << endl;
                    }else{
                        Vector3d differenza = coordinateFrattura.col(endIdx) - coordinateFrattura.col(startIdx);
                        double quadratoDistanza = differenza.squaredNorm();
                        fileVTK2 << setprecision(4) << quadratoDistanza << endl;
                    }
                }
            }


            fileVTK2 << "SCALARS SegmentId int 1" << endl;
            fileVTK2 << "LOOKUP_TABLE default" << endl;
            int segmentId = 0;
            for(const auto& frattura : Fract.coordinateFratture){
                for(unsigned int i = 0; i < frattura.second.cols(); i++){
                    fileVTK2 << segmentId++ << endl;
                }
            }

            fileVTK2.close();
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

