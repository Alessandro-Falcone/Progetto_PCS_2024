#include "Paraview.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

namespace DFNLibrary{

bool stampaDatiSulFileFrattureParaview(const string &percorsoFileFrattureParaview, DFN &Fract){

    ofstream fileFrattureParaview;
    fileFrattureParaview.open(percorsoFileFrattureParaview);

    if(fileFrattureParaview.fail()){
        cerr << "errore: impossibile aprire il file " << percorsoFileFrattureParaview << endl;
        return false;
    }else{
        // Scrivo l'intestazione
        fileFrattureParaview << "# vtk DataFile Version 3.0" << endl;
        fileFrattureParaview << "VTK file for point and line visualization" << endl;
        fileFrattureParaview << "ASCII" << endl;
        fileFrattureParaview << "DATASET POLYDATA" << endl;

        unsigned int puntiTotali = 0;
        unsigned int segmentiTotali = 0;
        unsigned int indiciTotali = 0;

        // Calcolo il numero totale di punti e segmenti
        for(const auto& frattura : Fract.coordinateFratture){
            puntiTotali += frattura.second.cols();
            segmentiTotali += frattura.second.cols();
            indiciTotali += 3 * frattura.second.cols();
        }

        // Scrivo i punti
        fileFrattureParaview << "POINTS " << puntiTotali << " double" << endl;
        for(const auto& frattura : Fract.coordinateFratture){
            const MatrixXd& coordinateFrattura = frattura.second;
            for(unsigned int i = 0; i < coordinateFrattura.cols(); i++){
                fileFrattureParaview << setprecision(4) << coordinateFrattura(0, i) << " "
                                     << coordinateFrattura(1, i) << " " << coordinateFrattura(2, i) << endl;
            }
        }

        // Scrivo i segmenti
        fileFrattureParaview << "LINES " << segmentiTotali << " " << indiciTotali << endl;
        unsigned int indicePunto = 0;
        for(const auto& frattura : Fract.coordinateFratture){
            const MatrixXd& coordinateFrattura = frattura.second;
            for(unsigned int i = 0; i < coordinateFrattura.cols(); i++){
                unsigned int indiceDiInizio = (indicePunto + i) % puntiTotali; // Indice inizio del segmento
                unsigned int indiceDiFine = (indicePunto + (i + 1)) % puntiTotali; // Indice fine del segmento
                if(indiceDiFine == (coordinateFrattura.cols() + indicePunto) || indiceDiFine == 0){
                    unsigned int end = 0;
                    fileFrattureParaview << "2 " << indiceDiInizio << " " << end + indicePunto << endl;
                }else{
                    fileFrattureParaview << "2 " << indiceDiInizio << " " << indiceDiFine << endl;
                }
            }
            indicePunto += coordinateFrattura.cols();
        }

        // Aggiungo dati scalari ai punti
        fileFrattureParaview << "POINT_DATA " << puntiTotali << endl;
        fileFrattureParaview << "SCALARS FractureId int 1" << endl;
        fileFrattureParaview << "LOOKUP_TABLE default" << endl;
        for(const auto& frattura : Fract.coordinateFratture){
            for(unsigned int i = 0; i < frattura.second.cols(); i++){
                fileFrattureParaview << frattura.first << endl;
            }
        }

        indicePunto = 0;
        // Aggiungo vettori ai punti (direzioni dei segmenti)
        fileFrattureParaview << "VECTORS SegmentDirection double" << endl;
        for(const auto& frattura : Fract.coordinateFratture){
            const MatrixXd& coordinateFrattura = frattura.second;
            for(unsigned int i = 0; i < coordinateFrattura.cols(); i++){
                unsigned int indiceDiInizio = (i) % puntiTotali; // Indice inizio del segmento
                unsigned int indiceDiFine = (i + 1) % puntiTotali; // Indice fine del segmento
                if(indiceDiFine == coordinateFrattura.cols()){
                    unsigned int end = 0;
                    Vector3d segmenti = coordinateFrattura.col(end) - coordinateFrattura.col(indiceDiInizio);
                    fileFrattureParaview << setprecision(4) << segmenti(0) << " " << segmenti(1) << " " << segmenti(2) << endl;
                }else{
                    Vector3d segmenti = coordinateFrattura.col(indiceDiFine) - coordinateFrattura.col(indiceDiInizio);
                    fileFrattureParaview << setprecision(4) << segmenti(0) << " " << segmenti(1) << " " << segmenti(2) << endl;
                }
            }
            indicePunto += coordinateFrattura.cols();
        }

        // Aggiungo dati scalari ai segmenti
        fileFrattureParaview << "CELL_DATA " << segmentiTotali << endl;
        fileFrattureParaview << "SCALARS SegmentLength double 1" << endl;
        fileFrattureParaview << "LOOKUP_TABLE default" << endl;
        for(const auto& frattura : Fract.coordinateFratture){
            const MatrixXd& coordinateFrattura = frattura.second;
            for(unsigned int i = 0; i < coordinateFrattura.cols(); i++){
                unsigned int startIdx = (i) % puntiTotali; // Indice inizio del segmento
                unsigned int endIdx = (i + 1) % puntiTotali; // Indice fine del segmento
                if(endIdx == coordinateFrattura.cols()){
                    unsigned int end = 0;
                    Vector3d differenza = coordinateFrattura.col(end) - coordinateFrattura.col(startIdx);
                    double quadratoDistanza = differenza.squaredNorm();
                    fileFrattureParaview << setprecision(4) << quadratoDistanza << endl;
                }else{
                    Vector3d differenza = coordinateFrattura.col(endIdx) - coordinateFrattura.col(startIdx);
                    double quadratoDistanza = differenza.squaredNorm();
                    fileFrattureParaview << setprecision(4) << quadratoDistanza << endl;
                }
            }
        }

        fileFrattureParaview << "SCALARS SegmentId int 1" << endl;
        fileFrattureParaview << "LOOKUP_TABLE default" << endl;
        int idSegmento = 0;
        for(const auto& frattura : Fract.coordinateFratture){
            for(unsigned int i = 0; i < frattura.second.cols(); i++){
                fileFrattureParaview << idSegmento++ << endl;
            }
        }

        fileFrattureParaview.close();
        return true;
        }
    }

bool stampaDatiSulFileTracceParaview(const string &percorsoFileTracceParaview, DFN &Fract){

    ofstream fileTracceParaview;
    fileTracceParaview.open(percorsoFileTracceParaview);

    if(fileTracceParaview.fail()){
        cerr << "errore: impossibile aprire il file " << percorsoFileTracceParaview << endl;
        return false;
    }else{

        // Scrivo l'intestazione
        fileTracceParaview << "# vtk DataFile Version 3.0" << endl;
        fileTracceParaview << "VTK file for point and line visualization" << endl;
        fileTracceParaview << "ASCII" << endl;
        fileTracceParaview << "DATASET POLYDATA" << endl;

        // Scrivo i punti
        fileTracceParaview << "POINTS " << 2 * Fract.coordinateIntersezioniTracce.size() << " double" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            fileTracceParaview << fixed << scientific << setprecision(9) << coords(0, 0) << " " << coords(0, 1) << " " << coords(0, 2) << endl;
            fileTracceParaview << fixed << scientific << setprecision(9) << coords(1, 0) << " " << coords(1, 1) << " " << coords(1, 2) << endl;
        }

        // Scrivo i segmenti
        fileTracceParaview << "LINES " << Fract.coordinateIntersezioniTracce.size() << " " << 3 * Fract.coordinateIntersezioniTracce.size() << endl;
        for(unsigned int i = 0; i < Fract.coordinateIntersezioniTracce.size(); i++){
            fileTracceParaview << "2 " << 2 * i << " " << 2 * i + 1 << endl;
        }

        // Aggiungo dati scalari ai punti
        fileTracceParaview << "POINT_DATA " << 2 * Fract.coordinateIntersezioniTracce.size() << endl;

        fileTracceParaview << "SCALARS TraceId int 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(unsigned int i = 0; i < Fract.coordinateIntersezioniTracce.size(); i++){
            fileTracceParaview << i << endl;
            fileTracceParaview << i << endl;
        }

        fileTracceParaview << "SCALARS FractureId1 int 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            fileTracceParaview << intersezioni.first.first << endl;
            fileTracceParaview << intersezioni.first.first << endl;
        }

        fileTracceParaview << "SCALARS FractureId2 int 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            fileTracceParaview << intersezioni.first.second << endl;
            fileTracceParaview << intersezioni.first.second << endl;
        }

        // Aggiungo vettori ai punti (direzioni dei segmenti)
        fileTracceParaview << "VECTORS SegmentDirection double" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            double dx = coords(1, 0) - coords(0, 0);
            double dy = coords(1, 1) - coords(0, 1);
            double dz = coords(1, 2) - coords(0, 2);
            fileTracceParaview << fixed << scientific << setprecision(9) << dx << " " << dy << " " << dz << endl;
            fileTracceParaview << fixed << scientific << setprecision(9) << dx << " " << dy << " " << dz << endl;
        }

        // Aggiungo dati scalari ai segmenti
        fileTracceParaview << "CELL_DATA " << Fract.coordinateIntersezioniTracce.size() << endl;

        fileTracceParaview << "SCALARS SegmentLength double 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(const auto& intersezioni : Fract.coordinateIntersezioniTracce){
            const MatrixXd& coords = intersezioni.second;
            Vector3d punto1(coords(0, 0), coords(0, 1), coords(0, 2));
            Vector3d punto2(coords(1, 0), coords(1, 1), coords(1, 2));
            double quadratoDistanza = (punto2 - punto1).squaredNorm();
            fileTracceParaview << fixed << scientific << setprecision(9) << quadratoDistanza << endl;
        }

        fileTracceParaview << "SCALARS SegmentId int 1" << endl;
        fileTracceParaview << "LOOKUP_TABLE default" << endl;
        for(unsigned int i = 0; i < Fract.coordinateIntersezioniTracce.size(); i++){
            fileTracceParaview << i << endl;
        }

        fileTracceParaview.close();
        return true;
        }
    }
}
