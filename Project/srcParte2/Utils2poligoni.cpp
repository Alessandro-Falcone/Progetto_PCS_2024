#include "Utils2poligoni.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace DFNLibrary;
using namespace std;

namespace PolygonalLibrary{

    bool letturaMesh(Frattura &Fratt, Traccia &Trac, PolygonalMesh &mesh){

        if(!letturaDatiCell0Ds(Fratt, Trac, mesh)){
            return false;
        }

        return true;
    }

    bool letturaDatiCell0Ds(Frattura &Fratt, Traccia &Trac, PolygonalMesh &mesh){

        unsigned int numeroTracceTotali = 0;

        vector<unsigned int> idFract1;
        vector<unsigned int> idFract2;
        vector<bool> passanteONonPassante1;
        vector<bool> passanteONonPassante2;
        vector<double> lunghezzaTracce;

        idFract1.reserve(Trac.coordinateIntersezioniTracce.size());
        idFract2.reserve(Trac.coordinateIntersezioniTracce.size());
        passanteONonPassante1.reserve(Trac.coordinateIntersezioniTracce.size());
        passanteONonPassante2.reserve(Trac.coordinateIntersezioniTracce.size());
        lunghezzaTracce.reserve(Trac.coordinateIntersezioniTracce.size());

        for(const auto &frattura : Trac.coordinateIntersezioniTracce){

            unsigned int idFratt1 = frattura.first.first;
            unsigned int idFratt2 = frattura.first.second;
            bool passONonPass1 = Trac.traccePassantiONonPassanti1[make_pair(idFratt1, idFratt2)];
            bool passONonPass2 = Trac.traccePassantiONonPassanti2[make_pair(idFratt1, idFratt2)];
            double lunghezzaTraccia = Trac.lunghezzaTracce[make_pair(idFratt1, idFratt2)];

            idFract1.push_back(idFratt1); // salvo nel vettore idFract1 associato alla prima frattura l'id della frattura 1
            idFract2.push_back(idFratt2); // salvo nel vettore idFract2 associato alla seconda frattura l'id della frattura 2
            passanteONonPassante1.push_back(passONonPass1); // salvo nel vettore passanteONonPassante1 il booleano che dice se la prima frattura è passante o non passante
            passanteONonPassante2.push_back(passONonPass2); // salvo nel vettore passanteONonPassante2 il booleano che dice se la seconda frattura è passante o non passante
            lunghezzaTracce.push_back(lunghezzaTraccia); // salvo nel vettore lunghezzaTracce la lunghezza della traccia
        }

        numeroTracceTotali = Trac.coordinateIntersezioniTracce.size(); // numero di tracce totali per il file scelto

        // unsigned int idCell0Ds = 0; // identificativo celle 0Ds
        // unsigned int idCell1Ds = 0; // identificativo celle 1Ds
        // unsigned int idCell2Ds = 0; // identificativo celle 2Ds
        // unsigned int Cell0Ds = 0, Cell1Ds = 0, Cell2Ds = 0; // contatori celle 0Ds, 1Ds, 2Ds
        for(unsigned int idFrattura = 0; idFrattura < Fratt.coordinateFratture.size(); idFrattura++){

            unsigned int numeroTraccePerFrattura = 0; // contatore numero delle tracce per ciascuna frattura
            unsigned int numeroTraccePassantiPerFrattura = 0; // contatore numero delle tracce passanti per ciascuna frattura

            vector<unsigned int> idTraccia; // vettore id di tutte le tracce
            vector<unsigned int> idTracciaPassante; // vettore id traccia passante
            vector<unsigned int> idTracciaNonPassante; // vettore id traccia non passante
            vector<double> lunghezzaTracciaPassante; // vettore lunghezza traccia passante
            vector<double> lunghezzaTracciaNonPassante; // vettore lunghezza traccia non passante

            for(unsigned int i = 0;  i < numeroTracceTotali; i++){
                if(idFrattura == idFract1[i] ||  idFrattura == idFract2[i]){
                    numeroTraccePerFrattura++;
                    idTraccia.push_back(i);
                }
            }

            if(numeroTraccePerFrattura > 0){
                cout << idFrattura << "; " << numeroTraccePerFrattura << endl;
                for(unsigned int i = 0; i < numeroTraccePerFrattura; i++){

                    if((idFrattura == idFract1[idTraccia[i]]) && (passanteONonPassante1[idTraccia[i]] == false)){
                        unsigned int idTr = idTraccia[i];
                        double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
                        idTracciaPassante.push_back(idTr); // salvo l'id della traccia perchè è passante
                        lunghezzaTracciaPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è passante
                        numeroTraccePassantiPerFrattura++;
                    }else if((idFrattura == idFract2[idTraccia[i]]) && (passanteONonPassante2[idTraccia[i]] == false)){
                        unsigned int idTr = idTraccia[i];
                        double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
                        idTracciaPassante.push_back(idTr); // salvo l'id della traccia perchè è passante
                        lunghezzaTracciaPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è passante
                        numeroTraccePassantiPerFrattura++;
                    }

                    if((idFrattura == idFract1[idTraccia[i]]) && (passanteONonPassante1[idTraccia[i]] == true)){
                        unsigned int idTr = idTraccia[i];
                        double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
                        idTracciaNonPassante.push_back(idTr); // salvo l'id della traccia perchè è non passante
                        lunghezzaTracciaNonPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è non passante
                    }else if((idFrattura == idFract2[idTraccia[i]]) && (passanteONonPassante2[idTraccia[i]] == true)){
                        unsigned int idTr = idTraccia[i];
                        double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
                        idTracciaNonPassante.push_back(idTr); // salvo l'id della traccia perchè è non passante
                        lunghezzaTracciaNonPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è non passante
                    }
                }

                if(idTracciaPassante.size() > 1){
                    for(unsigned int i = 0; i < idTracciaPassante.size() - 1; i++){
                        for(unsigned int j = i + 1; j < idTracciaPassante.size(); j++){
                            if(lunghezzaTracciaPassante[i] < lunghezzaTracciaPassante[j]){
                                swap(idTracciaPassante[i], idTracciaPassante[j]);
                                swap(lunghezzaTracciaPassante[i], lunghezzaTracciaPassante[j]);
                            }
                        }
                    }
                }

                if(idTracciaNonPassante.size() > 1){
                    for(unsigned int i = 0; i < idTracciaNonPassante.size() - 1; i++){
                        for(unsigned int j = i + 1; j < idTracciaNonPassante.size(); j++){
                            if(lunghezzaTracciaNonPassante[i] < lunghezzaTracciaNonPassante[j]){
                                swap(idTracciaNonPassante[i], idTracciaNonPassante[j]);
                                swap(lunghezzaTracciaNonPassante[i], lunghezzaTracciaNonPassante[j]);
                            }
                        }
                    }
                }

                bool passante = false;
                MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];
                MatrixXd matrPuntiDiIntersezionePrecTracciaPassante = MatrixXd::Zero(2,3);
                Vector3d diffPuntiDiIntersezioneTracciaPassante = Vector3d::Zero();
                MatrixXd matrPuntiDiIntersezioneSuccTracciaNonPassante = MatrixXd::Zero(2,3);

                if(idTracciaNonPassante.size() > 0){
                    matrPuntiDiIntersezioneSuccTracciaNonPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[0]],idFract2[idTracciaNonPassante[0]])];
                }

                for(unsigned int t = 0; t < idTracciaPassante.size(); t++){
                    cout << idTracciaPassante[t] << "; " << boolalpha << passante << endl;
                    MatrixXd matrPuntiDiIntersezioneTracciaPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaPassante[t]],idFract2[idTracciaPassante[t]])];
                    Vector3d diffPuntiDiIntersezioneTracciaPass = matrPuntiDiIntersezioneTracciaPassante.row(1) - matrPuntiDiIntersezioneTracciaPassante.row(0);

                    Vector3d diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0) - matrPuntiDiIntersezioneTracciaPassante.row(0);
                    // prodotto vettoriale 1 tra AB e AP1
                    Vector3d prodVett1 = diffPuntiDiIntersezioneTracciaPass.cross(diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1);
                    // z positiva a sinistra z negativa a destra
                    // differenza secondo dei due punti di intersezione vett AP2 = P2 - A
                    Vector3d diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1) - matrPuntiDiIntersezioneTracciaPassante.row(0);
                    // prodotto vettoriale 2 tra AB e AP2
                    Vector3d prodVett2 = diffPuntiDiIntersezioneTracciaPass.cross(diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2);

                    Vector3d absVec1 = prodVett1.cwiseAbs();

                    // Trova il massimo valore e la sua posizione
                    Index maxIndex1; // Definisce l'indice massimo usando Eigen::Index
                    absVec1.maxCoeff(&maxIndex1);

                    Vector3d absVec2 = prodVett2.cwiseAbs();

                    // Trova il massimo valore e la sua posizione
                    Index maxIndex2; // Definisce l'indice massimo usando Eigen::Index
                    absVec2.maxCoeff(&maxIndex2);

                    for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                        // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
                        Vector3d diffCoordinataFratturaPuntoDiIntersezioneTracciaPassante = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaPassante.row(0);
                        // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
                        Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTracciaPassante);

                        cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
                        Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

                        // Trova il massimo valore e la sua posizione
                        Index maxIndex; // Definisce l'indice massimo usando Index
                        absVec.maxCoeff(&maxIndex);

                        if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && prodVett1(maxIndex1) > 0 && prodVett2(maxIndex2) > 0){
                            cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                        }else if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && idTracciaNonPassante.size() == 0 && idTracciaPassante.size() == 1){
                            cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                        }
                    }

                    for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                        // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
                        Vector3d diffCoordinataFratturaPuntoDiIntersezioneTracciaPassante = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaPassante.row(0);
                        // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
                        Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTracciaPassante);

                        cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
                        Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

                        // Trova il massimo valore e la sua posizione
                        Index maxIndex; // Definisce l'indice massimo usando Index
                        absVec.maxCoeff(&maxIndex);

                        if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && prodVett1(maxIndex1) < 0 && prodVett2(maxIndex2) < 0){
                            cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                        }else if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && idTracciaNonPassante.size() == 0 && idTracciaPassante.size() == 1){
                            cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                        }
                    }
                    // matrice punti di intersezione precedente
                    matrPuntiDiIntersezionePrecTracciaPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaPassante[idTracciaPassante.size() -1]],idFract2[idTracciaPassante[idTracciaPassante.size() -1]])];
                    // differenza tra questi punti vett AB = B - A
                    diffPuntiDiIntersezioneTracciaPassante = matrPuntiDiIntersezionePrecTracciaPassante.row(1) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
                }

                bool nonPassante = true;
                for(unsigned int t = 0; t < idTracciaNonPassante.size(); t++){
                    cout << idTracciaNonPassante[t] << "; " << boolalpha << nonPassante << endl;

                    if((t < 1) && (numeroTraccePassantiPerFrattura > 0) && idTracciaNonPassante.size() > 0){
                        MatrixXd matrPuntiDiIntersezioneSuccTracciaNonPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[0]],idFract2[idTracciaNonPassante[0]])];
                        // differenza primo dei due punti di intersezione vett AP1 = P1 - A
                        Vector3d diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
                        // prodotto vettoriale 1 tra AB e AP1
                        Vector3d prodVett1 = diffPuntiDiIntersezioneTracciaPassante.cross(diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1);
                        // z positiva a sinistra z negativa a destra
                        // differenza secondo dei due punti di intersezione vett AP2 = P2 - A
                        Vector3d diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
                        // prodotto vettoriale 2 tra AB e AP2
                        Vector3d prodVett2 = diffPuntiDiIntersezioneTracciaPassante.cross(diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2);
                        // z positiva a sinistra z negativa a destra

                        Vector3d absVec1 = prodVett1.cwiseAbs();

                        // Trova il massimo valore e la sua posizione
                        Index maxIndex1; // Definisce l'indice massimo usando Eigen::Index
                        absVec1.maxCoeff(&maxIndex1);

                        Vector3d absVec2 = prodVett2.cwiseAbs();

                        // Trova il massimo valore e la sua posizione
                        Index maxIndex2; // Definisce l'indice massimo usando Eigen::Index
                        absVec2.maxCoeff(&maxIndex2);
                        // Calcola i vettori AB già calcolato è diffPuntiDiIntersezioneTracciaPassante e AP già calcolato diffSecondoPuntoDiIntersezioneSucc1

                        // Proiezione di AP su AB
                        double t = diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2.dot(diffPuntiDiIntersezioneTracciaPassante) / diffPuntiDiIntersezioneTracciaPassante.dot(diffPuntiDiIntersezioneTracciaPassante);

                        // Punto di proiezione punto di intersezione traccia non passante su traccia passante
                        Vector3d proiezionePuntoDiInterSuTracciaPass = matrPuntiDiIntersezionePrecTracciaPassante.row(0) + t * diffPuntiDiIntersezioneTracciaPassante.transpose();

                        Vector3d vettProiezionePuntoDiIntersezioneTracciaNonPass = proiezionePuntoDiInterSuTracciaPass.transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);

                        if(prodVett1(maxIndex1) > 0 && prodVett2(maxIndex2) > 0){

                            // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
                            // Per ogni cella 0D: un identificativo e le coordinate 3D
                            MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];

                            for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                                // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
                                Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
                                // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
                                Vector3d prodVettCoordFrattPunDiInterTraccia = vettProiezionePuntoDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
                                cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
                                Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

                                // Trova il massimo valore e la sua posizione
                                Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
                                absVec.maxCoeff(&maxIndex);

                                if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(0).squaredNorm()){
                                    cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                                }

                                if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(1).squaredNorm()){
                                    cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                                }
                            }
                        }

                        if(prodVett1(maxIndex1) < 0 && prodVett2(maxIndex2) < 0){

                            // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
                            // Per ogni cella 0D: un identificativo e le coordinate 3D
                            MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];

                            for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                                // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
                                Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
                                // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
                                Vector3d prodVettCoordFrattPunDiInterTraccia = vettProiezionePuntoDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
                                cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
                                Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

                                // Trova il massimo valore e la sua posizione
                                Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
                                absVec.maxCoeff(&maxIndex);

                                if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(0).squaredNorm()){
                                    cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                                }

                                if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(1).squaredNorm()){
                                    cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                                }
                            }
                        }

                        if(prodVett1(maxIndex1) > 0 && prodVett2(maxIndex2) < 0){

                            // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
                            // Per ogni cella 0D: un identificativo e le coordinate 3D
                            MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];

                            for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                                // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
                                Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
                                // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
                                Vector3d prodVettCoordFrattPunDiInterTraccia = vettProiezionePuntoDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
                                cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
                                Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

                                // Trova il massimo valore e la sua posizione
                                Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
                                absVec.maxCoeff(&maxIndex);

                                if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(0).squaredNorm()){
                                    cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                                }

                                if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(1).squaredNorm()){
                                    cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                                }
                            }
                        }

                        if(prodVett1(maxIndex1) < 0 && prodVett2(maxIndex2) > 0){

                            // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
                            // Per ogni cella 0D: un identificativo e le coordinate 3D
                            MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];

                            for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                                // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
                                Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
                                // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
                                Vector3d prodVettCoordFrattPunDiInterTraccia = vettProiezionePuntoDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
                                cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
                                Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

                                // Trova il massimo valore e la sua posizione
                                Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
                                absVec.maxCoeff(&maxIndex);

                                if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(0).squaredNorm()){
                                    cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                                }

                                if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(1).squaredNorm()){
                                    cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                                }
                            }
                        }
                    }

                    if(idTracciaNonPassante.size() > 0 && idTracciaPassante.size() == 0){

                        MatrixXd matrPuntiDiIntersezioneTracciaNonPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[t]],idFract2[idTracciaNonPassante[t]])];
                        Vector3d diffPuntiDiIntersezioneTracciaNonPass = matrPuntiDiIntersezioneTracciaNonPassante.row(1) - matrPuntiDiIntersezioneTracciaNonPassante.row(0);

                        for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                            Vector3d diffPuntoDiIntersezioneCoordinataFrattura = matrPuntiDiIntersezioneTracciaNonPassante.row(1) - coordinataFrattura.col(i).transpose();
                            cout << "diffPuntoDiIntersezioneCoordinataFrattura: " << diffPuntoDiIntersezioneCoordinataFrattura.transpose() << endl;
                            Vector3d diffCoordinateFratture = coordinataFrattura.col((i+1) % coordinataFrattura.cols()) - coordinataFrattura.col(i);
                            cout << "diffCoordinateFratture: " << diffCoordinateFratture.transpose() << endl;
                            double num = diffPuntoDiIntersezioneCoordinataFrattura.dot(diffCoordinateFratture);
                            double den = diffCoordinateFratture.dot(diffCoordinateFratture);
                            double t = num/den;
                            cout << "t: " << t << endl;
                            Vector3d proiezionePuntoDiIntersezioneSuDiffCoordinateFratture = coordinataFrattura.col(i).transpose() + t * diffCoordinateFratture.transpose();
                            proiezionePuntoDiIntersezioneSuDiffCoordinateFratture = (proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.array().abs() < 1e-15).select(0, proiezionePuntoDiIntersezioneSuDiffCoordinateFratture);
                            // cout << "proiezione non passante: " << proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() << endl;
                            if((proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() != coordinataFrattura.col((i+1) % coordinataFrattura.cols()).transpose() && proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() != coordinataFrattura.col(i).transpose()) && proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() != matrPuntiDiIntersezioneTracciaNonPassante.row(1)){
                                cout << "proiezione non passante: " <<proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() << endl;
                            }
                        }

                        for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                            // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
                            Vector3d diffCoordinataFratturaPuntoDiIntersezioneTracciaNonPassante = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaNonPassante.row(0);
                            // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
                            Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTracciaNonPassante);

                            cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
                            Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

                            // Trova il massimo valore e la sua posizione
                            Index maxIndex; // Definisce l'indice massimo usando Index
                            absVec.maxCoeff(&maxIndex);

                            if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && idTracciaNonPassante.size() == 1 && idTracciaPassante.size() == 0){
                                cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                            }
                        }

                        for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
                            // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
                            Vector3d diffCoordinataFratturaPuntoDiIntersezioneTracciaNonPassante = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaNonPassante.row(0);
                            // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
                            Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTracciaNonPassante);

                            cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
                            Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

                            // Trova il massimo valore e la sua posizione
                            Index maxIndex; // Definisce l'indice massimo usando Index
                            absVec.maxCoeff(&maxIndex);

                            if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && idTracciaNonPassante.size() == 1 && idTracciaPassante.size() == 0){
                                cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
                            }
                        }
                    }
                }
            }
        }
        return true;
    }
}

// #include "Utils2poligoni.hpp"
// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <iomanip>

// using namespace std;

// namespace PolygonalLibrary{

// bool letturaMesh(DFN &Fract, PolygonalMesh &mesh){

//     if(!letturaDatiCell0Ds(Fract, mesh)){
//         return false;
//     }

//     return true;
//     }

// bool letturaDatiCell0Ds(DFN &Fract, PolygonalMesh &mesh){

//     unsigned int numeroTracceTotali = 0;

//     vector<unsigned int> idFract1;
//     vector<unsigned int> idFract2;
//     vector<bool> passanteONonPassante1;
//     vector<bool> passanteONonPassante2;
//     vector<double> lunghezzaTracce;

//     idFract1.reserve(Fract.coordinateIntersezioniTracce.size());
//     idFract2.reserve(Fract.coordinateIntersezioniTracce.size());
//     passanteONonPassante1.reserve(Fract.coordinateIntersezioniTracce.size());
//     passanteONonPassante2.reserve(Fract.coordinateIntersezioniTracce.size());
//     lunghezzaTracce.reserve(Fract.coordinateIntersezioniTracce.size());

//     for(const auto &frattura : Fract.coordinateIntersezioniTracce){

//         unsigned int idFratt1 = frattura.first.first;
//         unsigned int idFratt2 = frattura.first.second;
//         bool passONonPass1 = Fract.traccePassantiONonPassanti1[make_pair(idFratt1, idFratt2)];
//         bool passONonPass2 = Fract.traccePassantiONonPassanti2[make_pair(idFratt1, idFratt2)];
//         double lunghezzaTraccia = Fract.lunghezzaTracce[make_pair(idFratt1, idFratt2)];

//         idFract1.push_back(idFratt1); // salvo nel vettore idFract1 associato alla prima frattura l'id della frattura 1
//         idFract2.push_back(idFratt2); // salvo nel vettore idFract2 associato alla seconda frattura l'id della frattura 2
//         passanteONonPassante1.push_back(passONonPass1); // salvo nel vettore passanteONonPassante1 il booleano che dice se la prima frattura è passante o non passante
//         passanteONonPassante2.push_back(passONonPass2); // salvo nel vettore passanteONonPassante2 il booleano che dice se la seconda frattura è passante o non passante
//         lunghezzaTracce.push_back(lunghezzaTraccia); // salvo nel vettore lunghezzaTracce la lunghezza della traccia
//     }

//     numeroTracceTotali = Fract.coordinateIntersezioniTracce.size(); // numero di tracce totali per il file scelto

//     unsigned int idCell0Ds = 0; // identificativo celle 0Ds
//     unsigned int idCell1Ds = 0; // identificativo celle 1Ds
//     unsigned int idCell2Ds = 0; // identificativo celle 2Ds
//     unsigned int Cell0Ds = 0, Cell1Ds = 0, Cell2Ds = 0; // contatori celle 0Ds, 1Ds, 2Ds
//     for(const auto& frattura : Fract.coordinateFratture){

//         unsigned int idFrattura = frattura.first; // estraggo l'id della frattura
//         unsigned int numeroTraccePerFrattura = 0; // contatore numero delle tracce per ciascuna frattura
//         unsigned int numeroTraccePassantiPerFrattura = 0; // contatore numero delle tracce passanti per ciascuna frattura

//         vector<unsigned int> idTraccia; // vettore id di tutte le tracce
//         vector<unsigned int> idTracciaPassante; // vettore id traccia passante
//         vector<unsigned int> idTracciaNonPassante; // vettore id traccia non passante
//         vector<double> lunghezzaTracciaPassante; // vettore lunghezza traccia passante
//         vector<double> lunghezzaTracciaNonPassante; // vettore lunghezza traccia non passante

//         for(unsigned int i = 0;  i < numeroTracceTotali; i++){
//             if(idFrattura == idFract1[i] ||  idFrattura == idFract2[i]){
//                 numeroTraccePerFrattura++;
//                 idTraccia.push_back(i);
//             }
//         }

//         if(numeroTraccePerFrattura > 0){
//             cout << idFrattura << "; " << numeroTraccePerFrattura << endl;
//             for(unsigned int i = 0; i < numeroTraccePerFrattura; i++){

//                 if((idFrattura == idFract1[idTraccia[i]]) && (passanteONonPassante1[idTraccia[i]] == false)){
//                     unsigned int idTr = idTraccia[i];
//                     double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
//                     idTracciaPassante.push_back(idTr); // salvo l'id della traccia perchè è passante
//                     lunghezzaTracciaPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è passante
//                     numeroTraccePassantiPerFrattura++;
//                 }else if((idFrattura == idFract2[idTraccia[i]]) && (passanteONonPassante2[idTraccia[i]] == false)){
//                     unsigned int idTr = idTraccia[i];
//                     double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
//                     idTracciaPassante.push_back(idTr); // salvo l'id della traccia perchè è passante
//                     lunghezzaTracciaPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è passante
//                     numeroTraccePassantiPerFrattura++;
//                 }

//                 if((idFrattura == idFract1[idTraccia[i]]) && (passanteONonPassante1[idTraccia[i]] == true)){
//                     unsigned int idTr = idTraccia[i];
//                     double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
//                     idTracciaNonPassante.push_back(idTr); // salvo l'id della traccia perchè è non passante
//                     lunghezzaTracciaNonPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è non passante
//                 }else if((idFrattura == idFract2[idTraccia[i]]) && (passanteONonPassante2[idTraccia[i]] == true)){
//                     unsigned int idTr = idTraccia[i];
//                     double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
//                     idTracciaNonPassante.push_back(idTr); // salvo l'id della traccia perchè è non passante
//                     lunghezzaTracciaNonPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è non passante
//                 }
//             }

//             if(idTracciaPassante.size() > 1){
//                 for(unsigned int i = 0; i < idTracciaPassante.size() - 1; i++){
//                     for(unsigned int j = i + 1; j < idTracciaPassante.size(); j++){
//                         if(lunghezzaTracciaPassante[i] < lunghezzaTracciaPassante[j]){
//                             swap(idTracciaPassante[i], idTracciaPassante[j]);
//                             swap(lunghezzaTracciaPassante[i], lunghezzaTracciaPassante[j]);
//                         }
//                     }
//                 }
//             }

//             if(idTracciaNonPassante.size() > 1){
//                 for(unsigned int i = 0; i < idTracciaNonPassante.size() - 1; i++){
//                     for(unsigned int j = i + 1; j < idTracciaNonPassante.size(); j++){
//                         if(lunghezzaTracciaNonPassante[i] < lunghezzaTracciaNonPassante[j]){
//                             swap(idTracciaNonPassante[i], idTracciaNonPassante[j]);
//                             swap(lunghezzaTracciaNonPassante[i], lunghezzaTracciaNonPassante[j]);
//                         }
//                     }
//                 }
//             }

//             bool passante = false;
//             MatrixXd matrPuntiDiIntersezionePrecTracciaPassante = MatrixXd::Zero(2,3);
//             Vector3d diffPuntiDiIntersezioneTracciaPassante = Vector3d::Zero();
//             bool puntoDiIntersezione1aSinistra = false;
//             bool puntoDiIntersezione2aSinistra = false;
//             bool puntoDiIntersezione1aDestra = false;
//             bool puntoDiIntersezione2aDestra = false;
//             for(unsigned int t = 0; t < idTracciaPassante.size(); t++){
//                 cout << idTracciaPassante[t] << "; " << boolalpha << passante << endl;

//                 MatrixXd matrPuntiDiIntersezioneTracciaPassante = MatrixXd::Zero(2,3);
//                 Vector3d diffPuntiDiIntersezioneTracciaPass = Vector3d::Zero();
//                 // se la frattura considerata ha una sola traccia passante e nessuna traccia non passante entra in questo if
//                 if((idTracciaPassante.size() > 0 && idTracciaPassante.size() < 2) && (idTracciaNonPassante.size() < 1)){
//                     cout << "if (idTracciaPassante.size() > 0 && idTracciaPassante.size() < 2) && (idTracciaNonPassante.size() < 1)" << endl;
//                     // matrice punti di intersezione precedente
//                     matrPuntiDiIntersezioneTracciaPassante = Fract.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaPassante[t]],idFract2[idTracciaPassante[t]])];
//                     // differenza tra questi punti vett AB = B - A
//                     diffPuntiDiIntersezioneTracciaPass = matrPuntiDiIntersezioneTracciaPassante.row(1) - matrPuntiDiIntersezioneTracciaPassante.row(0);
//                     // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
//                     // Per ogni cella 0D: un identificativo e le coordinate 3D
//                     MatrixXd coordinataFrattura = Fract.coordinateFratture[idFrattura];
//                     unsigned int num01 = 0;
//                     unsigned int num2 = 0;
//                     unsigned int nVertices = 0, nEdges = 0;
//                     // vector<Vector3d> coordFratt;
//                     Vector2i verticiCell1Ds;
//                     vector<unsigned int> vertices, edges;

//                     for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){

//                         // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                         Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaPassante.row(0);
//                         // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                         Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);

//                         cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                         Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

//                         // Trova il massimo valore e la sua posizione
//                         Index maxIndex; // Definisce l'indice massimo usando Index
//                         absVec.maxCoeff(&maxIndex);

//                         if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0){
//                             cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             idCell0Ds++;
//                             idCell1Ds++;
//                             num01++;
//                             Cell0Ds++;
//                             Cell1Ds++;
//                         }
//                     }
//                     mesh.NumberCell0Ds = num01;
//                     mesh.NumberCell1Ds = num01;
//                     cout << "NumberCell0Ds: " << mesh.NumberCell0Ds << endl;
//                     cout << "NumberCell1Ds: " << mesh.NumberCell1Ds << endl;

//                     // mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneTracciaPassante.row(1));
//                     // mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneTracciaPassante.row(0));

//                     mesh.VerticesCell1Ds.reserve(mesh.NumberCell1Ds);
//                     for(unsigned int i = 0; i < mesh.NumberCell1Ds; i++){
//                         verticiCell1Ds(0) = mesh.IdCell0Ds[i];
//                         verticiCell1Ds(1) = mesh.IdCell0Ds[(i+1) % mesh.NumberCell1Ds];
//                         mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                     }

//                     if(mesh.NumberCell0Ds >= 2){
//                         mesh.IdCell2Ds.push_back(idCell2Ds);
//                         idCell2Ds++;
//                         num2++;
//                         Cell2Ds++;
//                         mesh.NumberCell2Ds = num2;
//                     }
//                     cout << "NumberCell2Ds: " << mesh.NumberCell2Ds << endl;

//                     for(unsigned int i = 0; i < mesh.NumberCell0Ds; i++){
//                         nVertices++;
//                         vertices.push_back(mesh.IdCell0Ds[i]);
//                     }
//                     mesh.VerticesCell2Ds.push_back(vertices);

//                     for(unsigned int i = 0; i < mesh.NumberCell1Ds; i++){
//                         nEdges++;
//                         edges.push_back(mesh.IdCell1Ds[i]);
//                     }
//                     mesh.EdgesCell2Ds.push_back(edges);

//                     mesh.NumVertices.push_back(nVertices);
//                     mesh.NumEdges.push_back(nEdges);
//                     // cout << "numVertices: " << mesh.NumVertices << endl;
//                     // cout << "numEdges: " << mesh.NumEdges << endl;

//                     num01 = 0;
//                     num2 = 0;
//                     nVertices = 0;
//                     nEdges = 0;
//                     // mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneTracciaPassante.row(0));
//                     // mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneTracciaPassante.row(1));
//                     for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){

//                         // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                         Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaPassante.row(0);
//                         // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                         Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);

//                         cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                         Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

//                         // Trova il massimo valore e la sua posizione
//                         Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
//                         absVec.maxCoeff(&maxIndex);

//                         if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0){
//                             cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             // coordFratt.push_back(coordinataFrattura.col(i).transpose());
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             idCell0Ds++;
//                             idCell1Ds++;
//                             num01++;
//                             mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                             Cell0Ds++;
//                             Cell1Ds++;
//                             // mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                         }
//                     }
//                     mesh.NumberCell0Ds = num01;
//                     mesh.NumberCell1Ds = num01;
//                     cout << "NumberCell0Ds: " << mesh.NumberCell0Ds << endl;
//                     cout << "NumberCell1Ds: " << mesh.NumberCell1Ds << endl;
//                     // mesh.CoordinatesCell0Ds.push_back(coordFratt[1].transpose());
//                     // mesh.CoordinatesCell0Ds.push_back(coordFratt[0].transpose());

//                     mesh.VerticesCell1Ds.reserve(mesh.NumberCell1Ds);
//                     for(unsigned int i = 0; i < mesh.NumberCell1Ds; i++){
//                         verticiCell1Ds(0) = mesh.IdCell0Ds[i];
//                         verticiCell1Ds(1) = mesh.IdCell0Ds[(i+1) % mesh.NumberCell1Ds];
//                         mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                     }

//                     if(mesh.NumberCell0Ds >= 2){
//                         mesh.IdCell2Ds.push_back(idCell2Ds);
//                         idCell2Ds++;
//                         num2++;
//                         Cell2Ds++;
//                         mesh.NumberCell2Ds = num2;
//                     }
//                     cout << "NumberCell2Ds: " << mesh.NumberCell2Ds << endl;

//                     for(unsigned int i = 0; i < mesh.NumberCell0Ds; i++){
//                         nVertices++;
//                         vertices.push_back(mesh.IdCell0Ds[i]);
//                     }
//                     mesh.VerticesCell2Ds.push_back(vertices);

//                     for(unsigned int i = 0; i < mesh.NumberCell1Ds; i++){
//                         nEdges++;
//                         edges.push_back(mesh.IdCell1Ds[i]);
//                     }
//                     mesh.EdgesCell2Ds.push_back(edges);

//                     mesh.NumVertices.push_back(nVertices);
//                     mesh.NumEdges.push_back(nEdges);

//                     for(const auto& id : mesh.IdCell0Ds){cout << "IdCell0Ds: " << id << endl;}
//                     for(const auto& vec : mesh.CoordinatesCell0Ds){cout << "coordinate Cell0Ds: " << vec.transpose() << endl;}
//                     for(const auto& id : mesh.IdCell1Ds){cout << "IdCell1Ds: " << id << endl;}
//                     for(const auto& vec : mesh.VerticesCell1Ds){cout << "vertici lati Cell1Ds: " << vec.transpose() << endl;}
//                     for(const auto& id : mesh.IdCell2Ds){cout << "IdCell2Ds: " << id << endl;}
//                     // cout << "numVertices: " << mesh.NumVertices << endl;
//                     // cout << "numEdges: " << mesh.NumEdges << endl;
//                 }

//                 puntoDiIntersezione1aSinistra = false;
//                 puntoDiIntersezione2aSinistra = false;
//                 puntoDiIntersezione1aDestra = false;
//                 puntoDiIntersezione2aDestra = false;
//                 matrPuntiDiIntersezioneTracciaPassante = MatrixXd::Zero(2,3);
//                 diffPuntiDiIntersezioneTracciaPass = Vector3d::Zero();
//                 // se la frattura considerata ha più di una traccia passante e 0 < t < idTracciaPassante.size()
//                 // in quanto andiamo a confrontare la posizione dei punti di intersezioni di traccia successiva rispetto alla traccia precedente
//                 // per capire se si trova a destra o a sinistra o se la attraversa
//                 if(idTracciaPassante.size() > 1 && (t > 0) && (t < idTracciaPassante.size())){

//                     // matrice punti di intersezione precedente
//                     MatrixXd matrPuntiDiIntersezionePrec = Fract.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaPassante[t-1]],idFract2[idTracciaPassante[t-1]])];
//                     // differenza tra questi punti vett AB = B - A
//                     Vector3d diffPuntiDiIntersezione = matrPuntiDiIntersezionePrec.row(1) - matrPuntiDiIntersezionePrec.row(0);
//                     // matrice punti di intersezione successiva
//                     MatrixXd matrPuntiDiIntersezioneSucc = Fract.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaPassante[t]],idFract2[idTracciaPassante[t]])];
//                     // differenza primo dei due punti di intersezione vett AP1 = P1 - A
//                     Vector3d diffPrimoPuntoDiIntersezioneSucc1 = matrPuntiDiIntersezioneSucc.row(0) - matrPuntiDiIntersezionePrec.row(0);
//                     // prodotto vettoriale 1 tra AB e AP1
//                     Vector3d prodVett1 = diffPuntiDiIntersezione.cross(diffPrimoPuntoDiIntersezioneSucc1);
//                     // z positiva a sinistra z negativa a destra
//                     // differenza secondo dei due punti di intersezione vett AP2 = P2 - A
//                     Vector3d diffSecondoPuntoDiIntersezioneSucc1 = matrPuntiDiIntersezioneSucc.row(1) - matrPuntiDiIntersezionePrec.row(0);
//                     // prodotto vettoriale 2 tra AB e AP2
//                     Vector3d prodVett2 = diffPuntiDiIntersezione.cross(diffSecondoPuntoDiIntersezioneSucc1);
//                     // z positiva a sinistra z negativa a destra

//                     if(prodVett1(0) > 0 || prodVett1(1) > 0 || prodVett1(2) > 0){
//                         puntoDiIntersezione1aSinistra = true;
//                         cout << "sinistra Passante 1: " << idTracciaPassante[t-1] << " " <<  puntoDiIntersezione1aSinistra << endl;
//                     }else{
//                         puntoDiIntersezione1aDestra = true;
//                         cout << "destra Passante 1: " << idTracciaPassante[t-1] << " " << puntoDiIntersezione1aDestra << endl;
//                     }

//                     if(prodVett2(0) > 0 || prodVett2(1) > 0 || prodVett2(2) > 0){
//                         puntoDiIntersezione2aSinistra = true;
//                         cout << "sinistra Passante 2: " << idTracciaPassante[t-1] << " " << puntoDiIntersezione2aSinistra << endl;
//                     }else{
//                         puntoDiIntersezione2aDestra = true;
//                         cout << "destra Passante 2: " << idTracciaPassante[t-1] << " " << puntoDiIntersezione2aDestra << endl;
//                     }

//                     if(puntoDiIntersezione1aSinistra && puntoDiIntersezione2aSinistra){
//                         cout << "sinistra Passante" << endl;
//                     }else if(puntoDiIntersezione1aDestra && puntoDiIntersezione2aDestra){
//                         cout << "destra Passante" << endl;
//                     }else if(puntoDiIntersezione1aSinistra && puntoDiIntersezione2aDestra){
//                         cout << "sinistra1-destra2 passante" << endl;
//                     }else if(puntoDiIntersezione1aDestra && puntoDiIntersezione2aSinistra){
//                         cout << "destra1-sinistra2 passante" << endl;
//                     }
//                 }
//                 // matrice punti di intersezione precedente
//                 matrPuntiDiIntersezionePrecTracciaPassante = Fract.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaPassante[idTracciaPassante.size() -1]],idFract2[idTracciaPassante[idTracciaPassante.size() -1]])];
//                 // differenza tra questi punti vett AB = B - A
//                 diffPuntiDiIntersezioneTracciaPassante = matrPuntiDiIntersezionePrecTracciaPassante.row(1) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
//             }

//             bool nonPassante = true;
//             for(unsigned int t = 0; t < idTracciaNonPassante.size(); t++){
//                 cout << idTracciaNonPassante[t] << "; " << boolalpha << nonPassante << endl;

//                 puntoDiIntersezione1aSinistra = false;
//                 puntoDiIntersezione2aSinistra = false;
//                 puntoDiIntersezione1aDestra = false;
//                 puntoDiIntersezione2aDestra = false;
//                 // mi serve questo if in quanto arrivo da una traccia precendente passante e devo verificare rispetto a questo dove si trova la prima traccia non passante che incontro
//                 if((t < 1) && (numeroTraccePassantiPerFrattura > 0) && idTracciaNonPassante.size() > 0){
//                     cout << "if (t < 1) && (numeroTraccePassantiPerFrattura > 0) && idTracciaNonPassante.size() > 0" << endl;
//                     MatrixXd matrPuntiDiIntersezioneSuccTracciaNonPassante = Fract.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[0]],idFract2[idTracciaNonPassante[0]])];
//                     // differenza primo dei due punti di intersezione vett AP1 = P1 - A
//                     Vector3d diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
//                     // prodotto vettoriale 1 tra AB e AP1
//                     Vector3d prodVett1 = diffPuntiDiIntersezioneTracciaPassante.cross(diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1);
//                     // z positiva a sinistra z negativa a destra
//                     // differenza secondo dei due punti di intersezione vett AP2 = P2 - A
//                     Vector3d diffSecondoPuntoDiIntersezioneSucc1 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
//                     // prodotto vettoriale 2 tra AB e AP2
//                     Vector3d prodVett2 = diffPuntiDiIntersezioneTracciaPassante.cross(diffSecondoPuntoDiIntersezioneSucc1);
//                     // z positiva a sinistra z negativa a destra

//                     if(prodVett1(0) > 0 || prodVett1(1) > 0 || prodVett1(2) > 0){
//                         puntoDiIntersezione1aSinistra = true;
//                         cout << "sinistra1 t < 1" << " " << puntoDiIntersezione1aSinistra << endl;
//                     }else{
//                         puntoDiIntersezione1aDestra = true;
//                         cout << "destra1 t < 1" << " " << puntoDiIntersezione1aDestra << endl;
//                     }

//                     if(prodVett2(0) > 0 || prodVett2(1) > 0 || prodVett2(2) > 0){
//                         puntoDiIntersezione2aSinistra = true;
//                         cout << "sinistra2 t < 1" << " " <<  puntoDiIntersezione2aSinistra << endl;
//                     }else{
//                         puntoDiIntersezione2aDestra = true;
//                         cout << "destra2 t < 1" << " " << puntoDiIntersezione2aDestra << endl;
//                     }

//                     // Calcola i vettori AB già calcolato è diffPuntiDiIntersezioneTracciaPassante e AP già calcolato diffSecondoPuntoDiIntersezioneSucc1

//                     // Proiezione di AP su AB
//                     double t = diffSecondoPuntoDiIntersezioneSucc1.dot(diffPuntiDiIntersezioneTracciaPassante) / diffPuntiDiIntersezioneTracciaPassante.dot(diffPuntiDiIntersezioneTracciaPassante);

//                     // Punto di proiezione punto di intersezione traccia non passante su traccia passante
//                     Vector3d proiezionePuntoDiInterSuTracciaPass = matrPuntiDiIntersezionePrecTracciaPassante.row(0) + t * diffPuntiDiIntersezioneTracciaPassante.transpose();

//                     if(puntoDiIntersezione1aSinistra && puntoDiIntersezione2aSinistra){

//                         // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
//                         // Per ogni cella 0D: un identificativo e le coordinate 3D
//                         MatrixXd coordinataFrattura = Fract.coordinateFratture[idFrattura];
//                         Vector2i verticiCell1Ds;
//                         unsigned int num01 = 0;
//                         unsigned int num2 = 0;
//                         unsigned int nVertices = 0, nEdges = 0;
//                         vector<unsigned int> vertices, edges;

//                         for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                             // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                             Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
//                             // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                             Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaPassante.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
//                             cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                             Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

//                             // Trova il massimo valore e la sua posizione
//                             Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
//                             absVec.maxCoeff(&maxIndex);

//                             if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0){
//                                 cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                 mesh.IdCell0Ds.push_back(idCell0Ds);
//                                 mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                                 mesh.IdCell1Ds.push_back(idCell1Ds);
//                                 idCell0Ds++;
//                                 idCell1Ds++;
//                                 num01++;
//                                 Cell0Ds++;
//                                 Cell1Ds++;
//                             }
//                         }
//                         mesh.NumberCell0Ds = num01 + matrPuntiDiIntersezionePrecTracciaPassante.rows();
//                         mesh.NumberCell1Ds = num01 + matrPuntiDiIntersezionePrecTracciaPassante.rows();
//                         cout << "NumberCell0Ds: " << mesh.NumberCell0Ds << endl;
//                         cout << "NumberCell1Ds: " << mesh.NumberCell1Ds << endl;

//                         mesh.IdCell0Ds.reserve(2);
//                         mesh.CoordinatesCell0Ds.reserve(2);
//                         for(unsigned int i = 0; i < 2; i++){
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             Cell0Ds++;
//                             Cell1Ds++;
//                             idCell0Ds++;
//                             idCell1Ds++;
//                         }
//                         mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezionePrecTracciaPassante.row(1));
//                         mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezionePrecTracciaPassante.row(0));

//                         mesh.VerticesCell1Ds.reserve(mesh.NumberCell1Ds);
//                         for(unsigned int i = 0; i < mesh.IdCell0Ds.size(); i++){
//                             verticiCell1Ds(0) = mesh.IdCell0Ds[i];
//                             verticiCell1Ds(1) = mesh.IdCell0Ds[(i+1) % mesh.IdCell0Ds.size()];
//                             mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                         }

//                         if(mesh.NumberCell0Ds > 2){
//                             mesh.IdCell2Ds.push_back(idCell2Ds);
//                             idCell2Ds++;
//                             Cell2Ds++;
//                             num2++;
//                             mesh.NumberCell2Ds = num2;
//                         }
//                         cout << "NumberCell2Ds: " << mesh.NumberCell2Ds << endl;

//                         for(unsigned int i = 0; i < mesh.IdCell0Ds.size(); i++){
//                             nVertices++;
//                             vertices.push_back(mesh.IdCell0Ds[i]);
//                         }
//                         mesh.VerticesCell2Ds.push_back(vertices);

//                         for(unsigned int i = 0; i < mesh.IdCell1Ds.size(); i++){
//                             nEdges++;
//                             edges.push_back(mesh.IdCell1Ds[i]);
//                         }
//                         mesh.EdgesCell2Ds.push_back(edges);

//                         mesh.NumVertices.push_back(nVertices);
//                         mesh.NumEdges.push_back(nEdges);

//                         num01 = 0;
//                         num2 = 0;
//                         nVertices = 0;
//                         nEdges = 0;

//                         // Calcolare i vettori
//                         Vector3d AB = proiezionePuntoDiInterSuTracciaPass.transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
//                         Vector3d AP = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1) - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
//                         Vector3d BP = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1) - proiezionePuntoDiInterSuTracciaPass.transpose();

//                         // Verificare che P si trovi tra A e B
//                         // Controlla che P non superi i limiti di A e B per tutte le dimensioni
//                         if(AP.dot(AB) >= 0 && BP.dot(AB) <= 0){
//                             mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0));
//                             mesh.CoordinatesCell0Ds.push_back(proiezionePuntoDiInterSuTracciaPass.transpose());
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             idCell0Ds++;
//                             idCell1Ds++;
//                             Cell0Ds++;
//                             Cell1Ds++;
//                             num01++;
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             idCell0Ds++;
//                             idCell1Ds++;
//                             Cell0Ds++;
//                             Cell1Ds++;
//                             num01++;
//                         }

//                         mesh.NumberCell0Ds = num01;
//                         mesh.NumberCell1Ds = num01;
//                         cout << "NumberCell0Ds: " << mesh.NumberCell0Ds << endl;
//                         cout << "NumberCell1Ds: " << mesh.NumberCell1Ds << endl;

//                         // mesh.VerticesCell1Ds.reserve(mesh.NumberCell1Ds);
//                         // for(unsigned int i = 0; i < mesh.IdCell0Ds.size(); i++){
//                         //     verticiCell1Ds(0) = mesh.IdCell0Ds[i];
//                         //     verticiCell1Ds(1) = mesh.IdCell0Ds[(i+1) % mesh.IdCell0Ds.size()];
//                         //     mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                         // }

//                         if(mesh.NumberCell0Ds > 2){
//                             mesh.IdCell2Ds.push_back(idCell2Ds);
//                             idCell2Ds++;
//                             Cell2Ds++;
//                             num2++;
//                             mesh.NumberCell2Ds = num2;
//                         }
//                         cout << "NumberCell2Ds: " << mesh.NumberCell2Ds << endl;

//                         for(unsigned int i = 0; i < mesh.IdCell0Ds.size(); i++){
//                             nVertices++;
//                             vertices.push_back(mesh.IdCell0Ds[i]);
//                         }
//                         mesh.VerticesCell2Ds.push_back(vertices);

//                         for(unsigned int i = 0; i < mesh.IdCell1Ds.size(); i++){
//                             nEdges++;
//                             edges.push_back(mesh.IdCell1Ds[i]);
//                         }
//                         mesh.EdgesCell2Ds.push_back(edges);

//                         for(const auto& id : mesh.IdCell0Ds){cout << "IdCell0Ds: " << id << endl;}
//                         for(const auto& vec : mesh.CoordinatesCell0Ds){cout << "coordinate Cell0Ds: " << vec.transpose() << endl;}
//                         for(const auto& id : mesh.IdCell1Ds){cout << "IdCell1Ds: " << id << endl;}
//                         for(const auto& vec : mesh.VerticesCell1Ds){cout << "vertici lati Cell1Ds: " << vec.transpose() << endl;}
//                         for(const auto& id : mesh.IdCell2Ds){cout << "IdCell2Ds: " << id << endl;}
//                         // cout << "numVertices: " << mesh.NumVertices << endl;
//                         // cout << "numEdges: " << mesh.NumEdges << endl;

//                         cout << "sinistra t < 1" << endl;

//                     }else if(puntoDiIntersezione1aDestra && puntoDiIntersezione2aDestra){
//                         // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
//                         // Per ogni cella 0D: un identificativo e le coordinate 3D
//                         MatrixXd coordinataFrattura = Fract.coordinateFratture[idFrattura];
//                         unsigned int idCell0Ds = 0; // identificativo celle 0Ds
//                         vector<Vector3d> coordFratt;
//                         mesh.IdCell0Ds.reserve(2);
//                         mesh.CoordinatesCell0Ds.reserve(2);
//                         for(unsigned int i = 0; i < 2; i++){
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezionePrecTracciaPassante.row(idCell0Ds));
//                             idCell0Ds++; // aumento di uno k in quanto i punti di intersezione delle tracce sono sempre 2
//                         }
//                         for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                             // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                             Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
//                             // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                             Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaPassante.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
//                             if(prodVettCoordFrattPunDiInterTraccia(2) > 0){
//                                 cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                 mesh.IdCell0Ds.push_back(idCell0Ds);
//                                 coordFratt.push_back(coordinataFrattura.col(i).transpose());
//                                 idCell0Ds++;
//                             }
//                         }
//                         mesh.NumberCell0Ds = idCell0Ds + matrPuntiDiIntersezionePrecTracciaPassante.rows();
//                         cout << "NumberCell0Ds: " << mesh.NumberCell0Ds << endl;
//                         mesh.CoordinatesCell0Ds.push_back(coordFratt[1].transpose());
//                         mesh.CoordinatesCell0Ds.push_back(coordFratt[0].transpose());
//                         for(const auto& id : mesh.IdCell0Ds){
//                             cout << "IdCell0Ds: " << id << endl;
//                         }
//                         for(const auto& vec : mesh.CoordinatesCell0Ds){
//                             cout << vec.transpose() << endl;
//                         }

//                         cout << "destra t < 1" << endl;

//                     }else if(puntoDiIntersezione1aSinistra && puntoDiIntersezione2aDestra){
//                         cout << "sinistra1-destra2 t < 1" << endl;
//                     }else if(puntoDiIntersezione1aDestra && puntoDiIntersezione2aSinistra){
//                         cout << "destra1-sinistra2 t < 1" << endl;
//                     }
//                 }

//                 puntoDiIntersezione1aSinistra = false;
//                 puntoDiIntersezione2aSinistra = false;
//                 puntoDiIntersezione1aDestra = false;
//                 puntoDiIntersezione2aDestra = false;
//                 if(idTracciaNonPassante.size() > 1 && t > 0 && t < idTracciaNonPassante.size()){
//                     cout << "if idTracciaNonPassante.size() > 1 && t > 0 && t < idTracciaNonPassante.size()" << endl;
//                     // matrice punti di intersezione precedente
//                     MatrixXd matrPuntiDiIntersezionePrec = Fract.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[t-1]],idFract2[idTracciaNonPassante[t-1]])];
//                     // differenza tra questi punti vett AB = B - A
//                     Vector3d diffPuntiDiIntersezione = matrPuntiDiIntersezionePrec.row(1) - matrPuntiDiIntersezionePrec.row(0);
//                     // matrice punti di intersezione successiva
//                     MatrixXd matrPuntiDiIntersezioneSucc = Fract.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[t]],idFract2[idTracciaNonPassante[t]])];
//                     // differenza primo dei due punti di intersezione vett AP1 = P1 - A
//                     Vector3d diffPrimoPuntoDiIntersezioneSucc1 = matrPuntiDiIntersezioneSucc.row(0) - matrPuntiDiIntersezionePrec.row(0);
//                     // prodotto vettoriale 1 tra AB e AP1
//                     Vector3d prodVett1 = diffPuntiDiIntersezione.cross(diffPrimoPuntoDiIntersezioneSucc1);
//                     // z positiva a sinistra z negativa a destra
//                     // differenza secondo dei due punti di intersezione vett AP2 = P2 - A
//                     Vector3d diffSecondoPuntoDiIntersezioneSucc1 = matrPuntiDiIntersezioneSucc.row(1) - matrPuntiDiIntersezionePrec.row(0);
//                     // prodotto vettoriale 2 tra AB e AP2
//                     Vector3d prodVett2 = diffPuntiDiIntersezione.cross(diffSecondoPuntoDiIntersezioneSucc1);
//                     // z positiva a sinistra z negativa a destra

//                     if(prodVett1(0) > 0 || prodVett1(1) > 0 || prodVett1(2) > 0){
//                         puntoDiIntersezione1aSinistra = true;
//                         cout << "sinistra Non Passante 1: " << idTracciaNonPassante[t-1] << " " <<  puntoDiIntersezione1aSinistra << endl;
//                     }else{
//                         puntoDiIntersezione1aDestra = true;
//                         cout << "destra Non Passante 1: " << idTracciaNonPassante[t-1] << " " << puntoDiIntersezione1aDestra << endl;
//                     }

//                     if(prodVett2(0) > 0 || prodVett2(1) > 0 || prodVett2(2) > 0){
//                         puntoDiIntersezione2aSinistra = true;
//                         cout << "sinistra Non Passante 2: " << idTracciaNonPassante[t-1] << " " <<  puntoDiIntersezione2aSinistra << endl;
//                     }else{
//                         puntoDiIntersezione2aDestra = true;
//                         cout << "destra Non Passante 2: " << idTracciaNonPassante[t-1] << " " << puntoDiIntersezione2aDestra << endl;
//                     }

//                     if(puntoDiIntersezione1aSinistra && puntoDiIntersezione2aSinistra){
//                         cout << "sinistra Non Passante" << endl;
//                     }else if(puntoDiIntersezione1aDestra && puntoDiIntersezione2aDestra){
//                         cout << "destra Non Passante" << endl;
//                     }else if(puntoDiIntersezione1aSinistra && puntoDiIntersezione2aDestra){
//                         cout << "sinistra1-destra2 non passante" << endl;
//                     }else if(puntoDiIntersezione1aDestra && puntoDiIntersezione2aSinistra){
//                         cout << "destra1-sinistra2 non passante" << endl;
//                     }
//                 }
//             }
//         }
//     }
//     cout << "Cell0Ds: " << Cell0Ds << endl;
//     cout << "Cell1Ds: " << Cell1Ds << endl;
//     cout << "Cell2Ds: " << Cell2Ds << endl;

//     return true;
//     }
// }
