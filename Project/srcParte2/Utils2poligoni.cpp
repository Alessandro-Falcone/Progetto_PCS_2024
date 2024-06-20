// #include "Utils2poligoni.hpp"
// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <iomanip>

// using namespace DFNLibrary;
// using namespace std;

// namespace PolygonalLibrary{

//     bool letturaMesh(Frattura &Fratt, Traccia &Trac, PolygonalMesh &mesh){

//         unsigned int numeroTracceTotali = 0;

//         vector<unsigned int> idFract1;
//         vector<unsigned int> idFract2;
//         vector<bool> passanteONonPassante1;
//         vector<bool> passanteONonPassante2;
//         vector<double> lunghezzaTracce;

//         idFract1.reserve(Trac.coordinateIntersezioniTracce.size());
//         idFract2.reserve(Trac.coordinateIntersezioniTracce.size());
//         passanteONonPassante1.reserve(Trac.coordinateIntersezioniTracce.size());
//         passanteONonPassante2.reserve(Trac.coordinateIntersezioniTracce.size());
//         lunghezzaTracce.reserve(Trac.coordinateIntersezioniTracce.size());

//         for(const auto &frattura : Trac.coordinateIntersezioniTracce){

//             unsigned int idFratt1 = frattura.first.first;
//             unsigned int idFratt2 = frattura.first.second;
//             bool passONonPass1 = Trac.traccePassantiONonPassanti1[make_pair(idFratt1, idFratt2)];
//             bool passONonPass2 = Trac.traccePassantiONonPassanti2[make_pair(idFratt1, idFratt2)];
//             double lunghezzaTraccia = Trac.lunghezzaTracce[make_pair(idFratt1, idFratt2)];

//             idFract1.push_back(idFratt1); // salvo nel vettore idFract1 associato alla prima frattura l'id della frattura 1
//             idFract2.push_back(idFratt2); // salvo nel vettore idFract2 associato alla seconda frattura l'id della frattura 2
//             passanteONonPassante1.push_back(passONonPass1); // salvo nel vettore passanteONonPassante1 il booleano che dice se la prima frattura è passante o non passante
//             passanteONonPassante2.push_back(passONonPass2); // salvo nel vettore passanteONonPassante2 il booleano che dice se la seconda frattura è passante o non passante
//             lunghezzaTracce.push_back(lunghezzaTraccia); // salvo nel vettore lunghezzaTracce la lunghezza della traccia
//         }

//         numeroTracceTotali = Trac.coordinateIntersezioniTracce.size(); // numero di tracce totali per il file scelto

//         unsigned int idCell0Ds = 0, idCell1Ds = 0, idCell2Ds = 0; // identificativo celle 0Ds, 1Ds, 2Ds
//         unsigned int Cell0Ds = 0, Cell1Ds = 0, Cell2Ds = 0; // contatori celle 0Ds, 1Ds, 2Ds
//         unsigned int contatoreTraccia = 0;
//         for(unsigned int idFrattura = 0; idFrattura < Fratt.coordinateFratture.size(); idFrattura++){

//             unsigned int numeroTraccePerFrattura = 0; // contatore numero delle tracce per ciascuna frattura
//             unsigned int numeroTraccePassantiPerFrattura = 0; // contatore numero delle tracce passanti per ciascuna frattura

//             vector<unsigned int> idTraccia; // vettore id di tutte le tracce
//             vector<unsigned int> idTracciaPassante; // vettore id traccia passante
//             vector<unsigned int> idTracciaNonPassante; // vettore id traccia non passante
//             vector<double> lunghezzaTracciaPassante; // vettore lunghezza traccia passante
//             vector<double> lunghezzaTracciaNonPassante; // vettore lunghezza traccia non passante

//             for(unsigned int i = 0;  i < numeroTracceTotali; i++){
//                 if(idFrattura == idFract1[i] ||  idFrattura == idFract2[i]){
//                     numeroTraccePerFrattura++;
//                     idTraccia.push_back(i);
//                 }
//             }

//             if(numeroTraccePerFrattura > 0){
//                 cout << idFrattura << "; " << numeroTraccePerFrattura << endl;
//                 for(unsigned int i = 0; i < numeroTraccePerFrattura; i++){

//                     if((idFrattura == idFract1[idTraccia[i]]) && (passanteONonPassante1[idTraccia[i]] == false)){
//                         unsigned int idTr = idTraccia[i];
//                         double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
//                         idTracciaPassante.push_back(idTr); // salvo l'id della traccia perchè è passante
//                         lunghezzaTracciaPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è passante
//                         numeroTraccePassantiPerFrattura++;
//                     }else if((idFrattura == idFract2[idTraccia[i]]) && (passanteONonPassante2[idTraccia[i]] == false)){
//                         unsigned int idTr = idTraccia[i];
//                         double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
//                         idTracciaPassante.push_back(idTr); // salvo l'id della traccia perchè è passante
//                         lunghezzaTracciaPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è passante
//                         numeroTraccePassantiPerFrattura++;
//                     }

//                     if((idFrattura == idFract1[idTraccia[i]]) && (passanteONonPassante1[idTraccia[i]] == true)){
//                         unsigned int idTr = idTraccia[i];
//                         double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
//                         idTracciaNonPassante.push_back(idTr); // salvo l'id della traccia perchè è non passante
//                         lunghezzaTracciaNonPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è non passante
//                     }else if((idFrattura == idFract2[idTraccia[i]]) && (passanteONonPassante2[idTraccia[i]] == true)){
//                         unsigned int idTr = idTraccia[i];
//                         double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
//                         idTracciaNonPassante.push_back(idTr); // salvo l'id della traccia perchè è non passante
//                         lunghezzaTracciaNonPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è non passante
//                     }
//                 }

//                 if(idTracciaPassante.size() > 1){
//                     for(unsigned int i = 0; i < idTracciaPassante.size() - 1; i++){
//                         for(unsigned int j = i + 1; j < idTracciaPassante.size(); j++){
//                             if(lunghezzaTracciaPassante[i] < lunghezzaTracciaPassante[j]){
//                                 swap(idTracciaPassante[i], idTracciaPassante[j]);
//                                 swap(lunghezzaTracciaPassante[i], lunghezzaTracciaPassante[j]);
//                             }
//                         }
//                     }
//                 }

//                 if(idTracciaNonPassante.size() > 1){
//                     for(unsigned int i = 0; i < idTracciaNonPassante.size() - 1; i++){
//                         for(unsigned int j = i + 1; j < idTracciaNonPassante.size(); j++){
//                             if(lunghezzaTracciaNonPassante[i] < lunghezzaTracciaNonPassante[j]){
//                                 swap(idTracciaNonPassante[i], idTracciaNonPassante[j]);
//                                 swap(lunghezzaTracciaNonPassante[i], lunghezzaTracciaNonPassante[j]);
//                             }
//                         }
//                     }
//                 }

//                 bool passante = false;
//                 MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];
//                 MatrixXd matrPuntiDiIntersezionePrecTracciaPassante = MatrixXd::Zero(2,3);
//                 Vector3d diffPuntiDiIntersezioneTracciaPassante = Vector3d::Zero();
//                 MatrixXd matrPuntiDiIntersezioneSuccTracciaNonPassante = MatrixXd::Zero(2,3);

//                 if(idTracciaNonPassante.size() > 0){
//                     matrPuntiDiIntersezioneSuccTracciaNonPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[0]],idFract2[idTracciaNonPassante[0]])];
//                 }

//                 for(unsigned int t = 0; t < idTracciaPassante.size(); t++){
//                     cout << idTracciaPassante[t] << "; " << boolalpha << passante << endl;

//                     unsigned int vertici = 0; // contatore vertici che definiscono un lato
//                     Vector2i verticiCell1Ds;
//                     Vector2i vertCell1Ds;
//                     vector<unsigned int> vertices;

//                     MatrixXd matrPuntiDiIntersezioneTracciaPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaPassante[t]],idFract2[idTracciaPassante[t]])];
//                     contatoreTraccia++;

//                     if(contatoreTraccia < 2){
//                         mesh.CoordinatesCell0Ds.reserve(matrPuntiDiIntersezioneTracciaPassante.rows());

//                         // cout << "punti di intersezione:" << endl;
//                         // cout << matrPuntiDiIntersezioneTracciaPassante << endl;
//                         mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneTracciaPassante.row(0));
//                         mesh.IdCell0Ds.push_back(idCell0Ds);
//                         mesh.IdCell1Ds.push_back(idCell1Ds);
//                         verticiCell1Ds(vertici) = idCell0Ds;

//                         // vertices.push_back(idCell0Ds);
//                         // cout << "vertice: " << idCell0Ds << endl;

//                         idCell0Ds++;
//                         idCell1Ds++;
//                         Cell0Ds++;
//                         Cell1Ds++;
//                         vertici++;

//                         mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneTracciaPassante.row(1));
//                         mesh.IdCell0Ds.push_back(idCell0Ds);
//                         mesh.IdCell1Ds.push_back(idCell1Ds);
//                         verticiCell1Ds(vertici) = idCell0Ds;

//                         // vertices.push_back(idCell0Ds);
//                         // cout << "vertice: " << idCell0Ds << endl;

//                         idCell0Ds++;
//                         idCell1Ds++;
//                         Cell0Ds++;
//                         Cell1Ds++;
//                         vertici++;

//                         idCell1Ds++;
//                         Cell1Ds++;
//                     }

//                     if(vertici == 2){
//                         mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                         mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                         // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                     }

//                     Vector3d diffPuntiDiIntersezioneTracciaPass = matrPuntiDiIntersezioneTracciaPassante.row(1) - matrPuntiDiIntersezioneTracciaPassante.row(0);

//                     Vector3d diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0) - matrPuntiDiIntersezioneTracciaPassante.row(0);
//                     // prodotto vettoriale 1 tra AB e AP1
//                     Vector3d prodVett1 = diffPuntiDiIntersezioneTracciaPass.cross(diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1);
//                     // z positiva a sinistra z negativa a destra
//                     // differenza secondo dei due punti di intersezione vett AP2 = P2 - A
//                     Vector3d diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1) - matrPuntiDiIntersezioneTracciaPassante.row(0);
//                     // prodotto vettoriale 2 tra AB e AP2
//                     Vector3d prodVett2 = diffPuntiDiIntersezioneTracciaPass.cross(diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2);

//                     Vector3d absVec1 = prodVett1.cwiseAbs();

//                     // Trova il massimo valore e la sua posizione
//                     Index maxIndex1; // Definisce l'indice massimo usando Eigen::Index
//                     absVec1.maxCoeff(&maxIndex1);

//                     Vector3d absVec2 = prodVett2.cwiseAbs();

//                     // Trova il massimo valore e la sua posizione
//                     Index maxIndex2; // Definisce l'indice massimo usando Eigen::Index
//                     absVec2.maxCoeff(&maxIndex2);

//                     vertici = 0;
//                     for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                         // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                         Vector3d diffCoordinataFratturaPuntoDiIntersezioneTracciaPassante = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaPassante.row(0);
//                         // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                         Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTracciaPassante);

//                         // cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                         Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

//                         // Trova il massimo valore e la sua posizione
//                         Index maxIndex; // Definisce l'indice massimo usando Index
//                         absVec.maxCoeff(&maxIndex);

//                         if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && prodVett1(maxIndex1) > 0 && prodVett2(maxIndex2) > 0 && idTracciaNonPassante.size() != 0 && vertici < 2){

//                             // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                             verticiCell1Ds(0) = idCell0Ds;
//                             verticiCell1Ds(1) = vertici;
//                             if(vertici < 1){
//                                 mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                 mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                 // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                             }
//                             vertCell1Ds(vertici) = idCell0Ds;

//                             // vertices.push_back(idCell0Ds);
//                             // cout << "vertice: " << idCell0Ds << endl;

//                             idCell0Ds++;
//                             idCell1Ds++;
//                             Cell0Ds++;
//                             Cell1Ds++;
//                             vertici++;

//                             // idCell1Ds++;
//                             // Cell1Ds++;
//                         }else if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && idTracciaNonPassante.size() == 0 && idTracciaPassante.size() == 1 && vertici < 2){

//                             // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                             verticiCell1Ds(0) = idCell0Ds;
//                             verticiCell1Ds(1) = vertici;
//                             if(vertici < 1){
//                                 mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                 mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                 // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                             }
//                             vertCell1Ds(vertici) = idCell0Ds;

//                             // vertices.push_back(idCell0Ds);
//                             // cout << "vertice: " << idCell0Ds << endl;

//                             idCell0Ds++;
//                             idCell1Ds++;
//                             Cell0Ds++;
//                             Cell1Ds++;
//                             vertici++;

//                             // idCell1Ds++;
//                             // Cell1Ds++;
//                         }
//                     }

//                     if(vertici == 2){
//                         mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                         mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                         // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                         mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                         mesh.VerticesCell1Ds.push_back(vertCell1Ds);
//                         // cout << "vertice id: " << vertCell1Ds.transpose() << endl;
//                     }

//                     // mesh.VerticesCell2Ds.push_back(vertices);
//                     // for (const auto& innerVector : mesh.VerticesCell2Ds) {
//                     //     for (unsigned int value : innerVector) {
//                     //         cout << value << " ";
//                     //     }
//                     //     cout << endl;
//                     // }

//                     vertici = 0;
//                     for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                         // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                         Vector3d diffCoordinataFratturaPuntoDiIntersezioneTracciaPassante = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaPassante.row(0);
//                         // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                         Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTracciaPassante);

//                         // cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                         Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

//                         // Trova il massimo valore e la sua posizione
//                         Index maxIndex; // Definisce l'indice massimo usando Index
//                         absVec.maxCoeff(&maxIndex);

//                         if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && prodVett1(maxIndex1) < 0 && prodVett2(maxIndex2) < 0 && idTracciaNonPassante.size() != 0 && vertici < 2){

//                             // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                             verticiCell1Ds(0) = idCell0Ds;
//                             verticiCell1Ds(1) = vertici;
//                             if(vertici < 1){
//                                 mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                 mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                 // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                             }
//                             vertCell1Ds(vertici) = idCell0Ds;
//                             // vertices.push_back(idCell0Ds);
//                             // cout << "vertice: " << idCell0Ds << endl;
//                             idCell0Ds++;
//                             idCell1Ds++;
//                             Cell0Ds++;
//                             Cell1Ds++;
//                             vertici++;
//                         }else if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && idTracciaNonPassante.size() == 0 && idTracciaPassante.size() == 1 && vertici < 2){

//                             // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                             mesh.IdCell0Ds.push_back(idCell0Ds);
//                             mesh.IdCell1Ds.push_back(idCell1Ds);
//                             mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                             verticiCell1Ds(0) = idCell0Ds;
//                             verticiCell1Ds(1) = vertici;
//                             if(vertici < 1){
//                                 mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                 mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                 // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                             }
//                             vertCell1Ds(vertici) = idCell0Ds;
//                             vertices.push_back(idCell0Ds);
//                             // cout << "vertice: " << idCell0Ds << endl;
//                             idCell0Ds++;
//                             idCell1Ds++;
//                             Cell0Ds++;
//                             Cell1Ds++;
//                             vertici++;
//                         }
//                     }

//                     idCell1Ds++;
//                     Cell1Ds++;
//                     if(vertici == 2){
//                         mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                         mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                         // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                         mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                         mesh.VerticesCell1Ds.push_back(vertCell1Ds);
//                         // cout << "vertice id: " << vertCell1Ds.transpose() << endl;
//                     }

//                     // matrice punti di intersezione precedente
//                     matrPuntiDiIntersezionePrecTracciaPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaPassante[idTracciaPassante.size() -1]],idFract2[idTracciaPassante[idTracciaPassante.size() -1]])];
//                     // differenza tra questi punti vett AB = B - A
//                     diffPuntiDiIntersezioneTracciaPassante = matrPuntiDiIntersezionePrecTracciaPassante.row(1) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
//                 }

//                 bool nonPassante = true;
//                 for(unsigned int t = 0; t < idTracciaNonPassante.size(); t++){
//                     cout << idTracciaNonPassante[t] << "; " << boolalpha << nonPassante << endl;

//                     unsigned int vertici = 0;
//                     Vector2i verticiCell1Ds;
//                     Vector2i vertCell1Ds;
//                     Vector2i vCell1Ds;
//                     Vector2i verticesCell1Ds;
//                     vector<unsigned int> vertices;

//                     if((t < 1) && (numeroTraccePassantiPerFrattura > 0) && idTracciaNonPassante.size() > 0){
//                         idCell1Ds++;
//                         Cell1Ds++;
//                         idCell1Ds++;
//                         Cell1Ds++;
//                         MatrixXd matrPuntiDiIntersezioneSuccTracciaNonPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[0]],idFract2[idTracciaNonPassante[0]])];
//                         mesh.IdCell0Ds.push_back(idCell0Ds);
//                         mesh.IdCell1Ds.push_back(idCell1Ds);
//                         mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0));
//                         verticiCell1Ds(vertici) = idCell0Ds;
//                         // cout << "id cell0Ds: " << idCell0Ds << endl;
//                         vertCell1Ds(0) = idCell0Ds;
//                         vCell1Ds(0) = idCell0Ds;
//                         vertices.push_back(idCell0Ds);
//                         // cout << "vertice: " << idCell0Ds << endl;
//                         idCell0Ds++;
//                         idCell1Ds++;
//                         Cell0Ds++;
//                         Cell1Ds++;
//                         vertici++;
//                         // cout << "punto di intersezione 0 traccia non passante: " << matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0) << endl;

//                         // differenza primo dei due punti di intersezione vett AP1 = P1 - A
//                         Vector3d diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
//                         // prodotto vettoriale 1 tra AB e AP1
//                         Vector3d prodVett1 = diffPuntiDiIntersezioneTracciaPassante.cross(diffPrimoPuntoDiIntersezioneSuccTracciaNonPassante1);
//                         // z positiva a sinistra z negativa a destra
//                         // differenza secondo dei due punti di intersezione vett AP2 = P2 - A
//                         Vector3d diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2 = matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1) - matrPuntiDiIntersezionePrecTracciaPassante.row(0);
//                         // prodotto vettoriale 2 tra AB e AP2
//                         Vector3d prodVett2 = diffPuntiDiIntersezioneTracciaPassante.cross(diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2);
//                         // z positiva a sinistra z negativa a destra

//                         Vector3d absVec1 = prodVett1.cwiseAbs();

//                         // Trova il massimo valore e la sua posizione
//                         Index maxIndex1; // Definisce l'indice massimo usando Eigen::Index
//                         absVec1.maxCoeff(&maxIndex1);

//                         Vector3d absVec2 = prodVett2.cwiseAbs();

//                         // Trova il massimo valore e la sua posizione
//                         Index maxIndex2; // Definisce l'indice massimo usando Eigen::Index
//                         absVec2.maxCoeff(&maxIndex2);
//                         // Calcola i vettori AB già calcolato è diffPuntiDiIntersezioneTracciaPassante e AP già calcolato diffSecondoPuntoDiIntersezioneSucc1

//                         // Proiezione di AP su AB
//                         double t = diffSecondoPuntoDiIntersezioneSuccTracciaNonPassante2.dot(diffPuntiDiIntersezioneTracciaPassante) / diffPuntiDiIntersezioneTracciaPassante.dot(diffPuntiDiIntersezioneTracciaPassante);

//                         // Punto di proiezione punto di intersezione traccia non passante su traccia passante
//                         Vector3d proiezionePuntoDiInterSuTracciaPass = matrPuntiDiIntersezionePrecTracciaPassante.row(0) + t * diffPuntiDiIntersezioneTracciaPassante.transpose();
//                         mesh.IdCell0Ds.push_back(idCell0Ds);
//                         mesh.IdCell1Ds.push_back(idCell1Ds);
//                         mesh.CoordinatesCell0Ds.push_back(proiezionePuntoDiInterSuTracciaPass.transpose());
//                         verticiCell1Ds(vertici) = idCell0Ds;
//                         vertici++;
//                         if(vertici == 2){
//                             mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                             mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                             // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                         }
//                         verticiCell1Ds(0) = idCell0Ds;
//                         verticiCell1Ds(1) = 0;
//                         vertices.push_back(idCell0Ds);
//                         // cout << "vertice: " << idCell0Ds << endl;
//                         mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                         mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                         // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;

//                         verticiCell1Ds(0) = idCell0Ds;
//                         verticiCell1Ds(1) = 1;
//                         mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                         mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                         // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;

//                         // cout << "vertCell1Ds 0: " << vertCell1Ds(0) << endl;
//                         idCell0Ds++;
//                         idCell1Ds++;
//                         Cell0Ds++;
//                         Cell1Ds++;
//                         // cout << "proiezione passante: " << proiezionePuntoDiInterSuTracciaPass.transpose() << endl;

//                         Vector3d vettProiezionePuntoDiIntersezioneTracciaNonPass = proiezionePuntoDiInterSuTracciaPass.transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);

//                         if(prodVett1(maxIndex1) > 0 && prodVett2(maxIndex2) > 0){

//                             // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
//                             // Per ogni cella 0D: un identificativo e le coordinate 3D
//                             MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];

//                             for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                                 // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                                 Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
//                                 // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                                 Vector3d prodVettCoordFrattPunDiInterTraccia = vettProiezionePuntoDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
//                                 // cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                                 Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

//                                 // Trova il massimo valore e la sua posizione
//                                 Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
//                                 absVec.maxCoeff(&maxIndex);

//                                 if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(0).squaredNorm()){
//                                     // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                     mesh.IdCell0Ds.push_back(idCell0Ds);
//                                     mesh.IdCell1Ds.push_back(idCell1Ds);
//                                     mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());

//                                     // vertices.push_back(idCell0Ds);
//                                     // cout << "vertice: " << idCell0Ds << endl;

//                                     vertCell1Ds(1) = idCell0Ds;
//                                     // cout << "vertCell1Ds 1: " << vertCell1Ds(1) << endl;
//                                     mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                     mesh.VerticesCell1Ds.push_back(vertCell1Ds);

//                                     verticiCell1Ds(0) = idCell0Ds;
//                                     verticiCell1Ds(1) = 0;
//                                     mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                     mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                     // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                                     idCell0Ds++;
//                                     idCell1Ds++;
//                                     Cell0Ds++;
//                                     Cell1Ds++;
//                                 }

//                                 if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(1).squaredNorm()){
//                                     // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                     mesh.IdCell0Ds.push_back(idCell0Ds);
//                                     mesh.IdCell1Ds.push_back(idCell1Ds);
//                                     mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());

//                                     // vertices.push_back(idCell0Ds);
//                                     // cout << "vertice: " << idCell0Ds << endl;

//                                     vCell1Ds(1) = idCell0Ds;
//                                     // cout << "vertCell1Ds 1: " << vCell1Ds(1) << endl;
//                                     mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                     mesh.VerticesCell1Ds.push_back(vCell1Ds);

//                                     verticiCell1Ds(0) = idCell0Ds;
//                                     verticiCell1Ds(1) = 1;
//                                     mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                     mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                     // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                                     idCell0Ds++;
//                                     idCell1Ds++;
//                                     Cell0Ds++;
//                                     Cell1Ds++;
//                                 }
//                             }
//                         }

//                         vertici = 0;
//                         if(prodVett1(maxIndex1) < 0 && prodVett2(maxIndex2) < 0){

//                             // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
//                             // Per ogni cella 0D: un identificativo e le coordinate 3D
//                             MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];

//                             for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                                 // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                                 Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
//                                 // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                                 Vector3d prodVettCoordFrattPunDiInterTraccia = vettProiezionePuntoDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
//                                 // cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                                 Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

//                                 // Trova il massimo valore e la sua posizione
//                                 Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
//                                 absVec.maxCoeff(&maxIndex);

//                                 if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(0).squaredNorm()){
//                                     // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                     mesh.IdCell0Ds.push_back(idCell0Ds);
//                                     mesh.IdCell1Ds.push_back(idCell1Ds);
//                                     mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());

//                                     // vertices.push_back(idCell0Ds);
//                                     // cout << "vertice: " << idCell0Ds << endl;

//                                     verticiCell1Ds(vertici) = idCell0Ds;
//                                     idCell0Ds++;
//                                     idCell1Ds++;
//                                     Cell0Ds++;
//                                     Cell1Ds++;
//                                     vertici++;
//                                 }

//                                 if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(1).squaredNorm()){
//                                     // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                     mesh.IdCell0Ds.push_back(idCell0Ds);
//                                     mesh.IdCell1Ds.push_back(idCell1Ds);
//                                     mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                                     verticiCell1Ds(vertici) = idCell0Ds;

//                                     // vertices.push_back(idCell0Ds);
//                                     // cout << "vertice: " << idCell0Ds << endl;

//                                     idCell0Ds++;
//                                     idCell1Ds++;
//                                     Cell0Ds++;
//                                     Cell1Ds++;
//                                     vertici++;
//                                 }
//                             }

//                             if(vertici == 2){
//                                 mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                 mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                 // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                             }

//                         }

//                         vertici = 0;
//                         if(prodVett1(maxIndex1) > 0 && prodVett2(maxIndex2) < 0){

//                             // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
//                             // Per ogni cella 0D: un identificativo e le coordinate 3D
//                             MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];

//                             for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                                 // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                                 Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
//                                 // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                                 Vector3d prodVettCoordFrattPunDiInterTraccia = vettProiezionePuntoDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
//                                 // cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                                 Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

//                                 // Trova il massimo valore e la sua posizione
//                                 Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
//                                 absVec.maxCoeff(&maxIndex);

//                                 if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(0).squaredNorm()){
//                                     // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                     mesh.IdCell0Ds.push_back(idCell0Ds);
//                                     mesh.IdCell1Ds.push_back(idCell1Ds);
//                                     mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                                     verticiCell1Ds(vertici) = idCell0Ds;

//                                     // vertices.push_back(idCell0Ds);
//                                     // cout << "vertice: " << idCell0Ds << endl;

//                                     idCell0Ds++;
//                                     idCell1Ds++;
//                                     Cell0Ds++;
//                                     Cell1Ds++;
//                                     vertici++;
//                                 }

//                                 if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(1).squaredNorm()){
//                                     // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                     mesh.IdCell0Ds.push_back(idCell0Ds);
//                                     mesh.IdCell1Ds.push_back(idCell1Ds);
//                                     mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                                     verticiCell1Ds(vertici) = idCell0Ds;

//                                     // vertices.push_back(idCell0Ds);
//                                     // cout << "vertice: " << idCell0Ds << endl;

//                                     idCell0Ds++;
//                                     idCell1Ds++;
//                                     Cell0Ds++;
//                                     Cell1Ds++;
//                                     vertici++;
//                                 }
//                             }

//                             if(vertici == 2){
//                                 mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                 mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                 // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                             }

//                         }

//                         vertici = 0;

//                         if(prodVett1(maxIndex1) < 0 && prodVett2(maxIndex2) > 0){

//                             // prendo la matrice delle coordinate della frattura associata all'idFrattura adesso devo salvarmi i vertici
//                             // Per ogni cella 0D: un identificativo e le coordinate 3D
//                             MatrixXd coordinataFrattura = Fratt.coordinateFratture[idFrattura];

//                             for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                                 // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                                 Vector3d diffCoordinataFratturaPuntoDiIntersezioneTraccia = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneSuccTracciaNonPassante.row(0);
//                                 // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                                 Vector3d prodVettCoordFrattPunDiInterTraccia = vettProiezionePuntoDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTraccia);
//                                 // cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                                 Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs();

//                                 // Trova il massimo valore e la sua posizione
//                                 Index maxIndex; // Definisce l'indice massimo usando Eigen::Index
//                                 absVec.maxCoeff(&maxIndex);

//                                 if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(0).squaredNorm()){
//                                     // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                     mesh.IdCell0Ds.push_back(idCell0Ds);
//                                     mesh.IdCell1Ds.push_back(idCell1Ds);
//                                     mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                                     verticiCell1Ds(vertici) = idCell0Ds;

//                                     // vertices.push_back(idCell0Ds);
//                                     // cout << "vertice: " << idCell0Ds << endl;

//                                     idCell0Ds++;
//                                     idCell1Ds++;
//                                     Cell0Ds++;
//                                     Cell1Ds++;
//                                     vertici++;
//                                 }

//                                 if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && coordinataFrattura.col(i).transpose().squaredNorm() < matrPuntiDiIntersezionePrecTracciaPassante.row(1).squaredNorm()){
//                                     // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                     mesh.IdCell0Ds.push_back(idCell0Ds);
//                                     mesh.IdCell1Ds.push_back(idCell1Ds);
//                                     mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                                     verticiCell1Ds(vertici) = idCell0Ds;

//                                     // vertices.push_back(idCell0Ds);
//                                     // cout << "vertice: " << idCell0Ds << endl;

//                                     idCell0Ds++;
//                                     idCell1Ds++;
//                                     Cell0Ds++;
//                                     Cell1Ds++;
//                                     vertici++;
//                                 }
//                             }

//                             if(vertici == 2){
//                                 mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                 mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                                 // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                             }

//                         }
//                     }

//                     vertici = 0;
//                     if(idTracciaNonPassante.size() > 0 && idTracciaPassante.size() == 0){

//                         idCell1Ds++;
//                         Cell1Ds++;
//                         MatrixXd matrPuntiDiIntersezioneTracciaNonPassante = Trac.coordinateIntersezioniTracce[make_pair(idFract1[idTracciaNonPassante[t]],idFract2[idTracciaNonPassante[t]])];
//                         Vector3d diffPuntiDiIntersezioneTracciaNonPass = matrPuntiDiIntersezioneTracciaNonPassante.row(1) - matrPuntiDiIntersezioneTracciaNonPassante.row(0);

//                         for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){

//                             Vector3d diffPuntoDiIntersezioneCoordinataFrattura = matrPuntiDiIntersezioneTracciaNonPassante.row(1) - coordinataFrattura.col(i).transpose();
//                             Vector3d diffCoordinateFratture = coordinataFrattura.col((i+1) % coordinataFrattura.cols()) - coordinataFrattura.col(i);
//                             double num = diffPuntoDiIntersezioneCoordinataFrattura.dot(diffCoordinateFratture);
//                             double den = diffCoordinateFratture.dot(diffCoordinateFratture);
//                             double t = num / den;
//                             Vector3d proiezionePuntoDiIntersezioneSuDiffCoordinateFratture = coordinataFrattura.col(i).transpose() + t * diffCoordinateFratture.transpose();
//                             proiezionePuntoDiIntersezioneSuDiffCoordinateFratture = (proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.array().abs() < 1e-15).select(0, proiezionePuntoDiIntersezioneSuDiffCoordinateFratture);

//                             if((proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() != coordinataFrattura.col((i+1) % coordinataFrattura.cols()).transpose() && proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() != coordinataFrattura.col(i).transpose()) && proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() != matrPuntiDiIntersezioneTracciaNonPassante.row(1)){
//                                 // cout << "punto di intersezione 1 traccia non passante: " << matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1) << endl;
//                                 mesh.IdCell0Ds.push_back(idCell0Ds);
//                                 mesh.IdCell1Ds.push_back(idCell1Ds);
//                                 mesh.CoordinatesCell0Ds.push_back(matrPuntiDiIntersezioneSuccTracciaNonPassante.row(1));
//                                 verticiCell1Ds(vertici) = idCell0Ds;
//                                 vertCell1Ds(1) = idCell0Ds;
//                                 vCell1Ds(0) = idCell0Ds;

//                                 // vertices.push_back(idCell0Ds);
//                                 // cout << "vertice: " << idCell0Ds << endl;

//                                 idCell0Ds++;
//                                 idCell1Ds++;
//                                 Cell0Ds++;
//                                 Cell1Ds++;
//                                 vertici++;
//                                 // cout << "proiezione non passante: " << proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose() << endl;
//                                 // vertices.push_back(idCell0Ds);
//                                 // cout << "vertice: " << idCell0Ds << endl;

//                                 mesh.IdCell0Ds.push_back(idCell0Ds);
//                                 mesh.IdCell1Ds.push_back(idCell1Ds);
//                                 mesh.CoordinatesCell0Ds.push_back(proiezionePuntoDiIntersezioneSuDiffCoordinateFratture.transpose());
//                                 verticesCell1Ds(1) = idCell0Ds;
//                                 verticiCell1Ds(vertici) = idCell0Ds;
//                                 idCell0Ds++;
//                                 idCell1Ds++;
//                                 vertCell1Ds(0) = idCell0Ds;
//                                 mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                 mesh.VerticesCell1Ds.push_back(vertCell1Ds);
//                                 Cell0Ds++;
//                                 Cell1Ds++;
//                                 vertici++;
//                             }
//                         }

//                         if(vertici == 2){
//                             mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                             mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                             // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                         }

//                         vertici = 0;
//                         for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                             // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                             Vector3d diffCoordinataFratturaPuntoDiIntersezioneTracciaNonPassante = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaNonPassante.row(0);
//                             // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                             Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTracciaNonPassante);

//                             // cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                             Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

//                             // Trova il massimo valore e la sua posizione
//                             Index maxIndex; // Definisce l'indice massimo usando Index
//                             absVec.maxCoeff(&maxIndex);

//                             if(prodVettCoordFrattPunDiInterTraccia(maxIndex) < 0 && idTracciaNonPassante.size() == 1 && idTracciaPassante.size() == 0){
//                                 // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                 mesh.IdCell0Ds.push_back(idCell0Ds);
//                                 mesh.IdCell1Ds.push_back(idCell1Ds);
//                                 mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                                 verticiCell1Ds(vertici) = idCell0Ds;
//                                 vertices.push_back(idCell0Ds);
//                                 // cout << "vertice: " << idCell0Ds << endl;

//                                 idCell0Ds++;
//                                 idCell1Ds++;
//                                 if(vertici < 1){
//                                     verticesCell1Ds(0) = idCell0Ds;
//                                     mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                     mesh.VerticesCell1Ds.push_back(verticesCell1Ds);
//                                 }
//                                 Cell0Ds++;
//                                 Cell1Ds++;
//                                 vertici++;
//                             }
//                         }

//                         if(vertici == 2){
//                             mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                             mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                             // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                         }
//                         vCell1Ds(1) = idCell0Ds;
//                         mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                         mesh.VerticesCell1Ds.push_back(vCell1Ds);

//                         vertici = 0;
//                         for(unsigned int i = 0; i < coordinataFrattura.cols(); i++){
//                             // differenza primo dei due punti di intersezione vett matrPuntiDiIntersezionePrecTracciaPassante.row(0)coordinataFrattura.col(i) = coordinataFrattura.col(i) - matrPuntiDiIntersezionePrecTracciaPassante.row(0)
//                             Vector3d diffCoordinataFratturaPuntoDiIntersezioneTracciaNonPassante = coordinataFrattura.col(i).transpose() - matrPuntiDiIntersezioneTracciaNonPassante.row(0);
//                             // prodotto vettoriale n1 tra diffPuntiDiIntersezioneTracciaPassante e diffCoordinataFratturaPuntoDiIntersezioneTraccia
//                             Vector3d prodVettCoordFrattPunDiInterTraccia = diffPuntiDiIntersezioneTracciaNonPass.cross(diffCoordinataFratturaPuntoDiIntersezioneTracciaNonPassante);

//                             // cout << "vettore: " << prodVettCoordFrattPunDiInterTraccia.transpose() << endl;
//                             Vector3d absVec = prodVettCoordFrattPunDiInterTraccia.cwiseAbs(); // calcolo il valore assoluto di tutti gli elementi del vettore

//                             // Trova il massimo valore e la sua posizione
//                             Index maxIndex; // Definisce l'indice massimo usando Index
//                             absVec.maxCoeff(&maxIndex);

//                             if(prodVettCoordFrattPunDiInterTraccia(maxIndex) > 0 && idTracciaNonPassante.size() == 1 && idTracciaPassante.size() == 0){
//                                 // cout << "coordinate frattura: " << coordinataFrattura.col(i).transpose() << endl;
//                                 mesh.IdCell0Ds.push_back(idCell0Ds);
//                                 mesh.IdCell1Ds.push_back(idCell1Ds);
//                                 mesh.CoordinatesCell0Ds.push_back(coordinataFrattura.col(i).transpose());
//                                 verticiCell1Ds(vertici) = idCell0Ds;

//                                 // vertices.push_back(idCell0Ds);
//                                 // cout << "vertice: " << idCell0Ds << endl;

//                                 idCell0Ds++;
//                                 idCell1Ds++;
//                                 if(vertici < 1){
//                                     verticesCell1Ds(0) = idCell0Ds;
//                                     mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                                     mesh.VerticesCell1Ds.push_back(verticesCell1Ds);
//                                 }
//                                 Cell0Ds++;
//                                 Cell1Ds++;
//                                 vertici++;
//                             }
//                         }

//                         if(vertici == 2){
//                             mesh.VerticesCell1Ds.reserve(mesh.IdCell1Ds.size());
//                             mesh.VerticesCell1Ds.push_back(verticiCell1Ds);
//                             // cout << "vertice id: " << verticiCell1Ds.transpose() << endl;
//                         }
//                     }
//                 }
//             }
//         }

//         mesh.NumberCell0Ds = Cell0Ds;
//         mesh.NumberCell1Ds = Cell1Ds;
//         cout << "Cell0Ds: " << mesh.NumberCell0Ds << endl;
//         cout << "Cell1Ds: " << mesh.NumberCell1Ds << endl;
//         cout << "idCell0Ds: " << idCell0Ds << endl;
//         cout << "idCell1Ds: " << idCell1Ds << endl;

//         cout << "Celle 0D" << endl;
//         cout << "Id X Y Z" << endl;
//         unsigned int numCell0Ds = 0;
//         for(const auto& CoordinatesCell0Ds : mesh.CoordinatesCell0Ds){
//             cout << numCell0Ds << " " << fixed << scientific << setprecision(16) << CoordinatesCell0Ds.transpose() << endl;
//             numCell0Ds++;
//         }

//         cout << "Celle 1D" << endl;
//         cout << "Id Origin End" << endl;
//         unsigned int numCell1Ds = 0;
//         // unsigned int contatore = 0;
//         for(unsigned int i = 0; i < mesh.VerticesCell1Ds.size(); i++){
//             Vector2i vertices = mesh.VerticesCell1Ds[i];
//             cout << numCell1Ds << " " << vertices.transpose() << endl;
//             // contatore++;
//             // cout << "idCell2Ds: " << idCell2Ds << endl;
//             // if(contatore > 3){
//             //     Cell2Ds++;
//             //     idCell2Ds++;
//             //     cout << "Cell2Ds: " << Cell2Ds << endl;
//             //     contatore = 0;
//             // }
//             numCell1Ds++;
//         }

//         return true;
//     }
// }
