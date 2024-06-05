#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

namespace DFNLibrary{

// funzione per la lettura dei dati da file, in input il percorso dove si trova il file, Fract che corrisponde al riferimento
// della struttura per memorizzare i dati e il numero di fratture presenti, passando inizialmente il valore 0
bool letturaDatiFileFR(const string &percorsoFileFR, DFN &Fract, unsigned int &numFract){

    ifstream fileFR;
    fileFR.open(percorsoFileFR);

    if(fileFR.fail()){
        return false;

        // verificato che l'apertura del file sia avvenuta in modo corretto
    }else{

        string line;
        unsigned int idFrattura = 0, numVertici = 0;
        bool rigaTrovata = false;

        getline(fileFR, line);
        getline(fileFR, line);
        numFract = stoi(line); // converte la stringa con il numero di fratture in un intero
        // vengono lette le prime righe e memorizzate il numero di fratture presenti nel file

        for(unsigned int i = 0; i < numFract; i++){

            while(!fileFR.eof() && !rigaTrovata){ // procede la lettura del file saltando le righe con #, appoggiandosi a una variabile booleana
                getline(fileFR, line);
                // Salta riga di commento
                if(line[0] != '#'){
                    rigaTrovata = true;
                }
            }
            rigaTrovata = false;

            istringstream converter(line);
            string riga;

            // memorizza la prima informazione in cui ci sono l'id della frattura e il numero di vertici convertendo le stringhe a interi
            getline(converter,riga,';');
            idFrattura = stoi(riga);
            converter >> numVertici;

            // procede la lettura riga per riga, sempre separando i valori tramite il ; e escludendo le righe con #
            // vengono definite vettori (uno per ogni coordinate di lunghezza numVertices) e la matrice 3xnumVertices
            VectorXd coordX = VectorXd::Zero(numVertici);
            VectorXd coordY = VectorXd::Zero(numVertici);
            VectorXd coordZ = VectorXd::Zero(numVertici);
            MatrixXd coordinateVertici = MatrixXd::Zero(3,numVertici);

            while(!fileFR.eof() && !rigaTrovata){
                getline(fileFR, line);
                // Salta riga di commento
                if(line[0] != '#'){
                    rigaTrovata = true;
                }
            }
            rigaTrovata = false;

            // se non si trova il # si procede con la lettura delle righe per le coordinate x y z
            istringstream convert(line);
            for(unsigned int j = 0; j < numVertici; j++){
                getline(convert, riga, ';');
                coordX(j) = stod(riga);
            }

            getline(fileFR, line);
            istringstream conv(line);
            for(unsigned int j = 0; j < numVertici; j++){
                getline(conv, riga, ';');
                coordY(j) = stod(riga);
            }

            getline(fileFR, line);
            istringstream con(line);
            for(unsigned int j = 0; j < numVertici; j++){
                getline(con, riga, ';');
                coordZ(j) = stod(riga);
            }
            // ottenuti i vettori riga si vanno a mettere all'interno della matrice
            coordinateVertici << coordX.transpose(), coordY.transpose(), coordZ.transpose();

            // si riempie la mappa per le coordinate dei vertici definita nella struttura
            Fract.coordinateFratture[idFrattura] = coordinateVertici;

            // cout << fixed << scientific << setprecision(16) << Fract.coordinateFratture[idFractures] << endl;
            // cout << endl;

            // cout << fixed << scientific << setprecision(16) << coordinateVertici << endl;
            // cout << endl;
        }
        fileFR.close();
        return true;
        }
    }

// questa funzione servirà per estrarre il baricentro da ogni frattura, passando la mappa dove si è memorizzata e la variabile
// con il numero di fratture ottenuta dalla lettura del file, utile come criterio di arresto per i nostri cicli
// funzione per la scremature delle fratture lontane
bool calcoloBaricentriEDistBaricentroVertici(DFN &Fract, unsigned int &numFract, unsigned int &numIntersezioniFratture){

    vector<unsigned int> idFrattura;
    vector<Vector3d> coordinateBaricentro; // vettore dove a ogni id della frattura viene associato un vettore con le coordinate del baricentro
    vector<double> distMaxBaricentro; // vettore che a ogni id, quindi frattura, associa il valore massimo di distanza tra baricentro e vertici

    idFrattura.reserve(Fract.coordinateFratture.size());
    coordinateBaricentro.reserve(Fract.coordinateFratture.size());
    distMaxBaricentro.reserve(Fract.coordinateFratture.size());

    for(const auto& frattura : Fract.coordinateFratture){ // ciclo per scorrere le fratture, l'indice idFrattura corrisponde all'id della frattura

        MatrixXd& matrCoordinateFratture = Fract.coordinateFratture[frattura.first]; // prende la matrice con le coordinate dei vertici della frattura associata ad idFrattura
        unsigned int numVertici = matrCoordinateFratture.cols(); // numero di vertici della frattura

        // con questo comando .rowwise().sum() si fa la somma sulle righe della matrice che corrispondono alle x, y e z,
        // si fa la divisione tra la somma e il numero di colonne che corrisponde al numero di vertici
        // e si mette dentro il vettore coordBaricentro // prende il vettore delle coordinate del baricentro della frattura associata ad idFrattura
        Vector3d coordBaricentro = (matrCoordinateFratture.rowwise().sum()) / (matrCoordinateFratture.cols()); // vettore tridimensionale con le coordinate x y z del baricentro di ogni frattura

        // vettore distanza delle distanze al quadrato baricentro-vertice di ogni frattura
        VectorXd distanza = VectorXd::Zero(numVertici);

        for(unsigned int v = 0; v < numVertici; v++){
            distanza(v) = (coordBaricentro - matrCoordinateFratture.col(v)).squaredNorm(); // più efficiente e meno costoso computazionalmente del doppio ciclo for
        }
        double valMaxDistanza = distanza.maxCoeff(); // valore della distanza massimo per ogni frattura, da utilizzare per la scrematura

        idFrattura.push_back(frattura.first);
        coordinateBaricentro.push_back(coordBaricentro);
        distMaxBaricentro.push_back(valMaxDistanza); // valore della distanza massimo per ogni frattura si incolla nel vettore appositamente creato per memorizzare la massima distanza tra baricentro vertice
    }

    for(unsigned int i = 0; i < numFract; i++){
        for(unsigned int j = i+1; j < numFract; j++){

            // sostituzione con questo comando per calcolare la distanza al quadrato tra i due baricentri
            // ChatGPT dice questo:
            // metodo squaredNorm() di Eigen è probabilmente più efficiente e computazionalmente meno costoso rispetto al calcolo manuale della distanza al quadrato,
            // soprattutto in termini di semplicità del codice e possibili ottimizzazioni di libreria
            // calcola la distanza al quadrato tra i baricentri delle due fratture
            double quadratoDistanzaBaricentri = (coordinateBaricentro[idFrattura[i]] - coordinateBaricentro[idFrattura[j]]).squaredNorm();

            // calcola la somma delle massime distanze tra baricentro e vertici
            double sommaDistanzaMaxBaricentri = distMaxBaricentro[idFrattura[i]] + distMaxBaricentro[idFrattura[j]];

            // ora con il confronto si vanno a identificare le fratture che si intersecano e quindi ad escludere quelle che non si intersecano
            if(quadratoDistanzaBaricentri < sommaDistanzaMaxBaricentri){
                Fract.idFrattureCheSiIntersecano[numIntersezioniFratture] = make_pair(idFrattura[i], idFrattura[j]);
                numIntersezioniFratture++; // si tiene traccia di quante fratture si intersecano
            }
        }
    }
    return true;
    }

// con questa funzione andiamo a calcolare le equazioni dei piani che contengono le fratture e le equazioni dei lati delle fratture
bool calcoloEqPianoEdEqRetteLati(DFN& Fract){

    for(const auto& frattura : Fract.coordinateFratture){

        // vado a prendere l'id della frattura e la rispettiva matrice delle coordinate e del numero dei vertici associata
        const MatrixXd& matrCoordinateFratture = frattura.second; // scorre nel ciclo la matrice delle coordinate delle fratture
        // vado a prendere la matrice con le coordinate dei vertici salvati nella mappa

        // Trova due vettori direttori non paralleli contenuti nel piano
        Vector3d AB = matrCoordinateFratture.col(1) - matrCoordinateFratture.col(0); // Considera la colonna 1 come il vertice B e la colonna 0 come il vertice A
        Vector3d AC = matrCoordinateFratture.col(2) - matrCoordinateFratture.col(0); // Considera la colonna 2 come il vertice C e la colonna 0 come il vertice A

        // Calcola il vettore normale al piano utilizzando il prodotto vettoriale
        Vector3d normalePiano = AB.cross(AC);

        // Calcola il termine costante dell'equazione del piano
        double dPiano = - matrCoordinateFratture.col(0).dot(normalePiano);

        // Stampa l'equazione del piano
        // cout << "Per la frattura " << pair.first << ", l'equazione del piano : " << normale.transpose() << " * (x, y, z) + " << D_term << " = 0" << std::endl;

        Vector3d abcPiano(normalePiano(0), normalePiano(1), normalePiano(2));
        double termineNoto = dPiano;
        Fract.coeffabcPiano[frattura.first] = abcPiano;
        Fract.coeffdPiano[frattura.first] = termineNoto;

        // calcolo rette dei lati
        // conto le righe e le colonne di tale matrice, ovvero quanti vertici ha e le coordinate associate a questi vertici
        unsigned int numVertici = matrCoordinateFratture.cols(); // vertici sulle colonne
        unsigned int coordinatexyzVertici = matrCoordinateFratture.rows(); // coordinate associate ad ogni vertice sulle righe

        // definisco la matrice dei coefficienti direttori e dei punti iniziali delle rispettive rette
        MatrixXd coeffDirettoriRettaLati = MatrixXd::Zero(numVertici, coordinatexyzVertici);
        MatrixXd puntoRetta = MatrixXd::Zero(numVertici, coordinatexyzVertici);

        for(unsigned int v = 0; v < numVertici; v++){ // ciclo for che scorre il numero di vertici

            // si può sostituire l'if che computazionalmente costa tanto in questo modo
            // questo comando gestisce in modo più efficiente la differenza tra il primo e l'ultimo punto perchè non serve l'if

            // fa la differenza delle coordinate x, y e z tra due punti, da cui si ricavano i coefficienti direttori della retta
            // passante per quei due punti
            coeffDirettoriRettaLati.row(v) = matrCoordinateFratture.col((v + 1) % numVertici) - matrCoordinateFratture.col(v);

            // il punto iniziale della retta che mi serve per esprimere la retta in forma parametrica
            // lo ho già salvato praticamente sulle colonne della matrice delle coordinate delle fratture
            // inserisce in una mappa Fract.coordinateFratture i punti iniziali corrispondenti ad ogni retta dei lati della frattura considerata
        }

        // inserisce in una mappa i coefficienti direttori delle rette dei lati della frattura considerata
        Fract.coeffDirettoriRettaLati[frattura.first] = coeffDirettoriRettaLati;
    }
    return true;
    }

bool calcoloIntersezionePiani(DFN &Fract, unsigned int &numIntersezioniFratture){

    // calcolo retta di intersezione tra i due piani

    MatrixXd coeffDirettoriRettaTraccia = MatrixXd::Zero(numIntersezioniFratture, 3);
    Matrix3d A = Matrix3d::Zero();
    Vector3d terminiNoti = Vector3d::Zero();
    double tol = 100 * numeric_limits<double>::epsilon();

    // calcolo direttamente la retta della traccia delle fratture che si intersecano
    unsigned int k = 0;
    unsigned int intersezioneP = 0;
    for(const auto& frattura : Fract.idFrattureCheSiIntersecano){

        unsigned int idFrattura1 = frattura.second.first;
        unsigned int idFrattura2 = frattura.second.second;

        // aggiunta vettori temporanei se no non mi fa fare il prodotto vettoriale
        Vector3d vett1 = Fract.coeffabcPiano[idFrattura1];
        Vector3d vett2 = Fract.coeffabcPiano[idFrattura2];
        coeffDirettoriRettaTraccia.row(k) = vett1.transpose().cross(vett2.transpose());
        coeffDirettoriRettaTraccia.row(k) = (coeffDirettoriRettaTraccia.row(k).array() == -0).select(0, coeffDirettoriRettaTraccia.row(k)); // può succedere che ho -0 allora lo sostituisco con 0 con questa riga di codice

        A << Fract.coeffabcPiano[idFrattura1].transpose(), Fract.coeffabcPiano[idFrattura2].transpose(), coeffDirettoriRettaTraccia.row(k);
        terminiNoti << Fract.coeffdPiano[idFrattura1], Fract.coeffdPiano[idFrattura2], 0;
        terminiNoti = (terminiNoti.array() == -0).select(0, terminiNoti);

        // risolvo il sistema lineare per ricavare le coordinate del punto P se la matrice A ha determinante maggiore di una certa tolleranza
        if(abs(A.determinant()) > tol){

            // risolvo il sistema lineare A coordinatePuntoP = -terminiNoti e ricavo le coordinate del punto P
            Vector3d coordinatePuntoP = A.fullPivLu().solve(-terminiNoti);
            // può succedere che ho -0 allora lo sostituisco con 0 con questa riga di codice
            coordinatePuntoP = (coordinatePuntoP.array() == -0).select(0, coordinatePuntoP);

            // mi salvo gli id delle fratture che soddisfano questo if, questo mi servirà nella prossima funzione
            // mi salvo le coordinate del punto P che soddisfano la condizione dell'if
            Fract.coordinatePuntoP[make_pair(idFrattura1, idFrattura2)] = coordinatePuntoP;            
            // mi salvo i coefficienti direttori della traccia che soddisfano la condizione dell'if
            Fract.coeffDirettoriRettaTraccia[make_pair(idFrattura1, idFrattura2)] = coeffDirettoriRettaTraccia.row(intersezioneP);
            intersezioneP++;

            cout << "traccia numero: " << k << fixed << scientific << setprecision(16) << " P: " << Fract.coordinatePuntoP[make_pair(idFrattura1, idFrattura2)].transpose() << endl;
            cout << " (" << coeffDirettoriRettaTraccia.row(k) << ") * t" << endl;
        }
        k++;
    }
    return true;
    }

bool calcoloIntersezioneRettaTracciaERettalati(DFN &Fract, unsigned int &numeroTracceTotali){

    double tol = 1e+4 * numeric_limits<double>::epsilon();

    for(const auto& frattura : Fract.coordinatePuntoP){

        // prendo gli id delle due fratture che sto confrontando
        unsigned int idFrattura1 = frattura.first.first;
        unsigned int idFrattura2 = frattura.first.second;

        // vado a prendere il punto P e i coefficienti direttori della retta della traccia associati ai 2 id delle 2 fratture considerate
        Vector3d& coordinatePuntoP = Fract.coordinatePuntoP[make_pair(idFrattura1, idFrattura2)];
        Vector3d& coeffDirettoriRettaTraccia = Fract.coeffDirettoriRettaTraccia[make_pair(idFrattura1, idFrattura2)];

        // prendo i coefficienti direttori delle rette dei lati e i punti iniziali delle rette associate al primo id della frattura
        MatrixXd& coeffDirettoriRettaLatiFrattura1 = Fract.coeffDirettoriRettaLati[idFrattura1];
        const MatrixXd& puntoInizialeRettaFrattura1 = Fract.coordinateFratture[idFrattura1].transpose();

        // numero di vertici della prima frattura
        unsigned int numVerticiFrattura1 = coeffDirettoriRettaLatiFrattura1.rows();

        // prendo i coefficienti direttori delle rette dei lati e i punti iniziali delle rette associate al secondo id della frattura
        MatrixXd& coeffDirettoriRettaLatiFrattura2 = Fract.coeffDirettoriRettaLati[idFrattura2];
        const MatrixXd& puntoInizialeRettaFrattura2 = Fract.coordinateFratture[idFrattura2].transpose();

        // numero di vertici della seconda frattura
        unsigned int numVerticiFrattura2 = coeffDirettoriRettaLatiFrattura2.rows();

        // salvo tutte le intersezioni che soddisfano determinate condizioni specificate più avanti
        vector<Vector3d> intersezione;
        // mi salvo anche gli id delle tracce associate alle intersezioni che prendo
        vector<unsigned int> idFract;

        intersezione.reserve(4); // so già che devo avere massimo 4 punti di intersezione per questo faccio la reserve con 4
        idFract.reserve(4); // so già che devo avere massimo 4 id delle fratture 2 per frattura per questo faccio la reserve con 4

        // dichiaro ed inizializzo un vettore in cui vado a salvare le coordinate del punto P associato ad idFrattura1 ed idFrattura2
        Vector3d coordinateP = Vector3d::Zero();

        // dichiaro ed inizializzo un vettore in cui vado a salvare i coefficienti direttori della traccia associato ad idFrattura1 ed idFrattura2
        Vector3d coeffDirettoriTraccia = Vector3d::Zero();

        // contatore dei punti di intersezione che servirà più avanti
        unsigned int numPuntiDiIntersezione = 0;

        // parte il primo for in cui vado a calcolare i punti di intersezione tra la traccia e le rette dei lati della prima frattura
        for(unsigned int i = 0;  i < numVerticiFrattura1; i++){

            Vector3d coeffDirettoriRettaLatiFract1 = coeffDirettoriRettaLatiFrattura1.row(i);

            Vector3d prodottoVettorialeTracciaLati = coeffDirettoriRettaTraccia.cross(coeffDirettoriRettaLatiFract1);

            if(!prodottoVettorialeTracciaLati.isZero()){

                MatrixXd coeffDirettoriRette = MatrixXd::Zero(3,2);

                // prendo i coefficienti direttori della traccia e della retta del lato preso in considerazione
                coeffDirettoriRette.col(0) = coeffDirettoriRettaTraccia;
                coeffDirettoriRette.col(1) = - coeffDirettoriRettaLatiFrattura1.row(i);
                coeffDirettoriRette.col(1) = (coeffDirettoriRette.col(1).array() == -0).select(0, coeffDirettoriRette.col(1));
                // cout << fixed << scientific << setprecision(2) << coeffDirettoriRette << endl;
                // cout << fixed << scientific << setprecision(2) << puntoInizialeRettaFrattura1.row(i) << " " << - coordinatePuntoP.transpose() << endl;
                VectorXd diffPuntiIniziali = puntoInizialeRettaFrattura1.row(i) - coordinatePuntoP.transpose();

                // risolvo il sistema lineare coeffDirettoriRette parametrits = diffPuntiIniziali e ricavo i parametri t ed s
                // e quindi il punto di intersezione tra la traccia e le rette dei lati
                // t parametro retta in forma parametrica della traccia
                // s parametro retta in forma parametrica della prima frattura
                Vector2d parametrits = coeffDirettoriRette.fullPivLu().solve(diffPuntiIniziali);

                // separo i parametri t e s
                double t = parametrits(0);
                // double s = parametrits(1);

                // calcolo il punto di intersezione tra la retta della traccia e la retta del lato
                // ottengo le coordinate del punto di intersezione ma che devono passare ancora un controllo
                // punto di intersezione retta traccia con la prima frattura
                Vector3d puntoDiIntersezione1 = coordinatePuntoP + t * coeffDirettoriRettaTraccia;
                // Vector3d puntoDiIntersezione2 = puntoInizialeRettaFrattura1.row(i) + s * coeffDirettoriRettaLatiFrattura1.row(i);

                // calcolo del vettore che rappresenta il segmento
                Vector3d AB = puntoInizialeRettaFrattura1.row((i + 1) % numVerticiFrattura1) - puntoInizialeRettaFrattura1.row(i);
                // calcolo del vettore che va dal puntoInizialeRettaFrattura1.row(i) al punto di intersezione
                Vector3d AP = puntoDiIntersezione1.transpose() - puntoInizialeRettaFrattura1.row(i);

                // numeratore prodotto scalare tra AB e AP
                double numeratore = AP.dot(AB);
                // denominatore norma di AB elevata al quadrato quindi calcolo del prodotto scalare tra AB e AB
                double denominatore = AB.dot(AB);
                // double prod = AP.dot(AB) / AB.squaredNorm();

                // controllo se il punto di intersezione è interno usando il metodo della proiezione del punto di intersezione sul segmento del lato
                if(numeratore >= 0 && numeratore <= denominatore){
                    // cout << "punto di intersezione frattura " << idFrattura1 << ": " << puntoDiIntersezione1.transpose() << endl;
                    // i punti di intersezione che soddisfano questo if li salvo nel vector<Vector3d> intersezione che mi serivirà più avanti
                    intersezione.push_back(puntoDiIntersezione1);
                    // mi salvo anche gli idFrattura1
                    idFract.push_back(idFrattura1);
                    // mi salvo le coordinate del punto P
                    // (volendo questo si può togliere? lo metto anche nel prossimo for le coordinate del punto P basta salvarle una sola volta)
                    coordinateP = coordinatePuntoP;
                    // mi salvo i coefficienti direttori della retta della traccia
                    // (volendo questo si può togliere? lo metto anche nel prossimo for i coefficienti direttori della retta traccia basta salvarle una sola volta)
                    coeffDirettoriTraccia = coeffDirettoriRettaTraccia;
                    numPuntiDiIntersezione++;
                }
            }
        }

        // parte il secondo for in cui vado a calcolare i punti di intersezione tra la traccia e i lati della seconda frattura
        for(unsigned int i = 0;  i < numVerticiFrattura2; i++){

            Vector3d coeffDirettoriRettaLatiFract2 = coeffDirettoriRettaLatiFrattura2.row(i);

            Vector3d prodottoVettorialeTracciaLati = coeffDirettoriRettaTraccia.cross(coeffDirettoriRettaLatiFract2);

            if(!prodottoVettorialeTracciaLati.isZero()){

                MatrixXd coeffDirettoriRette = MatrixXd::Zero(3,2);

                coeffDirettoriRette.col(0) = coeffDirettoriRettaTraccia;
                coeffDirettoriRette.col(1) = - coeffDirettoriRettaLatiFrattura2.row(i);
                coeffDirettoriRette.col(1) = (coeffDirettoriRette.col(1).array() == -0).select(0, coeffDirettoriRette.col(1));
                VectorXd diffPuntiIniziali = puntoInizialeRettaFrattura2.row(i) - coordinatePuntoP.transpose();

                // risolvo il sistema lineare coeffDirettoriRette parametrits = diffPuntiIniziali e ricavo i parametri t ed s
                // e quindi il punto di intersezione tra la traccia e le rette dei lati
                // t parametro retta in forma parametrica della traccia
                // s parametro retta in forma parametrica della prima frattura
                Vector2d parametrits = coeffDirettoriRette.fullPivLu().solve(diffPuntiIniziali);

                // separo i parametri t e s
                double t = parametrits(0);
                // double s = parametrits(1);

                // calcolo il punto di intersezione tra la retta della traccia e la retta del lato
                // ottengo le coordinate del punto di intersezione ma che devono passare ancora un controllo
                // punto di intersezione retta traccia con la seconda frattura
                Vector3d puntoDiIntersezione2 = coordinatePuntoP + t * coeffDirettoriRettaTraccia;
                // Vector3d puntoDiIntersezione2 = puntoInizialeRettaFrattura2.row(i) + s * coeffDirettoriRettaLatiFrattura2.row(i);

                // calcolo del vettore che rappresenta il segmento
                Vector3d AB = puntoInizialeRettaFrattura2.row((i + 1) % numVerticiFrattura2) - puntoInizialeRettaFrattura2.row(i);
                // calcolo del vettore che va dal puntoInizialeRettaFrattura2.row(i) al punto di intersezione
                Vector3d AP = puntoDiIntersezione2.transpose() - puntoInizialeRettaFrattura2.row(i);

                // numeratore prodotto scalare tra AB e AP
                double numeratore = AP.dot(AB);
                // denominatore norma di AB elevata al quadrato quindi calcolo del prodotto scalare tra AB e AB
                double denominatore = AB.dot(AB);
                // double prod = AP.dot(AB) / AB.squaredNorm();

                // controllo se il punto di intersezione è interno usando il metodo della proiezione del punto di intersezione sul segmento del lato
                if(numeratore >= 0 && numeratore <= denominatore){
                    // cout << "punto di intersezione frattura " << idFrattura2 << ": " << puntoDiIntersezione2.transpose() << endl;
                    // i punti di intersezione che soddisfano questo if li salvo nel vector<Vector3d> intersezione che mi serivirà più avanti
                    intersezione.push_back(puntoDiIntersezione2);
                    // mi salvo anche gli idFrattura2
                    idFract.push_back(idFrattura2);
                    // mi salvo le coordinate del punto P
                    coordinateP = coordinatePuntoP;
                    // mi salvo i coefficienti direttori della retta della traccia
                    coeffDirettoriTraccia = coeffDirettoriRettaTraccia;
                    numPuntiDiIntersezione++;
                }
            }
        }
        cout << endl;

        // // controllo sul numero di punti di intersezione,
        // // ho tolto il controllo su 2 punti di intersezioni perchè tanto quando ne ho 2 ho l'id della stessa frattura e quindi non ha senso come condizione????
        if(numPuntiDiIntersezione == 4){

            VectorXd coordinataCurvilinea = VectorXd::Zero(numPuntiDiIntersezione);
            // inizializza il vettore posCoordCurvilinea con gli indici originali (0 1 2 3)
            VectorXi posCoordCurvilinea = VectorXi::LinSpaced(numPuntiDiIntersezione, 0, numPuntiDiIntersezione - 1);

            for(unsigned int i = 0; i < intersezione.size(); i++){
                cout << "Questi sono i punti di intersezione della frattura " << idFract[i] << " ecco le coordinate ";
                for(unsigned int j = 0; j < 3; j++){
                    cout << fixed << scientific << setprecision(7) << intersezione[i][j] << " ";
                }
                cout << endl;
            }
            cout << endl;

            cout << "coordinate punto P, che serve per l'ascissa curvilinea: " << fixed << scientific << setprecision(7) << coordinateP.transpose() << endl;
            cout << endl;

            cout << "coefficienti direttori retta traccia, che servono per ascissa curvilinea: " << fixed << scientific << setprecision(7) << coeffDirettoriTraccia.transpose() << endl;
            cout << endl;

            // /////// Controllo se gli alfa della stessa retta coincidono??????? ////

            unsigned int k = 0;
            for(unsigned int i = 0; i < intersezione.size(); i++){
                cout << idFract[i] << " ";
                for(unsigned int j = 0; j < 3; j++){
                    if(abs(coeffDirettoriTraccia(j)) > tol){
                        coordinataCurvilinea(k) = (intersezione[i][j] - coordinateP(j)) / coeffDirettoriTraccia(j);
                        cout << "coordinata curvilinea: " << fixed << scientific << setprecision(7) << coordinataCurvilinea(k) << " ";
                        k++;
                        break;
                    }
                }
                cout << endl;
            }
            cout << endl;

            // ordinamento dei vettori degli indici delle fratture (idFract) e delle coordinate curvilinee (coordinataCurvilinea)
            // in contemporanea in modo da non perdere l'informazione sull'id delle fratture
            for(unsigned int i = 0; i < coordinataCurvilinea.size() - 1; i++){
                for(unsigned int j = i + 1; j < coordinataCurvilinea.size(); j++){
                    if(coordinataCurvilinea(i) > coordinataCurvilinea(j)){
                        // Scambio dei dati
                        swap(coordinataCurvilinea[i], coordinataCurvilinea[j]);
                        // Scambio degli id delle fratture
                        swap(idFract[i], idFract[j]);
                        // Scambio delle posizioni degli indici originali
                        swap(posCoordCurvilinea[i], posCoordCurvilinea[j]);
                    }
                }
            }

            cout << "indici associati alle coordinate curvilinee ordinati: ";
            for(unsigned int i = 0; i < idFract.size(); i++){
                cout << idFract[i] << " ";
            }
            cout << endl;

            cout << "vettore della posizione delle coordinate ordinato: " << posCoordCurvilinea.transpose() << endl;

            cout << "vettore delle coordinate curvilinee ordinato: " << fixed << scientific << setprecision(7) << coordinataCurvilinea.transpose() << endl;
            cout << endl;

            MatrixXd coordinateIntersezioniTraccia = MatrixXd::Zero(2,3);
            // in questo primo caso prendo in considerazione il vettore delle coordinate curvilinee perchè ho coordinate curvilinee uguali 2 a 2
            // e in questo caso la traccia è passante per entrambe le fratture
            if(abs(coordinataCurvilinea[0] - coordinataCurvilinea[1]) < 1e-4 && abs(coordinataCurvilinea[2] - coordinataCurvilinea[3]) < 1e-4){
                cout << "traccia passante per la frattura " << idFract[0] << " e passante per la frattura "  << idFract[1] << endl;
                cout << fixed << scientific << setprecision(7) << intersezione[posCoordCurvilinea[1]].transpose() << " " << intersezione[posCoordCurvilinea[2]].transpose() << endl;

                coordinateIntersezioniTraccia.row(0) = intersezione[posCoordCurvilinea[1]].transpose();
                coordinateIntersezioniTraccia.row(1) = intersezione[posCoordCurvilinea[2]].transpose();
                bool tracciaPassante1 = false;
                bool tracciaPassante2 = false;
                double lunghezzaTraccia = (intersezione[posCoordCurvilinea[2]].transpose() - intersezione[posCoordCurvilinea[1]].transpose()).squaredNorm();

                if(idFract[0] < idFract[1]){
                    Fract.coordinateIntersezioniTracce[make_pair(idFract[0], idFract[1])] = coordinateIntersezioniTraccia;
                    Fract.traccePassantiONonPassanti1[make_pair(idFract[0], idFract[1])] = tracciaPassante1;
                    Fract.traccePassantiONonPassanti2[make_pair(idFract[0], idFract[1])] = tracciaPassante2;
                    Fract.lunghezzaTracce[make_pair(idFract[0], idFract[1])] = lunghezzaTraccia;
                }else{
                    Fract.coordinateIntersezioniTracce[make_pair(idFract[1], idFract[0])] = coordinateIntersezioniTraccia;
                    Fract.traccePassantiONonPassanti1[make_pair(idFract[1], idFract[0])] = tracciaPassante1;
                    Fract.traccePassantiONonPassanti2[make_pair(idFract[1], idFract[0])] = tracciaPassante2;
                    Fract.lunghezzaTracce[make_pair(idFract[1], idFract[0])] = lunghezzaTraccia;
                }

                numeroTracceTotali++; // si conta il numero di tracce

            // se gli indici centrali sono uguali la frattura associata a quell'id ha traccia passante invece
            // l'altra frattura associata all'altro id ha traccia non passante, perchè per quest'ultima frattura la traccia è interna
            }else if(idFract[1] == idFract[2]){
                cout << "traccia non passante per la frattura " << idFract[0] << " e passante per la frattura "  << idFract[1] << endl;
                cout << fixed << scientific << setprecision(7) << intersezione[posCoordCurvilinea[1]].transpose() << " " << intersezione[posCoordCurvilinea[2]].transpose() << endl;

                coordinateIntersezioniTraccia.row(0) = intersezione[posCoordCurvilinea[1]].transpose();
                coordinateIntersezioniTraccia.row(1) = intersezione[posCoordCurvilinea[2]].transpose();
                bool tracciaPassante = false;
                bool tracciaNonPassante = true;
                double lunghezzaTraccia = (intersezione[posCoordCurvilinea[2]].transpose() - intersezione[posCoordCurvilinea[1]].transpose()).squaredNorm();

                if(idFract[0] < idFract[1]){
                    Fract.coordinateIntersezioniTracce[make_pair(idFract[0], idFract[1])] = coordinateIntersezioniTraccia;
                    Fract.traccePassantiONonPassanti1[make_pair(idFract[0], idFract[1])] = tracciaNonPassante;
                    Fract.traccePassantiONonPassanti2[make_pair(idFract[0], idFract[1])] = tracciaPassante;
                    Fract.lunghezzaTracce[make_pair(idFract[0], idFract[1])] = lunghezzaTraccia;
                }else{
                    Fract.coordinateIntersezioniTracce[make_pair(idFract[1], idFract[0])] = coordinateIntersezioniTraccia;
                    Fract.traccePassantiONonPassanti1[make_pair(idFract[1], idFract[0])] = tracciaPassante;
                    Fract.traccePassantiONonPassanti2[make_pair(idFract[1], idFract[0])] = tracciaNonPassante;
                    Fract.lunghezzaTracce[make_pair(idFract[1], idFract[0])] = lunghezzaTraccia;
                }

                numeroTracceTotali++; // si conta il numero di tracce

            // se gli indici centrali sono diversi la traccia per entrambe le fratture è non passante,
            // ma bisogna stare anche attenti ed escludere i casi in cui abbiamo (idFract1 idFract1 idFract2 idFract2) o (idFract2 idFract2 idFract1 idFract1)
            }else if(idFract[1] != idFract[2] && (idFract[0] != idFract[1]) && (idFract[2] != idFract[3])){
                cout << "traccia non passante per la frattura " << idFract[1] << " e non passante per la frattura "  << idFract[2] << endl;
                cout << fixed << scientific << setprecision(7) << intersezione[posCoordCurvilinea[1]].transpose() << " " << intersezione[posCoordCurvilinea[2]].transpose() << endl;

                coordinateIntersezioniTraccia.row(0) = intersezione[posCoordCurvilinea[1]].transpose();
                coordinateIntersezioniTraccia.row(1) = intersezione[posCoordCurvilinea[2]].transpose();
                bool tracciaNonPassante1 = true;
                bool tracciaNonPassante2 = true;
                double lunghezzaTraccia = (intersezione[posCoordCurvilinea[2]].transpose() - intersezione[posCoordCurvilinea[1]].transpose()).squaredNorm();

                if(idFract[1] < idFract[2]){
                    Fract.coordinateIntersezioniTracce[make_pair(idFract[1], idFract[2])] = coordinateIntersezioniTraccia;
                    Fract.traccePassantiONonPassanti1[make_pair(idFract[1], idFract[2])] = tracciaNonPassante1;
                    Fract.traccePassantiONonPassanti2[make_pair(idFract[1], idFract[2])] = tracciaNonPassante2;
                    Fract.lunghezzaTracce[make_pair(idFract[1], idFract[2])] = lunghezzaTraccia;
                }else{
                    Fract.coordinateIntersezioniTracce[make_pair(idFract[2], idFract[1])] = coordinateIntersezioniTraccia;
                    Fract.traccePassantiONonPassanti1[make_pair(idFract[2], idFract[1])] = tracciaNonPassante1;
                    Fract.traccePassantiONonPassanti2[make_pair(idFract[2], idFract[1])] = tracciaNonPassante2;
                    Fract.lunghezzaTracce[make_pair(idFract[2], idFract[1])] = lunghezzaTraccia;
                }

                numeroTracceTotali++; // si conta il numero di tracce
            }
        }else{
            cout << numPuntiDiIntersezione << " " << idFrattura1 << " " << idFrattura2 << " pochi punti di intersezione" << endl;
        }

        cout << "fine for piu' esterno" << endl;
    }    
    return true;
    }

bool stampaDatiSuiFileDiOutput(const string &percorsoFileOutputPuntiDiIntersezione, const string &percorsoFileOutputLunghezzaTracce, DFN &Fract, unsigned int &numeroTracceTotali){

    // scrivo sul file di output punti di intersezione
    ofstream filePuntiDiIntersezione;
    // apro il file di output punti di intersezione
    filePuntiDiIntersezione.open(percorsoFileOutputPuntiDiIntersezione);

    if(filePuntiDiIntersezione.fail()){
        return false;
    }

    unsigned int idTracce = 0; // contatore id Tracce

    // definisco già i prossimi vettori, perchè anche se sto lavorando sulla stampa dei dati sul primo file punti di intersezione
    // posso già estrarre tutta l'informazione/i dati che mi serve/servono per la stampa dei dati sul secondo file lunghezza tracce
    // per i vettori definiti qui sotto consideriamo make_pair(idFratt1, idFratt2)
    vector<unsigned int> idFract1; // vettore in cui salvo il primo id della frattura (idFratt1)
    vector<unsigned int> idFract2; // vettore in cui salvo il secondo id della frattura (idFratt2)

    // per i prossimi 2 vettori consideriamo le due mappe
    // map<pair<unsigned int, unsigned int>, bool> traccePassantiONonPassanti1 = {};
    // map<pair<unsigned int, unsigned int>, bool> traccePassantiONonPassanti2 = {};
    vector<bool> passanteONonPassante1; // qui salvo bool riferito a traccePassantiONonPassanti1 quindi l'informazione associata alla prima colonna riferita a idFratt1
    vector<bool> passanteONonPassante2; // qui salvo bool riferito a traccePassantiONonPassanti2 quindi l'informazione associata alla seconda colonna riferita a idFratt2

    // in questo vettore invece salvo la lunghezza delle tracce
    vector<double> lunghezzaTracce;

    idFract1.reserve(Fract.coordinateIntersezioniTracce.size());
    idFract2.reserve(Fract.coordinateIntersezioniTracce.size());
    passanteONonPassante1.reserve(Fract.coordinateIntersezioniTracce.size());
    passanteONonPassante2.reserve(Fract.coordinateIntersezioniTracce.size());
    lunghezzaTracce.reserve(Fract.coordinateIntersezioniTracce.size());

    // continuiamo a lavorare per la stampa dei dati sul primo file punti di intersezione
    // stampo sul file prima l'intestazione e poi il numero di tracce
    filePuntiDiIntersezione << "# Number of Traces" <<endl;
    filePuntiDiIntersezione << numeroTracceTotali << endl;

    // stampo sul file l'intestazione Id traccia, Id frattura 1, Id frattura 2, punti di intersezione (x1; y1; z1; x2, y2; z2)
    filePuntiDiIntersezione << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    // in questo for faccio 2 cose
    // 1. vado a stampare i valori riferiti a questa intestazione (# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2)
    // 2. mi salvo tutti i dati nei vettori definiti in precedenza, vettori che mi serviranno per la stampa dei dati sul secondo file lunghezza tracce
    for(const auto &frattura : Fract.coordinateIntersezioniTracce){

        // le prossime righe di codice stampano i dati sul primo file di output punti di intersezione
        unsigned int idFratt1 = frattura.first.first;
        unsigned int idFratt2 = frattura.first.second;

        // estraggo l'informazione dei punti di intersezione che definiscono le tracce associate alle fratture idFratt1 e idFratt2
        MatrixXd stampaPuntiDiIntersezione = Fract.coordinateIntersezioniTracce[make_pair(idFratt1,idFratt2)];
        // stampo sul file punti di intersezione la seguente informazione # TraceId; FractureId1; FractureId2;
        filePuntiDiIntersezione << idTracce << "; " << idFratt1 << "; " << idFratt2 << "; ";

        // stampo sul file punti di intersezione le coordinate dei due punti di intersezione che definiscono le tracce X1; Y1; Z1; X2; Y2; Z2
        for(unsigned int t = 0; t < 2; t++){
            for(unsigned int s = 0; s < 3; s++){
                filePuntiDiIntersezione << fixed << scientific << setprecision(9) << stampaPuntiDiIntersezione(t,s);
                // questo if serve per capire se sono a fine riga o no, se sono a fine riga non metto il punto e virgola
                if(t == 1 && s == 2){
                    filePuntiDiIntersezione << "";
                }else{
                    filePuntiDiIntersezione << "; ";
                }
            }
        }
        filePuntiDiIntersezione << endl;
        // incremento il contatore del numero di tracce che definisce quindi anche l'id associato alla traccia
        idTracce++;

        // a questo punto del for estraggo prima dalle mappe relative tutta l'informazione e poi
        // salvo tutti i dati, che servono per la stampa sul secondo file di output lunghezza tracce, nei vettori definiti in precendenza
        // nelle prossime 2 mappe estraggo
        // nella prima la prima colonna associata alla prima frattura, in cui ho l'informazione se per la prima frattura la traccia è passante o non passante
        // nella seconda la seconda colonna associata alla seconda frattura, in cui ho l'informazione se per la seconda frattura la traccia è passante o non passante
        bool passONonPass1 = Fract.traccePassantiONonPassanti1[make_pair(idFratt1, idFratt2)];
        bool passONonPass2 = Fract.traccePassantiONonPassanti2[make_pair(idFratt1, idFratt2)];
        // in questa mappa estraggo la lunghezza della traccia associata alle 2 fratture idFratt1 e idFratt2
        double lunghezzaTraccia = Fract.lunghezzaTracce[make_pair(idFratt1, idFratt2)];

        idFract1.push_back(idFratt1); // salvo nel vettore idFract1 associato alla prima frattura l'id della frattura 1
        idFract2.push_back(idFratt2); // salvo nel vettore idFract2 associato alla seconda frattura l'id della frattura 2
        passanteONonPassante1.push_back(passONonPass1); // salvo nel vettore passanteONonPassante1 il booleano che dice se la prima frattura è passante o non passante
        passanteONonPassante2.push_back(passONonPass2); // salvo nel vettore passanteONonPassante2 il booleano che dice se la seconda frattura è passante o non passante
        lunghezzaTracce.push_back(lunghezzaTraccia); // salvo nel vettore lunghezzaTracce la lunghezza della traccia
    }
    filePuntiDiIntersezione.close(); // chiudo il primo file di output punti di intersezione

    // scrivo sul file di output tracce passanti o non passanti/ lunghezza tracce
    ofstream fileLunghezzaTracce;
    // apro il file di output lunghezza tracce
    fileLunghezzaTracce.open(percorsoFileOutputLunghezzaTracce);

    if(fileLunghezzaTracce.fail()){
        return false;
    }

    for(const auto& frattura : Fract.coordinateFratture){

        unsigned int idFrattura = frattura.first; // estraggo l'id della frattura
        unsigned int numeroTraccePerFrattura = 0; // contatore numero delle tracce per ciascuna frattura

        // segue un'altra serie di vettori in cui questa volta nel primo vettore salvo l'id della traccia e negli altri vettori
        // separo l'informazione in modo da avere tutte le informazioni sulle tracce passanti e sulle tracce non passanti per ciscuna frattura
        vector<unsigned int> idTraccia; // vettore id di tutte le tracce
        vector<unsigned int> idTracciaPassante; // vettore id traccia passante
        vector<unsigned int> idTracciaNonPassante; // vettore id traccia non passante
        vector<double> lunghezzaTracciaPassante; // vettore lunghezza traccia passante
        vector<double> lunghezzaTracciaNonPassante; // vettore lunghezza traccia non passante

        // for per contare il numero di tracce per ciascunaq frattura e salvare l'id delle tracce associato alla tracce presenti in ciascuna frattura
        for(unsigned int i = 0;  i < numeroTracceTotali; i++){
            if(idFrattura == idFract1[i] ||  idFrattura == idFract2[i]){
                numeroTraccePerFrattura++;
                idTraccia.push_back(i);
            }
        }

        // se il numero di tracce per quella determinata frattura è maggiore di 0 stampo tutta l'informazione sul file di output lunghezza tracce,
        // se il numero di tracce è zero non ha senso stampare solo l'intestazione è che il numero di tracce per quella frattura è 0
        if(numeroTraccePerFrattura > 0){
            // stampo sul secondo file di output lunghezza tracce la prima intestazione
            fileLunghezzaTracce << "# FractureId; NumTraces" << endl; // prima riga del secondo file di output

            cout << idFrattura << "; " << numeroTraccePerFrattura << endl; // cout di verifica da togliere

            // stampo sul file l'id della frattura e il relativo numero di tracce
            fileLunghezzaTracce << idFrattura << "; " << numeroTraccePerFrattura << endl;

            // stampo sul secondo file di output lunghezza tracce la seconda intestazione
            fileLunghezzaTracce << "# TraceId; Tips; Length" << endl; // seconda riga del secondo file di output
            // for in cui vado a verificare per ogni singola frattura quali tracce sono passanti e quali tracce sono non passanti
            for(unsigned int i = 0; i < numeroTraccePerFrattura; i++){

                cout << "id Traccia: " << idTraccia[i] << endl; // cout di verifica da togliere
                cout << "passante o no: " << boolalpha << passanteONonPassante1[idTraccia[i]] << " passante o no: " << boolalpha << passanteONonPassante2[idTraccia[i]] << endl; // cout di verifica da togliere

                // in questo if ed else if lavoro sulla colonna dell'id frattura 1 (if) e sulla colonna dell'id frattura 2 (else if)
                // e verifico se idFract1 (if) e se idFract2 (else if) sulla riga idTraccia[i] (id traccia) è passante se è passante salvo i dati
                if((idFrattura == idFract1[idTraccia[i]]) && (passanteONonPassante1[idTraccia[i]] == false)){
                    unsigned int idTr = idTraccia[i];
                    double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
                    idTracciaPassante.push_back(idTr); // salvo l'id della traccia perchè è passante
                    lunghezzaTracciaPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è passante
                }else if((idFrattura == idFract2[idTraccia[i]]) && (passanteONonPassante2[idTraccia[i]] == false)){
                    unsigned int idTr = idTraccia[i];
                    double lunghezzaTraccia = lunghezzaTracce[idTraccia[i]];
                    idTracciaPassante.push_back(idTr); // salvo l'id della traccia perchè è passante
                    lunghezzaTracciaPassante.push_back(lunghezzaTraccia); // salvo la lunghezza della traccia perchè è passante
                }

                // in questo if ed else if lavoro sulla colonna dell'id frattura 1 (if) e sulla colonna dell'id frattura 2 (else if)
                // e verifico se idFract1 (if) e se idFract2 (else if) sulla riga idTraccia[i] (id traccia) è non passante se è non passante salvo i dati
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

            // se c'è più di una traccia passante vado a riordinare le tracce passanti in ordine di lunghezza decrescente
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
            bool passante = false;
            // for che mi serve per stampare sul file di output lunghezza tracce il seguente contenuto
            // # TraceId; Tips; Length ovvero id traccia passante, false e lunghezza traccia passante in ordine decrescente
            for(unsigned int t = 0; t < idTracciaPassante.size(); t++){
                fileLunghezzaTracce << idTracciaPassante[t] << "; " << boolalpha << passante << "; " << fixed << scientific << setprecision(9) << lunghezzaTracciaPassante[t] << endl;
            }

            // se c'è più di una traccia non passante vado a riordinare le tracce non passanti in ordine di lunghezza decrescente
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
            bool nonPassante = true;
            // for che mi serve per stampare sul file di output lunghezza tracce il seguente contenuto
            // # TraceId; Tips; Length ovvero id traccia non passante, true e lunghezza traccia non passante in ordine decrescente
            for(unsigned int t = 0; t < idTracciaNonPassante.size(); t++){
                fileLunghezzaTracce << idTracciaNonPassante[t] << "; " << boolalpha << nonPassante << "; " << fixed << scientific << setprecision(9) << lunghezzaTracciaNonPassante[t] << endl;
            }
        }
    }
    fileLunghezzaTracce.close(); // chiudo il secondo file di output lunghezza tracce

    return true;
    }
}
