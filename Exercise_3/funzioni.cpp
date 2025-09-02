#include "funzioni.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <numeric> 
#include <algorithm>
#include <iomanip> //setprecision 

using namespace std; 

double Media(const vector<double>& dati) {
    double somma = 0.0;
    for (double valore : dati) {
        somma += valore;
    }
    return somma / dati.size();
}

double DeviazioneStandard(const vector<double>& dati) {
    double media = Media(dati);
    double sommaQuadrati = 0.0;
    for (double valore : dati) {
        sommaQuadrati += pow(valore - media, 2);
    }
    return sqrt(sommaQuadrati / dati.size());
}


// Per questa funzione bisogna passare un vettore con i valori medi dei blocchi!
 //funzione indipendente dalla dimensione del vettore medie_blocchi
 //         COSE DA MODIFICARE:      
 // !!!! USARE FUNZIONE MEDIA PER FARE LE MEDIE IN QUESTA FUNZIONE !!!

 
vector<double> STD_media(const vector<double>& medie_blocchi) {
    vector<double> dev_std_blocchi;
    double somma = 0.0;
    double somma_quad = 0.0;

    for (size_t i = 0; i < medie_blocchi.size(); i++) {
        somma += medie_blocchi[i];
        somma_quad += medie_blocchi[i] * medie_blocchi[i];

        // Calcolare il valore usando la tua formula
        double N = i + 1;  // i parte da 0, quindi N = i + 1
        if (N > 1) {
            double media = somma / N;
            double media_quad = somma_quad / N;

            double sigma = sqrt((media_quad - media * media) / (N - 1));
            dev_std_blocchi.push_back(sigma);
        } else {
            dev_std_blocchi.push_back(0.0);  // Per il primo blocco, la deviazione è 0
        }
    }

    return dev_std_blocchi;
}

//funzione datablocking: accetta un vettore , lo suddivide in M blocchi, calcola la media su ogni blocco e la deviazione standard
// ritorna una matrice di due vettori con le componenti della media progressiva e la deviaizone standard a blocchi progressiva

//chat:
vector<vector<double>> Datablocking(const vector<double>& data, int M) {
    int N = data.size();
    int L = N / M;  // lunghezza blocco

    vector<double> ave(M, 0.0);       // medie dei singoli blocchi
    vector<double> prog_mean(M, 0.0); // media progressiva
    vector<double> errori(M, 0.0);    // deviazione standard progressiva

    // Calcola media di ogni blocco
    for (int i = 0; i < M; i++) {
        vector<double> blocco(data.begin() + i * L, data.begin() + (i + 1) * L);
        ave[i] = Media(blocco);
    }

    // Calcola media progressiva
    for (int i = 0; i < M; i++) {
        vector<double> subvector(ave.begin(), ave.begin() + i + 1);
        prog_mean[i] = Media(subvector);
    }

    // Calcola deviazione standard progressiva usando STD_media
    errori = STD_media(ave);

    return {prog_mean, errori};
}


/* mia versione sbagliata 22/03/2024
vector<vector<double>> Datablocking( const vector<double>& data, int M){
        
        int N=data.size();
        int L=N/M; //lunghezza blocco
        vector<double> prog_mean;
        vector<double> ave;
        vector<double> errori;

        double mean;
        
        for(int i=0; i<M; i++){
                vector <double> block;
            //ciclo sul blocco
            for(int j=0; j<L; j++){

                    block.push_back(data[j+L*i]);
            }
            //media sul blocck
            mean=Media(block);
            ave.push_back(mean);
            prog_mean.push_back(mean);
        }
            errori=STD_media(ave);

    vector<vector<double>> r={prog_mean, errori};

    return r;

}
*/

//SCRIVERE FUNZIONE CHI^2 GENERICA che mi restituisce un vettore di chi_quadro; 
// cioè: chi^2=sum(i=0; M)[chi^2_i]=sum(i=0; M) [(O_i-E_i)^2/var_i]
// mi resituisce gli addendi della sommatoria
// ACCETTA TRE VECTOR O_i (osservazioni) E_i (expected values), Var_i varianza sperimentale

vector<double> Chiquadro_components(const vector<double>& osservation, const vector<double>& expected, const vector<double>& var ) {
    vector<double> chi;
    
    for(size_t i=0; i<osservation.size(); i++){
        chi.push_back(pow((osservation[i]-expected[i]),2)/var[i]);
        
    }

    return chi;
}

double Chiquadro(const vector<double>& freq, const vector<double>& expected, const vector<double>& var) {
    double chi2 = 0.0;
    for (size_t i = 0; i < freq.size(); i++) {
        chi2 += pow(freq[i] - expected[i], 2) / var[i];
    }
    return chi2;
}



/*
OSS: size_t nel ciclo for è un tipo di dato senza segno,definito in <cstddef>, 
comunemente usato per rappresentare dimensioni e indici di array o vettori in C++.

Maggiore sicurezza per loop che coinvolgono .size(), poiché .size() restituisce size_t,
usare size_t per i rende il codice più coerente e riduce il rischio di bug legati a segni negativi.
*/


/*funzione stampa su file
 accetta in input:
 -1 una matrice di dati, le colonne sono le diverse quantità fisiche mentre le righe i rispettivi valori
 -2 vettori di stringhe con i nomi delle colonne
 -3 il nome con cui vogliamo chiamare file di output
 -4 il numero di cifre decimali 

 prima di chiamare questa funzione bisogna quindi mettere i dati una matrice
 vector <vector<double>> data={medie_prog, errori_blocchi, energia_media }

*/
void print_columns_to_file(const vector<vector<double> > & data,const vector<string>& nomi_colonne ,const string& filename,const int precisione) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Errore nell'apertura del file: " << filename << endl;
        return;
    }

    size_t N = data[0].size(); // dimensione dei vettori
    size_t cols = data.size(); // numero di colonne da stampare

    
    

    for(size_t i=0; i<cols; i++){
        outfile << nomi_colonne[i] << "\t";
    }

    outfile<< endl;

    outfile << fixed << setprecision(precisione);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < cols; j++) {
            outfile << data[j][i];
            if (j < cols - 1)
                outfile <<setprecision(precisione)<< "\t";  // separatore tra colonne
        }
        outfile << endl;
    }

    outfile.close();
}