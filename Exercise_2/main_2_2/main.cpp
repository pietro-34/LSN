/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <numeric> 
#include "random.h"
#include "funzioni.h"

using namespace std;
 
int main (int argc, char *argv[]){
//Passo da terminale il numero di elementi che voglio generare e la lunghezza del blocco
/*
if (argc != 3) {  // Verifica che siano passati 2 argomenti (oltre al nome del programma)
        cout << "Errore: devi passare due interi come argomenti, M randomwalks generati, N numero di blocchi." << endl;
        return 1;  // Termina il programma con codice di errore
    }

    // atoi Converti i parametri della riga di comando in interi
    int M = atoi(argv[1]);  // M numero di RN generati
    int N = atoi(argv[2]);  // N numero di blocchi

if (N == 0) {  // Evita la divisione per zero
        std::cout << "Errore: divisione per zero!" << std::endl;
        return 1;
    }
    // Usa i valori di M e N
  cout << "Numero di RW generati: " << M << endl;
  cout << "NUmero di blocchi: " << N << endl;

 int L= N/M; //elementi in un blocco
*/


//INIZIALIZZAZIONE GENERATORE NUMERI CASUALI
   //passo un file "Primes" con  i parametri del generatore
   Random rnd;
   rnd.setup_random_generator(rnd);

int M = 10000;           // numero totale di random walk
int Nstep = 100;         // passi per ogni random walk
int Nblk = 100;          // numero di blocchi
int L = M / Nblk;        // RW per blocco

vector<vector<double> > r2_blk(Nblk, vector<double>(Nstep, 0.0));  

// --- Loop sui blocchi
for (int iblk = 0; iblk < Nblk; iblk++) {
    // somma temporanea dei r^2(n) nel blocco
    vector<double> sum_r2(Nstep, 0.0);

    // genera L random walk per questo blocco
    for (int i = 0; i < L; i++) {
        vector<double> r2_walk(Nstep, 0.0);
        r2_walk=path_continuous(Nstep, rnd); //path_discrete(Nstep, rnd );

        

        // aggiungi i contributi del walk alla somma
        for (int n = 0; n < Nstep; n++) {
            sum_r2[n] += r2_walk[n];
        }
    }

    // calcola media di blocco
    for (int n = 0; n < Nstep; n++) {
        r2_blk[iblk][n] = sum_r2[n] / (double)L;
    }
}

// --- Ora hai r2_blk[iblk][n], cioè la media <r^2(n)> per ogni blocco

// Calcola media finale e incertezza sui blocchi
vector<double> mean(Nstep, 0.0), mean2(Nstep, 0.0), err(Nstep, 0.0);

for (int n = 0; n < Nstep; n++) {
    for (int iblk = 0; iblk < Nblk; iblk++) {
        mean[n]  += r2_blk[iblk][n];
        mean2[n] += r2_blk[iblk][n] * r2_blk[iblk][n];
    }
    mean[n]  /= (double)Nblk;
    mean2[n] /= (double)Nblk;

    if (Nblk > 1) {
        err[n] = sqrt((mean2[n] - mean[n]*mean[n]) / (Nblk - 1));
    } else {
        err[n] = 0.0;
    }
}

// --- Se vuoi lavorare con sqrt(<r^2>) invece che <r^2>
// se vuoi <r>
vector<double> R(Nstep, 0.0), R_err(Nstep, 0.0);
for (int n = 0; n < Nstep; n++) {
    R[n] = sqrt(mean[n]);
    if (mean[n] > 0) {
        R_err[n] = 0.5 * err[n] / sqrt(mean[n]);  // propagazione errore
    } else {
        R_err[n] = 0.0;
    }
}

//stampa:
/*
r= <r^2>
errore su R=sqrt(<r^2>) ---> \simga_R=[ 1/(2* sqrt(r))] * \sigma_r 
*/
//stampo su due fue file i risultati delle medie progressive e della deviazione standard a blocchi per importance sampling
ofstream outfile("mean_continuos.txt");  // crea (o sovrascrive) il file output.txt
if (!outfile) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
//medie_prog di r
for (int i = 0; i < R.size(); i++) {
    outfile << R[i] << endl;

}

ofstream outfile1("err_continuos.txt");  // crea (o sovrascrive) il file output.txt
if (!outfile1) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
//errori_blocchi di r
for (int i = 0; i < R_err.size(); i++) {
    outfile1 << R_err[i] << endl;

}


   
   return 0;
}


