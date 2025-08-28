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
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <numeric> 
#include "random.h"
#include "funzioni.h"

using namespace std;
 
int main (int argc, char *argv[]){
//Passo da terminale il numero di elementi che voglio generare e la lunghezza del blocco

if (argc != 3) {  // Verifica che siano passati 2 argomenti (oltre al nome del programma)
        cout << "Errore: devi passare due interi come argomenti, M numero di S_N random generati, N numero di rnd utilizzati per S_N." << endl;
        return 1;  // Termina il programma con codice di errore
    }

    // atoi Converti i parametri della riga di comando in interi
    int M = atoi(argv[1]);  // N numero di RN generati
    int N = atoi(argv[2]);  // M numero di blocchi

if (N == 0) {  // Evita la divisione per zero
        std::cout << "Errore: divisione per zero!" << std::endl;
        return 1;
    }
    // Usa i valori di M e N
  cout << "Numero di S_N generati: " << M << endl;
  cout << "Numero di elementi delle somme parziali: " << N << endl;

 
//INIZIALIZZAZIONE GENERATORE NUMERI CASUALI
 Random rnd;
 rnd.setup_random_generator(rnd);

vector<double> sum_uniform;
vector<double> sum_exp;
vector<double> sum_cl;



//RICHIESTE 1.2 generare 10^4 S_N= 1/N sum_{i=1}_^{N} x_i con N=1,2,10,100 per uniforme, exp, cauchy-lorentz

for(int i=0; i<M; i++){
        double somma_uniform=0;
        double somma_exp=0;
        double somma_cl=0;

        for(int j=0; j<N;j++){

                somma_uniform += rnd.Rannyu();
                somma_exp += rnd.Exponential();
                somma_cl += rnd.CauchyLorentz();
        }
        sum_uniform.push_back((1.0/N)*somma_uniform);
        sum_exp.push_back((1.0/N)*somma_exp);
        sum_cl.push_back((1.0/N)*somma_cl);
}
// Costruisco il nome file
    ostringstream filename1;
    filename1 << "distrib_uniforme_" << N << ".txt";

    // Apro il file
    ofstream outfile(filename1.str());
    if (!outfile) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
    for(int i=0; i<sum_uniform.size() ; i++){
            outfile << sum_uniform[i] <<endl;
    }

    outfile.close();

ostringstream filename2;
    filename2 << "distrib_exp_" << N << ".txt";

    // Apro il file
    ofstream outfile1(filename2.str());
    if (!outfile1) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
    for(int i=0; i<sum_exp.size() ; i++){
            outfile1 << sum_exp[i] <<endl;
    }

    outfile1.close();

ostringstream filename3;
    filename3 << "distrib_cl_" << N << ".txt";

    // Apro il file
    ofstream outfile2(filename3.str());
    if (!outfile2) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
    for(int i=0; i<sum_cl.size() ; i++){
            outfile2 << sum_cl[i] <<endl;
    }

    outfile2.close();



   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
