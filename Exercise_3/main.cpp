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
#include <algorithm>
#include "random.h"
#include "funzioni.h"
#include <iomanip> 

using namespace std;
 
int main (int argc, char *argv[]){
//Passo da terminale il numero di elementi che voglio generare e la lunghezza del blocco

if (argc != 3) {  // Verifica che siano passati 2 argomenti (oltre al nome del programma)
        cout << "Errore: devi passare due interi come argomenti, N numero di numeri random generati, M numero di blocchi." << endl;
        return 1;  // Termina il programma con codice di errore
    }

    // atoi Converti i parametri della riga di comando in interi
    int M = atoi(argv[1]);  // M numero di RN generati
    int N = atoi(argv[2]);  // N numero di blocchi (ovvero numero di espirimenti)

if (N == 0) {  // Evita la divisione per zero
        std::cout << "Errore: divisione per zero!" << std::endl;
        return 1;
    }
    // Usa i valori di M e N
  cout << "Numero di RN generati: " << M << endl;
  cout << "Numero di blocchi: " << N << endl;

 //int L= M/N; //elementi in un blocco



//INIZIALIZZAZIONE GENERATORE NUMERI CASUALI
Random rnd;

rnd.setup_random_generator(rnd);

//parametri:
double S_0=100; //asset price
double T=1 ;    //delivery time
double K=100;   // strike price
double r=0.1;   // risk_free interest rate
double sigma=0.25; // volatility
double dt=0.01 ; //intervallo discreto
double steps=100; //numeri di passi per comprire l'intervallo [0,1] con un passo di 0.01

vector<double> call_dir, put_dir;
vector<double> call_ind, put_ind;
double d;
double n;
double c_dir, p_dir;
double c_ind, p_ind;
//genero 10^4 processi GBM
for(int i=0; i<M; i++){

    d=rnd.GBM_diretto(S_0,r,sigma,T);
    n=rnd.GBM_ricorsivo(S_0,r,sigma,dt,steps);

    //call
    c_dir= exp(-r*T) * ( max(d-K , 0.0));
    c_ind=exp(-r*T) * ( max(n-K , 0.0));
    
    call_dir.push_back(c_dir);
    call_ind.push_back(c_ind);

    //put
    p_dir=exp(-r*T) * ( max(K-d , 0.0));
    p_ind=exp(-r*T) * ( max(K-n , 0.0));

    put_dir.push_back(p_dir);
    put_ind.push_back(p_ind);
}

//DATBLOCKING
vector <vector <double > > data_call_dir, data_call_ind, data_put_dir, data_put_ind;
data_call_dir=Datablocking(call_dir, N);
data_call_ind=Datablocking(call_ind, N);
data_put_dir=Datablocking(put_dir,N);
data_put_ind=Datablocking(put_ind,N);

 

//STAMPO SU FILE

vector <string> names;
string filename;
names={"call","err"};
print_columns_to_file(data_call_dir, names ,"call_dir.txt", 6);

//stampo dati call;
vector <vector <double> > data_call={data_call_dir[0], data_call_dir[1], data_call_ind[0],data_call_ind[1]};

names={"CALL_dir", "Err_CALL_dir", "CALL_ind", "Err_CALL_ind"}; 
filename="output_call.txt";

print_columns_to_file(data_call, names ,filename, 6);
//stampo dati put;
vector <vector <double> > data_put={data_put_dir[0],data_put_dir[1],data_put_ind[0], data_put_ind[1]};

names={"PUT_dir", "Err_PUT_dir", "PUT_ind", "Err_PUT_ind"}; 
filename="output_put.txt";

print_columns_to_file(data_put, names ,filename, 6);
   return 0;
}


