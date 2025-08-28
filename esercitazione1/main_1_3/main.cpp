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
        cout << "Errore: devi passare due interi come argomenti, M numero di lanci, N numero di blocchi." << endl;
        return 1;  // Termina il programma con codice di errore
    }

    // atoi Converti i parametri della riga di comando in interi
    int M = atoi(argv[1]);  // M numero di lanci
    int N = atoi(argv[2]);  // N numero di blocchi

if (N == 0) {  // Evita la divisione per zero
        std::cout << "Errore: divisione per zero!" << std::endl;
        return 1;
    }
    // Usa i valori di M e N
  cout << "Numero di lanci: " << M << endl;
  cout << "Numero di blocchi: " << N << endl;

int m=M/N; //numero di lanci per blocco

//INIZIALIZZAZIONE GENERATORE NUMERI CASUALI
 Random rnd;
 rnd.setup_random_generator(rnd);

double d=10.0; //distanza tra le linee parallele
double L=9.0;  //lunghezza ago


double r,x,y,Y; //distanza tra il punto medio dell'ago con la la linea parallela più vicina
vector <double> pi_prog, pi_gen;


for(int i=0; i<N; i++){
int hit=0;
double pi, prob;

    for(int j=0; j<m; j++){
            x= rnd.Gauss(0,1);
            y= rnd.Gauss(0,1);
            Y= abs(y)/sqrt(x*x+y*y);

            r=rnd.Rannyu(0, d/2.0);

            if(r <= L/2.0 * Y){ hit++; }

    }   
        prob=hit/double(m);

        pi=2*L/ (d*prob);
        cout<< pi <<endl;
        pi_gen.push_back(pi);
        pi_prog.push_back(Media(pi_gen));
}

vector<double> errori= DeviazioneStandardBlocchi(pi_gen);
cout<< "media progressiva "<<endl;
for(int i=0; i<pi_prog.size(); i++){

    cout<< pi_prog[i]<<"err"<< errori[i]<<endl;
}


ofstream outfile("medie_prog.txt");
if (!outfile) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
for (int i = 0; i < pi_prog.size(); i++) {
    outfile << pi_prog[i] << endl;

}

ofstream outfile1("errori.txt");
if (!outfile1) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
for (int i = 0; i < errori.size(); i++) {
    outfile1 << errori[i] << endl;

}



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
