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

if (argc != 3) {  // Verifica che siano passati 2 argomenti (oltre al nome del programma)
        cout << "Errore: devi passare due interi come argomenti, N numero di numeri random generati, M numero di blocchi." << endl;
        return 1;  // Termina il programma con codice di errore
    }

    // atoi Converti i parametri della riga di comando in interi
    int N = atoi(argv[1]);  // N numero di RN generati
    int M = atoi(argv[2]);  // M numero di blocchi

if (N == 0) {  // Evita la divisione per zero
        std::cout << "Errore: divisione per zero!" << std::endl;
        return 1;
    }
    // Usa i valori di M e N
  cout << "Numero di RN generati: " << N << endl;
  cout << "NUmero di blocchi: " << M << endl;

 int L= N/M; //elementi in un blocco



//INIZIALIZZAZIONE GENERATORE NUMERI CASUALI
   //passo un file "Primes" con  i parametri del generatore
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");  //primes è un file con i parametri del generatore lineare congruenziale
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;    //leggiamo la prima riga di primes; primes su righe diverse generano sequenze statisticamente indipendenti(OSS: meglio cambiare i primes che il seed; cambiando il primes abbiamo la  certezza di avere sequenze statisticamente indipendenti)
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;



//RICHIESTE 1.01 1.02 

//dichiarazione vettori per punto 1 es1
vector<double> medie; //dichiaro un vettore in cui inserisco le sui blocchi
vector<double> errori_blocchi; //dichiaro un vettore in cui inserisco la deviazione standard a blocchi
vector<double> medie_prog;// vettore con media sui blocchi calcolata progressivamente

//dichiarazione vettore per punto 2 es2 importance sampling
vector<double> medie_sampling; 
vector<double> errori_blocchi_sampling; 
vector<double> medie_prog_sampling;



double pi_2=M_PI/2 ;
double y=0;
double s=0;
double t=0;
double z=0;
//ciclo con cui genero i RN e calcolo media su ogni blocco e poi media progressiva per r e \sigma^2

for(int i=0; i<M; i++){
   vector<double> dati; 
   vector<double> dati_sampling;

   for(int j=0 ; j<L; j++ ){
       //campiono con distribuzione uniforme
       y=rnd.Rannyu();
       s= pi_2*cos(pi_2 * y);
       dati.push_back(s);

       //campiono con taylor 
       t=rnd.taylor();
       z=pi_2*cos(pi_2 * t)/(2*(1-t));
       dati_sampling.push_back(z);
      
   }
   
   //sul singolo blocco per campionamento uniforme
   medie.push_back(Media(dati));
   medie_prog.push_back(Media(medie));

   //sul singolo blocco per campionamento con taylor
   medie_sampling.push_back(Media(dati_sampling));
   medie_prog_sampling.push_back(Media(medie_sampling));
   
   
   
}

errori_blocchi=DeviazioneStandardBlocchi(medie);
errori_blocchi_sampling=DeviazioneStandardBlocchi(medie_sampling);

//stampo su due fue file i risultati delle medie progressive e della deviazione standard a blocchi per il campionamento unforme
ofstream outfile("output_medie.txt");  // crea (o sovrascrive) il file output.txt
if (!outfile) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
//medie_prog di r
for (int i = 0; i < medie_prog.size(); i++) {
    outfile << medie_prog[i] << endl;

}

ofstream outfile2("output_errori.txt");  // crea (o sovrascrive) il file output.txt
if (!outfile2) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
//errori_blocchi di r
for (int i = 0; i < errori_blocchi.size(); i++) {
    outfile2 << errori_blocchi[i] << endl;

}

//stampo su due fue file i risultati delle medie progressive e della deviazione standard a blocchi per importance sampling
ofstream outfile3("output_medie_sampling.txt");  // crea (o sovrascrive) il file output.txt
if (!outfile3) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
//medie_prog di r
for (int i = 0; i < medie_prog_sampling.size(); i++) {
    outfile3 << medie_prog_sampling[i] << endl;

}

ofstream outfile4("output_errori_sampling.txt");  // crea (o sovrascrive) il file output.txt
if (!outfile4) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
//errori_blocchi di r
for (int i = 0; i < errori_blocchi_sampling.size(); i++) {
    outfile4 << errori_blocchi_sampling[i] << endl;

}


outfile.close();  // chiudi il file


   //rnd.SaveSeed();
   
   return 0;
}


