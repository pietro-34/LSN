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

//dichiarazione vettori
vector<double> medie; //dichiaro un vettore in cui inserisco le sui blocchi
vector<double> errori_blocchi; //dichiaro un vettore in cui inserisco la deviazione standard a blocchi
vector<double> medie_prog;// vettore con media sui blocchi calcolata progressivamente
vector<double> var_medie;
vector<double> var_errori_blocchi;
vector<double> var_medie_prog;



//ciclo con cui genero i RN e calcolo media su ogni blocco e poi media progressiva per r e \sigma^2

for(int i=0; i<M; i++){
   vector<double> dati;  // dati riferiti al blocco i-esimo per la quantità <r>
   vector<double> dati_var; // dati riferiti al blocco i-esimo per la quantità sigma^2
   for(int j=0 ; j<L; j++ ){

      double y=rnd.Rannyu();
      double s= pow((y-0.5),2);
      dati.push_back(y);
      dati_var.push_back( s ) ;
   }
   
   //sul singolo blocco
   medie.push_back(Media(dati));
   medie_prog.push_back(Media(medie));

   var_medie.push_back(Media(dati_var));
   var_medie_prog.push_back(Media(var_medie));
   
   
}

errori_blocchi=DeviazioneStandardBlocchi(medie);
var_errori_blocchi=DeviazioneStandardBlocchi(var_medie);

//stampo su due fue file i risultati delle medie progressive e della deviazione standard a blocchi per r
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


//stampo su due fue file i risultati delle medie progressive e della deviazione standard a blocchi per var
//var_medie
ofstream outfile3("output_var_medie.txt");
if (!outfile3) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
for (int i = 0; i < var_medie_prog.size(); i++) {
    outfile3 << var_medie_prog[i] << endl;

}
//var_errori_blocchi
ofstream outfile4("output_errori_var.txt");
if (!outfile4) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
for (int i = 0; i < var_errori_blocchi.size(); i++) {
    outfile4 << var_errori_blocchi[i] << endl;

}

//RICHIESTA 1.03 ; test del chi^2 per valutare l'indipendenza statistica dei RN



int bin = 100;       // Numero di intervalli
int trials = 10000;  // Numero di numeri casuali per ogni test
int tests = 100;     // Numero di test chi-quadrato

  
vector<double> expected(bin, double(trials) / bin);  // Atteso n/M per ogni bin

vector<double> chi_values;  // Per salvare i valori di Chi^2

// Loop per 100 esperimenti
for (int test = 0; test < tests; test++) {
    vector<double> frequencies(bin, 0.0);  // Reset delle frequenze

    // Genera 10^4 numeri casuali e assegna ai bin
    for (int j = 0; j < trials; j++) {
        double rv = rnd.Rannyu();
        int bin_index = int(rv * bin);  // Trova il bin corrispondente
        frequencies[bin_index]++;
    }

    // Calcola Chi^2 per questo test
    double chi2 = Chiquadro(frequencies, expected, expected);
    chi_values.push_back(chi2);
}
//stampo chi su file
ofstream outfile5("chiquadro.txt");
if (!outfile5) {
        cerr << "Errore nell'apertura del file!" << endl;
        return 1;
    }
for (int i = 0; i < chi_values.size(); i++) {
    outfile5 << chi_values[i] << endl;

}




outfile5.close();  // chiudi il file












  //rnd.SaveSeed();

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
