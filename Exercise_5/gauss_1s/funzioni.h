#ifndef FUNZIONI_H
#define FUNZIONI_H

#include <vector>
#include "random.h"
using namespace std;

//print_file : stampa i dati su un file
void print_columns_to_file(const vector<vector<double> > & data, const vector<string>& nomi_colonne ,const string& filename, const int precisione);

// Calcola la media di un vettore di double
double Media(const vector<double>& dati);

// Calcola la deviazione standard di un vettore di double
double DeviazioneStandard(const vector<double>& dati);

//Calcola la deviazione standard a blocchi
vector<double> STD_media(const vector<double>& medie_blocchi);

//calcola la media progressiva e la std della media, sui blocchi
vector<vector <double> > Datablocking( const vector<double>& data, int N);


//calcola le componenti della sommatoria del chi_quadro
vector<double> Chiquadro_components(const vector<double>& osservation, const vector<double>& expected, const vector<double>& var ) ;

//chi_quadro
double Chiquadro(const vector<double>& freq, const vector<double>& expected, const vector<double>& var) ;


//funzione norma di un vettore
double norm(const vector<double>& v) ;

//funzione densità di prob GS atomo H
double p_1s (vector<double> v);

//funzione densità di probabilità n=2 l=1 m=0, orbitale 2p, atomo di H
double p_2p (vector<double> v);

//funzione metropolis, restituisce il punto dopo la mossa e l'accettanza di essa
pair <vector<double>, int> metropolis(const vector<double> & x_0, 
                                                     Random & rnd ,
                                                     double( *p) (vector<double> ) ,
                                                      double sigma ,
                                                      string flag="uniforme");
#endif