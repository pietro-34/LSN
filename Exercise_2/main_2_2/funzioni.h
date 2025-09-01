#ifndef FUNZIONI_H
#define FUNZIONI_H

#include <vector>
#include "random.h"
using namespace std;

// Calcola la media di un vettore di double
double Media(const vector<double>& dati);

// Calcola la deviazione standard di un vettore di double
double DeviazioneStandard(const vector<double>& dati);

//Calcola la deviazione standard a blocchi
vector<double> DeviazioneStandardBlocchi(const vector<double>& medie_blocchi);

//calcola le componenti della sommatoria del chi_quadro
vector<double> Chiquadro_components(const vector<double>& osservation, const vector<double>& expected, const vector<double>& var ) ;

//chi_quadro
double Chiquadro(const vector<double>& freq, const vector<double>& expected, const vector<double>& var) ;

//genera percorso lungo Nsteps in reticolo cubico discreto di passo a=1
vector<double> path_discrete(int Nstep, Random &rnd);

//genera un RW 3D continuo lungo Nstep
vector<double> path_continuous(int Nstep, Random &rnd);
#endif