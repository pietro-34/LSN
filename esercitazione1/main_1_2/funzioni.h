#ifndef FUNZIONI_H
#define FUNZIONI_H

#include <vector>
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

//distribuzione exp
double exp_distribution(double r, double lambda=1);

//distibuzione di cauchy-lorentz
double cl_distribution(double r, double Gamma=1, double mu=0 ) ;

#endif