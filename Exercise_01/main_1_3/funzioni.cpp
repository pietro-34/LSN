#include "funzioni.h"
#include <cmath>

using namespace std;  // Ora non è necessario scrivere std:: ogni volta

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
vector<double> DeviazioneStandardBlocchi(const vector<double>& medie_blocchi) {
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

//ritoerna chi quadro totale
double Chiquadro(const vector<double>& freq, const vector<double>& expected, const vector<double>& var) {
    double chi2 = 0.0;
    for (size_t i = 0; i < freq.size(); i++) {
        chi2 += pow(freq[i] - expected[i], 2) / var[i];
    }
    return chi2;
}

//ritorna numero distribuito secon do distribuzione exp
double exp_distribution(double r, double lambda){

    double x;
    x= - (1/lambda)*log(r);
    return x;
}
//ritorna numero distribuito secondo distribu cauchy lorentz p(x)=1/pi * Gamma/[(x-mu)^2+Gamma^2]
double cl_distribution(double r, double Gamma, double mu ){

    double x;
    x=mu+Gamma*(tan(M_PI*(r-0.5)));

    return x;
}
/*
OSS: size_t nel ciclo for è un tipo di dato senza segno,definito in <cstddef>, 
comunemente usato per rappresentare dimensioni e indici di array o vettori in C++.

Maggiore sicurezza per loop che coinvolgono .size(), poiché .size() restituisce size_t,
usare size_t per i rende il codice più coerente e riduce il rischio di bug legati a segni negativi.
*/
