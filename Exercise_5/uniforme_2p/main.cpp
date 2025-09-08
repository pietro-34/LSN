#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "random.h"
#include "funzioni.h"  // deve contenere p_1s e metropolis
using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: ./main.exe M N sigma\n";
        return 1;
    }

    int M = atoi(argv[1]);      // total steps
    int N = atoi(argv[2]);      // blocks (unused here)
    double sigma = atof(argv[3]);

    // Inizializza generatore random
    Random rnd;
    rnd.setup_random_generator(rnd);

    // Posizione iniziale (origine)
    vector<double> x = {2.0, 2.0, 2.0};

    // Output file punti campionati
    ofstream outfile("prob_2p.dat");
    outfile << "x\ty\tz\n";
    outfile << x[0] << "\t" << x[1] << "\t" << x[2] << "\n";

    // Metropolis parameters
    string metodo = "uniforme";
    int accepted = 0;
    
    vector<double> radialpoint_2p; //vettore dei moduli su cui poi applico datablocking

    //transitorio
    /*
    int K = 100000;  // burn-in di 100k passi //10 blocchi
    for (int i = 0; i < K; ++i) {
                auto m = metropolis(x, rnd, p_1s, sigma, "uniforme");
                    x = m.first;      // aggiorno lo stato ma non salvo
    }

    */

    for (int i = 0; i < M; i++) {
        auto m = metropolis(x, rnd, p_2p, sigma);
        x = m.first;     //
        radialpoint_2p.push_back(norm(x)); // se non uso tutti i punti che campiono per farci la media:
                                           // sottostimo gli stati da cui è difficile allontanarsi (bassa probabilità di accettazione), 
        accepted += m.second;              //e sovrastimo quelli con alta accettazione.

        //salvo i punti campionati per visualizzarli; salvo solo i punti generati da mosse accettate
        if (m.second == 1) {
            
            outfile << x[0] << "\t" << x[1] << "\t" << x[2] << "\n";
            
        }

        
    }


    outfile.close();
    
    cout << "Numero di mosse: " << M << endl;
    cout << "Accettate: " << accepted << endl;
    cout << "Accettanza: " << double(accepted) / double(M) << endl;


    //Datablocking
    vector<vector<double> > dati_2p;
    dati_2p=Datablocking(radialpoint_2p, N );

    //print to file:
    string filename="dati_2p.dat";
    vector<string> colonne={"Media blocco","Media prog", "errore"};

    print_columns_to_file( dati_2p , colonne,filename,10) ;
 
    return 0;
}
