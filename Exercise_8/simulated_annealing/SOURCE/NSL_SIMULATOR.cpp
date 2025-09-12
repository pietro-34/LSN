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
#include "particle.h"
#include "system.h"

using namespace std;


int main(int argc, char *argv[]) {
    // Inizializza il contatore per le configurazioni salvate in formato XYZ (usato per output periodico)
    int nconf = 1;
    
    // Crea un oggetto System che rappresenta il sistema fisico da simulare e ne gestisce lo stato
    System SYS;
    
    // Inizializza il sistema leggendo i file di input e impostando le condizioni iniziali.
    
    SYS.initialize();
    
    // Inizializza le proprietà da misurare durante la simulazione, leggendo l'elenco da "../INPUT/properties.dat".
    
    SYS.initialize_properties();
    
    // Esegue un reset degli accumulatori del blocco iniziale (blocco 0) prima di iniziare la simulazione a blocchi.
    
    SYS.block_reset(0);
    
    
    for(double t=SYS.get_temp(); t >0.002 ; t -= SYS.get_inverse_dt() ){

                        

                    SYS.simulated_annealing(t);
    }
    
    SYS.finalize();
    
    // Termina il programma restituendo 0 al sistema operativo (indicando che l'esecuzione è avvenuta senza errori).
    return 0;
} // fine della funzione main


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
