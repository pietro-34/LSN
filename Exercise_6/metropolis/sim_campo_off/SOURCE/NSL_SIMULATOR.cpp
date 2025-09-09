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
    
    // Ciclo esterno sui blocchi di simulazione:
    
    for(int i = 0; i < SYS.get_nbl(); i++) { // loop sui blocchi
    
        // Ciclo interno sui passi di simulazione all'interno del blocco i-esimo:
        
        for(int j = 0; j < SYS.get_nsteps(); j++) { // loop sugli step nel blocco
            // Esegue un passo di simulazione.
    
            //       - Nel caso Ising (sim_type 2 o 3), System::move() seleziona uno spin casuale e prova a flipparlo. 
            //       -Se sim_type=2 (Metropolis), la variazione di energia ΔE è calcolata e si accetta il flip con probabilità exp(-βΔE)
            //       -Se sim_type=3 (Gibbs sampling), il flip può essere deciso direttamente in base alla probabilità condizionale. Anche qui _nattempts aumenta di 1 per tentativo e _naccepted di 1 per ogni flip accettato.
            //   * Alla fine di System::step(), viene aggiornato _nattempts += _npart (vengono conteggiati _npart tentativi per step, sia in MD che MC, per uniformità di definizione dello step). 
            SYS.step();
            
            // Esegue la misurazione delle proprietà fisiche dopo il passo.
            
            SYS.measure();
            
            // Ogni 50 passi circa, salva la configurazione attuale su file (opzione per analisi visuale).
            //if(j % 50 == 0) {
                // SYS.write_XYZ(nconf); 
                //nconf++;  // Incrementa il contatore di configurazioni anche se la scrittura è disabilitata, mantenendo la numerazione.
            //}
        } // fine del loop sui passi all'interno del blocco i
        
        
        
        SYS.averages(i+1);
        
        
        SYS.block_reset(i+1);
    } // fine del loop sui blocchi
    
    
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
