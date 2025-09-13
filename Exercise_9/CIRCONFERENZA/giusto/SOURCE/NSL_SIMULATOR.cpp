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
#include "system.h"

using namespace std;


int main(int argc, char *argv[]) {
    
    
    // Crea un oggetto System che rappresenta il sistema fisico da simulare e ne gestisce lo stato
    System SYS;
    
    // Inizializza il sistema leggendo i file di input e impostando le condizioni iniziali.
    
    SYS.initialize();//stampa population.dat

    
    
    SYS.compute_loss();//stampa pop_loss
    
    SYS.sort_population();

    
/*
   for(int i=0; i< SYS.get_ngen(); i++){ 
    SYS.crossover();
    SYS.mutation();//qui chromo_son viene asseganto a _chromo e rifatto il sort;

   }
    SYS.print_population("final_gen.dat");

*/
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
