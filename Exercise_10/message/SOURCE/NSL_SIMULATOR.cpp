/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <mpi.h>
#include <iostream>
#include "system.h"

using namespace std;


int main(int argc, char *argv[]) {
    
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Crea un oggetto System che rappresenta il sistema fisico da simulare e ne gestisce lo stato
    System SYS;
    
    // Inizializza il sistema leggendo i file di input e impostando le condizioni iniziali.

    
    
    SYS.initialize(rank);//passo rank in modo tale da inizializzare il random generator con primes diversi per ogni processo  
    

    SYS.initialize_properties(rank);//passo rank in modo tale da creare i file di output per ogni rank 
    
    SYS.compute_loss();
    
    SYS.sort_population();
    

    

   for(int i=0; i< SYS.get_ngen(); i++){ 
    SYS.crossover();
    SYS.mutation();//qui chromo_son viene asseganto a _chromo e rifatto il sort;
    
    
    if ((i+1) % 100 == 0) {
            SYS.exchange_best_chromosomes(rank, size);
        }
    

   SYS.measure(i,rank); // stampo la best loss e la loss mediata sulla migliore metà della popolazione; passo rank per la stampa dei risulati sui giusti file


   }
    

    MPI_Finalize();
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
