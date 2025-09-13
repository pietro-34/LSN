/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Particle__
#define __Particle__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

class Chromosome {

private:
  int _ndim ;           // Dimensionality of the chromosome
  ivec _g;              // genetic code vector
  double _loss=0.0;          // loss 

public: // Function declarations
  void initialize(int n);                      // Initialize particle properties
  void check();                                // check the constrains of the TSP probem; the first gene is always 1; no gene repetition;  
  void permutation(int i, int j);              //Permutate two genes
  void shift(int i);                            //shifts genes to i-position
  void mix();                                  // mix the genes with suffle() respecting the bound;
  void set_loss(mat L);                      // compute and set the loss, by passing a mat;
  int get_gene(int k) const;                        //get the gene in position k
  double get_loss() const;                          //get the loss of chromosome
  Chromosome& operator=(const Chromosome& other); //overload = operator 
  bool operator<(const Chromosome& other) const;  //overload < operator; It confronts the loss of the chromos
  void print(int i);
  ivec get_rest_chromo(int m);         //get the _ndim-k genes of the chromo
  void set_rest_chromo(ivec sub, mat L);     //set the sub chromo after a crossover, and recompute the loss
};

#endif // __Particle__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
