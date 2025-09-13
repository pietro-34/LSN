/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <armadillo>
#include <stdlib.h> //exit
#include "chromosome.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {

private:

  int _n_city;              //number of city
  int _population;          //population---> total number of cromosomas
  string _spatial_config;   //spazial city configuration 
  int _ngen;                //number of generation 
  double p_c=0.7;           //probability of crossover
  double p_m_swap=0.08;         //probability of swap of two genes
  double p_m_shift=0.08;       //probability of shift
  double p_m_exchange=0.08;    //probability of exchange two blocks
  double p_m_block=0.08;       //probability of inversio ordered of a block
  mat city;              //city è una mat con 2 colonne ed _n_city righe; colonna 0 x e colonna 1 y
  mat distance2;         //matrice _n_city x _n_city contenente i valori delle distanze al quadrato tra le città
  Random _rnd;          // Random number generator
  field <Chromosome> _chromo; // Field of particle objects representing the system-- !!DIVENTERA' field <Cromosoma> _cromo;
  field <Chromosome> _chromo_son; //filed of the new generation of chromo
  field <Chromosome> _chromo_elite; // field of chromosome inherited by the new generation
  Chromosome _best_chromo; // auxiliar variable to store the best chromo ef each generation
  int elite=200; //number of chromo inherited by the new generation
  
public: // Function declarations
  int get_ngen();               //get the number of generation 
  int get_ncity();              // Get the number of blocks
  int get_population();           // Get the number of steps in each block
  void initialize();          // Initialize system properties
  void initialize_properties(); //read and initialize the file's property; (best_loss e avareged_block )
  void generate_city(bool C);  //genera le città su una circonferenza o dentro ad un quadrato
  void compute_distance2();    //calcola le distanze al quadrato tra le città e le assegna alla mat distance2
  void finalize();            // Finalize system and clean up
  void write_configuration(); // Write final system configuration to XYZ file
  void population();    //generate the chromosome population
  void compute_loss();  //compute loss for each particle of the field using set_loss() method by chromosome.h
  mat get_distance2(); //get distance squared matrix
  void sort_population(); //sort by loss with the overcast operator
  void print_population(string filename); //print the population on file
  void print_new_gen(string filename1); //print the _chormo_son on file
  void print_terminal(); //print loss of each chromo to terminal  
  ivec selection();   //selection operator, select a couple of chormos (THE FIELD MUST BE ORDERED BY LOSS)  retuns a couple of index
  void crossover();   //crossover operator over the population 
  void crossover_2chromo(Chromosome& _chromo1, Chromosome& _chromo2, int k); //cross two chromos
  ivec search_sub_in_chromo(ivec& sub, Chromosome & _chromo); //method used by crossover_2chromo
  void mutation(); //method that invoc permutation(i,j) by chromosome with probability P_m for all _chromo_son
  void measure(int i);  //mesure and print best Loss and loss averaged on best half population;
  
};

#endif // __System__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
