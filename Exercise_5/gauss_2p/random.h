/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

#include <string>


// This class contains functions for generating random numbers using the RANNYU algorithm
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // Default constructor
  Random();
  // Destructor
  ~Random();
  // Method to set the seed for the RNG
  void SetRandom(int * , int, int);
  // Method to save the seed to a file
  void SaveSeed(); 
  // Method to generate a random number in the range [0,1)
  double Rannyu(void);
  // Method to generate a random number in the range [min,max)
  double Rannyu(double min, double max);
  // Method to generate a random number with a Gaussian distribution
  double Gauss(double mean, double sigma);

  void setup_random_generator(Random& rnd,
                            const std::string& primes_file = "Primes",
                            const std::string& seed_file = "seed.in");
                            

 
  //Metodo per generare un numero compreso tra [0,1] distribuito come d(x)=2(x-1) 
  double taylor(void) ;

  //metedo per generare un GBM in modo diretto
  double GBM_diretto(double S_0, double mu, double sigma, double T_fin); 

  //metdod per generare un GBM in modo tra t_i+1 e t_i
  double GBM_indiretto(double S_0, double mu, double sigma, double T_fin, double T_in);

  //meteodo per generare un GBM in modo ricorsivo
  double GBM_ricorsivo(double S_t, double mu, double sigma, double dt, int steps) ;
  
};

 
#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
