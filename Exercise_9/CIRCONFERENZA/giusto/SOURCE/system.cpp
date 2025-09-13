/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <cstdlib>
#include <string>
#include "system.h"

using namespace std;
using namespace arma;


void System :: initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  //questa parte rimane uguale- LETTURA PRIMES
  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  
  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;                      
  double delta;
  while ( !input.eof() ){
    input >> property;
    if( property == "N_CITY" ){
      input >> _n_city;
      coutf << "N CITY= " << _n_city << endl;    

    } else if( property == "POPULATION" ){
      input >> _population;
      coutf << "POPULATION= " << _population << endl;
      _chromo.set_size(_population);
      _chromo_son.set_size(_population);
    } else if( property == "SPAZIAL_CONFIGURATION" ){
      input >> _spatial_config;
      
      coutf << "SPAZIAL CONFIGURATION= " << _spatial_config << endl;
    } else if( property == "N_GENERATION" ){
      input >> _ngen;
      coutf << "NUMBER OF GENERATION= " << _ngen << endl;
    } else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  
  //chiudo l'input e creo le città; metodo city() che prende in argomento un bool ; se 0--> circumference se 1--> square
  
  input.close();
  bool C;
  if (_spatial_config=="CIRCUMFERENCE"){ C=false;}
    else{C=true;}
  this->generate_city(C); // read_configuration diventa generate_city
  this->compute_distance2();// calcola distanze al quadrato tra le città e le assegna alla matrice
  this->population(); //genera la popolazione
  coutf << "Population generated!" << endl;
  coutf.close();
  return;
}

//genera le coordinate della città e le assegna alla matrice city; stampa su file con .save()
void System :: generate_city(bool C){

    city.set_size(_n_city,2); //matrice  _n_city x 2=  34 x 2 le colonne corrispondono alle coordinate delle città; colonna zero x, colonna 1 y

  if(C==false){    //genero le città su di una circonferenza di raggio 10

                  double r=10.0;
                  double theta,x,y;
                  for(int i=0; i<city.n_rows; i++){
                        theta= _rnd.Rannyu(0, 2*M_PI);
                        x=r*cos(theta);
                        y=r*sin(theta);
                        
                        city(i,0)=x;
                        city(i,1)=y;}
                      
    
          }else{            double x,y;  //genero le città all'interno di un quadrato di lato L=10
                            for(int i=0; i<city.n_rows; i++){
                                  
                                  x=_rnd.Rannyu(0, 10);
                                  y=_rnd.Rannyu(0, 10);
                                  
                                  city(i,0)=x;
                                  city(i,1)=y;}
  }

          //salvo su file con il metodo .save()
          city.save("../OUTPUT/coordinate_città.csv",  csv_ascii);

  return;
}

//calcola distanze al quadrato e le assegna alla mat distance2
void System :: compute_distance2(){

    distance2.set_size(_n_city,_n_city);
    distance2.zeros();
     
    double dx,dy, d;

    for(int i=0; i<_n_city;i++){
      for(int j=0;j<_n_city;j++){
        
                      dy=city(i,1)-city(j,1);
                      dx=city(i,0)-city(j,0);
                      d= dx*dx + dy*dy;
                      distance2(i,j)=d;
                      }
                  }
      distance2.save("../OUTPUT/distance2.csv",  csv_ascii);
  return;
}

void System :: population(){

  

  ofstream coutp;
  coutp.open("../OUTPUT/population.dat");

    for (int i=0; i<_population; i++){

      _chromo(i).initialize(_n_city);
      _chromo(i).mix();
      _chromo(i).check();


      for (int k=0; k<_n_city; k++){

            coutp << _chromo(i).get_gene(k)<<"\t";
      }
      coutp<<endl;
    }
  coutp.close();

  //initialize the chromo son field
  for (int i=0; i<_population; i++){

      _chromo_son(i).initialize(_n_city);
      
      _chromo_son(i).check();
  }

    

    return;
}

void System :: print_population(string filename ){

  ofstream coutps;
  coutps.open("../OUTPUT/" + filename);

  for (int i=0; i<_population; i++){
                for (int k=0; k<_n_city; k++){

            coutps << _chromo(i).get_gene(k)<<"\t";
      }
      coutps<<" loss "<< _chromo(i).get_loss() <<endl;
    }
  coutps.close();
    

  return;
}

void System :: print_new_gen(string filename ){

  ofstream coutng;
  coutng.open("../OUTPUT/" + filename);

  for (int i=0; i<_population; i++){
                for (int k=0; k<_n_city; k++){

            coutng << _chromo_son(i).get_gene(k)<<"\t";
      }
      coutng<<" loss "<< _chromo_son(i).get_loss() <<endl;
    }
  coutng.close();
    

  return;
}

void System :: compute_loss(){

    //stampo il cromosoma con la rispettiva loss
    ofstream coutl;
    coutl.open("../OUTPUT/pop_loss.dat");
    
    for(int i=0; i<_population; i++){

        _chromo(i).set_loss(distance2);

        //stampa di debug:
        for (int k=0; k<_n_city; k++){

            coutl << _chromo(i).get_gene(k)<<"\t";
      }
      coutl<< "loss "<< _chromo(i).get_loss() <<endl;

    }

    return;
}

void System::sort_population() {
    // bubble sort su _population elementi
    for(int pass = 0; pass < _population - 1; ++pass) {
        // dopo pass iterazioni, gli ultimi pass elementi sono già al loro posto
        for(int j = 0; j < _population - pass - 1; ++j) {
            // se il cromosoma j+1 ha loss minore di j, scambiali
            if (_chromo(j+1) < _chromo(j)) {
                swap(_chromo(j), _chromo(j+1));
            }
        }
    }
}

//implemento j=int(M * r^p)+1 ; il +1 non selezionerebbe mai il _chromo(0) quindi non lo aggiungo
//ritona la coppia di indici 
ivec System :: selection(){
  
  ivec index(2);
  
   index(0)=int(_population * pow(_rnd.Rannyu(),2));
   index(1)=int(_population * pow(_rnd.Rannyu(),2));
   return index;
}

void System ::  mutation(){
  int index1, index2;
  //int mutation=0;
  for(int i=0; i< _population ; i++){

    if(_rnd.Rannyu() < p_m){

        //mutation++;
        index1= int(_rnd.Rannyu()*(_n_city-1))+1; //in questo modo non può essere zero, quindi non permuto il primo gene(rispetto il bound)
        index2= int(_rnd.Rannyu()*(_n_city-1))+1; 

          //cout<<"cromosoma "<< i <<"  indici di gene "<<index1<<"  "<< index2<<endl; 

        _chromo_son(i).permutation(index1, index2);
        _chromo_son(i).set_loss(distance2);
        _chromo_son(i).check();

    }
    if(_rnd.Rannyu()< p_m){
          int k=(_n_city * _rnd.Rannyu())+1;
        _chromo_son(i).shift(k);
        _chromo_son(i).check();
        _chromo_son(i).set_loss(distance2);

    }
  }
  //cout<< "cromosomi mutati "<< double(mutation)/double(_population)<< endl;
  _chromo=_chromo_son;
  sort_population();

  return;
}

void System :: crossover(){

  ivec index(2);
  double n_c;
  n_c=_population/2;
  int coppie_cros=0;
//ciclo sulle coppie
  for(int i=0; i< n_c ; i++){
          
            //seleziono la coppia e assegno gli indici all'ivec index
            index = selection();

            //se la _rnd < probabilità di fare crossover allora faccio il crossover 
            if(_rnd.Rannyu()<p_c){
                
                //cout<<"coppia "<< i <<" indici " << index(0)<<"\t"<< index(1) <<endl;
              //ora devo selezionare un indice di gene a caso tra quelli disponibili; e poi fare il cut 
              //in questo modo non considero mai l'indice 0
              int k=int( _n_city * _rnd.Rannyu())+1;


                crossover_2chromo( _chromo(index(0)), _chromo(index(1)),  k);
                
                  //i genitori crossoverizzati diventano i figli
                _chromo_son(i)=_chromo(index(0)); 
                _chromo_son(n_c+i)=_chromo(index(1));
                _chromo_son(i).check();   //after crossover check the bounds
                _chromo_son(n_c+1).check(); //
                coppie_cros++;
              
              
              } else{   _chromo_son(i)=_chromo(index(0)); //altrimenti copia genitori in figli
                        _chromo_son(n_c+i)=_chromo(index(1));
                      
                    }
  }
    

}

void System :: crossover_2chromo(Chromosome& _chromo1, Chromosome& _chromo2, int k){

            int dimsub= _n_city-k;
        //se k=0 e k=_n_city -1 madre e padre non cambiano, in realtà non può essere k=0 nella chiamata in crossover() 
        if (k != (_n_city-1) ){
            ivec sub_c1(dimsub), sub_c2(dimsub), sub_new1(dimsub), sub_new2(dimsub);

            sub_c1=_chromo1.get_rest_chromo(k); //metodo che ritorna gli ultimi _n_city-k elementi in un vettore ivec
            sub_c2=_chromo2.get_rest_chromo(k);

            
            sub_new1=search_sub_in_chromo(sub_c1, _chromo2);
            sub_new2=search_sub_in_chromo(sub_c2, _chromo1);

            _chromo1.set_rest_chromo(sub_new1, distance2); //sostituisce sub_new1 in chromo1 e RICALCOLA
            _chromo2.set_rest_chromo(sub_new2, distance2);

            
            
        } 


}
ivec System :: search_sub_in_chromo(ivec& sub, Chromosome & _chromo){

        ivec sub_new(sub.size());
        int j=0;

        //facendolo partire da 1 evitiamo l'inutile ricerca del chromo(0)=1=sub(i)
        for(int k=1; k < _n_city; k++ ){
          for(int i=0; i< sub.size(); i++){

                if( _chromo.get_gene(k) == sub(i) ){      
                                                      sub_new(j)=sub(i);
                                                      j++;
                }   
          }

        }

        return sub_new;

}



/*
void System :: finalize(){
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}
*/

int System :: get_ncity(){
  return _n_city;
}

int System :: get_population(){
  return _population;
}

mat System :: get_distance2(){

  return distance2;
}

int System :: get_ngen(){

    return _ngen;
}

void System :: print_terminal(){

  for(int i=0; i< _population ; i++){

      _chromo(i).print(i);

  }
  return;
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
