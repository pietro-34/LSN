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
#include <math.h>
#include "chromosome.h"

using namespace std;

//sets the chromo's dimension and assigns the identical permutation  on vector _g
void Chromosome :: initialize(int n){
  
    _ndim=n;
   _g.resize(_ndim);
  for(int i=0; i<_ndim; i++){ 
   _g(i)=i+1;
   }
   
   return;
}

Chromosome& Chromosome::operator=(const Chromosome& other) {
    if (this != &other) {
        // Copio i membri base
        _ndim = other._ndim;
        _g     = other._g;      // ivec supporta copy assignment
        _loss  = other._loss;   // se hai un membro loss

        // Se avessi puntatori o risorse dinamiche, le gestiresti qui
    }
    return *this;
}

bool Chromosome::operator<(const Chromosome& other) const {
    return this->get_loss() < other.get_loss();
}

void Chromosome :: check(){
    
    //il metodo unique crea un ivec con solo gli elementi unici del vettore 
    ivec v = unique(_g);

    if(v.size() != _ndim ){

        cerr << "Errore: il vettore _g contiene elementi duplicati.\n";
        exit(EXIT_FAILURE);
    }else if(_g(0) !=1 ){

        cerr << "Errore: il vettore _g non ha 1 come primo gene.\n";
        exit(EXIT_FAILURE);

    }else {return;}

}

void Chromosome :: permutation(int i, int j){

    if(i == 0 || j == 0 ){ cout<<"Non si può permutarè il primo gene "<<endl;
                            return;}
    int n=_g(i);
    int m=_g(j);

    _g(i)=m;
    _g(j)=n;

    return;

}

void Chromosome :: shift(int i){   //i è lo shift [1,2,3,4,5]-- shift 1-- [1,5,2,3,4]
    if(i <= 0 || _ndim <= 1) return;

    int len = _ndim - 1;       // parti da posizione 1 in poi
    i = i % len;               // evita overflow
    ivec temp = _g.subvec(1, _ndim - 1);  // copia parte da ruotare

    // shift circolare a destra
    ivec shifted = join_vert(temp.tail(i), temp.head(len - i));
    

    
    // aggiorna il vettore originale
    _g(span(1, _ndim - 1)) = shifted;
}


void Chromosome :: mix(){

        //assegno alla variabile tail, il sottovettore di _g ; non cambio primo gene
        ivec tail = _g.subvec(1, _g.n_elem - 1);

        // lo mescolo
        tail = shuffle(tail);

        // lo rimetto in v (sempre posizioni 1..n_elem-1)
        _g.subvec(1, _g.n_elem - 1) = tail;

return;
}

void Chromosome :: set_loss(mat L){

        int indice_riga, indice_colonna, next;
        double l=0.0;
        //controllo dimensione matrice
        if (L.size() != _ndim *_ndim ){cout << "La matrice non è della dimensione corretta \n"
                                            << "Dimensione attuale "<< sqrt(L.size()) <<"x"<< sqrt(L.size())<<"\n"
                                            << "dimensione corretta "<< _ndim<<"x"<<_ndim<<endl;}

        //devo contare gli elementi della matrice in base alla stringa di codice
        // IMPORTANTE AGGIUNGERE elemento di matrice tra l'ultima componente del cromo e la prima
        for (int k=0; k<_ndim; k++){

            next= (k+1) % _ndim;  //in questo modo conto pure il percorso dall'ultima città alla prima

            indice_riga=_g(k)-1;      //gli indici della matrice partono da 0 ; mentre gli alleli di da 1
                               
            indice_colonna=_g(next)-1;

            if(next == 0){  indice_riga=0; indice_colonna=_g(k)-1;     } //la matrice distance2 è triangolare strettamente superiore


            l += L( indice_riga, indice_colonna);

        }

        _loss=l;
    
return;
}

void Chromosome :: set_rest_chromo(ivec sub, mat L){

    int dimsub= sub.size();
    int cut= _g.size()-dimsub; //questo è l'indice dove è avvenuto il crossing; 
                               //da questo indice fino a _g.size() devo mettere le componenti di sub in quelle di _g
    int j=0;
    for(int i=cut; i<_ndim; i++){
        _g(i)=sub(j);
        j++;
    } 
    this->set_loss(L);
}

ivec Chromosome :: get_rest_chromo(int m){

        int index=_ndim -m;
        ivec rest(index);
        int j=0;
        for(int i=m; i<_ndim; i++){
            rest(j)=_g(i);
            j++;
        }
    return rest;
}

int Chromosome :: get_gene(int k) const{
    return _g(k);
}

double Chromosome :: get_loss() const{
    return _loss;
}

void Chromosome :: print(int i){

    cout<<"Chromosome "<< i <<" ";
    
    for(int k=0; k < _ndim; k++){

            cout<< _g(k) << "\t";
    } 

            cout<< "loss "<< _loss <<endl;




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
