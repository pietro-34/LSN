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
#include <sstream> //libreria ottieni accesso a tre classi fondamentali della Standard Library C++ che trattano le stringhe come flussi (stream), 
                  //esattamente come faresti con cin/cout ma agendo su std::string anziché sul terminale o sui file; lo usiamo per -->istringstream.
using namespace std;
using namespace arma;

void System :: step(){ // Perform a simulation step
  if(_sim_type == 0) this->Verlet();  // Perform a MD step
  else if(_sim_type == 4) {this->move_vmc();} //Add by me; move vmc non accetta argomenti in input perchè agisce solo la particle(0)
  else for(int i=0; i<_npart; i++) this->move(int(_rnd.Rannyu()*_npart)); // Perform a MC step on a randomly choosen particle
  _nattempts += _npart; //update number of attempts performed on the system
  return;
}

void System :: Verlet(){
  double xnew, ynew, znew;
  for(int i=0; i<_npart; i++){ //Force acting on particle i
    _fx(i) = this->Force(i,0);
    _fy(i) = this->Force(i,1);
    _fz(i) = this->Force(i,2);
  }
  for(int i=0; i<_npart; i++){ //Verlet integration scheme
    xnew = this->pbc( 2.0 * _particle(i).getposition(0,true) - _particle(i).getposition(0,false) + _fx(i) * pow(_delta,2), 0);
    ynew = this->pbc( 2.0 * _particle(i).getposition(1,true) - _particle(i).getposition(1,false) + _fy(i) * pow(_delta,2), 1);
    znew = this->pbc( 2.0 * _particle(i).getposition(2,true) - _particle(i).getposition(2,false) + _fz(i) * pow(_delta,2), 2);
    _particle(i).setvelocity(0, this->pbc(xnew - _particle(i).getposition(0,false), 0)/(2.0 * _delta));
    _particle(i).setvelocity(1, this->pbc(ynew - _particle(i).getposition(1,false), 1)/(2.0 * _delta));
    _particle(i).setvelocity(2, this->pbc(znew - _particle(i).getposition(2,false), 2)/(2.0 * _delta));
    _particle(i).acceptmove(); // xold = xnew
    _particle(i).setposition(0, xnew);
    _particle(i).setposition(1, ynew);
    _particle(i).setposition(2, znew);
  }
  _naccepted += _npart;
  return;
}

double System :: Force(int i, int dim){
  double f=0.0, dr;
  vec distance;
  distance.resize(_ndim);
  for (int j=0; j<_npart; j++){
    if(i != j){
      distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
      distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
      distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
      dr = sqrt( dot(distance,distance) );
      if(dr < _r_cut){
        f += distance(dim) * (48.0/pow(dr,14) - 24.0/pow(dr,8));
      }
    }
  }
  return f;
}

void System :: move(int i){ // Propose a MC move for particle i
  if(_sim_type == 3){ //Gibbs sampler for Ising
    int s_left  = _particle(pbc(i-1)).getspin();
    int s_right = _particle(pbc(i+1)).getspin();

    double h_i = _J * (s_left + s_right) + _H;
    double p_up = 1.0 / (1.0 + exp(-2.0 * _beta * h_i));

    if (_rnd.Rannyu() < p_up)
        _particle(i).setspin(+1);
    else
        _particle(i).setspin(-1);
  } else {           // M(RT)^2
    if(_sim_type == 1){       // LJ system
      vec shift(_ndim);       // Will store the proposed translation
      for(int j=0; j<_ndim; j++){
        shift(j) = _rnd.Rannyu(-1.0,1.0) * _delta; // uniform distribution in [-_delta;_delta)
      }
      _particle(i).translate(shift, _side);  //Call the function Particle::translate
      if(this->metro(i)){ //Metropolis acceptance evaluation
        _particle(i).acceptmove();
        _naccepted++;
      } else _particle(i).moveback(); //If translation is rejected, restore the old configuration
    } else {                  // Ising 1D
      if(this->metro(i)){     //Metropolis acceptance evaluation for a spin flip involving spin i
        _particle(i).flip();  //If accepted, the spin i is flipped
        _naccepted++;
      }
    }
  }
  return;
}

//add method move_vmc--> muove la particella che campiona la distrubuzione di probabilità; in questo metodo è come se ci fosse un metro_vmc 
//                        nella logica del simulatore una volta fatta la mossa chiamo metro(i) che mi acceta o meno la kossa, in questo caso non ho implementato un metodo metro_vmc. 

void System::move_vmc() {
    // Proponi una nuova posizione per la particella
    double x_trial = _particle(0).getposition(0, true) + _rnd.Rannyu(-_delta, _delta);

    // Imposta la posizione trial
    _particle(0).setposition(0, x_trial);

    // Calcola il rapporto di probabilità
    double psi2_old = psi2_trial(_particle(0).getposition(0, false)); // x_old
    double psi2_new = psi2_trial(_particle(0).getposition(0, true));  // x_trial

    double acceptance = psi2_new / psi2_old;

    // Accetta o rifiuta la mossa
    if (_rnd.Rannyu() < acceptance) {
        _particle(0).acceptmove();  // xold = x
        _naccepted++;
    } else {
        _particle(0).moveback();    // x = xold
    }
}
bool System :: metro(int i){ // Metropolis algorithm
  bool decision = false;
  double delta_E, acceptance;
  if(_sim_type == 1) delta_E = this->Boltzmann(i,true) - this->Boltzmann(i,false);
  else delta_E = 2.0 * _particle(i).getspin() * 
                 ( _J * (_particle(this->pbc(i-1)).getspin() + _particle(this->pbc(i+1)).getspin() ) + _H );
  acceptance = exp(-_beta*delta_E);
  if(_rnd.Rannyu() < acceptance ) decision = true; //Metropolis acceptance step
  return decision;
}

double System :: Boltzmann(int i, bool xnew){
  double energy_i=0.0;
  double dx, dy, dz, dr;
  for (int j=0; j<_npart; j++){
    if(j != i){
      dx = this->pbc(_particle(i).getposition(0,xnew) - _particle(j).getposition(0,1), 0);
      dy = this->pbc(_particle(i).getposition(1,xnew) - _particle(j).getposition(1,1), 1);
      dz = this->pbc(_particle(i).getposition(2,xnew) - _particle(j).getposition(2,1), 2);
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      if(dr < _r_cut){
        energy_i += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }
  return 4.0 * energy_i;
}

double System :: pbc(double position, int i){ // Enforce periodic boundary conditions
  return position - _side(i) * rint(position / _side(i));
}

int System :: pbc(int i){ // Enforce periodic boundary conditions for spins
  if(i >= _npart) i = i - _npart;
  else if(i < 0)  i = i + _npart;
  return i;
} 

void System :: initialize(){ // Initialize the System object according to the content of the input files in the ../INPUT/ directory

  int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
  ifstream Primes("../INPUT/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();
  int seed[4]; // Read the seed of the RNG
  ifstream Seed("../INPUT/seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  _rnd.SetRandom(seed,p1,p2);

  ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
  couta << "#   N_BLOCK:  ACCEPTANCE:" << endl;
  couta.close();

  ifstream input("../INPUT/input.dat"); // Start reading ../INPUT/input.dat
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat");
  string property;
  double delta;
  // !input.eof crea problemi con l'aggiunta di sim tipe =4 si può usare while(input >> property)
  while (input >> property ){
    //input >> property; da usare dopo while se la condizione è !input.eof
    if( property == "SIMULATION_TYPE" ){
      input >> _sim_type;
      if(_sim_type == 2 || _sim_type == 3){
        input >> _J;
        input >> _H;
      }
      // 1. Dopo aver letto con `input >> _sim_type;` il valore numerico di SIMULATION_TYPE,
      //    ci troviamo ancora sullo stesso rigo, subito dopo il numero, e non abbiamo ancora
      //    consumato il carattere di 'newline'.
      //    Se usassimo un secondo `input >> ...`, rischieremmo di leggere il token successivo
      //    (ad es. “RESTART”) invece di restOfLine.

      // 2. Vogliamo “spazzare via” tutto ciò che rimane sulla riga corrente, fino al '\n',
      //    senza intaccare il contenuto della riga successiva.
      //    Per farlo, usiamo `getline`, che legge fino al prossimo '\n' escluso.
      string restOfLine;
      getline(input, restOfLine);  
      //    — restOfLine ora contiene tutti i caratteri tra la fine del numero e il newline.
      //    — il cursore del file si sposta sulla riga successiva, pronto per la prossima `input >>`.

      // 3. Per poter trattare il contenuto di restOfLine come un piccolo stream di parole,
      //    creiamo un `istringstream`, che espone la stessa interfaccia di un `istream`.
      istringstream iss(restOfLine);

      // 4. Controlliamo se siamo nel caso MD (sim_type==0): solo lì ci interessa vedere se dopo il numero c’è la parola MAXBOLTZ o DIRAC
      if (_sim_type == 0) {
                      string tok;

                      // 5. Proviamo a estrarre il primo token da `iss`.
                      //    — Se `iss >> tok` restituisce true, significa che restOfLine NON era vuoto:
                      //      c’era almeno una parola dopo il numero.
                      //    — Se restituisce false, la riga era vuota (nessun token) e siamo nel caso puro di restart.
                      if (iss >> tok) {
                          // 6a. Se il token è esattamente "MAXBOLTZ" o "DIRAC", allora l’utente
                          //     ha specificato la distribuzione: la salviamo e la scriviamo su output.dat.
                          if (tok == "MAXBOLTZ" || tok == "DIRAC") {
                              _start_velocities = tok;
                              coutf << "STARTING VELOCITIES DISTRIBUTION=  "
                                    << _start_velocities << endl;
                          }
                          // 6b. Se invece c’era un token qualsiasi diverso (ad es. "FOOBAR"), lo interpretiamo come "errore"
                          //     
                          else {
                              cout << "PROBLEM: Unknown starting distribution; MAXBOLTZ and DIRAC are permissed, if RESTART=1 NO starting ditribution" << endl;
                              exit(EXIT_FAILURE);
                          }
                      }
                      // 7. Caso `iss >> tok` è false: non c’era **nessuna** parola dopo il numero,
                      //    perciò scriviamo il messaggio di restart.
                      else {
                          coutf << "INITIALIZE VELOCITIES FROM PREVIOUS SIMULATION" << endl;
                      }
                  }
    
      
      
      if(_sim_type > 4){
        cerr << "PROBLEM: unknown simulation type" << endl;
        exit(EXIT_FAILURE);
      }
      if(_sim_type == 0) coutf << "LJ MOLECULAR DYNAMICS (NVE) SIMULATION"  << endl; 
      else if(_sim_type == 1) coutf << "LJ MONTE CARLO (NVT) SIMULATION"         << endl;
      else if(_sim_type == 2) coutf << "ISING 1D MONTE CARLO (MRT^2) SIMULATION" << endl;
      else if(_sim_type == 3) coutf << "ISING 1D MONTE CARLO (GIBBS) SIMULATION" << endl;
      else if(_sim_type == 4) coutf << "VARIATIONAL QUANTUM MONTECARLO"  << endl; 
      
      
     } //property ora diventa diverso da SIMULATION TIME
     else if( property == "RESTART" ){
      input >> _restart;
    } else if( property == "TEMP" ){
      input >> _temp;
      _beta = 1.0/_temp;
      coutf << "TEMPERATURE= " << _temp << endl;
    } else if( property == "NPART" ){
      input >> _npart;
      _fx.resize(_npart);
      _fy.resize(_npart);
      _fz.resize(_npart);
      _particle.set_size(_npart);
      for(int i=0; i<_npart; i++){ 
        _particle(i).initialize();
       if(_rnd.Rannyu() > 0.5) _particle(i).flip(); // to randomize the spin configuration
      }
      coutf << "NPART= " << _npart << endl;
    } else if( property == "RHO" ){
      input >> _rho;
      _volume = _npart/_rho;
      _side.resize(_ndim);
      _halfside.resize(_ndim);
      double side = pow(_volume, 1.0/3.0);
      for(int i=0; i<_ndim; i++) _side(i) = side;
      _halfside=0.5*_side;
      coutf << "SIDE= ";
      for(int i=0; i<_ndim; i++){
        coutf << setw(12) << _side[i];
      }
      coutf << endl;
    } else if( property == "R_CUT" ){
      input >> _r_cut;
      coutf << "R_CUT= " << _r_cut << endl;
    } else if( property == "DELTA" ){
      input >> delta;
      coutf << "DELTA= " << delta << endl;
      _delta = delta;
    } 
    //add property by me for sim_type =4 (VMC)
      else if(property == "X0"){   //starting position for vmc
        input >> _x0;
        coutf << "X0 = " << _x0 << endl;
        _particle(0).setposition(0, _x0); //inizializzo la posizione della particella (la quale si muove nella distribuzione di probabilità)

      }else if( property == "MU"){ //mean value gaussian autofunction 
        input >> _mu;
        coutf << "MU = " << _mu << endl;

      }else if( property == "SIGMA"){
        input >> _sigma;
        coutf << "SIGMA = " << _sigma << endl;
      
    //
      }else if( property == "NBLOCKS" ){
      input >> _nblocks;
      coutf << "NBLOCKS= " << _nblocks << endl;
    } else if( property == "NSTEPS" ){
      input >> _nsteps;
      coutf << "NSTEPS= " << _nsteps << endl;
    } //SA INPUT
      else if(property == "INVERSE_DT"){
      input >> _inverse_dt;
      coutf << "INVERSE_DT= " << _inverse_dt << endl;
    } else if(property == "DELTA_MU"){
      input >> _delta_mu;
      coutf << "DELTA_MU= " << _delta_mu << endl;
    }else if(property == "DELTA_SIGMA"){
      input >> _delta_sigma;
      coutf << "DELTA_SIGMA= " << _delta_sigma << endl;
    }else if(property == "SA_STEPS"){
      input >> _sa_steps;
      coutf << "SA STEPS= " << _sa_steps << endl;
    }else if( property == "ENDINPUT" ){
      coutf << "Reading input completed!" << endl;
      break;
    } else cerr << "PROBLEM: unknown input" << endl;
  }
  
  input.close();

  if(_sim_type != 4){
  this->read_configuration();}
  if(_sim_type==0) this->initialize_velocities();
  coutf << "System initialized!" << endl;
  coutf.close();
  return;
}

void System :: initialize_velocities(){
  double xold, yold, zold;
  if(_restart){
    ifstream cinf;
    cinf.open("../INPUT/CONFIG/conf-1.xyz");
    if(cinf.is_open()){
      string comment;
      string particle;
      int ncoord;
      cinf >> ncoord;
      if (ncoord != _npart){
        cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
        exit(EXIT_FAILURE);
      }
      cinf >> comment;
      for(int i=0; i<_npart; i++){
        cinf >> particle >> xold >> yold >> zold; // units of coordinates in conf.xyz is _side
        _particle(i).setpositold(0, this->pbc(_side(0)*xold, 0));
        _particle(i).setpositold(1, this->pbc(_side(1)*yold, 1));
        _particle(i).setpositold(2, this->pbc(_side(2)*zold, 2));
      }
    } else cerr << "PROBLEM: Unable to open INPUT file conf-1.xyz"<< endl;
    cinf.close();
  } else if(_start_velocities == "MAXBOLTZ"){   
   // Inizializza le velocità delle particelle secondo la distribuzione di Maxwell-Boltzmann e imposta le posizioni precedenti per l'algoritmo di Verlet

  vec vx(_npart), vy(_npart), vz(_npart); // Vettori per le componenti delle velocità
  vec sumv(_ndim);                         // Somma delle componenti per rimuovere il centro di massa
  sumv.zeros();                            // Inizializza a zero

  // 1. Estrazione delle velocità gaussiane e somma per il drift
  for (int i = 0; i < _npart; i++) {
      // Estrae le tre componenti della velocità da una gaussiana di media 0 e varianza proporzionale alla temperatura
      vx(i) = _rnd.Gauss(0., sqrt(_temp));
      vy(i) = _rnd.Gauss(0., sqrt(_temp));
      vz(i) = _rnd.Gauss(0., sqrt(_temp));

      // Somma delle componenti per calcolare il momento totale
      sumv(0) += vx(i);
      sumv(1) += vy(i);
      sumv(2) += vz(i);
  }

  // 2. Calcolo della velocità media (drift del centro di massa)
  for (int idim = 0; idim < _ndim; idim++)
      sumv(idim) /= double(_npart);

  // 3. Rimozione del drift e calcolo dell'energia cinetica media
  double sumv2 = 0.0, scalef;
  for (int i = 0; i < _npart; i++) {
      // Sottrae la velocità media per rimuovere il drift del centro di massa
      vx(i) -= sumv(0);
      vy(i) -= sumv(1);
      vz(i) -= sumv(2);

      // Somma delle velocità quadratiche per calcolare la temperatura effettiva
      sumv2 += vx(i) * vx(i) + vy(i) * vy(i) + vz(i) * vz(i);
  }
  sumv2 /= double(_npart); // Energia cinetica media per particella

  // Calcola il fattore di scala per correggere la temperatura
  scalef = sqrt(3.0 * _temp / sumv2);

  // 4. Applica il fattore di scala e assegna le velocità alle particelle
  for (int i = 0; i < _npart; i++) {
      _particle(i).setvelocity(0, vx(i) * scalef);
      _particle(i).setvelocity(1, vy(i) * scalef);
      _particle(i).setvelocity(2, vz(i) * scalef);
  }

  // 5. Calcola le posizioni precedenti per il primo passo dell'algoritmo di Verlet
  for (int i = 0; i < _npart; i++) {
      xold = this->pbc(_particle(i).getposition(0, true) - _particle(i).getvelocity(0) * _delta, 0);
      yold = this->pbc(_particle(i).getposition(1, true) - _particle(i).getvelocity(1) * _delta, 1);
      zold = this->pbc(_particle(i).getposition(2, true) - _particle(i).getvelocity(2) * _delta, 2);

      _particle(i).setpositold(0, xold);
      _particle(i).setpositold(1, yold);
      _particle(i).setpositold(2, zold);
  }

    }
  else if( _start_velocities== "DIRAC"){

    vec vx(_npart), vy(_npart), vz(_npart);
vec sumv(_ndim);
sumv.zeros();

// 1. Imposta modulo costante e direzioni casuali uniformi
double v0 = sqrt(3.0 * _temp);  // modulo fisso delle velocità per dare la temperatura corretta

for (int i = 0; i < _npart; i++) {
    // Genera direzione uniforme su sfera
    double theta = acos(1.0 - 2.0 * _rnd.Rannyu());     // θ ∈ [0, π], distribuzione corretta per uniformità
    double phi   = 2.0 * M_PI * _rnd.Rannyu();          // φ ∈ [0, 2π]

    // Conversione in componenti cartesiane
    vx(i) = v0 * sin(theta) * cos(phi);
    vy(i) = v0 * sin(theta) * sin(phi);
    vz(i) = v0 * cos(theta);

    // Somma per il drift
    sumv(0) += vx(i);
    sumv(1) += vy(i);
    sumv(2) += vz(i);
}

// 2. Calcola la velocità media
for (int idim = 0; idim < _ndim; idim++)
    sumv(idim) /= double(_npart);

// 3. Rimuove il momento totale
double sumv2 = 0.0, scalef;
for (int i = 0; i < _npart; i++) {
    vx(i) -= sumv(0);
    vy(i) -= sumv(1);
    vz(i) -= sumv(2);

    sumv2 += vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i);
}
sumv2 /= double(_npart);

// 4. Rescale per ottenere la temperatura desiderata
scalef = sqrt(3.0 * _temp / sumv2);
for (int i = 0; i < _npart; i++) {
    _particle(i).setvelocity(0, vx(i) * scalef);
    _particle(i).setvelocity(1, vy(i) * scalef);
    _particle(i).setvelocity(2, vz(i) * scalef);
}

// 5. Calcolo delle posizioni precedenti per Verlet
for (int i = 0; i < _npart; i++) {
    xold = this->pbc(_particle(i).getposition(0,true) - _particle(i).getvelocity(0)*_delta, 0);
    yold = this->pbc(_particle(i).getposition(1,true) - _particle(i).getvelocity(1)*_delta, 1);
    zold = this->pbc(_particle(i).getposition(2,true) - _particle(i).getvelocity(2)*_delta, 2);

    _particle(i).setpositold(0, xold);
    _particle(i).setpositold(1, yold);
    _particle(i).setpositold(2, zold);
}



  } else{ cerr << "PROBLEM: unknown starting  velocities distribution"; //forse da problemi con il restart;
          exit(EXIT_FAILURE);
        }
  //posso scrivere un pezzo di codice che mi stampa la distribuzione iniziale di velocità appena creata

  return;
}

void System :: initialize_properties(){ // Initialize data members used for measurement of properties

  string property;
  int index_property = 0;
  _nprop = 0;

  _measure_penergy  = false; //Defining which properties will be measured
  _measure_kenergy  = false;
  _measure_tenergy  = false;
  _measure_pressure = false;
  _measure_gofr     = false;
  _measure_pofv     = false; 
  _measure_magnet   = false;
  _measure_cv       = false;
  _measure_chi      = false;
  _measure_instant_temp = false; //ADD BY ME; per misurare la temperatura istantanea;
  _measure_instant_magnet= false; //ADD BY ME; per misurare la magnetizazzione istantanea;
  _measure_instant_pot= false;    //ADD BY ME; 
  _measure_local_energy= false;   //ADD BY ME;
  
  ifstream input("../INPUT/properties.dat");
  if (input.is_open()){
    while ( !input.eof() ){
      input >> property;
      if( property == "POTENTIAL_ENERGY" ){
        _vtail = (8.0/3.0)*M_PI*_rho*((1.0/(3.0*pow(_r_cut,9)))-(1.0/pow(_r_cut,3))); //FIXED IN EXERCISE 7  (8/3 *Pi* rho*epsilon)(sigma^12/3 r_c^9 - sigma^6/r_c^3)--natural unit (8/3)*PI*rho(1/3*rc^9-1/r_c^3) 
        //ho spostato _vtail all'inzio di questo if in modo da stamparlo sulla prima riga dell'output pressione
        ofstream coutp("../OUTPUT/potential_energy.dat");
        coutp << "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:      V_TAIL="<< _vtail << endl;
        coutp.close();
        _nprop++;
        _index_penergy = index_property;
        _measure_penergy = true;
        index_property++;
         
      }else if( property== "ISTANT_POTENTIAL"){ //non aggiungo l'index così non recupera le misure da _measure e stampa subito
        ofstream coutipe("../OUTPUT/istant_potential.dat");
        coutipe << "#     ISTANT_POTENTIAL:" << endl;
        coutipe.close();
        _nprop++;
        _measure_instant_pot = true;

        }else if( property == "KINETIC_ENERGY" ){
        ofstream coutk("../OUTPUT/kinetic_energy.dat");
        coutk << "#     BLOCK:   ACTUAL_KE:    KE_AVE:      ERROR:" << endl;
        coutk.close();
        _nprop++;
        _measure_kenergy = true;
        _index_kenergy = index_property;
        index_property++;
      } else if( property == "TOTAL_ENERGY" ){
        ofstream coutt("../OUTPUT/total_energy.dat");
        coutt << "#     BLOCK:   ACTUAL_TE:    TE_AVE:      ERROR:" << endl;
        coutt.close();
        _nprop++;
        _measure_tenergy = true;
        _index_tenergy = index_property;
        index_property++;
      } else if( property == "TEMPERATURE" ){
        ofstream coutte("../OUTPUT/temperature.dat");
        coutte << "#     BLOCK:   ACTUAL_T:     T_AVE:       ERROR:" << endl;
        coutte.close();
        _nprop++;
        _measure_temp = true;
        _index_temp = index_property;
        index_property++;
      }else if( property== "ISTANT_TEMP"){ //non aggiungo l'index così non recupera le misure da _measure e stampa subito
        ofstream coutite("../OUTPUT/istant_temperature.dat");
        coutite << "#     ISTANT_TEMP:" << endl;
        coutite.close();
        _nprop++;
        _measure_instant_temp = true;
        //_index_instant_temp = index_property;
        //index_property++;
    
       } else if( property == "PRESSURE" ){
         _ptail = (32.0/3.0)*M_PI*_rho*((1.0/(3.0*pow(_r_cut,9)))-(1.0/(2*pow(_r_cut,3)))); //FIXED IN EXERCISE 7---> natural units  32/3*Pi*rho(1/3*r_c^9 - 1/2*r_c^3)
        //stesso commento del _vatil;
        ofstream coutpr("../OUTPUT/pressure.dat");
        coutpr << "#     BLOCK:   ACTUAL_P:     P_AVE:       ERROR:       P_TAIL="<<_ptail << endl;
        coutpr.close();
        _nprop++;
        _measure_pressure = true;
        _index_pressure = index_property;
        index_property++;
        
      } else if( property == "GOFR" ){
        ofstream coutgr("../OUTPUT/gofr.dat");
        coutgr << "# BLOCK:    DISTANCE:    ACTUAL_GOFR:    AVE_GOFR:    ERROR:" << endl;
        coutgr.close();
        input>>_n_bins;
        _nprop+=_n_bins;
        _bin_size = (_halfside.min() )/(double)_n_bins;
        _measure_gofr = true;
        _index_gofr = index_property;
        index_property+= _n_bins;
      } else if( property == "POFV" ){
        if(_sim_type > 0){
          cerr << "PROBLEM: DOES NOT MAKE SENSE COMPUTING POFV FOR MC" << endl;
          exit(EXIT_FAILURE);
        }
        ofstream coutpv("../OUTPUT/pofv.dat");
        coutpv << "# BLOCK:    VELOCITY:    BLOCK_POFV:    AVE_POFV:    ERROR:"  << endl;//modificato per una stampa più completa
        coutpv.close();
        input>>_n_bins_v;
        _nprop += _n_bins_v;
        _bin_size_v = 4.0*_temp/(double)_n_bins_v; // (FIXED IN EXERCISE 4)
        _measure_pofv = true;
        _index_pofv = index_property;
        index_property += _n_bins_v;
      } else if( property == "MAGNETIZATION" ){
        ofstream coutpr("../OUTPUT/magnetization.dat");
        coutpr << "#     BLOCK:   ACTUAL_M:     M_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_magnet = true;
        _index_magnet = index_property;
        index_property++;
      } else if( property== "ISTANT_MAGNETIZATION"){ //non aggiungo l'index così non recupera le misure da _measure e stampa subito con measure()
        ofstream coutim("../OUTPUT/istant_magnetization.dat");
        coutim << "#     ISTANT_MAGNET:" << endl;
        coutim.close();
        _nprop++;
        _measure_instant_magnet = true;
        
      } else if( property == "SPECIFIC_HEAT" ){
        ofstream coutpr("../OUTPUT/specific_heat.dat");
        coutpr << "#     BLOCK:   ACTUAL_CV:    CV_AVE:      ERROR:" << endl;
        coutpr.close();
        
        _measure_cv = true;
        _index_cv = index_property; //indice per il cv
        index_property++;

        _index_cvH2= index_property; //indice per H^2
         index_property++;
        _nprop += 2;
      } else if( property == "SUSCEPTIBILITY" ){
        ofstream coutpr("../OUTPUT/susceptibility.dat");
        coutpr << "#     BLOCK:   ACTUAL_X:     X_AVE:       ERROR:" << endl;
        coutpr.close();
        _nprop++;
        _measure_chi = true;
        _index_chi = index_property;
        index_property++;
      } else if(property== "LOCAL_ENERGY"){
        ofstream coutle("../OUTPUT/local_energy.dat");
        coutle << "#     BLOCK:   ACTUAL_LE:     LE_AVE:       ERROR:" << endl;
        coutle.close();
        _nprop++;
        _measure_local_energy = true;
        _index_local_energy = index_property;
        index_property++;

      }else if( property == "ENDPROPERTIES" ){
        ofstream coutf;
        coutf.open("../OUTPUT/output.dat",ios::app);
        coutf << "Reading properties completed!" << endl;
        coutf.close();
        break;
      } else cerr << "PROBLEM: unknown property" << endl;
    }
    input.close();
  } else cerr << "PROBLEM: Unable to open properties.dat" << endl;

  // according to the number of properties, resize the vectors _measurement,_average,_block_av,_global_av,_global_av2
  _measurement.resize(_nprop);
  _average.resize(_nprop);
  _block_av.resize(_nprop);
  _global_av.resize(_nprop);
  _global_av2.resize(_nprop);
  _average.zeros();
  _global_av.zeros();
  _global_av2.zeros();
  _nattempts = 0;
  _naccepted = 0;
  return;
}


//scrivo funzione per stampare distribuzione di velocità iniziale; la scrivo qui perchè nel main 
//uso initialize() e initialize_properties()--> qui leggo il bin di velocità e se voglio effettivamente misurare la pofv

void System::write_starting_velocities_distribution() {
    // Calcola e salva la distribuzione iniziale dei moduli di velocità

    // Definizione parametri dell'istogramma
    double v_max = 4.0 * sqrt(_temp); // velocità massima prevista (copre la distribuzione)
    _bin_size_v = v_max / _n_bins_v;  // larghezza del bin

    vector<int> hist(_n_bins_v, 0); // istogramma inizializzato a 0

    // Costruzione istogramma dei moduli delle velocità
    for (int i = 0; i < _npart; i++) {
        double vx = _particle(i).getvelocity(0);
        double vy = _particle(i).getvelocity(1);
        double vz = _particle(i).getvelocity(2);

        double v_mod = sqrt(vx * vx + vy * vy + vz * vz); // modulo velocità

        int bin = int(v_mod / _bin_size_v);
        if (bin < _n_bins_v) {
            hist[bin]++;
        }
    }

    // Scrittura su file
    ofstream out("../OUTPUT/starting_distribution_velocities.dat");
    out << "# v_bin_center    density" << std::endl;

    for (int i = 0; i < _n_bins_v; i++) {
        double v_center = (i + 0.5) * _bin_size_v;
        double density = double(hist[i]) / (_npart * _bin_size_v); // normalizzazione come densità
        out << v_center << "\t" << density << std::endl;
    }

    out.close();
}





void System :: finalize(){
  this->write_configuration();
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
  return;
}

// Write final configurations as .xyz files in the directory ../OUTPUT/CONFIG/
void System :: write_configuration(){
  ofstream coutf;
  if(_sim_type < 2){
    coutf.open("../OUTPUT/CONFIG/config.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  " 
              << setprecision(17) << _particle(i).getposition(0,true)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,true)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,true)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
    coutf.close();
    coutf.open("../OUTPUT/CONFIG/conf-1.xyz");
    if(coutf.is_open()){
      coutf << _npart << endl;
      coutf << "#Comment!" << endl;
      for(int i=0; i<_npart; i++){
        coutf << "LJ" << "  "
              << setprecision(17) << _particle(i).getposition(0,false)/_side(0) << "   " // x
              << setprecision(17) << _particle(i).getposition(1,false)/_side(1) << "   " // y
              << setprecision(17) << _particle(i).getposition(2,false)/_side(2) << endl; // z
      }
    } else cerr << "PROBLEM: Unable to open conf-1.xyz" << endl;
    coutf.close();
  } else {
    coutf.open("../OUTPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++) coutf << _particle(i).getspin() << " ";
    coutf.close();
  }
  return;
}

// Write configuration nconf as a .xyz file in directory ../OUTPUT/CONFIG/
void System :: write_XYZ(int nconf){
  ofstream coutf;
  coutf.open("../OUTPUT/CONFIG/config_" + to_string(nconf) + ".xyz");
  if(coutf.is_open()){
    coutf << _npart << endl;
    coutf << "#Comment!" << endl;
    for(int i=0; i<_npart; i++){
      coutf << "LJ" << "  " 
            << setw(16) << _particle(i).getposition(0,true)          // x
            << setw(16) << _particle(i).getposition(1,true)          // y
            << setw(16) << _particle(i).getposition(2,true) << endl; // z
    }
  } else cerr << "PROBLEM: Unable to open config.xyz" << endl;
  coutf.close();
  return;
}

// Read configuration from a .xyz file in directory ../OUTPUT/CONFIG/
void System :: read_configuration(){
  ifstream cinf;
  cinf.open("../INPUT/CONFIG/config.xyz");
  if(cinf.is_open()){
    string comment;
    string particle;
    double x, y, z;
    int ncoord;
    cinf >> ncoord;
    if (ncoord != _npart){
      cerr << "PROBLEM: conflicting number of coordinates in input.dat & config.xyz not match!" << endl;
      exit(EXIT_FAILURE);
    }
    cinf >> comment;
    for(int i=0; i<_npart; i++){
      cinf >> particle >> x >> y >> z; // units of coordinates in conf.xyz is _side
      _particle(i).setposition(0, this->pbc(_side(0)*x, 0));
      _particle(i).setposition(1, this->pbc(_side(1)*y, 1));
      _particle(i).setposition(2, this->pbc(_side(2)*z, 2));
      _particle(i).acceptmove(); // _x_old = _x_new
    }
  } else cerr << "PROBLEM: Unable to open INPUT file config.xyz"<< endl;
  cinf.close();
  if(_restart and _sim_type > 1){
    int spin;
    cinf.open("../INPUT/CONFIG/config.spin");
    for(int i=0; i<_npart; i++){
      cinf >> spin;
      _particle(i).setspin(spin);
    }
    cinf.close();
  }
  return;
}

void System :: block_reset(int blk){ // Reset block accumulators to zero
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_av.zeros();
  return;
}

void System :: measure(){ // Measure properties
  _measurement.zeros();
  // POTENTIAL ENERGY, VIRIAL, GOFR ///////////////////////////////////////////
  int bin;
  vec distance;
  distance.resize(_ndim);
  double penergy_temp=0.0, dr; // temporary accumulator for potential energy
  double kenergy_temp=0.0; // temporary accumulator for kinetic energy
  double tenergy_temp=0.0;
  double magnetization=0.0;
  double chi=0.0; //add by me
  double cv=0.0;  //add by me
  double virial=0.0;
  if (_measure_penergy or _measure_pressure or _measure_gofr) {
    for (int i=0; i<_npart-1; i++){
      for (int j=i+1; j<_npart; j++){
        distance(0) = this->pbc( _particle(i).getposition(0,true) - _particle(j).getposition(0,true), 0);
        distance(1) = this->pbc( _particle(i).getposition(1,true) - _particle(j).getposition(1,true), 1);
        distance(2) = this->pbc( _particle(i).getposition(2,true) - _particle(j).getposition(2,true), 2);
        dr = sqrt( dot(distance,distance) );
        // GOFR ...  (FIXED IN EXERCISE 7)
        if(_measure_gofr) {
                              int bin = int(dr/_bin_size);
                                  if(bin < _n_bins) {
                                                // sommo 2.0 perché conto sia i→j che j→i
                                                _measurement[_index_gofr + bin] += 2.0;
                                            }
                                        }

            
          
        
        if(dr < _r_cut){
          if(_measure_penergy)  penergy_temp += 1.0/pow(dr,12) - 1.0/pow(dr,6); // POTENTIAL ENERGY
          if(_measure_pressure) virial       += 1.0/pow(dr,12) - 0.5/pow(dr,6); // PRESSURE
        }
      }
    }
  }
  // POFV (FIXED IN EXERCISE 4)
if(_measure_pofv){ 
  // Azzera distribuzione istantanea
 /* for(int i=0; i < _n_bins_v; i++)
      _measurement[_index_pofv + i] = 0.0; inutile perchè c'è measure.zeros() all'inizio di measure()*/

  double v_mod;
  for(int i=0; i < _npart; i++){
      v_mod = sqrt( pow(_particle[i].getvelocity(0),2) +
                    pow(_particle[i].getvelocity(1),2) +
                    pow(_particle[i].getvelocity(2),2) );

      // Determina il bin di velocità corrispondente a v_mod
      int bin = int(v_mod / _bin_size_v);
      
      // Incrementa il conteggio del bin corrispondente, se bin < _n_bins_v cioè se bin è più piccolo del  numero di bin in cui suddivido [0, vmax]
      if(bin < _n_bins_v)
          _measurement[_index_pofv + bin] += 1.0; // incrementa bin corretto
  }

  // Normalizzazione della distribuzione istantanea
  for(int i=0; i < _n_bins_v; i++)
      _measurement[_index_pofv + i] /= (_npart * _bin_size_v);


}

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    penergy_temp = _vtail + 4.0 * penergy_temp / double(_npart);
    _measurement(_index_penergy) = penergy_temp;
  }
  // ISTANT POTENTIAL ENERGY /////////////////////////////////////////////////// add by me 
  if(_measure_instant_pot) {

      ofstream coutipe("../OUTPUT/istant_potential.dat",ios::app);   
      coutipe <<  _measurement(_index_penergy)<< endl;
      coutipe.close();
  }

  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    for (int i=0; i<_npart; i++) kenergy_temp += 0.5 * dot( _particle(i).getvelocity() , _particle(i).getvelocity() ); 
    kenergy_temp /= double(_npart);
    _measurement(_index_kenergy) = kenergy_temp;
  }
  // TOTAL ENERGY (kinetic+potential) ////////////////////////////////////////// (total energy per particle)
  if (_measure_tenergy){
    if (_sim_type < 2) _measurement(_index_tenergy) = kenergy_temp + penergy_temp;
    else { 
      double s_i, s_j;
      for (int i=0; i<_npart; i++){
        s_i = double(_particle(i).getspin());
        s_j = double(_particle(this->pbc(i+1)).getspin());
        tenergy_temp += - _J * s_i * s_j - 0.5 * _H * (s_i + s_j);
      }
      tenergy_temp /= double(_npart);
      _measurement(_index_tenergy) = tenergy_temp;
    }
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp and _measure_kenergy) _measurement(_index_temp) = (2.0/3.0) * kenergy_temp;
  // ISTANT TEMPERATURE ////////////////////////////////////////////////////////
  // se _measure_instant_temp=true allora faccio stampare direttamente la temperatura istantanea
  if(_measure_instant_temp) {

      ofstream coutf("../OUTPUT/istant_temperature.dat",ios::app);   
      coutf <<  _measurement(_index_temp)<< endl;
      coutf.close();


  }

  // PRESSURE ////////////////////////////////////////////////////////////////// 
  //ATTENZIONE: Le simulazioni MC non inizializzano le velocità quindi non viene calcolato il termine di gas ideale della pressione
  //            quindi sostituisco _rho * (2.0/3.0) * kenergy_temp con _rho * _temp
  if (_measure_pressure){  if(_sim_type==1){_measurement[_index_pressure] = _rho *_temp + (_ptail*_npart + 48.0*virial/3.0)/_volume;}
                            else{
                            _measurement[_index_pressure] = _rho * (2.0/3.0) * kenergy_temp + (_ptail*_npart + 48.0*virial/3.0)/_volume;}
  }
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  /*
        M=  < sum_i=0^npart s_i >]
  */
  if(_measure_magnet){
        for(int i=0; i<_npart; i++){   magnetization += _particle(i).getspin(); }

      _measurement(_index_magnet)=  double(magnetization/_npart) ;

  }

 

  // INSTANT MAGNETIZATION /////////////////////////////////////////////////////
  if(_measure_instant_magnet) {

      ofstream coutf("../OUTPUT/istant_magnetization.dat",ios::app);   
      coutf <<  _measurement(_index_magnet)<< endl;
      coutf.close();


  }



  // SPECIFIC HEAT /////////////////////////////////////////////////////////////
  // Cv= k_b * beta^2 (<H^2> - <H>^2)= k_b * beta^2 (sigma_H)^2
  // posso calcolare la deviazione standard 
  if(_measure_cv){
                    _measurement(_index_cvH2)= tenergy_temp * tenergy_temp ;

  }



  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  // chi= beta* [<(sum s_i)^2>] 
  if(_measure_chi){
      
      for(int i=0; i<_npart; i++){   chi += _particle(i).getspin(); }

          chi=chi*chi;
        //beta lo moltiplico qui e non in avarage, beta posso portarlo dentro e fuori la media essendo costante;
      _measurement(_index_chi)=   _beta * chi /double(_npart) ;
    }
  //LOCAL ENERGY ////////////////////////////////////////////////////////////

  if(_measure_local_energy){
                                double z =_particle(0).getposition(0, true);
                               _measurement(_index_local_energy)= local_energy(z);
  }
  
  _block_av += _measurement; //Update block accumulators

  return;
}

void System :: averages(int blk){

  ofstream coutf;
  double average, sum_average, sum_ave2;

  _average     = _block_av / double(_nsteps);
  _global_av  += _average;                        //accumula le medie di blocco
  _global_av2 += _average % _average; // % -> element-wise multiplication      //accumula i valori medi al quadrato

  // POTENTIAL ENERGY //////////////////////////////////////////////////////////
  if (_measure_penergy){
    coutf.open("../OUTPUT/potential_energy.dat",ios::app);
    average  = _average(_index_penergy);
    sum_average = _global_av(_index_penergy);
    sum_ave2 = _global_av2(_index_penergy);
    coutf << setw(12) << blk 
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // KINETIC ENERGY ////////////////////////////////////////////////////////////
  if (_measure_kenergy){
    coutf.open("../OUTPUT/kinetic_energy.dat",ios::app);
    average  = _average(_index_kenergy);
    sum_average = _global_av(_index_kenergy);
    sum_ave2 = _global_av2(_index_kenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TOTAL ENERGY //////////////////////////////////////////////////////////////
  if (_measure_tenergy){
    coutf.open("../OUTPUT/total_energy.dat",ios::app);
    average  = _average(_index_tenergy);
    sum_average = _global_av(_index_tenergy);
    sum_ave2 = _global_av2(_index_tenergy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // TEMPERATURE ///////////////////////////////////////////////////////////////
  if (_measure_temp){
    coutf.open("../OUTPUT/temperature.dat",ios::app);
    average  = _average(_index_temp);
    sum_average = _global_av(_index_temp);
    sum_ave2 = _global_av2(_index_temp);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // PRESSURE //////////////////////////////////////////////////////////////////
  if (_measure_pressure){
    coutf.open("../OUTPUT/pressure.dat",ios::app);
    average  = _average(_index_pressure);
    sum_average = _global_av(_index_pressure);
    sum_ave2 = _global_av2(_index_pressure);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  // POFV //////////////////////////////////////////////////////////////////////
if(_measure_pofv){

   
    ofstream coutf;
    coutf.open("../OUTPUT/pofv.dat", ios::app);

    // Cicla su tutti i bin della distribuzione di velocità
    for(int i = 0; i < _n_bins_v; i++){

        // Calcola l'indice assoluto nel vettore delle misure "_measurement"
        // _index_pofv è l'indice iniziale che indica il punto da cui partono i bin della distribuzione velocità nel vettore _measurement
        int idx = _index_pofv + i;

        // Estrae la media del singolo blocco per questo bin (già calcolata prima come _block_av[idx]/_nsteps)
        average = _average(idx);

        // Estrae le somme cumulative (globali) che sono state aggiornate poco sopra tramite:
        // _global_av  += _average;
        // _global_av2 += _average % _average;
        sum_average = _global_av(idx);
        sum_ave2    = _global_av2(idx);

        // Scrittura risultati su file:
        coutf << setw(12) << blk                             // numero del blocco attuale
              << setw(12) << (i + 0.5) * _bin_size_v         // centro del bin della velocità
              << setw(12) << average                         // media del blocco corrente
              << setw(12) << sum_average / double(blk)       // media progressiva fino al blocco corrente
              << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;  // errore statistico
    }

    // Lascia una riga vuota tra i blocchi per chiarezza del file di output
    coutf << endl;

    // Chiude il file aperto per la distribuzione
    coutf.close();
}

  

  // GOFR //////////////////////////////////////////////////////////////////////
  if(_measure_gofr){
    coutf.open("../OUTPUT/gofr.dat", ios::app);
    
  if(_nblocks-blk==0){ //stampo solo l'ultimo blocco
    double rho = double(_npart)/_volume;    
    for(int i=0; i<_n_bins; ++i){
        int idx = _index_gofr + i;

        // prendo i valori calcolati sopra
        average    = _average(idx);
        sum_average= _global_av(idx);
        sum_ave2   = _global_av2(idx);

        double r   = (i+0.5)*_bin_size;
        double dr  = _bin_size;
        double dV  = (4.0/3.0)*M_PI*(pow(r+dr,3)-pow(r,3));
        double norm= rho*_npart*dV;

        double g_block = average    / norm;
        double g_prog  = (sum_average/blk) / norm;
        double err     = this->error(sum_average, sum_ave2, blk) / norm;
        
        coutf << setw(12)<< blk
              << setw(12)<< r
              << setw(12)<< g_block
              << setw(12)<< g_prog
              << setw(12)<< err  << "\n";
    }
    coutf << "\n";
    coutf.close();
  }
}

  
  // MAGNETIZATION /////////////////////////////////////////////////////////////
  if (_measure_magnet){
    coutf.open("../OUTPUT/magnetization.dat",ios::app);
    average  = _average(_index_magnet);
    sum_average = _global_av(_index_magnet);
    sum_ave2 = _global_av2(_index_magnet);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  
  

   // SPECIFIC HEAT PER PARTICLE ////////////////////////////////////////////////

   if (_measure_cv){
    // Calcola le medie nel blocco
    
    double H2_blk = _average(_index_cvH2);         // <H^2> blocco
    double H_blk  = _average(_index_tenergy);    // <H> blocco
    _average(_index_cv) = _npart * _beta * _beta * (H2_blk - H_blk * H_blk);  //c_v del blocco


     _global_av(_index_cv) += _average(_index_cv);
     _global_av2(_index_cv) += _average(_index_cv) * _average(_index_cv) ;
    
    ofstream coutf("../OUTPUT/specific_heat.dat", ios::app);
    coutf << setw(12) << blk
          << setw(12) << _average(_index_cv)
          << setw(12) << _global_av(_index_cv)/double(blk)
          << setw(12) <<  this->error(_global_av(_index_cv), _global_av2(_index_cv), blk) << endl;
    coutf.close();
}




  
  // SUSCEPTIBILITY ////////////////////////////////////////////////////////////
  if(_measure_chi){
    coutf.open("../OUTPUT/susceptibility.dat",ios::app);
    average  = _average(_index_chi);
    sum_average = _global_av(_index_chi);
    sum_ave2 = _global_av2(_index_chi);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  //LOCAL ENERGY ////////////////////////////////////////////////////////////
  if(_measure_local_energy){
    coutf.open("../OUTPUT/local_energy.dat",ios::app);
    average  = _average(_index_local_energy);
    sum_average = _global_av(_index_local_energy);
    sum_ave2 = _global_av2(_index_local_energy);
    coutf << setw(12) << blk
          << setw(12) << average
          << setw(12) << sum_average/double(blk)
          << setw(12) << this->error(sum_average, sum_ave2, blk) << endl;
    coutf.close();
  }
  
  // ACCEPTANCE ////////////////////////////////////////////////////////////////
  double fraction;
  coutf.open("../OUTPUT/acceptance.dat",ios::app);
  if(_nattempts > 0) fraction = double(_naccepted)/double(_nattempts);
  else fraction = 0.0; 
  coutf << setw(12) << blk << setw(12) << fraction << endl;
  coutf.close();
  
  return;
}

double System :: error(double acc, double acc2, int blk){
  if(blk <= 1) return 0.0;
  else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
}

int System :: get_nbl(){
  return _nblocks;
}

int System :: get_nsteps(){
  return _nsteps;
}

void System :: spin(){
      for(int i=0; i< _npart  ; i++ ){
            cout<< _particle(i).getspin() <<endl;
    };

}

//implemento funzione psi_trial //modulo quadro di psi_t=|exp(-(x+mu)^2/2sigma^2) + exp(-(x-mu)^2/2sigma^2|^2

double System :: psi2_trial(double x){

                double g1=exp((-1.0) *(x - _mu)*(x - _mu)/(2* _sigma * _sigma));
                double g2=exp((-1.0) * (x + _mu)*(x + _mu)/(2* _sigma * _sigma));
                
                return pow((g1+g2), 2);
}
//implemnto funzione E_local; sarebbe H[psi_t]/psi_t= V(x)+ d^2/dx^2 psi_t= x^4 - 5/2 x^2+

double System :: local_energy(double x){

                double g1=exp((-1.0) *(x - _mu)*(x - _mu)/(2* _sigma * _sigma));
                double g2=exp((-1.0) * (x + _mu)*(x + _mu)/(2* _sigma * _sigma));
                double pot=pow(x,4)- 2.5 * pow(x,2);

                double local_energy=pot-0.5 *( ( (x - _mu)*(x - _mu) * g1 +(x + _mu)*(x + _mu) * g2 )/( pow( _sigma,4) *(g1+g2)) - 1.0/pow(_sigma, 2 )  );
                return local_energy;
}
//metodo get temperature;
double System :: get_temp(){
                              return _temp;
} 
//metodo get inverse_dt:
double System :: get_inverse_dt(){
                                    return _inverse_dt;
}
int System :: get_sa_steps(){
                                return _sa_steps;
}
//FUNZIONE DI DEBUG:
void System::print_position(int blk) {
    double x = _particle(0).getposition(0, true);
    cout << "Blocco " << blk << ": x = " << x << endl;
}
//esegue calcolo vmc per SA


 double System::run_vmc(double mu, double sigma) {
    block_reset(0);
   _global_av.zeros();
    _global_av2.zeros();
    _block_av.zeros();
    _nattempts = 0;
    _naccepted = 0;

    // Salvo le posizioni correnti
    double x_saved    = _particle(0).getposition(0, true);
    double xold_saved = _particle(0).getposition(0, false);

    // Setto i parametri della trial wavefunction
    _mu = mu;
    _sigma = sigma;

    // Simulazione VMC
    for (int i = 0; i < _nblocks; i++) {
        for (int j = 0; j < _nsteps; j++) {
            step();
            measure();
        }
        averages(i + 1);
        block_reset(i + 1);
    }

    // Ripristino la configurazione iniziale
    _particle(0).setposition(0, x_saved);
    _particle(0).setpositold(0, xold_saved);

    

    return _global_av(_index_local_energy)/double(_nblocks);
}
//esegue SA su mu e sigma
void System::simulated_annealing(double T) {
    // Parametri di SA
    const int n_steps = get_sa_steps();
    const double beta = 1.0 / T;

    // Punto di partenza
    double best_mu     = _mu;
    double best_sigma  = _sigma;
    double best_energy = run_vmc(best_mu, best_sigma);
    // Calcolo dell'errore statistico dopo il run VMC iniziale
    double best_err    = error(_global_av(_index_local_energy),
                               _global_av2(_index_local_energy),
                               _nblocks);

    int accepted = 0;

    // Ciclo di Metropolis per SA
    for (int i = 0; i < n_steps; ++i) {
        // Proposta simmetrica per mu e sigma
        double mu_trial    = best_mu    + _rnd.Rannyu(-_delta_mu,    _delta_mu);
        double sigma_trial = best_sigma + _rnd.Rannyu(-_delta_sigma, _delta_sigma);
        // vincola sigma > 0
        if (sigma_trial <= 1e-6) sigma_trial = 1e-6;

        // Calcolo energia trial (run_vmc resetta e riempie gli accumulatori a blocchi)
        double E_trial = run_vmc(mu_trial, sigma_trial);
        double err_trial = error(_global_av(_index_local_energy),
                                 _global_av2(_index_local_energy),
                                 _nblocks);

        double dE    = E_trial - best_energy;
        double alpha = exp(-beta * dE);

        // Criterio di accettazione
        if (_rnd.Rannyu() < alpha) {
            best_mu     = mu_trial;
            best_sigma  = sigma_trial;
            best_energy = E_trial;
            best_err    = err_trial;
            accepted   += 1;
        }
    }

    // Aggiorna i parametri del sistema
    _mu    = best_mu;
    _sigma = best_sigma;

    // Registra i risultati: aggiungi header se il file è vuoto
    std::ofstream f("../OUTPUT/energy_SA_new.dat", std::ios::app);
    if (f.tellp() == std::streampos(0)) {
        f << "# T\tmu\tsigma\tE\tEerr\tacc\n";
    }
    f << T << "\t" << _mu << "\t" << _sigma << "\t"
      << best_energy << "\t" << best_err << "\t"
      << double(accepted) / double(n_steps) << "\n";
    f.close();

    // Stampa a terminale (denominatore corretto: n_steps)
    std::cout << T << "\t" << _mu << "\t" << _sigma << "\t"
              << best_energy << "\t" << best_err << "\t"
              << double(accepted) / double(n_steps) << "\n";
}



/*
 QUESTA VERSIONE NON HA L'ERRORE IN OUTPUT!!!
void System::simulated_annealing(double T) {
    // Parametri di SA
    const int n_steps = get_sa_steps();
    const double beta = 1.0 / T;
    
    // Punto di partenza
    double best_mu     = _mu;
    double best_sigma  = _sigma;
    double best_energy = run_vmc(best_mu, best_sigma);

    int accepted = 0;
    

    // Ciclo di Metropolis per SA
    for (int i = 0; i < n_steps; ++i) {
        // Proposta simmetrica per mu e sigma
        double mu_trial = best_mu    + _rnd.Rannyu(-_delta_mu,    _delta_mu);
        double sigma_trial = best_sigma + _rnd.Rannyu(- _delta_sigma , _delta_sigma);

        
        // Calcolo energia trial (run_vmc si occupa di resettare block_accumulators)
        double E_trial = run_vmc(mu_trial, sigma_trial);
        double ΔE      = E_trial - best_energy;
        double alpha   = exp(-beta * ΔE);

        // Criterio di accettazione
        if (_rnd.Rannyu() < alpha) {
            best_mu     = mu_trial;
            best_sigma  = sigma_trial;
            best_energy = E_trial;
            accepted +=1;
        }

        
    }
    
    // Aggiorna i parametri del sistema
    _mu    = best_mu;
    _sigma = best_sigma;
    // Registra i risultati in ENERGY_SA.dat (sempre 4 colonne)
  std::ofstream f("../OUTPUT/energy_SA.dat", std::ios::app);
  f << T << "\t" << _mu << "\t" << _sigma << "\t" << best_energy << "\t" << double(accepted)/double(_nsteps)<< "\n";

    //stampa a terminale:
    cout<< T << "\t" << _mu << "\t" << _sigma << "\t" << best_energy << "\t" << double(accepted)/double(_nsteps)<< "\n";
}

*/

    /*
    //cout<< "TEMP " << T <<" mu "<< _mu << " sigma "<< _sigma << " accettanza sa "<< accettanza/double(get_sa_steps()) << " E_new "<< run_vmc(_mu, _sigma) <<endl; 
    //STAMPO DIRETTAMEBNTE mu e sigma SU FILE dopo sa_steps
    
    ofstream mu_file("../OUTPUT/mu.dat", ios::app);
    ofstream sigma_file("../OUTPUT/sigma.dat", ios::app);

    mu_file << setw(12) << T
            << setw(12) << _mu << endl;
    sigma_file << setw(12) << T
              << setw(12) << _sigma << endl;

    mu_file.close();
    sigma_file.close(); 
    
}
*/
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
