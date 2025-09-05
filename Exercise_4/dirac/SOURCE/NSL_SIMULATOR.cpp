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

// Funzione principale: esegue la simulazione suddivisa in blocchi e gestisce input/output
int main(int argc, char *argv[]) {
    // Inizializza il contatore per le configurazioni salvate in formato XYZ (usato per output periodico)
    int nconf = 1;
    
    // Crea un oggetto System che rappresenta il sistema fisico da simulare e ne gestisce lo stato
    System SYS;
    
    // Inizializza il sistema leggendo i file di input e impostando le condizioni iniziali.
    //  - Legge i semi per il generatore di numeri casuali dai file ../INPUT/Primes e ../INPUT/seed.in.
    //  - Apre il file di output "../OUTPUT/acceptance.dat" e vi scrive l'intestazione delle colonne (es. "#   N_BLOCK:  ACCEPTANCE:").
    //  - Legge i parametri di simulazione da "../INPUT/input.dat" (tipo di simulazione SIMULATION_TYPE, eventuali parametri aggiuntivi _J e _H per Ising, flag RESTART, temperatura TEMP, numero particelle NPART, densità RHO, raggio di cutoff R_CUT, passo di Monte Carlo DELTA, numero di blocchi NBLOCKS, numero di step per blocco NSTEPS, etc.).
    //  - Calcola e stampa su "../OUTPUT/output.dat" un riepilogo dei parametri letti (ad es. tipo di simulazione scelto, temperatura, numero di particelle, lunghezza del lato del volume simulato, R_CUT, DELTA, NBLOCKS, NSTEPS) e conferma la fine della lettura input.
    //  - Imposta le posizioni iniziali delle particelle leggendo il file di configurazione iniziale "../INPUT/CONFIG/config.xyz". In caso di simulazione di Ising e RESTART attivo, legge anche "../INPUT/CONFIG/config.spin" per impostare gli spin iniziali.
    //  - Se la simulazione è di dinamica molecolare (MD) e non è un restart, inizializza le velocità delle particelle con valori casuali gaussiani coerenti con la temperatura, regolandole per avere momento totale nullo e la giusta energia cinetica media. Se invece è un restart MD, legge dal file "../INPUT/CONFIG/conf-1.xyz" le posizioni "old" precedenti per avviare il Verlet in continuità.
    //  - Infine, scrive su output.dat la conferma "System initialized!" indicando che il sistema è stato inizializzato correttamente.
    SYS.initialize();
    
    // Inizializza le proprietà da misurare durante la simulazione, leggendo l'elenco da "../INPUT/properties.dat".
    //  - Per ogni proprietà richiesta (ad es. POTENTIAL_ENERGY, KINETIC_ENERGY, TOTAL_ENERGY, TEMPERATURE, PRESSURE, ecc.), apre il corrispondente file di output in "../OUTPUT/" (es. "potential_energy.dat", "kinetic_energy.dat", ...) e scrive la riga di intestazione con i titoli delle colonne (es. "#     BLOCK:  ACTUAL_PE:     PE_AVE:      ERROR:" per l'energia potenziale).
    //  - Attiva i flag interni che indicano quali grandezze misurare (_measure_penergy, _measure_kenergy, etc.) e assegna un indice a ciascuna proprietà per l'archiviazione nei vettori.
    //  - Se la proprietà richiede parametri aggiuntivi, li legge (ad esempio, per GOFR legge il numero di bin successivo e calcola _bin_size per l'istogramma g(r); per POFV legge il numero di bin di velocità).
    //  - Aggiunge su "../OUTPUT/output.dat" (in modalità append) un messaggio "Reading properties completed!" a conferma che la configurazione delle proprietà da misurare è completata.
    //  - Alloca e azzera i vettori interni per memorizzare i valori istantanei (_measurement), le medie di blocco (_average), i cumulativi di blocco (_block_av) e le somme globali (_global_av e _global_av2) delle proprietà misurate.
    //  - Resetta i contatori _nattempts e _naccepted a 0 per iniziare il conteggio delle mosse Monte Carlo.
    SYS.initialize_properties();

    //Stampa la distribuzione di velocità iniziale, creata con initialize()
    SYS.write_starting_velocities_distribution(); 

    // Esegue un reset degli accumulatori del blocco iniziale (blocco 0) prima di iniziare la simulazione a blocchi.
    //  - La funzione System::block_reset(int blk) azzera il vettore _block_av che accumula la somma delle misure nel blocco corrente.
    //  - Poiché blk=0 (primo avvio), non viene scritto alcun messaggio su file (la scrittura "Block completed" avviene solo per blk>0). Questo chiamata serve dunque per inizializzare correttamente gli accumulatori per il primo blocco.
    SYS.block_reset(0);
    
    // Ciclo esterno sui blocchi di simulazione:
    //   SYS.get_nbl() restituisce il numero totale di blocchi (NBLOCKS) specificato nell'input. 
    //   Ad ogni iterazione `i` corrisponde un intero blocco di simulazione in cui verranno raccolti dati per le statistiche.
    for(int i = 0; i < SYS.get_nbl(); i++) { // loop sui blocchi
    
        // Ciclo interno sui passi di simulazione all'interno del blocco i-esimo:
        //   SYS.get_nsteps() restituisce il numero di step (iterazioni) da eseguire per ciascun blocco.
        //   In ogni iterazione `j` si esegue un singolo passo di simulazione seguito da una misurazione.
        for(int j = 0; j < SYS.get_nsteps(); j++) { // loop sugli step nel blocco
            // Esegue un passo di simulazione.
            //   * Se la simulazione è di **dinamica molecolare NVE** (SIMULATION_TYPE = 0), la funzione System::step() chiamerà SYS.Verlet(), eseguendo un passo di integrazione di Verlet su tutte le particelle (calcola le forze attuali su ogni particella, poi aggiorna le posizioni a t+Δt e le velocità usando posizioni attuali e precedenti). Dopo l'integrazione, _naccepted viene incrementato di _npart (in MD consideriamo tutti i "movimenti" accettati perché deterministici).
            //   * Se la simulazione è di **Monte Carlo** (SIMULATION_TYPE = 1 per Lennard-Jones NVT, 2 o 3 per Ising 1D), System::step() esegue _npart tentativi di mossa MC:
            //       - Nel caso Lennard-Jones (sim_type 1), ad ogni chiamata System::move() seleziona una particella a caso e prova a spostarla di un piccolo delta (passo MAX definito da _delta) in una direzione casuale. La mossa viene accettata o rifiutata secondo il criterio di Metropolis (confrontando energia prima e dopo, calcolata dalla funzione SYS.Boltzmann). Ogni tentativo incrementa _nattempts e, se accettato, _naccepted.
            //       - Nel caso Ising (sim_type 2 o 3), System::move() seleziona uno spin casuale e prova a flipparlo. Se sim_type=2 (Metropolis), la variazione di energia ΔE è calcolata e si accetta il flip con probabilità exp(-βΔE); se sim_type=3 (Gibbs sampling), il flip può essere deciso direttamente in base alla probabilità condizionale. Anche qui _nattempts aumenta di 1 per tentativo e _naccepted di 1 per ogni flip accettato.
            //   * Alla fine di System::step(), viene aggiornato _nattempts += _npart (vengono conteggiati _npart tentativi per step, sia in MD che MC, per uniformità di definizione dello step). 
            SYS.step();
            
            // Esegue la misurazione delle proprietà fisiche dopo il passo.
            //   La funzione System::measure() calcola i valori istantanei di tutte le grandezze richieste:
            //   - **Energia potenziale** (se _measure_penergy attivo): somma i contributi del potenziale tra tutte le coppie di particelle (ad esempio potenziale di Lennard-Jones 4ε[(σ/r)^12 - (σ/r)^6] se simulazione LJ).
            //   - **Energia cinetica** (se _measure_kenergy attivo, in MD): calcola ½m v^2 per ogni particella (con m=1 convenzionale) e fa la media su tutte le particelle.
            //   - **Energia totale** (se _measure_tenergy attivo): somma di energia potenziale e cinetica per sistemi continui, oppure energia interna dell'Ising (calcolata da -J Σ_s s_i s_j - μH Σ_i s_i).
            //   - **Temperatura** (se _measure_temp attivo, in MD): calcolata da 2/3 * (energia cinetica media) per sistema di particelle.
            //   - **Pressione** (se _measure_pressure attivo, per LJ): calcola il viriale tramite la forza fra particelle (termine 48[(σ/r)^12 - 0.5(σ/r)^6]) e usa l'equazione di stato microscopica P = ρ k_B T + (viriale/volume).
            //   - **Magnetizzazione** (se _measure_magnet attivo, per Ising): calcola la magnetizzazione media Σ_i s_i / N.
            //   - **Altre proprietà** come distribuzione g(r) (se richiesta con GOFR) o distribuzione delle velocità P(v) (POFV) o specific heat (CV), suscettività (CHI) vengono anch'esse preparate, anche se nel codice base alcune sono segnaposto per esercizi successivi.
            //   La funzione azzera un vettore _measurement e vi inserisce i valori calcolati delle proprietà attive. 
            //   Infine, aggiunge (_measurement) ai cumulativi di blocco: _block_av += _measurement, per sommare il contributo di questo step alle somme del blocco corrente (servirà per la media di blocco).
            SYS.measure();
            
            // Ogni 50 passi circa, salva la configurazione attuale su file (opzione per analisi visuale).
            if(j % 50 == 0) {
                // SYS.write_XYZ(nconf); // Scrive la configurazione corrente in formato XYZ nel file "../OUTPUT/CONFIG/config_<nconf>.xyz"
                // (Questa riga è stata commentata per evitare la creazione di troppi file e il riempimento del filesystem)
                // Se fosse attiva, per ogni 50° step verrebbe creato un file con le coordinate di tutte le particelle,
                // in formato XYZ (prima riga numero particelle, seconda riga commento, poi coordinate x,y,z normalizzate al lato della scatola).
                // Questo output serve, ad esempio, per visualizzare l'evoluzione del sistema tramite strumenti grafici esterni.
                nconf++;  // Incrementa il contatore di configurazioni anche se la scrittura è disabilitata, mantenendo la numerazione.
            }
        } // fine del loop sui passi all'interno del blocco i
        
        // Al termine di ciascun blocco di simulazione (i-esimo blocco completato), calcola le medie e scrive i risultati del blocco.
        // La funzione System::averages(int blk) elabora i dati accumulati nel blocco appena concluso (blk = i+1, numerazione dei blocchi da 1):
        //   - Calcola la **media di blocco** di ogni proprietà: _average = _block_av / Nsteps, dove Nsteps è il numero di step per blocco.
        //   - Aggiorna i valori globali: accumula _global_av += _average e _global_av2 += (_average ^ 2) elemento per elemento, per utilizzare questi valori nel calcolo delle incertezze statistiche.
        //   - Per **ogni proprietà misurata**, apre il rispettivo file di output (in append) e aggiunge una riga con i risultati del blocco:
        //         numero del blocco, valore medio nel blocco, media cumulativa fino a quel blocco, e errore statistico (calcolato con la formula dell'errore standard usando _global_av e _global_av2 dei blocchi precedenti).
        //     Ad esempio, per l'energia potenziale (se attiva) scrive su "../OUTPUT/potential_energy.dat" il numero di blocco, la media di E_pot del blocco, la media progressiva E_pot fino a blocco blk, ed errore associato.
        //     Analogamente scrive su "kinetic_energy.dat", "total_energy.dat", "temperature.dat", "pressure.dat", etc. per le altre grandezze attive.
        //   - Calcola inoltre la **frazione di accettazione** delle mosse Monte Carlo fino a questo blocco: fraction = _naccepted/_nattempts (se _nattempts > 0), e aggiunge una riga nel file "../OUTPUT/acceptance.dat" con il numero di blocco e tale frazione. Questo permette di monitorare l'efficienza dell'algoritmo MC (ad esempio, per regolare Δ).
        //   - La funzione poi termina (non azzera _nattempts/_naccepted, che nel codice rimangono cumulativi sull'intera simulazione).
        SYS.averages(i+1);
        
        // Reset degli accumulatori per iniziare il blocco successivo.
        //   La funzione System::block_reset(int blk) azzera il vettore _block_av (somma dei valori per il nuovo blocco).
        //   Inoltre, se blk > 0, scrive un messaggio di log su "../OUTPUT/output.dat": nello specifico, qui blk = i+1 rappresenta il blocco appena concluso, quindi verrà aggiunta la riga "Block completed: i+1".
        //   In questo modo il file output.dat conterrà una traccia testuale del completamento di ogni blocco. Dopo il reset, _block_av è pronto a riaccumulare dati nel prossimo blocco.
        SYS.block_reset(i+1);
    } // fine del loop sui blocchi
    
    // Finalizzazione della simulazione una volta completati tutti i blocchi.
    // La funzione System::finalize() chiude correttamente la simulazione svolgendo diverse operazioni di output finale:
    //   - Chiama SYS.write_configuration() per salvare la configurazione finale del sistema:
    //         * Per simulazioni di Lennard-Jones (SIMULATION_TYPE 0 o 1): scrive il file "../OUTPUT/CONFIG/config.xyz" con le coordinate finali delle particelle (in formato XYZ, coordinate normalizzate alla scatola unitaria), e anche "../OUTPUT/CONFIG/conf-1.xyz" con le posizioni "old" (dell'ultimo step precedente) che possono servire per riavviare una simulazione MD in seguito.
    //         * Per simulazioni di Ising (SIMULATION_TYPE 2 o 3): scrive il file "../OUTPUT/CONFIG/config.spin" contenente la configurazione finale degli spin.
    //   - Chiama _rnd.SaveSeed() per salvare lo stato finale del generatore di numeri casuali su file (tipicamente il file "seed.out"), così da poter riutilizzare il seed finale per continuare la simulazione o per riproducibilità.
    //   - Apre il file "../OUTPUT/output.dat" in append e scrive "Simulation completed!" per segnalare nel log che la simulazione è terminata con successo.
    //   - (Eventuali file di output aperti vengono chiusi all'interno delle funzioni sopra, garantendo che tutti i dati siano scritti sul disco).
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
