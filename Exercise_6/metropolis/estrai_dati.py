import os
import numpy as np

#  Lista delle temperature in formato stringa
temp_list = ['2_0', '1_8', '1_6', '1_4', '1_2', '1_0', '0_8', '0_6', '0_5']

# ✅ Percorsi base
base_path = os.getcwd()  # cartella metropolis/
sim_off_dir = os.path.join(base_path, "sim_campo_off")
sim_mag_dir = os.path.join(base_path, "sim_campo")
output_dir = os.path.join(base_path, "output_metropolis")
os.makedirs(output_dir, exist_ok=True)  # crea output_metropolis se non esiste

#  Struttura dei dati da leggere da sim_campo_off
results = {
    "energy": {
        "filename": "total_energy_T_{}.dat",
        "output": "energy.dat",
        "label": "Energy"
    },
    "susceptibility": {
        "filename": "susceptibility_T_{}.dat",
        "output": "susceptibility.dat",
        "label": "Susceptibility"
    },
    "specific_heat": {
        "filename": "specific_heat_T_{}.dat",
        "output": "specific_heat.dat",
        "label": "SpecificHeat"
    }
}

#  Estrazione energia, suscettibilità, calore specifico da sim_campo_off/
for key in results:
    Temp = []
    Values = []
    Errors = []

    print(f"\n🔍 Elaborazione: {results[key]['label']}")

    for t in temp_list:
        folder = f"T_{t}"
        filename = results[key]["filename"].format(t)
        filepath = os.path.join(sim_off_dir, folder, filename)

        if os.path.exists(filepath):
            try:
                val, err = np.loadtxt(filepath, usecols=(2, 3), unpack=True)
                T_float = float(t.replace('_', '.'))
                Temp.append(T_float)
                Values.append(val[-1])
                Errors.append(err[-1])
                print(f"✅ {folder}/{filename}")
            except Exception as e:
                print(f"⚠️ Errore leggendo {folder}/{filename}: {e}")
        else:
            print(f"❌ File non trovato: {folder}/{filename}")

    # Scrittura file in output_metropolis
    output_file = os.path.join(output_dir, results[key]["output"])
    with open(output_file, "w") as f:
        f.write(f"# Temp\t{results[key]['label']}\tError\n")
        for T, V, E in zip(Temp, Values, Errors):
            f.write(f"{T:.2f}\t{V:.6f}\t{E:.6f}\n")

    print(f"📄 Salvato: output_metropolis/{results[key]['output']}")

#  Estrazione magnetizzazione da sim_campo/
print("\n🔍 Elaborazione: Magnetization")

Temp = []
Mag = []
ErrMag = []

for t in temp_list:
    folder = f"T_{t}"
    filename = f"magnet_T_{t}.dat"
    filepath = os.path.join(sim_mag_dir, folder, filename)

    if os.path.exists(filepath):
        try:
            val, err = np.loadtxt(filepath, usecols=(2, 3), unpack=True)
            T_float = float(t.replace('_', '.'))
            Temp.append(T_float)
            Mag.append(val[-1])
            ErrMag.append(err[-1])
            print(f"✅ {folder}/{filename}")
        except Exception as e:
            print(f"⚠️ Errore leggendo {folder}/{filename}: {e}")
    else:
        print(f"❌ File non trovato: {folder}/{filename}")

# Scrittura file magnetization.dat
output_file = os.path.join(output_dir, "magnetization.dat")
with open(output_file, "w") as f:
    f.write("# Temp\tMagnetization\tError\n")
    for T, M, EM in zip(Temp, Mag, ErrMag):
        f.write(f"{T:.2f}\t{M:.6f}\t{EM:.6f}\n")

print(f"📄 Salvato: output_metropolis/magnetization.dat")
