import os
import shutil
import subprocess

#  Temperature in ordine decrescente
temperature_list = [2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.5]

#  Percorsi base (lo script si trova dentro sim_campo_off)
sim_path = os.getcwd()
input_dir = os.path.join(sim_path, "INPUT")
input_config_dir = os.path.join(input_dir, "CONFIG")
output_dir = os.path.join(sim_path, "OUTPUT")
output_config_dir = os.path.join(output_dir, "CONFIG")
source_dir = os.path.join(sim_path, "SOURCE")
simulator_path = os.path.join(source_dir, "simulator.exe")  # nome corretto

#  Crea le cartelle T_2, T_1_8, ..., T_0_2
for temp in temperature_list:
    folder_T = f"T_{str(temp).replace('.', '_')}"
    os.makedirs(os.path.join(sim_path, folder_T), exist_ok=True)

def edit_input_file(file_path, temperature, restart):
    with open(file_path, "r") as file:
        lines = file.readlines()

    new_lines = []
    for line in lines:
        if "TEMP" in line:
            new_lines.append(f"TEMP          {temperature}\n")
        elif "RESTART" in line:
            new_lines.append(f"RESTART       {restart}\n")
        else:
            new_lines.append(line)

    with open(file_path, "w") as file:
        file.writelines(new_lines)


#  Loop sulle temperature
for i, temp in enumerate(temperature_list):
    folder_T = f"T_{str(temp).replace('.', '_')}"
    folder_path = os.path.join(sim_path, folder_T)

    print(f"\n=== Simulazione a T = {temp} ===")

    seed_out = os.path.join(output_dir, "seed.out")

    if i > 0:
        # 1. Copia config.spin da OUTPUT/CONFIG → INPUT/CONFIG/
        os.makedirs(input_config_dir, exist_ok=True)
        config_spin_src = os.path.join(output_config_dir, "config.spin")
        config_spin_dst = os.path.join(input_config_dir, "config.spin")
        if os.path.exists(config_spin_src):
            shutil.copy(config_spin_src, config_spin_dst)
        else:
            print(f"⚠️ config.spin non trovato per T = {temperature_list[i-1]}")

        # 2. Rimuove seed.in se esiste
        seed_in = os.path.join(input_dir, "seed.in")
        if os.path.exists(seed_in):
            os.remove(seed_in)

        # 3. Copia seed.out → seed.in
        if os.path.exists(seed_out):
            shutil.copy(seed_out, seed_in)
        else:
            print(f"⚠️ seed.out non trovato per T = {temperature_list[i-1]}")

    # 4. Modifica input.dat con la nuova temperatura e il flag RESTART
    restart_flag = 0 if i == 0 else 1
    input_file_path = os.path.join(input_dir, "input.dat")
    edit_input_file(input_file_path, temp, restart_flag)

    # Debug: stampa il contenuto di input.dat
    print("\n📄 Contenuto input.dat prima del run:")
    with open(input_file_path, "r") as f:
        print(f.read())

    # 4b. Subito dopo la modifica, salviamo l'input.dat nella cartella T_x_y
    input_archive = os.path.join(folder_path, f"input_{folder_T}.dat")
    shutil.copy(input_file_path, input_archive)

    # 5. Lancia il simulatore
    os.chmod(simulator_path, 0o755)
    result = subprocess.run(["./simulator.exe"], cwd=source_dir)

    if result.returncode == 0:
        print(f"✅ Simulazione a T = {temp} completata.")
    else:
        print(f"❌ Errore nella simulazione a T = {temp} (codice {result.returncode})")
        continue

    # 6. Archiviazione dei file di output

    # total_energy.dat
    energy_src = os.path.join(output_dir, "total_energy.dat")
    energy_dest = os.path.join(folder_path, f"total_energy_{folder_T}.dat")
    if os.path.exists(energy_src):
        shutil.move(energy_src, energy_dest)
        print(f"📁 Spostato {energy_src} → {energy_dest}")
    else:
        print(f"⚠️ total_energy.dat non trovato per T = {temp}")

     # susceptibility.dat
    susce_src = os.path.join(output_dir, "susceptibility.dat")
    susce_dest = os.path.join(folder_path, f"susceptibility_{folder_T}.dat")
    if os.path.exists(susce_src):
        shutil.move(susce_src, susce_dest)
        print(f"📁 Spostato {susce_src} → {susce_dest}")
    else:
        print(f"⚠️ susceptibility.dat non trovato per T = {temp}")

    #specific_heat.dat
    sheat_src = os.path.join(output_dir, "specific_heat.dat")
    sheat_dest = os.path.join(folder_path, f"specific_heat_{folder_T}.dat")
    if os.path.exists(sheat_src):
        shutil.move(sheat_src, sheat_dest)
        print(f"📁 Spostato {sheat_src} → {sheat_dest}")
    else:
        print(f"⚠️ specific_heat.dat non trovato per T = {temp}")


    # output.dat
    output_dat_src = os.path.join(output_dir, "output.dat")
    output_dat_dest = os.path.join(folder_path, f"output_{folder_T}.dat")
    if os.path.exists(output_dat_src):
        shutil.copy(output_dat_src, output_dat_dest)
        print(f"📄 Copiato {output_dat_src} → {output_dat_dest}")
    else:
        print(f"⚠️ output.dat non trovato per T = {temp}")

    # seed.out → archivio
    seed_out_dest = os.path.join(folder_path, f"seed_{folder_T}.out")
    if os.path.exists(seed_out):
        shutil.copy(seed_out, seed_out_dest)

    # config.spin → archivio
    config_spin_src = os.path.join(output_config_dir, "config.spin")
    config_spin_dest = os.path.join(folder_path, f"config_{folder_T}.spin")
    if os.path.exists(config_spin_src):
        shutil.copy(config_spin_src, config_spin_dest)

# 7. Scrive la riga fissa in seed.in
final_seed_path = os.path.join(input_dir, "seed.in")
with open(final_seed_path, "w") as f:
    f.write("0000 0000 0000 0001\n")

print(f"\n📝 seed.in sovrascritto con '0000 0000 0000 0001'")