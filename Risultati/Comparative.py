import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Parametri
# -------------------------
temps = [100, 200,400]
specie = ['CH4', 'H2']
coppie = ['CH4 - H2']
masse = [16.04, 2.016]  # masse in u
colori = ['orange', 'green', 'purple']
kB = 1.380649e-23  # [kg m²/s² K]
kg_to_g = 1e3    
g_to_uma = 6.02214076e23 
m_to_am = 1e10
s_to_fs = 1e15
kB = (kB * kg_to_g * g_to_uma * m_to_am**2) / (s_to_fs**2)  # kB in uma Å²/fs² K
c = 2.99792458e-5  # cm/fs

# -------------------------
# Dizionari per memorizzare i dati
# -------------------------
risultati_gr_intra = {}
risultati_gr_inter = {}
risultati_vacf = {}
risultati_vacf_dft = {}
risultati_vacf_cross = {}
risultati_vacf_dft_cross = {}
risultati_fft_cross = {}
risultati_fft_singolo = {}
temperatura = {}
vel_cm = {}
mean_cm = {}
etot = {}
E_kin = {}
e_pot = {}
mom = {}
step = {}
x_plot_En = {}
TempMedia = {}  
std_temp = {}
E_kin_media = {}    
E_pot_media = {}
E_tot_media = {}
mom_medio = {}
FluttazioneEn_tot = {}
FluttazioneEn_kin = {}
FluttazioneEn_pot = {}
fluttuazione_temp = {}
Fluttuazione_cm = {}
max_flutt_E = {}
flutt_media_E = {}
E_tot_std = {}      
Flutt_std_E = {}
Err_rel ={}

# -------------------------
# Lettura dati
# -------------------------
for T in temps:
    # Carico file
    dati_gr = np.loadtxt(f'risultati_gr_{T}.dat', unpack=True)
    dati_vacf = np.loadtxt(f'VACF_singolo_{T}.dat', unpack=True)
    dati_vacf_dft = np.loadtxt(f'VACF_singolo_{T}_DFT.dat', unpack=True)
    dati_vacf_cross = np.loadtxt(f'VACF_cross_{T}.dat', unpack=True)
    dati_vacf_dft_cross = np.loadtxt(f'VACF_cross_{T}_DFT.dat', unpack=True)
    
    r = dati_gr[0]
    tempi = dati_vacf[0]
    freq = dati_vacf_dft[0] / c  # Converti in cm^-1
    mask = (freq >= 0) & (freq < 5000)
    mask_gr = (r >= 0) & (r < 13)
    
    freq = freq[mask]
    dati_vacf_dft = dati_vacf_dft[:, mask]
    dati_vacf_dft_cross = dati_vacf_dft_cross[:, mask]
    r = r[mask_gr]
    dati_gr = dati_gr[:, mask_gr]

    # -------------------------
    # Dati intramolecolari
    # -------------------------
    for i, chem in enumerate(specie):
        risultati_gr_intra.setdefault(chem, {})[T] = dati_gr[i+1]
        risultati_vacf.setdefault(chem, {})[T] = dati_vacf[i+1] / dati_vacf[i+1][0]  # Normalizzazione
        risultati_vacf_dft.setdefault(chem, {})[T] = dati_vacf_dft[i+1]
        
        # FFT
        N = len(tempi)
        dt = tempi[1] - tempi[0]
        fft_vals = np.fft.fft(risultati_vacf[chem][T] - np.mean(risultati_vacf[chem][T]))
        fft_freq = np.fft.fftfreq(N, dt) / c
        mask_fft = (fft_freq >= 0) & (fft_freq < 3000)
        risultati_fft_singolo.setdefault(chem, {})[T] = {
            "freq_THz": fft_freq[mask_fft],
            "amp": np.abs(fft_vals[mask_fft]) * 2 
        }

    # -------------------------
    # Dati intermolecolari / coppie
    # -------------------------
    for j, chem in enumerate(coppie):
        idx = len(coppie)
        risultati_gr_inter.setdefault(chem, {})[T] = dati_gr[idx]
        risultati_vacf_cross.setdefault(chem, {})[T] = dati_vacf_cross[-1]/dati_vacf_cross[-1][0]
        risultati_vacf_dft_cross.setdefault(chem, {})[T] = dati_vacf_dft_cross[-1]
        
        fft_vals = np.fft.fft(risultati_vacf_cross[chem][T] - np.mean(risultati_vacf_cross[chem][T]))
        fft_freq = np.fft.fftfreq(N, dt) / c
        mask_fft = (fft_freq >= 0) & (fft_freq < 3000)
        risultati_fft_cross.setdefault(chem, {})[T] = {
            "freq_THz": fft_freq[mask_fft],
            "amp": np.abs(fft_vals[mask_fft]) * 2 
        }

    # -------------------------
    # Lettura energia e temperatura
    # -------------------------

    step, temperatura[T], etot[T], E_kin[T], e_pot[T] = np.loadtxt(
    f'Energia_pressione_{T}.dat',
    unpack=True,
    usecols=(0, 1, 2, 3, 4),  # controlla che il file abbia queste colonne
    skiprows=1      
    )
    steps, vel_cm[T] = np.loadtxt(f'Vel_cm_{T}.dat',    unpack=True,
    usecols=(0, 1))
    
    etot[T] = etot[T]/200  # Energia totale in unità sim.
    E_kin[T] = E_kin[T]/200  # Energia cinetica in unità sim.
    e_pot[T] = e_pot[T]/200  # Energia potenziale in unità sim. 
    #mom[T] = mom[T]/200  # Momento in unità sim.    
# Array di indici per plottare
    x_plot_En[T] = np.arange(len(step))

# Calcoli statistici
    TempMedia[T] = np.mean(temperatura[T])
    std_temp[T] = np.std(temperatura[T])

    E_tot_media[T] = np.mean(etot[T])
    E_pot_media[T] = np.mean(e_pot[T])
    E_kin_media[T] = np.mean(E_kin[T])
    mean_cm[T] = np.mean(vel_cm[T])

   
    FluttazioneEn_tot[T] = etot[T] - E_tot_media[T]
    FluttazioneEn_tot[T] = FluttazioneEn_tot[T]  # per una migliore visualizzazione
    FluttazioneEn_kin[T] = E_kin[T] - E_kin_media[T]
    FluttazioneEn_pot[T] = e_pot[T] - E_pot_media[T]
   # print (f"T={T}K: E_tot={E_tot_media[T]:.2f}±{np.std(etot[T]):.2f}, E_kin={E_kin_media[T]:.2f}±{np.std(E_kin[T]):.2f}, E_pot={E_pot_media[T]:.2f}±{np.std(e_pot[T]):.2f}, T={TempMedia[T]:.2f}±{std_temp[T]:.2f}, Mom={mom_medio[T]:.2f}±{np.std(mom[T]):.2f}")
    fluttuazione_temp[T] = temperatura[T] - TempMedia[T]
    Fluttuazione_cm[T] = vel_cm[T] - mean_cm[T]
    
    max_flutt_E[T] = np.max(np.abs(FluttazioneEn_tot[T]))
    flutt_media_E[T] = np.mean(FluttazioneEn_tot[T])
    E_tot_std[T] = np.std(etot[T])
    Err_rel[T] = E_tot_std[T]/E_tot_media[T] * 100
# -------------------------
# Plot g(r) specie singole
# -------------------------
for chem in specie:
    plt.figure(figsize=(7,5))
    for k, T in enumerate(temps):
        plt.plot(r, risultati_gr_intra[chem][T], label=f'T={T}K', color=colori[k])
    plt.axhline(1,linestyle='--', color = 'black')    
    plt.title(f'g(r) per {chem}')
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.legend()
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f'{chem}_gr.png')
    plt.close()

# -------------------------
# Plot VACF + FFT specie singole
# -------------------------
for chem in specie:
    fig, axs = plt.subplots(2,1, figsize=(8,10))
    for k, T in enumerate(temps):
        color = colori[k]
        axs[0].plot(tempi, risultati_vacf[chem][T], label=f'T={T}K', color=color)
        axs[0].axhline(0, linestyle = '--')
        axs[1].plot(risultati_fft_singolo[chem][T]["freq_THz"],
                    risultati_fft_singolo[chem][T]["amp"], 'o', label=f'T={T}K', color=color)
        
    axs[0].set_title(f'VACF per {chem}')
    axs[0].set_xlabel('Tempo (fs)')
    axs[0].set_ylabel('VACF (norm.)')
    axs[0].legend()
    axs[0].grid(False)
    axs[1].set_title(f'Spettro  {chem}')
    axs[1].set_xlabel('Frequenza (cm^-1)')
    axs[1].set_yticks([])
    axs[1].legend()
    axs[1].grid(False)
    plt.tight_layout()
    plt.savefig(f'{chem}_vacf_fft.png')
    plt.close()

# -------------------------
# Plot g(r) coppie
# -------------------------
for chem in coppie:
    plt.figure(figsize=(7,5))
    for k, T in enumerate(temps):
        plt.plot(r, risultati_gr_inter[chem][T], label=f'T={T}K', color=colori[k])
    plt.title(f'g(r) intermolecolare per {chem}')
    plt.xlabel('r (Å)')
    plt.ylabel('g(r)')
    plt.legend()
    plt.grid(False)
    plt.axhline(1,linestyle = '--', color = 'black')
    plt.tight_layout()
    plt.savefig(f'{chem}_gr_cross.png')
    plt.close()

# -------------------------
# Plot VACF + FFT coppie
# -------------------------
for chem in coppie:
    plt.figure(figsize=(8,10))
    for k, T in enumerate(temps):
        color = colori[k]
        plt.plot(tempi, risultati_vacf_cross[chem][T], label=f'T={T}K', color=color)
        plt.axhline(0,linestyle = '--')
    plt.title(f'VACF cross per {chem}')
    plt.xlabel('Tempo (fs)')
    plt.ylabel('VACF cross (norm.)')
    #plt.yticks([])
    plt.legend()
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f'{chem}_vacf_fft_cross.png')
    plt.close()

# -------------------------
# Plot temperatura
# -------------------------
plt.figure(figsize=(12,6))

for k, T in enumerate(temps):
    color = colori[k % len(colori)]
    plt.plot(steps, vel_cm[T], label=f'T={T}K', color=color, alpha=0.7)
    #plt.axhline(y=mean_cm[T], color=color, linestyle='--', 
    #            label=f'Fluttuazione media Momento Temp. {T}K = {mean_cm[T]:.3e}K')

plt.xlabel('Numero di step')
plt.ylabel('Fluttuazione Momento in Å/fs')
plt.title('Fluttuazioni Momento')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('Confronto_mom_tutte.png', dpi=300)
#plt.show()
plt.close()

plt.figure(figsize=(12,6))
for k, T in enumerate(temps):
    color = colori[k]
    plt.plot(x_plot_En[T], temperatura[T], label=f'T={T}K', color=color, alpha=0.7)
    plt.axhline(y=TempMedia[T], color=color, linestyle='--', label=f'T media {T}K = {TempMedia[T]:.1f}K')
plt.xlabel('Numero di step')
plt.ylabel('Temperatura (K)')
plt.title('Confronto temperature nel tempo')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('Confronto_temperature_tutte.png', dpi=300)
plt.close()

# -------------------------
# Plot energie
# -------------------------

for k, T in enumerate(temps):
    plt.figure(figsize=(12,6))
    color = colori[k]
    plt.plot(x_plot_En[T], E_kin[T], label=f'E_kin {T}K', color='green')
    plt.axhline(E_kin_media[T], color='green', linestyle='--', label=f'E_kin media {T}K = {E_kin_media[T]:.6f} eV')
    plt.plot(x_plot_En[T], e_pot[T], label=f'E_pot {T}K', color='orange')
    plt.axhline(E_pot_media[T], color='orange', linestyle='--', label=f'E_pot media {T}K = {E_pot_media[T]:.6f} eV')
    plt.axhline(E_tot_media[T], color='red', linestyle='--', label=f'E_tot_media{T} K ={E_tot_media[T]:.6f} eV')
    print(f'A temperatura {T} K la fluttuazione media % è', Err_rel[T] )
    plt.title(f'Energia  per particella - Temperatura {temps[k]} K')
    plt.xlabel('Numero di step')
    plt.ylabel('Energia (eV)')
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'Confronto_energie_tutte_temp_{T}.png', dpi=300)
    plt.close()

# -------------------------
# Fluttuazioni energia totale
# -------------------------
fig, axs = plt.subplots(3,1, figsize=(12,10))
for k, T in enumerate(temps):
    ax = axs[k]
    ax.plot(x_plot_En[T], FluttazioneEn_tot[T]/2, label=f'T={T}K', color=colori[k], linewidth=1)
   # ax.axhline(flutt_media_E[T], color=colori[k], linestyle='--', label=f'Media = {flutt_media_E[T]:.2e}')
   # ax.axhline(max_flutt_E[T]/2, color=colori[k], linestyle=':', label=f'Max Fluttuazione En T= {T} K')
    ax.set_ylabel('ΔE_tot (eV)')
    ax.set_title(f'Fluttuazioni Energia Totale per particella - T={T}K')
    ax.legend()
    ax.grid(alpha=0.3)
axs[-1].set_xlabel('Step')
plt.tight_layout()
plt.savefig('Fluttuazioni_subplots_method.png', dpi=300)
plt.close()
