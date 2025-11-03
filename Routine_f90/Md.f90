PROGRAM My_first_MD_Program
  USE Input           ! legge input.txt e alloca le variabili principali, def. dimensioni del box, cutoff ..
  USE Inizializza     ! inizializza posizioni random senza overlap e vel. secondo Boltz 
  USE Forze           ! Calcola forze e g(r)
  USE Integrazione    ! Verlet 
  USE Vacf            ! Calcola Velocity autocorr. function (prova openmp)
  USE SPECTRUM        ! Calcola trasformata di Fourier di dati contenuti in un file di input
  USE Scaling         ! Adatta le dimensioni del box ed equilibra, se la densità input è troppo alta 
  USE NVT             ! Termostato Bussi-Parinello-Donadio per ensemble NVT

  ! unità misura: fs, Amstr, uma, K

  IMPLICIT NONE

  INTEGER, PARAMETER             :: dp = KIND(1.0D0)
  INTEGER                        :: n_steps, n_particles, numero_specie, step, s, p, k, i, s1, ios, s2, p1,p2, q
  INTEGER                        :: end_idx_s1, end_idx_s2, start_idx_s1, start_idx_s2 ! per i loop su più specie 
  INTEGER                        :: nbins, fact_nspecie, step_sal, Max_delay, delay, gr_count, r_max ! Per g(r)
  INTEGER, ALLOCATABLE           :: n_particles_specie(:)
  REAL(KIND=dp)                  :: deltat, Temperatura, bin_size, r, en_pot, Temperatura_calc, hn, density 
  REAL(KIND=dp)                  :: E_kin, e_tot, viriale, pressione, sumsqt, r_max_gen_casuale, rcut2, tau, vel_cm
  REAL(KIND=dp)                  :: Volume_bin 
  CHARACTER(LEN=20), ALLOCATABLE :: specie(:)
  REAL(KIND=dp), ALLOCATABLE     :: density_intra(:), norm_intra(:), norm_inter(:), density_inter(:)
  REAL(KIND=dp), ALLOCATABLE     :: positions(:,:), velocities(:,:), forces(:,:)
  REAL(KIND=dp), ALLOCATABLE     :: massa_specie(:), old_positions(:,:), new_positions(:,:)
  REAL(KIND=dp), ALLOCATABLE     :: Sigma_Lennard_Jones(:,:), Epsilon_Lennard_Jones(:,:)
  REAL(KIND=dp)                  :: Box_size(3), min_val
  REAL(KIND=dp), ALLOCATABLE     :: gr_inter(:,:), gr_intra(:,:)
  REAL(KIND=dp), ALLOCATABLE     :: gr_inter_accum(:,:), gr_intra_accum(:,:) ! accumulare g(r) da forze.f90
  REAL, ALLOCATABLE              :: VACF_singolo(:,:), VACF_cross(:,:), D_cross(:)
  COMPLEX(KIND=dp), ALLOCATABLE  :: VACF_singolo_spectra(:,:), VACF_cross_spectra(:,:)
  REAL, ALLOCATABLE              :: D_singolo(:)
  LOGICAL                        :: do_scaling, Convergenza, Thermostat, Termalizzazione
  ! Locali
  REAL(KIND=dp), ALLOCATABLE     :: VACF_cross_0(:)      
  REAL(KIND = dp), ALLOCATABLE   :: H_inter(:,:), H_intra(:,:)                           
  REAL(KIND=dp), PARAMETER       :: pi = ACOS(-1.0d0)
  REAL(KIND=dp)                  :: kB = 1.380649D-23, umaAmfs2_to_kgms2, umadivAfs2_to_bar, uma_to_g, g_to_kg, &
                                    am_to_m, fs_to_s, umaA2fs2_to_kj, N_avogadro, j_to_ev, T_accumulato
  REAL(KIND=dp) :: sum_v2_dir(3),  Box_size_target(3), Delta, q_value 

  REAL    :: start_time, end_time, elapsed_time
  INTEGER :: count_rate, count_start, count_end
  
  N_avogadro = 6.02214076D23
  uma_to_g = 1/N_avogadro
  am_to_m = 1D-10
  g_to_kg = 1D-3
  fs_to_s = 1D-15

  umaAmfs2_to_kgms2 = uma_to_g*g_to_kg/(am_to_m*(fs_to_s)**2) ! 1/N_avogadro per passare da uma a g, 10-3 per g/kg, 10-30 per fs2/s2, 10-10 A/m. così sono pascal
  umadivAfs2_to_bar = umaAmfs2_to_kgms2 * (1D-5)
  j_to_ev = 1/(1.60217D-19)
  umaA2fs2_to_kj = uma_to_g*g_to_kg*am_to_m**2/(fs_to_s**2) * j_to_ev ! Energia in ev

                               ! ======= INPUT =======
  CALL leggi_input_alloca(n_steps, deltat, numero_specie, n_particles_specie, &
                            Temperatura, specie, n_particles, massa_specie, &
                            Sigma_Lennard_Jones, Epsilon_Lennard_Jones, Box_size, Temperatura_calc, &
                            e_tot, en_pot, E_kin, old_positions, new_positions, positions, &
                            forces, velocities, density,  &
                            rcut2, r_max_gen_casuale, min_val, &
                            gr_inter, gr_intra, viriale, nbins, fact_nspecie, Thermostat)
  bin_size = 0.1d0  ! dimensione del bin g(r)
  ALLOCATE(gr_inter_accum( fact_nspecie,nbins), gr_intra_accum( numero_specie,nbins))
  gr_inter_accum = 0.0_dp ! Per accumulare g(r) dalla chiamata della funz. Forze
  gr_intra_accum = 0.0_dp ! come sopra
  gr_count = 0
 

                ! ======= INIZIALIZZAZIONE POSIZIONI E VELOCITA' =======
  CALL inizializza_posizioni_velocita(n_particles, positions, velocities, numero_specie, &
                                    n_particles_specie, massa_specie, Temperatura, Box_size, &
                                    r_max_gen_casuale, sumsqt, do_scaling)
            ! ======= SCALING BOX E POS. (DENSITÀ TROPPO ALTA IN INPUT) ======== 

  IF (do_scaling .EQV. .TRUE.) THEN ! Scaling box size a densità target
    CALL Scaling_pos(positions, velocities, n_particles, Box_size, Temperatura, density, &
                      Sigma_Lennard_Jones, Epsilon_Lennard_Jones, n_particles_specie, &
                      numero_specie, massa_specie, r_max_gen_casuale, E_kin, En_pot, E_tot, &
                        bin_size, nbins, rcut2, viriale, gr_intra, gr_inter, Thermostat, Termalizzazione, Temperatura_calc, vel_cm)
    PRINT *, "Termalizzazione completata, temperatura finale pre produzione:", Temperatura_calc, "K" 
    !Finita la termalizzazione prendo la temperatura media della termalizzazione come target per la fase di produzione, vedi Scaling.f90 
    Temperatura = Temperatura_calc
    PRINT *, "Inizio produzione a temperatura:", Temperatura, "K"  
  ENDIF 

                                  ! ===== PRODUZIONE =====                              
  old_positions = positions - velocities * deltat ! Ottengo posizioni precedenti all'istante iniziale con relazione lineare

!  OPEN(unit=20, file='Posizioni_400.dat', status='replace', action='write') ! File per risultati 
!  OPEN(unit=50, file = 'Vel_cm_400.dat', status = 'replace', action= 'write')
  OPEN(unit=40, file='Velocità.dat', status='replace', action='write')
!  OPEN(unit=30, file='Energia_pressione_400.dat', status='replace', action='write') 

!  WRITE(30, '(A10,4A20)') 'Step', 'E_tot', 'E_pot', 'E_kin', 'Pressione'

  step_sal  = 10  ! Ogni 10 step salvo frame
  Max_delay = 1600 ! Massimo delay per calcolo vacf

  ALLOCATE(VACF_singolo(Max_delay, numero_specie), VACF_cross(Max_delay, fact_nspecie)) ! Per calcolare diffusione specie in "solvente" alloco
                                                                                        ! Vacf-cross con il numero di coppie
  CALL SYSTEM_CLOCK(count_rate=count_rate)
  CALL SYSTEM_CLOCK(count=count_start)
  
  DO step = 1, n_steps ! Ciclo di MD 

    CALL Calcola_forze(positions, forces, n_particles, n_particles_specie, numero_specie, &
                    Sigma_Lennard_Jones, Epsilon_Lennard_Jones, massa_specie, en_pot, & 
                    Box_size, rcut2, viriale, gr_intra, gr_inter, bin_size, nbins)
    CALL Integrazione_equazioni_moto(positions, velocities, forces, n_particles, n_particles_specie, &
                       numero_specie, deltat, E_kin, Box_size, old_positions, new_positions, &
                       massa_specie, Temperatura_calc, En_pot, E_tot, step, vel_cm)
    pressione = (2.0d0/(3.0d0*PRODUCT(Box_size))*(E_kin - 0.5d0*viriale)) 
    IF (Thermostat) THEN ! Chiamo il termostato fuori dall'integrazione, per generalità
      CALL BPD_Thermostat(velocities, n_particles, n_particles_specie, numero_specie, &
                                      massa_specie, deltat, E_kin, Temperatura_calc, En_pot, E_tot, &
                                      Temperatura, Termalizzazione, step)
    ENDIF

   ! IF (MOD(step, step_sal) == 0) THEN
     ! WRITE(20, '(3F15.8)') positions
   !  gr_inter_accum = gr_inter_accum + gr_inter
   !  gr_intra_accum = gr_intra_accum + gr_intra ! accumulo la gr calcolata nelle forze
   ! ENDIF
    
   ! IF(MOD(step, 10)  == 0) THEN ! Su file riporto dati, ogni 200 step
   !   WRITE(30, '(I10,5F15.8)') step, Temperatura_calc, e_tot, E_kin, En_pot
   !   WRITE(50, *) step, vel_cm
   ! ENDIF

    IF(MOD(step, step_sal) == 0 ) THEN ! Salvo le velocità per VACF e diffusione
    DO p = 1,n_particles 
      WRITE(40, '(3F15.8)') velocities(1,p), velocities(2,p), velocities(3,p)
!      WRITE(20,*) positions(1,p), positions(2,p), positions(3,p)
    ENDDO  
    ENDIF

    IF (MOD(step, 50000) == 0) THEN
       PRINT *, "Step Prod.:", step, "E_tot,(eV) :", E_tot*umaA2fs2_to_kj, "E_pot (eV):", &
                 en_pot*umaA2fs2_to_kj , "E_kin (eV):", E_kin*umaA2fs2_to_kj, &
               "Temperatura (K):", Temperatura_calc 
    END IF
  ENDDO
  
  CALL SYSTEM_CLOCK(count=count_end)
  elapsed_time = REAL(count_end - count_start) / REAL(count_rate)
  PRINT '("WALL TIME Produzione= ", F10.3, " secondi")', elapsed_time
  !CLOSE(20); CLOSE(30);
   CLOSE(40)

                      ! =========== FINE PRODUZIONE ==============
  
                      ! ===== Normalizzazione della g(r) ======
                  
  ALLOCATE(density_intra(numero_specie), norm_intra(numero_specie), &
  density_inter(fact_nspecie), norm_inter(fact_nspecie))

  ! Normalizzo la g(r)
 ! density_intra(:) = n_particles_specie(:)/PRODUCT(Box_size) 
 ! norm_intra(:) = density_intra(:)*n_particles_specie(:) * n_steps/step_sal
 ! DO s1 = 1, numero_specie - 1
 !   DO s2 = s1 + 1, numero_specie
  !     k = (s1-1)*(2*numero_specie-s1)/2 + (s2-s1)
  !     density_inter(k) = (n_particles_specie(s2))/PRODUCT(Box_size) ! Densità specie 2 a distanza r e r+dr da specie 1
  !     norm_inter(k) = density_inter(k)*n_particles_specie(s1)* n_steps/step_sal ! numero di origini della specie 1
  !  END DO
  !ENDDO

!  DO i = 1, nbins
!     r = bin_size*(i-0.5d0)  ! Centro del bin
!     Volume_bin = 4.0d0 * pi * r*r * bin_size
!     gr_intra_accum(:,i) = gr_intra_accum(:,i)/(Volume_bin * norm_intra(:))
!     gr_inter_accum(:,i) =  gr_inter_accum(:,i)/(Volume_bin * norm_inter(:))
!  ENDDO
  


!  OPEN(unit=70, file='risultati_gr_400.dat', status='replace')
!  DO i = 1, nbins
!    r = bin_size * (i - 0.5_dp)
!    WRITE(70, '(100F15.8)') r ,(gr_inter_accum(s1,i), s1=1,fact_nspecie), &
    !                         (gr_intra_accum(s1,i), s1=1,numero_specie)                       
!  ENDDO

!  CLOSE(70) 

                          ! =======================
                       
                 ! ====== Velocity autocorrelation function ======     

  !CALL  Calcola_Vacf(n_particles, Max_delay, numero_specie, &
  !                                n_particles_specie, VACF_singolo, VACF_cross, &
  !                                n_steps,Step_sal, deltat, D_cross, D_singolo, Temperatura, fact_nspecie)
                                 

  !OPEN(unit=50, file='VACF_singolo_400.dat', status='replace')
  !OPEN(unit=60, file='VACF_cross_400.dat', status='replace')
  !DO delay = 1, Max_delay
   ! WRITE(50, '(F12.4, 5X, 100(E25.8, 5X))') delay*deltat, (VACF_singolo(delay, s), s = 1, numero_specie)
  !ENDDO

!  DO delay = 1, Max_delay
!    WRITE(60, '(F12.4, 5X, 100(E25.8, 5X))') delay*deltat, (VACF_cross(delay, s), s = 1, fact_nspecie)
!  ENDDO
!  CLOSE(50); CLOSE(60)
!  DO s = 1, numero_specie
!    PRINT *, "D_singolo specie", s, ":", D_singolo(s)*1D-1 , "cm**2/s"
!  ENDDO
  ! ============ SPETTRI CON DISCRETE FOURIER TRANS. ========== 
!  PRINT *, "Calcolo spettro di Fourier della VACF..."

 ! Chiamata alla subroutine DFT per calcolare gli spettri

!  CALL  DFT('VACF_singolo_400.dat')
!  CALL  DFT('VACF_cross_400.dat') 
  ! ===========================================================
  PRINT *, "Simulazione completata con successo!"


END PROGRAM My_first_MD_Program
