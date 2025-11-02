MODULE Scaling
USE Forze
USE Integrazione
USE NVT 

IMPLICIT NONE
CONTAINS

! Recap. Il codice prende in input le box size e le riduce per avvicinarsi alla densità target. Ad ogni step di "Riduzione", genera le velocità secondo boltzmann perchè il volume si è riadattato.
! Fa un ciclo interno di equilibrio con termostato di berendsen ad ogni scaling delle dimensioni. Sono 15000 step per scaling, con tau 20deltat.
! Adesso può passare all'NVT vero con Bussi-Parrinello-Donadio, cercando di stabilizzarsi alla temperatura ultima in uscita dal ciclo interno.

SUBROUTINE Scaling_pos(positions, velocities, n_particles, Box_size, Temperatura, density, &
                      Sigma_Lennard_Jones, Epsilon_Lennard_Jones, n_particles_specie, &
                      numero_specie, massa_specie, r_max_gen_casuale, E_kin, En_pot, E_tot, &
                      bin_size, nbins, rcut2, viriale, gr_intra, gr_inter, Thermostat, Termalizzazione, Temperatura_calc, vel_cm)
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    INTEGER, INTENT(IN) :: n_particles, numero_specie, n_particles_specie(:), nbins
    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: velocities, positions
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)      ::  massa_specie
    REAL(KIND=dp), DIMENSION(3), INTENT(INOUT)   :: Box_size
    REAL(KIND=dp), INTENT(INOUT)                 :: Temperatura, density, Sigma_Lennard_Jones(:,:), Epsilon_Lennard_Jones(:,:)
    REAL(KIND=dp), INTENT(IN)                    :: r_max_gen_casuale
    REAL(KIND = dp), INTENT(INOUT)                  :: vel_cm  
    REAL(KIND=dp), INTENT(OUT)                   :: E_kin, En_pot, E_tot 
    REAL(KIND=dp), INTENT(INOUT)                 :: bin_size, viriale, rcut2
    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: gr_intra, gr_inter
    REAL(KIND=dp), INTENT(OUT)                   :: Temperatura_calc
    LOGICAL, INTENT(IN)                          :: Thermostat
    LOGICAL, INTENT(OUT)                         :: Termalizzazione

    ! Variabili locali
    REAL(KIND=dp) :: DeltaBox, Box_target, hn, scaling_temp,  & 
                     alpha, sumsq(3), sumsqt, sumv(3), Mtot
    REAL(KIND=dp) :: tau_T, deltat, tau_berends, DeltaTemp
    
    REAL(KIND=dp), DIMENSION(SIZE(positions,1), SIZE(positions,2)) :: forces, old_positions, new_positions
    INTEGER :: step, start_idx_s, end_idx_s, p, s, term_step
    LOGICAL :: Equilibrio, Equilibrio_termico
    REAL(KIND= dp) :: kB, kg_to_g, g_to_uma, m_to_am, s_to_fs, T_accumulato = 0.0_dp
    kB = 1.380649D-23 ![kg m**2/s**2*K]
    kg_to_g = 1D3    
    g_to_uma = 6.02214076D23 
    m_to_am = 1D10
    s_to_fs = 1D15
    kB = (kB*kg_to_g*g_to_uma*m_to_am**2d0)/((s_to_fs)**2d0)

    ! Parametri iniziali
    Box_target = (n_particles / density)**(1.0_dp/3.0_dp)
    Equilibrio = .FALSE.
    Equilibrio_termico = .FALSE.
    Termalizzazione = .TRUE.
 
    hn = 0.1_dp
    deltat = 1 ! fs
  
    PRINT*, "Temperatura target:", Temperatura
    
    DO WHILE(.NOT. Equilibrio)
      ! Scaling del box e delle posizioni
      Box_size = Box_size * (1.0_dp - hn*(Box_size(1)/Box_target - 1.0_dp))
      positions = positions * (1.0_dp - hn*(Box_size(1)/Box_target - 1.0_dp))

      ! Ad ogni scaling rigenero la distribuzione di maxwell boltzmann, ho ridotto il volume ! GOTO 136 per il codice
      CALL Termalizza_Velocita(velocities, Temperatura, massa_specie, n_particles_specie, numero_specie, n_particles)

      old_positions = positions - velocities * deltat  ! lo posso fare se le velocità sono generate ogni volta
    
      PRINT *, "Scaling Box Size ", Box_size(1), Box_size(2), Box_size(3)       

      Equilibrio_termico = .FALSE.
      term_step = 0
      
      ! Ciclo di equilibrazione interna
      tau_berends = 150.0d0

      DO step =1,20000 
        CALL  Calcola_forze(positions, forces, n_particles, n_particles_specie, numero_specie, &
              Sigma_Lennard_Jones, Epsilon_Lennard_Jones, massa_specie, en_pot, box_size, rcut2, viriale, &
              gr_intra, gr_inter, bin_size, nbins)      
        CALL Integrazione_equazioni_moto(positions, velocities, forces, n_particles, &
             n_particles_specie, numero_specie, deltat, E_kin, Box_size, &
             old_positions, new_positions, massa_specie, Temperatura_calc, &
             En_pot, E_tot, step, vel_cm)   
        velocities = velocities *SQRT(1 + (deltat/tau_berends)*((Temperatura/Temperatura_calc) -1)) ! scaling berendesen, con tau = 25*deltat
      ENDDO  

      PRINT*, "Post scaling interno, Temperatura :", Temperatura_calc

      DeltaBox = Box_size(1) - Box_target
      DeltaBox = Box_size(1) - Box_target
      IF (ABS(DeltaBox) < 0.01_dp) THEN
         Equilibrio = .TRUE.
      ELSE 
        IF (DeltaBox < 0.0_dp) THEN
            hn = -ABS(hn) * 0.2_dp
        ELSE      
            hn = ABS(hn) * 1.2_dp
        END IF
      END IF
    END DO  ! fine end do while, convergenza box raggiunta
    term_step = 0
    PRINT *, "Fase di equilibrio con termostato NVT"
    
    CALL Termalizza_Velocita(velocities, Temperatura, massa_specie, n_particles_specie, numero_specie, n_particles)

    old_positions = positions - velocities * deltat ! Stima posizioni precedenti con rel. lineare

    OPEN(Unit = 10, file='Momento_100.dat', status='replace', action='write'  )


    DO step = 1, 5000000 
      CALL  Calcola_forze(positions, forces, n_particles, n_particles_specie, numero_specie, &
                        Sigma_Lennard_Jones, Epsilon_Lennard_Jones, massa_specie, en_pot, box_size, rcut2, viriale, &
                        gr_intra, gr_inter, bin_size, nbins)
            
      CALL Integrazione_equazioni_moto(positions, velocities, forces, n_particles, &
                          n_particles_specie, numero_specie, deltat, E_kin, Box_size, &
                          old_positions, new_positions, massa_specie, Temperatura_calc, &
                          En_pot, E_tot, step, vel_cm)

      CALL BPD_Thermostat(velocities, n_particles, n_particles_specie, numero_specie, &
                                      massa_specie, deltat,  E_kin, Temperatura_calc, En_pot, E_tot, &
                                      Temperatura, Termalizzazione, step)  
      T_accumulato = T_accumulato + Temperatura_calc
      term_step = term_step + 1                                  
     IF (MOD(step, 20) == 0) THEN ! Per controllare il momento 
      WRITE(10, *) step, vel_cm, E_tot
     ENDIF                             
     IF (MOD(step, 5000) == 0) THEN
        PRINT *, "Step equilibrazione:", step, "E_tot:", E_tot, "E_pot:", En_pot, "E_kin:", E_kin, &
        "Temperatura:", Temperatura_calc
      ENDIF         
    ENDDO
    PRINT *, "Box size after scaling:", Box_size(1), Box_size(2), Box_size(3)   

    Temperatura_calc = T_accumulato / term_step
    PRINT *, "Temperatura media dopo termostato:", Temperatura_calc

  CLOSE(10)
END SUBROUTINE Scaling_pos

SUBROUTINE Termalizza_Velocita(velocities, Temperatura, massa_specie, n_particles_specie, numero_specie, n_particles)
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    INTEGER, INTENT(IN) :: n_particles, numero_specie, n_particles_specie(:)
    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: velocities
    REAL(KIND=dp), DIMENSION(:), INTENT(IN) :: massa_specie
    REAL(KIND=dp), INTENT(IN) :: Temperatura
    
    ! Variabili locali
    REAL(KIND=dp) :: kB, sigma_specie(numero_specie), u1, u2, Gaussian_random
    REAL(KIND=dp) :: sumv(3), sumsq(3), sumsqt, scaling_factor, Mtot
    INTEGER :: s, p, start_idx_s, end_idx_s, i
    
    ! Costante di Boltzmann in unità appropriate
    kB = 1.380649D-23
    kB = (kB*1D3*6.02214076D23*1D20)/(1D30)  ! Conversione a [uma*Å²/fs²*K]
    
    ! Inizializzazione
    sumv = 0.0_dp
    sumsq = 0.0_dp
    sumsqt = 0.0_dp
    Mtot = 0.0_dp
    
    ! Generazione velocità Maxwell-Boltzmann
    start_idx_s = 1
    DO s = 1, numero_specie
        sigma_specie(s) = SQRT(kB * Temperatura / massa_specie(s))
        end_idx_s = start_idx_s + n_particles_specie(s) - 1
        
        DO p = start_idx_s, end_idx_s
          
          DO i = 1,3  
           
            CALL RANDOM_NUMBER(u1)
            CALL RANDOM_NUMBER(u2)
            Gaussian_random = SQRT(-2.0d0 * LOG(u1)) * COS(2.0 * 3.141592653589793 * u2)
            velocities(i,p) = Gaussian_random * sigma_specie(s)
          
          ENDDO
            sumv(:) = sumv(:) + massa_specie(s) * velocities(:, p)
            sumsq(:) = sumsq(:) + massa_specie(s) * velocities(:, p)**2
        END DO
        
        Mtot = Mtot + massa_specie(s) * n_particles_specie(s)
        start_idx_s = end_idx_s + 1
    END DO
    
    ! Rimozione velocità del centro di massa
    IF (Mtot > 0.0_dp) THEN
        sumv = sumv / Mtot
        velocities = velocities - SPREAD(sumv, DIM=2, NCOPIES=n_particles)
    END IF
    
    ! Riscalamento alla temperatura esatta
    sumsq = 0.0_dp
    start_idx_s = 1
    DO s = 1, numero_specie
        end_idx_s = start_idx_s + n_particles_specie(s) - 1
        DO p = start_idx_s, end_idx_s
            sumsq(:) = sumsq(:) + massa_specie(s) * velocities(:, p)**2
        END DO
        start_idx_s = end_idx_s + 1
    END DO
    
    sumsqt = SUM(sumsq)
    scaling_factor = SQRT(3.0_dp * kB * n_particles * Temperatura / sumsqt)
    velocities = velocities * scaling_factor
    
END SUBROUTINE Termalizza_Velocita

END MODULE Scaling