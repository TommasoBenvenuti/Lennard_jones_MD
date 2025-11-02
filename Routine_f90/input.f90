MODULE Input
IMPLICIT NONE 
! fai un box grande e comprimo prima degli step 
CONTAINS
! ==================================================
! SUBROUTINE PER LETTURA INPUT E ALLOCAZIONE MEMORIA
! ==================================================
SUBROUTINE leggi_input_alloca(n_steps, deltat, numero_specie, n_particles_specie, &
                            Temperatura, specie, n_particles, massa_specie, &
                            Sigma_Lennard_Jones, Epsilon_Lennard_Jones, Box_size, Temperatura_calc, &
                            e_tot, en_pot, E_kin, old_positions, new_positions, positions, &
                            forces, velocities, density,  &
                            rcut2, r_max_gen_casuale, min_val, &
                            gr_inter, gr_intra, viriale, nbins, fact_nspecie, Thermostat)
IMPLICIT NONE
! --- Dichiarazione argomenti ---
INTEGER, PARAMETER                          :: dp = KIND(1.0D0)
INTEGER, INTENT(OUT)                        :: n_steps, numero_specie, n_particles, nbins, fact_nspecie
REAL(KIND=dp), INTENT(OUT)                  :: deltat, Temperatura, density
CHARACTER(LEN=20), ALLOCATABLE, INTENT(OUT) :: specie(:)
INTEGER, ALLOCATABLE, INTENT(OUT)           :: n_particles_specie(:)
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)     :: massa_specie(:)
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)     :: Sigma_Lennard_Jones(:,:), Epsilon_Lennard_Jones(:,:)
REAL(KIND=dp), DIMENSION(3), INTENT(OUT)    :: Box_size
REAL(KIND=dp),   INTENT(OUT)                :: en_pot, E_kin, e_tot, Temperatura_calc
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)     :: old_positions(:,:), positions(:,:), new_positions(:,:)
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)     :: velocities(:,:), forces(:,:)
REAL(KIND=dp), INTENT(OUT)                  :: viriale
REAL(KIND=dp), INTENT(OUT)                  :: rcut2, r_max_gen_casuale, min_val
REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)     :: gr_inter(:,:), gr_intra(:,:)
LOGICAL, INTENT(OUT)                        :: Thermostat
! --- Variabili locali ---
INTEGER :: i, j, k, dim_gr, ios, stat_var, max_part, n_switch, s
CHARACTER(LEN=40) :: coppia, solvente, switch, ensemble
REAL  :: sigma_weighted, max_sigma, bin_size = 0.1d0, rcut
REAL(KIND=dp), PARAMETER :: Pi = 3.14159265358979323846_dp
REAL(KIND=dp)            ::  kb, kg_to_g, g_to_uma, m_to_am, s_to_fs

kB = 1.380649D-23 ![kg m**2/s**2*K]
kg_to_g = 1D3    
g_to_uma = 6.02214076D23 
m_to_am = 1D10
s_to_fs = 1D15
kB = (kB*kg_to_g*g_to_uma*m_to_am**2d0)/((s_to_fs)**2d0)

!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    ! Apertura file di input
    OPEN(UNIT=10, FILE='input.txt', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
      PRINT *, "!!! ERRORE: Impossibile aprire il file input.txt !!!"
    END IF
     
     READ(10, *, IOSTAT=ios) 
     IF (ios /= 0) STOP "Prima riga di intestazione"
    ! Lettura parametri principali
    READ(10, *, IOSTAT=ios) n_steps, deltat, numero_specie, Temperatura, density, ensemble
    IF (ios /= 0) THEN
      PRINT *, "!!! ERRORE: Formato parametri non valido !!!"
      PRINT *, "Atteso: numero di steps, timestep (fs), numero di specie nella miscela, Temperatura (K), Ensemble"
      STOP
    END IF
    ! Controlli sui valori
    IF (numero_specie <= 0) THEN
      PRINT *, "!!! ERRORE: Numero di specie deve essere > 0 !!!"
      STOP
    END IF
    IF (n_steps <= 0) THEN
      PRINT *, "!!! ERRORE: Numero di passi deve essere > 0 !!!"
      STOP
    END IF
    IF (deltat <= 0.0) THEN
      PRINT *, "!!! ERRORE: Deltat deve essere > 0 !!!"
      STOP
    END IF
    IF (density > 1.04_dp) THEN
      PRINT *, "!!! ERRORE: Densità troppo alta, deve essere <= 1.04, unità ridotte !!!"
      STOP
    END IF
    IF (TRIM(ensemble) /= 'NVT' .AND. TRIM(ensemble) /= 'NVE') THEN
      PRINT *, "!!! ERRORE: Ensemble non valido, atteso 'NVT' o 'NVE' !!!"
      STOP
    END IF
    IF (TRIM(ensemble) == 'NVT') THEN
      Thermostat = .TRUE.
      ELSE 
      Thermostat = .FALSE.  
    END IF
    ! Allocazione array basata sul numero di specie
    ALLOCATE(specie(numero_specie), n_particles_specie(numero_specie), massa_specie(numero_specie), STAT=stat_var) 
    IF (stat_var /= 0) THEN
        PRINT *, "!!! ERRORE: Allocazione fallita per specie/n_particles_specie/massa_specie !!!"
        STOP
    ENDIF   
    ! Lettura specie e numero particelle per specie
    DO i = 1, numero_specie
      READ(10, *, IOSTAT=ios) specie(i), n_particles_specie(i)
    END DO
    CLOSE(10)

    solvente = specie(1) ! Per comodità preferisco avere la più numerosa in prima posizione
    max_part = n_particles_specie(1)
    DO s = 2, numero_specie
      IF (n_particles_specie(s)>max_part) THEN
        switch = specie(1) ! la specie da cambiare diventa la prima
        n_switch = n_particles_specie(1)

        solvente = specie(s) ! La specie abbondante è la s-ma
        max_part = n_particles_specie(s)

        specie(1) = solvente ! A questo punto scambio le particelle e le specie negli array 
        n_particles_specie(1) = max_part
        
        specie(s) = switch ! la s-ma è la meno abbondante fra le due
        n_particles_specie(s) = n_switch

      ENDIF
    ENDDO  

    n_particles =SUM(n_particles_specie)
  
    ! Allocazione matrici Lennard-Jones e array posizioni/velocità/forze
    ALLOCATE(Sigma_Lennard_Jones(numero_specie, numero_specie), &
             Epsilon_Lennard_Jones(numero_specie, numero_specie), & 
             old_positions(3, n_particles), new_positions(3, n_particles), &
             velocities(3, n_particles), & 
             forces(3, n_particles), positions( 3, n_particles), &
             STAT=stat_var)
    IF (stat_var /= 0) THEN
      PRINT *, "!!! ERRORE: Allocazione fallita per matrici LJ/posizioni !!!"
      STOP
    ENDIF       
 
    ! ==============================================
    ! 1. CALCOLO DELLE MASSE
    ! ==============================================
    DO i = 1, numero_specie
        SELECT CASE(TRIM(specie(i))) ! uma
        CASE('H2') 
          massa_specie(i) = 2.016 
        CASE('O2') 
          massa_specie(i) = 32.00    
        CASE('CH4')
          massa_specie(i) = 16.04     
        CASE('CO2')
          massa_specie(i) = 44.01     
        CASE('N2') 
          massa_specie(i) = 28.02     
        CASE DEFAULT
          PRINT *, "!!! ERRORE: Specie sconosciuta '", TRIM(specie(i)), "' !!!"
          PRINT *, "Specie disponibili: H2, O2, CH4, CO2, N2"
          STOP
      END SELECT
    END DO
    ! ==============================================
    ! 2. ASSEGNAZIONE PARAMETRI DI LENNARD-JONES
    ! ==============================================
    DO i = 1, numero_specie
      DO j = i, numero_specie ! Triangolo superiore
        IF (i == j) THEN
          ! Parametri per singola specie
          SELECT CASE(TRIM(specie(i)))
            CASE('H2') 
              Sigma_Lennard_Jones(i,i) = 2.827d0 
              Epsilon_Lennard_Jones(i,i) = 59.7d0
            CASE('O2') 
              Sigma_Lennard_Jones(i,i) = 3.467d0 
              Epsilon_Lennard_Jones(i,i) = 106.7d0 
            CASE('CH4') 
              Sigma_Lennard_Jones(i,i) = 3.678d0 
              Epsilon_Lennard_Jones(i,i) = 168.78d0 
            CASE('CO2') 
              Sigma_Lennard_Jones(i,i) = 3.832d0 
              Epsilon_Lennard_Jones(i,i) = 230.50d0 
            CASE('N2')
              Sigma_Lennard_Jones(i,i) = 3.663d0  
              Epsilon_Lennard_Jones(i,i) = 96.92d0 
          END SELECT
        ELSE
        ! Parametri per coppie miste
        coppia = TRIM(specie(i)) // "-" // TRIM(specie(j))
        SELECT CASE(coppia)
          CASE('CH4-H2', 'H2-CH4')
            Sigma_Lennard_Jones(i,j) = 3.05d0 
            Epsilon_Lennard_Jones(i,j) = 92.4d0
          CASE('CH4-O2', 'O2-CH4')
            Sigma_Lennard_Jones(i,j) = 3.35d0
            Epsilon_Lennard_Jones(i,j) = 126.0d0
          CASE('CH4-N2', 'N2-CH4')
            Sigma_Lennard_Jones(i,j) = 3.59d0
            Epsilon_Lennard_Jones(i,j) = 122.0d0
          CASE('CO2-H2', 'H2-CO2')
            Sigma_Lennard_Jones(i,j) = 3.30d0 
            Epsilon_Lennard_Jones(i,j) = 105.0d0
          CASE('CO2-O2', 'O2-CO2')
            Sigma_Lennard_Jones(i,j) = 3.65d0
            Epsilon_Lennard_Jones(i,j) = 140.2d0
          CASE('N2-O2', 'O2-N2')
            Sigma_Lennard_Jones(i,j) = 3.60d0
            Epsilon_Lennard_Jones(i,j) = 88.5d0
          CASE('N2-H2', 'H2-N2')
            Sigma_Lennard_Jones(i,j) = 3.18d0
            Epsilon_Lennard_Jones(i,j) = 60.1d0
          CASE('O2-H2', 'H2-O2')
            Sigma_Lennard_Jones(i,j) = 3.28d0
            Epsilon_Lennard_Jones(i,j) = 75.2d0
          CASE DEFAULT
            PRINT*, "Coppia", TRIM(coppia), "non disponibile "
          STOP
        END SELECT   
          ! Riempimento parte simmetrica
          Sigma_Lennard_Jones(j,i) = Sigma_Lennard_Jones(i,j)
          Epsilon_Lennard_Jones(j,i) = Epsilon_Lennard_Jones(i,j)
        END IF
      END DO
    END DO

    Epsilon_Lennard_Jones = Epsilon_Lennard_Jones * kB ! Conversione in unità di energia interna (K -> uma Å²/fs²)    

    ! Calcolo la minima distanza fra le particelle come media pesata dei valori di sigma
    sigma_weighted = 0.0_dp
    DO j = 1, numero_specie
      Sigma_weighted= sigma_weighted + n_particles_specie(j)*Sigma_Lennard_Jones(j,j)
    ENDDO 
    Sigma_weighted = Sigma_weighted/n_particles
    density = density/(Sigma_weighted**3) ! per comodità in input la chiedo in unità ridotte, convertiamola
    Box_size(:) = (n_particles/(density))**(1.0/3.0) 

    r_max_gen_casuale = (1.12_dp*MAXVAL(Sigma_Lennard_Jones)) 
    rcut =    (4.0d0*MAXVAL(Sigma_Lennard_Jones)) ! Cutoff interazioni LJ
    rcut2 = rcut * rcut

    nbins = INT(rcut/(bin_size)) + 1 
    ! Alloco g(r) in base al numero di bin e numero di specie per la intra specie e combinazioni di numero specie in classi di due per inter specie
    fact_nspecie = 1
    DO i = 1, numero_specie
      fact_nspecie = fact_nspecie*i
    ENDDO  
    DO i = 1,numero_specie -2
      fact_nspecie = fact_nspecie/i
    ENDDO   ! lo chiamo fact_nspecie ma sono le combinazioni del numero di specie in classi 2 
    fact_nspecie = fact_nspecie/2
    ALLOCATE(gr_inter(fact_nspecie, nbins), gr_intra(numero_specie, nbins), & 
             STAT=stat_var)
    IF (stat_var /= 0) THEN
      PRINT *, "!!! ERRORE: Allocazione fallita per funzione g(r)!!!"
      STOP
    ENDIF       


  END SUBROUTINE leggi_input_alloca
END MODULE Input

