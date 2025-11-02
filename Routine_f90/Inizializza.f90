MODULE Inizializza

IMPLICIT NONE
CONTAINS
  ! =================================================================a
  ! SUBROUTINE PER INIZIALIZZAZIONE POSIZIONI E VELOCITÀ
  ! =================================================================
  SUBROUTINE inizializza_posizioni_velocita(n_particles, positions, velocities, numero_specie, &
                                          n_particles_specie, massa_specie, &
                                           Temperatura, Box_size, r_max_gen_casuale, sumsqt, &
                                           do_scaling) 
    IMPLICIT NONE
    ! --- Dichiarazione argomenti ---
    INTEGER, PARAMETER                         :: dp = KIND(1.0D0)
    INTEGER, INTENT(IN)                        :: n_particles, numero_specie
    INTEGER, INTENT(IN)                        :: n_particles_specie(:)
    REAL(KIND=dp), DIMENSION(:,:), INTENT(OUT) :: velocities 
    REAL(KIND=dp), DIMENSION(:,:), INTENT(OUT) :: positions
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)    :: massa_specie
    REAL(KIND=dp), INTENT(IN)                  :: Temperatura
    REAL(KIND=dp), DIMENSION(3), INTENT(INOUT) :: Box_size
    REAL(KIND=dp), INTENT(IN)                  :: r_max_gen_casuale ! Distanza minima tra particelle 
    ! --- Variabili locali ---
    REAL(KIND=dp) :: x, y, z, dist, vel_casuale, dist_x, dist_y, dist_z
    REAL(KIND=dp) :: sumv(3) = 0.0, sumsq(3) = 0.0     ! Somma velocità e quadrati
    REAL(KIND=dp) :: scaling_factor
    REAL(KIND=dp), INTENT(INOUT) :: sumsqt ! Energia cinetica totale
    REAL(KIND=dp)                :: min_val_box, Mtot
    REAL(KIND = dp), ALLOCATABLE :: Sigma_specie(:)
    INTEGER                      :: p, s, p1, p2, k, i, j ! p1, p2 contatori processi dove istruzione esegue controllo su due part. p e s, part. e specie processi singoli
    INTEGER                      :: numero_vel, counter_pos, counter_box
    INTEGER                      :: end_idx_s, start_idx_s
    LOGICAL                      :: posizione_valida
    LOGICAL, INTENT(INOUT)       :: do_scaling ! Indica se il box è stato scalato
    ! ==== Costanti ====
    REAL(KIND= dp) :: kB, kg_to_g, g_to_uma, m_to_am, s_to_fs, scale_factor, u1, u2, Gaussian_random, T_in
    kB = 1.380649D-23 ![kg m**2/s**2*K]
    kg_to_g = 1.0D3    
    g_to_uma = 6.02214076D23 
    m_to_am = 1.0D10
    s_to_fs = 1.0D15
    kB = (kB*kg_to_g*g_to_uma*m_to_am**2.0d0)/((s_to_fs)**2.0d0)
    sumsqt = 0.0d0
    
    ALLOCATE(Sigma_specie(numero_specie))
    ! ==============================================
    ! 1. INIZIALIZZAZIONE POSIZIONI CASUALI
    ! ==============================================
    PRINT *, "Inizializzazione posizioni in corso..."
    counter_pos= 0
    scale_factor = 1.02d0
    do_scaling = .FALSE.

    CALL RANDOM_SEED()
    DO p1 = 1, n_particles     
      posizione_valida = .FALSE. ! per entrare nel loop       
      DO WHILE (.NOT. posizione_valida)

        IF (counter_pos > 50000) THEN ! Nel caso in cui non trovi 50000 posizioni valide, prova ad aumentare il box
          PRINT *, "Attenzione: impossibile trovare posizioni valide in box dim." , Box_size
          Box_size = Box_size * scale_factor ! Aumenta le dimensioni del box
          do_scaling = .TRUE.
          counter_pos = 0
        END IF

        ! Genera coordinate casuali
        CALL RANDOM_NUMBER(x)
        CALL RANDOM_NUMBER(y)
        CALL RANDOM_NUMBER(z)
        ! Scala le posizioni in base al box
        x = x * Box_size(1)  
        y = y * Box_size(2)
        z = z * Box_size(3)      

        IF (p1 == 1) THEN
          ! Prima particella: posizione sempre valida
          posizione_valida = .TRUE.
        ELSE
          ! Particelle successive: controllo overlap
          posizione_valida = .TRUE. ! assumo sia verificata
          DO p2 = 1, p1-1
            dist_x = (positions(1, p2) - x) - Box_size(1) * NINT((positions(1, p2) - x) / Box_size(1))
            dist_y = (positions(2, p2) - y) - Box_size(2) * NINT((positions(2, p2) - y) / Box_size(2))
            dist_z = (positions(3, p2) - z) - Box_size(3) * NINT((positions(3, p2) - z) / Box_size(3))
            dist = SQRT(dist_x**2 + dist_y**2 + dist_z**2)
            IF (dist < r_max_gen_casuale) THEN
              posizione_valida = .FALSE.
              counter_pos= counter_pos + 1
              EXIT
            END IF
          END DO
        END IF
        ! Assegnazione coordinate se valide
        IF (posizione_valida) THEN
          positions(1, p1) = x
          positions(2, p1) = y
          positions(3, p1) = z
        END IF
      END DO
    END DO
 
  OPEN(unit=20, file='Posizioni_in.dat', status='replace', action='write')

  WRITE(20, '(3F15.8)')  positions(1, :), positions(2, :), positions(3, :) 
  CLOSE(20)
    
  
  PRINT*, "Dimensioni del box", (Box_size(i), i=1,3)

   ! 2. INIZIALIZZAZIONE VELOCITÀ (MAXWELL-BOLTZMANN)
  Mtot = 0.0d0
  start_idx_s = 1
  sumv = 0.0d0
  sumsq = 0.0d0
  DO s = 1, numero_specie
    sigma_specie(s) = SQRT(kB*Temperatura/massa_specie(s))
    end_idx_s = start_idx_s + n_particles_specie(s) - 1
    DO p = start_idx_s, end_idx_s

      DO i = 1,3  
         CALL RANDOM_NUMBER(u1)
         CALL RANDOM_NUMBER(u2)
         Gaussian_random = SQRT(-2.0d0 * LOG(u1)) * COS(2.0 * 3.141592653589793 * u2)
         velocities(i,p) = Gaussian_random * sigma_specie(s)
      ENDDO

      sumv(:) = sumv(:) + massa_specie(s)*velocities(:,p)
      sumsq(:) = sumsq(:) + massa_specie(s)*velocities(:,p)**2
    END DO    
    Mtot= Mtot + massa_specie(s)*n_particles_specie(s)
    start_idx_s = end_idx_s + 1
  END DO

  ! Calcolo velocità media 
  sumv = sumv / Mtot    
  ! Rimozione vel. centro di massa
  velocities = velocities - SPREAD(sumv, DIM = 2, ncopies = n_particles) ! crea un array 3*n_particles con vcm poi lo sottrae
  ! Calcolo energia cinetica (per riscalare alla temperatura target)
  sumsq = 0.0d0
  start_idx_s = 1
  DO s = 1, numero_specie
    end_idx_s = start_idx_s + n_particles_specie(s) - 1
    DO p = start_idx_s, end_idx_s
      sumsq(:) = sumsq(:) + massa_specie(s)*velocities(:,p)**2
    END DO
    start_idx_s = end_idx_s + 1
  END DO
  ! Calcolo fattore di scala per temperatura desiderata
  sumsqt = sum(sumsq)

  scaling_factor = SQRT(kB*3*n_particles*Temperatura/sumsqt)
  ! Applico il fattore di scala
  velocities = velocities * scaling_factor

  END SUBROUTINE inizializza_posizioni_velocita

END MODULE  
