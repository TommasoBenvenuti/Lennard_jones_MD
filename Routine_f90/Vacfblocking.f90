MODULE Vacf

  IMPLICIT NONE
  CONTAINS

  SUBROUTINE Calcola_Vacf(n_particles, Max_delay, numero_specie, &
                                  n_particles_specie, VACF_singolo, VACF_cross, &
                                  n_steps,Step_sal, deltat, D_cross, D_singolo, Temperatura, fact_nspecie)

    USE OMP_LIB                           
    INTEGER, PARAMETER :: dp = KIND(1.0D0)

    INTEGER, INTENT(IN) :: n_particles, Max_delay, numero_specie, &
                          n_steps, Step_sal
    INTEGER, INTENT(IN) :: n_particles_specie(numero_specie), fact_nspecie
    REAL(KIND =dp), INTENT(IN) ::  deltat, Temperatura
    ! Output
    REAL, INTENT(INOUT) :: VACF_singolo(:,:)       ! Autocorrelazione singolo
    REAL, INTENT(OUT), ALLOCATABLE:: D_singolo(:)
    REAL, INTENT(INOUT) :: VACF_cross(:,:)  ! Correlazioni 'incrociate'
    REAL, INTENT(OUT), ALLOCATABLE :: D_cross(:)

    ! Variabili locali
    REAL, ALLOCATABLE :: save_velocities(:,:,:)
    REAL :: norm_singolo
    REAL, ALLOCATABLE :: Norm_cross(:,:)
    REAL :: sum_v2_dir(3) 
    INTEGER :: n_saved_steps, delay, time, time1, time2, dir, ios, step, p, s, &
              start_idx_s1, start_idx_s2, end_idx_s1, end_idx_s2, s1, s2, dim, k, &
              p_altrap, timet, pp, Block_size ! per blocking

    REAL    :: v1x, v1y, v1z, v1x_d, v1y_d, v1z_d
    REAL    :: v2x, v2y, v2z, v2x_d, v2y_d, v2z_d
    INTEGER :: Block_size_t    
    INTEGER ::  start_idx_s, end_idx_s, p_altra, p_solvente
    INTEGER :: start_solvente, end_solvente, start_altra, end_altra
    REAL    :: start_time, end_time, elapsed_time
    INTEGER :: count_rate, count_start, count_end
    ! Calcola il numero di frame salvati
    n_saved_steps = n_steps/Step_sal

    ALLOCATE(D_Cross(fact_nspecie))
    ALLOCATE(Norm_cross(Max_delay, fact_nspecie))
    ALLOCATE(D_singolo(numero_specie))

    ! Verifica che Max_delay sia ragionevole
    IF (Max_delay >= n_saved_steps) THEN
      PRINT *, "ERRORE: Max_delay troppo grande rispetto ai frame salvati"
      PRINT *, "Max_delay = ", Max_delay, " n_saved_steps = ", n_saved_steps
      STOP
    END IF

    ALLOCATE(save_velocities(3,n_particles, n_saved_steps), STAT=ios)
    IF (ios /= 0) THEN
      PRINT *, "Errore nell'allocazione di save_velocities"
      STOP
    END IF

    ! Lettura delle velocità
    OPEN(UNIT=20, FILE='Velocità.dat', STATUS='OLD', ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
      PRINT *, "Errore nell'aprire il file Velocità.dat"
      STOP
    END IF

    DO step = 1, n_saved_steps
      DO p = 1, n_particles
        READ(20, *, IOSTAT=ios) save_velocities(1,P,step), &
                               save_velocities(2,p, step), &
                               save_velocities(3,p, step)
        IF (ios /= 0) THEN
          PRINT *, "Errore lettura al passo ", step, " particella ", p
          CLOSE(20)
          STOP
        END IF
      END DO
    END DO
    CLOSE(20)

 ! Inizializzazione output
  VACF_singolo = 0.0_dp
  VACF_cross = 0.0_dp
  D_singolo = 0.0_dp
  D_Cross = 0.0_dp
  norm_cross = 0.0_dp

  ! 1. Calcolo autocorrelazione solvente con normalizzazione 
  PRINT *, "Calcolo VACF solvente..."
  start_idx_s = 1 
  DO s = 1, numero_specie   
    end_idx_s = start_idx_s + n_particles_specie(s) - 1
    DO delay = 1, Max_delay ! calcolo come somma su tutti i ritardi
    sum_v2_dir = 0.0_dp
      DO time = 1, n_saved_steps - delay ! medio anche su più origini
        DO P = start_idx_s, end_idx_s
          sum_v2_dir(:) = sum_v2_dir(:) + &
          (save_velocities(:, P, time) * &
          save_velocities(:, P, time + delay))
        END DO
      END DO
      norm_singolo = 3.0_dp * n_particles_specie(s) * (n_saved_steps - delay)
      VACF_singolo(delay, s) = SUM(sum_v2_dir)/norm_singolo
    END DO
    start_idx_s = end_idx_s + 1 
  END DO
  ! 2. Calcolo correlazioni incrociate

  ! Prima calcola tutte le norme
  start_idx_s1 = 1  
  DO s1 = 1, numero_specie - 1
    end_idx_s1 = start_idx_s1 + n_particles_specie(s1) - 1
    start_idx_s2 = end_idx_s1 + 1  
    DO s2 = s1 + 1, numero_specie
      end_idx_s2 = start_idx_s2 + n_particles_specie(s2) - 1
      k = (s1 - 1) * (2 * numero_specie - s1) / 2 + (s2 - s1)       
      DO delay = 1, Max_delay  
        Norm_cross(delay,k) = 2.0_dp * 3.0_dp * n_particles_specie(s1) * &
                           n_particles_specie(s2) * (n_saved_steps - delay)
      ENDDO 
    start_idx_s2 = end_idx_s2 + 1 
    END DO
  start_idx_s1 = end_idx_s1 + 1  
  END DO

  CALL SYSTEM_CLOCK(count_rate=count_rate)
  CALL SYSTEM_CLOCK(count=count_start)

  VACF_cross = 0.0d0
  start_idx_s1 = 1  
  Block_size_t = 128

  DO s1 = 1, numero_specie - 1
    end_idx_s1 = start_idx_s1 + n_particles_specie(s1) - 1
    start_idx_s2 = end_idx_s1 + 1
    DO s2 = s1 + 1, numero_specie
      end_idx_s2 = start_idx_s2 + n_particles_specie(s2) - 1
      k = (s1 - 1) * (2 * numero_specie - s1) / 2 + (s2 - s1)   
      PRINT*,"Processo la coppia", k, "s1:", s1, "s2:", s2
      !save velocities la lascio shared perchè fuori dal  registro openmp è shared
      !$OMP PARALLEL DO PRIVATE(delay, timet, time, p, p_altra, &
      !$OMP v1x, v1y, v1z, v1x_d, v1y_d, v1z_d, v2x, v2y, v2z, v2x_d, v2y_d, v2z_d) &
      !$OMP SHARED(start_idx_s1, end_idx_s1, start_idx_s2, end_idx_s2, &
      !$OMP save_velocities, k, n_saved_steps, Block_size_t) &
      !$OMP REDUCTION(+:VACF_cross) SCHEDULE(dynamic)
      DO delay = 1, Max_delay  
      DO timet = 1, n_saved_steps - delay, Block_size_t
        DO time = timet, MIN(timet + Block_size_t -1, n_saved_steps - delay)  
           DO p = start_idx_s1, end_idx_s1
            v1x = save_velocities(1, p, time)
            v1y = save_velocities(2, p, time)
            v1z = save_velocities(3, p, time)
            v1x_d = save_velocities(1, p, time + delay)
            v1y_d = save_velocities(2, p, time + delay)
            v1z_d = save_velocities(3, p, time + delay)        
            DO p_altra = start_idx_s2, end_idx_s2
              v2x = save_velocities(1, p_altra, time)
              v2y = save_velocities(2, p_altra, time)
              v2z = save_velocities(3, p_altra, time)
              v2x_d = save_velocities(1, p_altra, time + delay)
              v2y_d = save_velocities(2, p_altra, time + delay)
              v2z_d = save_velocities(3, p_altra, time + delay)

                  VACF_cross(delay,k) = VACF_cross(delay,k) + &
                              (v1x * v2x_d) + (v1x_d * v2x) + & ! Componente X
                              (v1y * v2y_d) + (v1y_d * v2y) + & ! Componente Y
                              (v1z * v2z_d) + (v1z_d * v2z)    ! Componente Z
            ENDDO
            ENDDO
        ENDDO
      ENDDO  
      ENDDO   
    !$OMP END PARALLEL DO   
    start_idx_s2 = end_idx_s2 + 1
    END DO
  start_idx_s1 = end_idx_s1 + 1
  END DO
  VACF_cross = VACF_cross/norm_cross
  CALL SYSTEM_CLOCK(count=count_end)
  elapsed_time = REAL(count_end - count_start) / REAL(count_rate)
  PRINT '("WALL TIME VACF= ", F10.3, " secondi")', elapsed_time

  D_singolo = 0.0_dp
  D_Cross = 0.0_dp


  DO s = 1, numero_specie  
    DO delay = 1,Max_delay -1
      D_singolo(s) = D_singolo(s) + ((VACF_singolo(delay+1, s) + VACF_singolo(delay, s))/2)*deltat
    ENDDO  
  ENDDO

  D_singolo = D_singolo/3
  DO delay = 1, Max_delay -1
    D_cross(:) = D_cross(:) + ((VACF_cross(delay+1, :) + VACF_cross(delay,:))/2)*deltat
  ENDDO 
  D_cross(:) = D_cross(:)/3


  !Deallocazione
  DEALLOCATE(save_velocities)
  PRINT *, "Calcolo VACF completato con successo!"
    
  END SUBROUTINE Calcola_Vacf
END MODULE Vacf
