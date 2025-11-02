  Mtot = 0.0d0
  start_idx_s = 1
  DO s = 1, numero_specie
    sigma_specie(s) = SQRT(kB*Temperatura/massa_specie(s))
    end_idx_s = start_idx_s + n_particles_specie(s) - 1
    DO p = start_idx_s, end_idx_s 
      CALL RANDOM_NUMBER(u1)
      CALL RANDOM_NUMBER(u2)
      Gaussian_random = SQRT(-2.0d0 * LOG(u1)) * COS(2.0 * 3.141592653589793 * u2)
      velocities(:,p) = Gaussian_random * sigma_specie(s)
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
  scaling_factor = SQRT(3.0d0 * kB * n_particles*Temperatura / sumsqt)
  ! Applico il fattore di scala
  velocities = velocities * scaling_factor
