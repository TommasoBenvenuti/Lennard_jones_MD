MODULE Integrazione

IMPLICIT NONE
CONTAINS
  ! =================================================================
  ! SUBROUTINE PER INTEGRAZIONE (ALGORITMO DI VERLET)
  ! =================================================================
  SUBROUTINE  Integrazione_equazioni_moto(positions, velocities, forces, n_particles, n_particles_specie, &
                       numero_specie, deltat, E_kin, Box_size, old_positions, new_positions, &
                       massa_specie, Temperatura_calc, En_pot, E_tot, step, vel_cm)
  IMPLICIT NONE
  ! --- Dichiarazione argomenti ---
  INTEGER, PARAMETER :: dp = KIND(1.0D0)
  REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: positions, velocities, forces
  REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: old_positions, new_positions 
  INTEGER, INTENT(IN)                          :: n_particles, numero_specie, step
  INTEGER, DIMENSION(:), INTENT(IN)            :: n_particles_specie
  REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: massa_specie
  REAL(KIND=dp), INTENT(IN)                    :: deltat
  REAL(KIND=dp), INTENT(INOUT)                 :: E_kin, Temperatura_calc
  REAL(KIND=dp), INTENT(IN)                    :: En_pot
  REAL(KIND=dp), INTENT(OUT)                   :: E_tot, vel_cm
  REAL(KIND=dp), DIMENSION(3), INTENT(IN)      :: Box_size
  ! --- Variabili locali ---
  REAL(KIND=dp)                                :: sumv(3), sumsq(3), sumsqt
  INTEGER                                      :: p, s, k, start_idx_s, end_idx_s, DOF
  REAL(KIND=dp)                                :: u1, u2, theta, r
  REAL(KIND=dp)                                :: E_kin_target, alpha, decay, & ! Bussi-Parinello thermostat
                                                  R1, sumR, Tau_t, Gaussian, factor, Mtot
  ! ==== Costanti ====
  REAL(KIND= dp) :: kB, kg_to_g, g_to_uma, m_to_am, s_to_fs
  kB = 1.380649D-23 ![kg m**2/s**2*K]
  kg_to_g = 1D3    
  g_to_uma = 6.02214076D23 
  m_to_am = 1D10
  s_to_fs = 1D15
  kB = (kB*kg_to_g*g_to_uma*m_to_am**2d0)/((s_to_fs)**2d0)
  sumv             = 0.0
  sumsq            = 0.0
  sumsqt           = 0.0
  E_kin            = 0.0d0
  Temperatura_calc = 0.0d0
  E_tot            = 0.0d0

  ! Verlet Algorithm
  start_idx_s = 1
  DO s = 1, numero_specie
    end_idx_s = start_idx_s + n_particles_specie(s) - 1
    DO p = start_idx_s, end_idx_s
      ! Algoritmo di Verlet
      new_positions(:,p) = 2*positions(:,p) - old_positions(:,p) + (forces(:,p) * deltat**2 / massa_specie(s))
      ! Calcolo velocità
      velocities(:,p) = (new_positions(:,p) - old_positions(:,p)) / (2.0d0 * deltat) 
      ! Aggiornamento
      sumsq(:) = sumsq(:) + (velocities(:,p)**2)*massa_specie(s)
      old_positions(:,p) = positions(:,p)
      positions(:,p) = new_positions(:,p)
      sumv(:) = sumv(:) + velocities(:,p)*massa_specie(s) ! per velocità centro di massa
    END DO    
    start_idx_s = end_idx_s +1
  END DO 

  Mtot = SUM(massa_specie(:) * n_particles_specie(:))
  sumv = sumv / Mtot   ! vel. cm

  vel_cm = SQRT(SUM(sumv**2))

  sumsqt = SUM(sumsq) ! calcolo energia totale cinetica per applicare termostato
  E_kin = 0.5d0 * sumsqt
  E_tot = E_kin + en_pot ! energia totale per particella

  Temperatura_calc = (2.0d0*E_kin) / (3.0d0 * n_particles*kB)


  END SUBROUTINE Integrazione_equazioni_moto

END MODULE Integrazione
