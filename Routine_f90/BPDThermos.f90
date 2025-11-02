MODULE NVT

IMPLICIT NONE
CONTAINS
SUBROUTINE BPD_Thermostat(velocities, n_particles, n_particles_specie, numero_specie, &
                                      massa_specie, deltat,  E_kin, Temperatura_calc, En_pot, E_tot, &
                                      Temperatura, Termalizzazione, step)
  ! REF. "Canonical sampling through velocity rescaling", 
  !G.Bussi, D. Donadio, M. Parrinello, J. Chem. Phys, 126, 2007
  IMPLICIT NONE 
  INTEGER, PARAMETER :: dp = KIND(1.0D0)
  REAL(KIND = dp), INTENT(INOUT) :: Velocities(:,:)
  REAL(KIND = dp), INTENT(IN)    :: Temperatura, deltat ! Per energia cinetica target
  REAL(KIND = dp), INTENT(INOUT) :: Temperatura_calc
  REAL(KIND = dp), INTENT(OUT)   :: E_kin
  REAL(KIND = dp), INTENT(INOUT) :: En_pot, E_tot ! Per calcolare l'energia totale in uscita 
  REAL(KIND = dp), INTENT(IN)    :: massa_specie(:)
  INTEGER, INTENT(IN)            :: n_particles, numero_specie, n_particles_specie(:), step
  LOGICAL                        :: Termalizzazione

  ! .-.-.-. Locali -.-.-.-.-.
  REAL(KIND = dp) :: E_kin_target, Tau_t, sumR, &
                     decay, alpha, factor, R1, Gaussian_random, &
                     sumsq(3), sumsqt, u1, u2, sumv(3), Mtot! Variabili reali per Termos. 
  INTEGER         :: DOF  ! Gradi di libertà 
  INTEGER         :: p, s, end_idx_s, start_idx_s ! contatori
   
  ! .-.-.- Costanti .-.-.-
  REAL(KIND= dp) :: kB, kg_to_g, g_to_uma, m_to_am, s_to_fs
  REAL(KIND=dp), PARAMETER       :: pi = ACOS(-1.0d0)
  kB       = 1.380649D-23 ![kg m**2/s**2*K]
  kg_to_g  = 1D3    
  g_to_uma = 6.02214076D23 
  m_to_am  = 1D10
  s_to_fs  = 1D15
  kB       = (kB*kg_to_g*g_to_uma*m_to_am**2d0)/((s_to_fs)**2d0)

  DOF = 3*n_particles - 3

  IF (Termalizzazione) THEN
    Tau_t = deltat * 250.0d0 ! Termostato "rapido" per equilibrazione, in fase di equilibrazione 
  ELSE 
    Tau_t = deltat*275.0d0  ! più lento per la fase di produzione 
  ENDIF    

  CALL RANDOM_NUMBER(u1)
  CALL RANDOM_NUMBER(u2)
  Gaussian_random = SQRT(-2.0d0 * LOG(u1)) * COS(2.0 * pi * u2)

  R1 = (Gaussian_random)
  sumR = 0.0_dp
  sumsq = 0.0_dp
  sumsqt = 0.0_dp
  alpha = 0.0_dp

  DO p = 1, DOF - 1 ! Somma di numeri distribuiti secondo Gaussiana, come riportato nel paper
    CALL RANDOM_NUMBER(u1)
    CALL RANDOM_NUMBER(u2)
    Gaussian_random = SQRT(-2.0d0 * LOG(u1)) * COS(2.0 * 3.141592653589793 * u2)
    sumR = sumR + (Gaussian_random)**2
  END DO

  decay = EXP(-deltat/Tau_t)
  factor = Temperatura / (Temperatura_calc*DOF)
  alpha = SQRT(decay + factor * (1.0_dp - decay) * ((R1)**2 + sumR) + &
          2.0_dp * SQRT(decay * factor * (1.0_dp - decay) )*R1)

  velocities = velocities * alpha
  start_idx_s = 1
  DO s = 1, numero_specie
    end_idx_s = start_idx_s + n_particles_specie(s) - 1
    DO p = start_idx_s, end_idx_s
      sumsq(:) = sumsq(:) + (velocities(:,p)**2)*massa_specie(s)
    END DO    
    start_idx_s = end_idx_s +1
  END DO 

  sumsqt = SUM(sumsq) ! ricalcolo energia totale cinetica dopo termostato, temperatura ecc..
  E_kin = 0.5d0 * sumsqt
  E_tot = En_pot + E_kin
  Temperatura_calc = (2.0d0*E_kin) / (3.0d0 * n_particles*kB)

  END SUBROUTINE BPD_Thermostat

END MODULE NVT
