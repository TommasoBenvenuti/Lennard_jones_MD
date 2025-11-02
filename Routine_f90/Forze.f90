MODULE Forze
IMPLICIT NONE
CONTAINS

SUBROUTINE Calcola_forze(positions, forces, n_particles, n_particles_specie, numero_specie, &
                        Sigma_Lennard_Jones, Epsilon_Lennard_Jones, massa_specie, en_pot, box_size, rcut2, viriale, &
                        gr_intra, gr_inter, bin_size, nbins)
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)
  
  ! Dichiarazione argomenti
  REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)    :: positions 
  REAL(KIND=dp), DIMENSION(:,:), INTENT(OUT)   :: forces
  REAL(KIND=dp), INTENT(INOUT)                 :: en_pot
  INTEGER, INTENT(IN)                          :: n_particles, numero_specie, nbins
  INTEGER, DIMENSION(:), INTENT(IN)            :: n_particles_specie
  REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)    ::  Sigma_Lennard_Jones, Epsilon_Lennard_Jones
  REAL(KIND=dp), DIMENSION(:), INTENT(IN)      :: massa_specie
  REAL(KIND=dp), DIMENSION(3), INTENT(IN)      :: box_size
  REAL(KIND=dp), INTENT(OUT)                   :: viriale
  REAL(KIND=dp), INTENT(IN)                    :: rcut2, bin_size
  REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: gr_intra, gr_inter

  ! Variabili locali
  INTEGER :: s1, s2, p1, p2, k, l, i 
  INTEGER :: start_idx_s1, end_idx_s1, start_idx_s2, end_idx_s2, counter
  REAL(KIND=dp) :: sigma, epsilon, r2, r2i, r6i, f_magnitude, r, density, norm, Volume_bin
  REAL(KIND=dp) :: sigma6, shift, rcut6
  REAL(KIND=dp), DIMENSION(3) :: dk
  REAL(KIND=dp), PARAMETER :: kB = 1.380649D-23 * 6.02214076D23*1D3 * 1D20 / 1D30  ! Conversione K -> uma Å²/fs²
  REAL(KIND = dp ) :: pi = 3.141592653589793d0

  ! Inizializzazione
  forces = 0.0_dp
  en_pot = 0.0_dp
  viriale = 0.0_dp
  rcut6 = rcut2**3
  gr_intra = 0.0_dp
  gr_inter = 0.0_dp

  ! ==============================================
  ! CALCOLO FORZE E RDF INTRA-SPECIE
  ! ==============================================
  start_idx_s1 = 1
  DO s1 = 1, numero_specie
    end_idx_s1 = start_idx_s1 + n_particles_specie(s1) - 1
    sigma = Sigma_Lennard_Jones(s1, s1)
    epsilon = Epsilon_Lennard_Jones(s1, s1)
    sigma6 = sigma**6
    shift = 4.0_dp * epsilon * (sigma6/rcut6) * (sigma6/rcut6 - 1.0_dp)
      DO p1 = start_idx_s1, end_idx_s1
      DO p2 = p1 + 1, end_idx_s1
        ! Calcolo distanza con condizioni periodiche
        dk(:) = positions(:, p1) - positions(:, p2)
        dk(:) = dk(:) - box_size(:) * NINT(dk(:)/box_size(:))
        r2 = dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3)
        IF (r2 < rcut2) THEN
          r = SQRT(r2)
          l = INT(r / bin_size) ! +1
          gr_intra(s1, l) = gr_intra(s1, l) + 2.0_dp  ! +2 per coppie simmetriche

          r2i = 1.0_dp/r2
          r6i = r2i*r2i*r2i
          f_magnitude = 48.0_dp * epsilon * r2i * sigma6 * r6i * (sigma6 * r6i - 0.5_dp)

          forces(:, p1) = forces(:, p1) + f_magnitude * dk(:)
          forces(:, p2) = forces(:, p2) - f_magnitude * dk(:)
          en_pot = en_pot + 4.0_dp * epsilon * sigma6 * r6i * (sigma6 * r6i - 1.0_dp) - shift
        !  viriale = viriale + f_magnitude * r2
        END IF
      END DO
    END DO
    start_idx_s1 = end_idx_s1 + 1
  END DO
  counter= 0
  ! ==============================================
  ! CALCOLO FORZE E RDF INTER-SPECIE
  ! ==============================================
  start_idx_s1 = 1
  DO s1 = 1, numero_specie - 1
    end_idx_s1 = start_idx_s1 + n_particles_specie(s1) - 1
    start_idx_s2 = end_idx_s1 + 1
    DO s2 = s1 + 1, numero_specie
      end_idx_s2 = start_idx_s2 + n_particles_specie(s2) - 1
      k = (s1 - 1) * (2 * numero_specie - s1) / 2 + (s2 - s1)  ! Indice combinazioni
      sigma = Sigma_Lennard_Jones(s1, s2)
      epsilon = Epsilon_Lennard_Jones(s1, s2)
      sigma6 = sigma*sigma*sigma*sigma*sigma*sigma 
      shift = 4.0_dp * epsilon * (sigma6/rcut6) * (sigma6/rcut6 - 1.0_dp)
      DO p1 = start_idx_s1, end_idx_s1
        DO p2 = start_idx_s2, end_idx_s2
          dk(:) = positions(:, p1) - positions(:, p2)
          dk(:) = dk(:) - box_size(:) * NINT(dk(:)/box_size(:))
          r2 = SUM(dk(:)*dk(:))
          IF (r2 < rcut2) THEN
          counter = counter + 1
            r = SQRT(r2)
            l = INT(r / bin_size) ! + 1
            gr_inter(k, l) = gr_inter(k, l) + 1.0_dp ! conto +1 perchè A-B è diverso da B-A
            r2i = 1.0_dp/r2
            r6i = r2i*r2i*r2i 
            f_magnitude = 48.0_dp * epsilon * r2i * sigma6 * r6i * (sigma6 * r6i - 0.5_dp)
            forces(:, p1) = forces(:, p1) + f_magnitude * dk(:)
            forces(:, p2) = forces(:, p2) - f_magnitude * dk(:)
            en_pot = en_pot + 4.0_dp * epsilon * sigma6 * r6i * (sigma6 * r6i - 1.0_dp) - shift
         
          END IF
        END DO
      END DO
      start_idx_s2 = end_idx_s2 + 1
    END DO
    start_idx_s1 = end_idx_s1 + 1
  END DO

  viriale = viriale * 0.5_dp 

END SUBROUTINE Calcola_forze
END MODULE Forze
