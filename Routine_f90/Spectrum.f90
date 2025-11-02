MODULE SPECTRUM 
  IMPLICIT NONE

CONTAINS

  SUBROUTINE DFT(filename)
 
  INTEGER, PARAMETER :: dp = KIND(1.0d0)
  CHARACTER (LEN=*), INTENT(IN) :: filename

  !.-.-.-. LOCALI -.-.-.-.
  CHARACTER(LEN=1000) ::  outputfilename, line
  INTEGER :: unit, ios, i, j, n_colonne, n_x, pos, omega, colonna, x, n, k
  COMPLEX(dp), ALLOCATABLE :: spectrum(:,:)
  REAL(dp), ALLOCATABLE :: Dati_in(:,:), Dati_out(:,:), &
                           power_spectrum(:,:), window(:), ind_var(:)
  REAL(dp) :: freq, dt, freq_camp 
  COMPLEX(dp) :: imaginary_unit
  REAL(KIND = dp), PARAMETER :: pi = ACOS(-1.0d0)
  


  OPEN(unit=20, file= TRIM(filename), status='OLD', action='READ', iostat =ios) ! File per risultati 
  n_x = 0
  DO
    READ(20, '(A)', iostat=ios) line
      IF (ios < 0) THEN  ! Fine del file raggiunta
      EXIT
        ELSE IF (ios > 0) THEN  ! Errore di lettura
          PRINT *, 'Errore durante la lettura del file'
          STOP 
      ENDIF
    n_x = n_x + 1
  ENDDO

  n_colonne = 0
  ! Conta le colonne, assumo che i dati siano scritti con dei punti a separare i decimali
  ! Non è evidentemente la migliore delle opzioni ma non mi veniva in mente altro
  DO i = 1, LEN_TRIM(line)
     IF (line(i:i) == '.') THEN
       n_colonne = n_colonne + 1
     END IF
  END DO

  rewind(20)
  
  ALLOCATE(ind_var(n_x), Dati_in(n_x, n_colonne-1))
  ALLOCATE(Dati_out(n_x, n_colonne-1), power_spectrum(n_x, n_colonne-1), &
           spectrum(n_x, n_colonne -1), window(n_x))

  DO x = 1,n_x
    READ(20,*,iostat = ios) ind_var(x), (Dati_in(x, colonna), colonna = 1, n_colonne-1)  
  ENDDO 

  CLOSE(20) 

    spectrum = (0.0_dp, 0.0_dp)
    power_spectrum = 0.0_dp
    imaginary_unit = (0.0d0, 1.0d0)
    dt = ind_var(2) - ind_var(1)
    freq_camp = 1.0_dp / dt
    
    ! Rimozione media  
    DO colonna = 1, n_colonne-1
      Dati_in(:, colonna) = Dati_in(:, colonna) - SUM(Dati_in(:, colonna)) / REAL(n_x, dp)
    ENDDO

    ! DFT 
    DO colonna = 1, n_colonne-1
      DO omega = 1, n_x  ! omega = indice frequenza 
        k = omega - 1    ! k = indice matematico
        DO x = 1, n_x    ! x = indice tempo 
          n = x - 1      ! n = indice matematico (0-based)  
          spectrum(omega, colonna) = spectrum(omega, colonna) + &
          Dati_in(x, colonna) * &
          EXP(-imaginary_unit * 2.0_dp * pi * &
          REAL(k, dp) * REAL(n, dp) / REAL(n_x, dp))
        END DO
        power_spectrum(omega, colonna) = ABS(spectrum(omega, colonna))**2
      END DO
    ENDDO

 
 ! Trova la posizione del punto
  pos = INDEX(filename, '.', BACK=.TRUE.)  ! Cerca dall'indietro
  outputfilename = filename(1:pos-1) // "_DFT.dat"
  
  ! Scrittura file spettro, freq = indice * freq_campionamento/ n_tot_punti
  OPEN(unit=30, file=outputfilename, status='replace', action='write')
  DO omega = 1, n_x/2 +1
    freq = REAL(omega -1, dp) *freq_camp / (REAL(n_x, dp))
    WRITE(30, '(F15.12, 5X, 4(5X, E25.8))') freq, (power_spectrum(omega, colonna), colonna = 1, n_colonne-1)
  END DO
  CLOSE(30)
 

  DEALLOCATE(power_spectrum, window, Dati_in, Dati_out)
END SUBROUTINE DFT

END MODULE SPECTRUM