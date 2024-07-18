SUBROUTINE ZBESJ_F77(FNU,ZR,ZI,CYR,CYI)
IMPLICIT NONE
REAL(8),INTENT(IN) :: FNU
REAL(8),INTENT(IN) :: ZR, ZI
INTEGER :: KODE=1
REAL(8), DIMENSION(1),INTENT(OUT) :: CYR, CYI 
INTEGER :: NZ, IERR
CALL ZBESJ(ZR, ZI, FNU, KODE, 1, CYR, CYI, NZ, IERR)
END SUBROUTINE ZBESJ_F77

SUBROUTINE ZBESY_F77(FNU,ZR,ZI,CYR,CYI)
IMPLICIT NONE
REAL(8),INTENT(IN) :: FNU
REAL(8),INTENT(IN) :: ZR, ZI
INTEGER :: KODE=1
REAL(8), DIMENSION(1),INTENT(OUT) :: CYR, CYI 
INTEGER :: NZ, IERR
COMPLEX(8), DIMENSION(1) :: CWRKR, CWRKI 
CALL ZBESY(ZR, ZI, FNU, KODE, 1, CYR, CYI, NZ, CWRKR, CWRKI, IERR)
END SUBROUTINE ZBESY_F77

SUBROUTINE ZBESI_F77(FNU,ZR,ZI,CYR,CYI)
IMPLICIT NONE
REAL(8),INTENT(IN) :: FNU
REAL(8),INTENT(IN) :: ZR, ZI
INTEGER :: KODE=1
REAL(8), DIMENSION(1),INTENT(OUT) :: CYR, CYI 
INTEGER :: NZ, IERR
CALL ZBESI(ZR, ZI, FNU, KODE, 1, CYR, CYI, NZ, IERR)
END SUBROUTINE ZBESI_F77

SUBROUTINE ZBESK_F77(FNU,ZR,ZI,CYR,CYI)
IMPLICIT NONE
REAL(8),INTENT(IN) :: FNU
REAL(8),INTENT(IN) :: ZR, ZI
INTEGER :: KODE=1
REAL(8), DIMENSION(1),INTENT(OUT) :: CYR, CYI 
INTEGER :: NZ, IERR
CALL ZBESK(ZR, ZI, FNU, KODE, 1, CYR, CYI, NZ, IERR)
END SUBROUTINE ZBESK_F77

SUBROUTINE ZBESH_F77(K,FNU,ZR,ZI,CYR,CYI)
IMPLICIT NONE
REAL(8),INTENT(IN) :: FNU
INTEGER, INTENT(IN) :: K
REAL(8),INTENT(IN) :: ZR, ZI
INTEGER :: KODE=1
REAL(8), DIMENSION(1),INTENT(OUT) :: CYR, CYI 
INTEGER :: NZ, IERR
CALL ZBESH(ZR, ZI, FNU, KODE, K, 1, CYR, CYI, NZ, IERR)
END SUBROUTINE ZBESH_F77



