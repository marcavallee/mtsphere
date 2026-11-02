SUBROUTINE csphjy(n, z, nm, csj, cdj, csy, cdy)
    ! ===================================================================
    ! Objectif : Calculer les fonctions de Bessel sphériques j_n(z) et
    !            y_n(z) ainsi que leurs dérivées pour un argument complexe z.
    !
    ! Arguments :
    !    n   (In)  : Ordre maximal souhaité (Entier)
    !    z   (In)  : Argument de la fonction (Complexe double précision)
    !    nm  (Out) : Ordre maximal effectivement calculé avec succès
    !    csj (Out) : Tableau (0:n) contenant j_m(z)
    !    cdj (Out) : Tableau (0:n) contenant j'_m(z)
    !    csy (Out) : Tableau (0:n) contenant y_m(z)
    !    cdy (Out) : Tableau (0:n) contenant y'_m(z)
    ! ===================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    COMPLEX(8), INTENT(IN) :: z
    INTEGER, INTENT(OUT) :: nm
    COMPLEX(8), DIMENSION(0:n), INTENT(OUT) :: csj, cdj, csy, cdy

    COMPLEX(8) :: cj0, cj1, cy0, cy1, cs, cs0, cs1, cf, cf0, cf1
    REAL(8) :: a0
    INTEGER :: k, m, mtop
    
    INTEGER, EXTERNAL :: msta2

    nm = n
    a0 = ABS(z)
    
    ! --- Cas d'erreur ou de valeur limite proche de 0 ---
    IF (a0 < 1.0d-100) THEN
        csj(0) = (1.0d0, 0.0d0)
        cdj(0) = (0.0d0, 0.0d0)
        csy(0) = (-1.0d0 / 0.0d0, 0.0d0) ! Représente l'infini négatif
        cdy(0) = (1.0d0 / 0.0d0, 0.0d0)
        DO k = 1, n
            csj(k) = (0.0d0, 0.0d0)
            cdj(k) = (0.0d0, 0.0d0)
            csy(k) = (-1.0d0 / 0.0d0, 0.0d0)
            cdy(k) = (1.0d0 / 0.0d0, 0.0d0)
        END DO
        RETURN
    END IF

    ! --- Étape 1 : Calcul de j_n(z) par récurrence descendante ---
    csj(0) = SIN(z) / z
    csj(1) = (csj(0) - COS(z)) / z
    
    mtop = msta2(a0, n, 15) ! Recherche du point de départ optimal via msta2
    IF (mtop < n) THEN
        nm = mtop
    ELSE
        mtop = n
    END IF

    ! Initialisation de la récurrence descendante (Algorithme de Miller)
    cf1 = (0.0d0, 0.0d0)
    cf0 = (1.0d-100, 0.0d0)
    cf = (0.0d0, 0.0d0)
    
    DO k = mtop, 0, -1
        cf = (2.0d0 * k + 3.0d0) * cf0 / z - cf1
        IF (k <= nm) csj(k) = cf0
        cf1 = cf0
        cf0 = cf
    END DO

    ! Normalisation des valeurs de j_n(z) obtenues
    a0 = ABS(csj(0))
    IF (a0 >= ABS(csj(1))) THEN
        cs = SIN(z) / (z * csj(0))
    ELSE
        cs = (SIN(z) / z - COS(z)) / (z * csj(1))
    END IF
    
    DO k = 0, nm
        csj(k) = csj(k) * cs
    END DO

    ! --- Étape 2 : Calcul de y_n(z) par récurrence ascendante ---
    cy0 = -COS(z) / z
    cy1 = (cy0 - SIN(z)) / z
    csy(0) = cy0
    csy(1) = cy1
    
    DO k = 2, nm
        csy(k) = (2.0d0 * k - 1.0d0) * csy(k-1) / z - csy(k-2)
    END DO

    ! --- Étape 3 : Calcul des dérivées premières j'_n(z) et y'_n(z) ---
    cdj(0) = (COS(z) - SIN(z) / z) / z
    DO k = 1, nm
        cdj(k) = csj(k-1) - (k + 1.0d0) * csj(k) / z
    END DO
    
    cdy(0) = (SIN(z) + COS(z) / z) / z
    DO k = 1, nm
        cdy(k) = csy(k-1) - (k + 1.0d0) * csy(k) / z
    END DO

    RETURN
END SUBROUTINE csphjy
