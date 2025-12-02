!=====================================================================
!     HILL48 PLASTICITY UMAT WITH VOCE ISOTROPIC HARDENING
!     IMPLEMENTATION: GENERALIZED RADIAL-RETURN IN EIGENSPACE
!===============================================================================
!
!  AUTHOR: 
!    Mohammad Hasaninia
!    Computational Advanced Manufacturing and Materials Laboratory (CAMML)
!    Department of Mechanical Engineering
!    University of Wyoming
!
!  REFERENCE:
!    Versino, D. and Bennett, K.C. (2018). Generalized radial-return mapping 
!    algorithm for anisotropic von Mises plasticity framed in material 
!    eigenspace. Int. J. Numer. Meth. Engng, 116: 202-222.
!    https://doi.org/10.1002/nme.5921
!
!-------------------------------------------------------------------------------
!  INPUT PROPERTIES (PROPS) DEFINITION:
!    PROPS(1)  : Young's Modulus (E)
!    PROPS(2)  : Poisson's Ratio (nu)
!    PROPS(3)  : Initial Yield Stress (Sigma_Y0)
!    PROPS(4)  : Saturation Stress (Q_inf) - Voce Parameter
!    PROPS(5)  : Hardening Rate (b)        - Voce Parameter
!    PROPS(6)  : Hill F
!    PROPS(7)  : Hill G
!    PROPS(8)  : Hill H
!    PROPS(9)  : Hill L
!    PROPS(10) : Hill M
!    PROPS(11) : Hill N
!
!  STATE VARIABLES (STATEV) OUTPUT:
!    STATEV(1)   : Equivalent Plastic Strain (alpha)
!    STATEV(2)   : Current Yield Stress
!    STATEV(3-8) : Plastic Strain Tensor (11, 22, 33, 12, 13, 23)
!    STATEV(9)   : Convergence Flag (0=Converged, 1=Failed)
!=====================================================================

SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,           &
                DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME,   &
                TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, &
                NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT,          &
                CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT,       &
                KSTEP, KINC)

    IMPLICIT NONE
    
    ! Abaqus UMAT interface variables
    CHARACTER*80 CMNAME
    INTEGER NDI, NSHR, NTENS, NSTATV, NPROPS
    INTEGER NOEL, NPT, LAYER, KSPT, KSTEP, KINC
    
    REAL(8) STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS)
    REAL(8) SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, DTEMP, PNEWDT, CELENT
    REAL(8) STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1)
    REAL(8) PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
    REAL(8) DDSDDT(NTENS), DRPLDE(NTENS)
    
    ! Local variables for plasticity computation
    REAL(8) :: alpha_tr, alpha, StressYield, s_y, s_y_D, s_y_tr, f_gamma_tr
    REAL(8) :: Sigma_Y0, Q_inf, b_voce
    REAL(8) :: MatM(6,6), MatD(6,6), MatI(6,6)
    REAL(8) :: MatM_copy(6,6), Mat_inv(6,6)
    REAL(8) :: W(6), Q(6,6), Gamma(6,6), K(6,6), dK_dDelta_gamma_matrix(6,6)
    REAL(8) :: s_tr(6), strain_p(6), Stress_New(6), s_tilde_tr(6)
    REAL(8) :: Vec(6), dev_corrected(6)
    REAL(8) :: E, XNUE, EBULK3, EG2, EG, ELAM
    REAL(8) :: HillF, HillG, HillH, HillL, HillM, HillN
    REAL(8) :: pressure, d_ga, omega, d_al, f_gamma_delta, delta_alpha
    REAL(8) :: Df_gamma_d_d_gamma, d_omega_d_d_gamma, d_sy_d_d_gamma
    REAL(8) :: TOLERANCE, exp_term
    INTEGER :: iteration, K3, K5, K6, K7, K8, INFO
    LOGICAL :: converged
    
    ! Extract material properties
    E = PROPS(1)              ! Young's modulus
    XNUE = PROPS(2)           ! Poisson's ratio
    Sigma_Y0 = PROPS(3)       ! Initial yield stress
    Q_inf = PROPS(4)          ! Saturation stress (Q∞)
    b_voce = PROPS(5)         ! Voce hardening rate parameter
    HillF = PROPS(6)          ! Hill parameter F
    HillG = PROPS(7)          ! Hill parameter G
    HillH = PROPS(8)          ! Hill parameter H
    HillL = PROPS(9)          ! Hill parameter L
    HillM = PROPS(10)         ! Hill parameter M
    HillN = PROPS(11)         ! Hill parameter N
    
    ! Initialize state variables
    IF (TIME(2) < 1.0D-6) THEN
        alpha_tr = 0.0D0
        ! Voce: σ_y0 = σ_Y0 + Q∞(1 - e^0) = σ_Y0
        StressYield = Sigma_Y0
        strain_p = 0.0D0
        STATEV(1:9) = 0.0D0
    ELSE
        alpha_tr = STATEV(1)
        StressYield = STATEV(2)
        strain_p(1:6) = STATEV(3:8)
    END IF
    
    ! Create elasticity matrix
    EBULK3 = E/(1.0D0-2.0D0*XNUE)
    EG2 = E/(1.0D0+XNUE)
    EG = EG2/2.0D0
    ELAM = (EBULK3-EG2)/3.0D0
    
    MatD = 0.0D0
    MatD(1:3,1:3) = ELAM
    DO K3 = 1, 3
        MatD(K3,K3) = EG2 + ELAM
    END DO
    DO K3 = 4, 6
    MatD(K3, K3) = EG
    END DO
    
    ! Create Hill anisotropy matrix
    MatM = 0.0D0
    MatM(1,1) = HillG + HillH
    MatM(2,2) = HillF + HillH
    MatM(3,3) = HillF + HillG
    MatM(1,2) = -HillH
    MatM(1,3) = -HillG
    MatM(2,3) = -HillF
    MatM(2,1) = MatM(1,2)
    MatM(3,1) = MatM(1,3)
    MatM(3,2) = MatM(2,3)
    MatM(4,4) = 2.0D0*HillN
    MatM(5,5) = 2.0D0*HillM
    MatM(6,6) = 2.0D0*HillL

    ! Identity matrix
    MatI = 0.0D0
    DO K5 = 1, 6
        MatI(K5,K5) = 1.0D0
    END DO
    
    ! Compute trial stress (elastic predictor)
    Vec = MATMUL(MatD, DSTRAN)
    Stress_New = Vec + STRESS
    s_tr = Stress_New
    
    ! Extract deviatoric part
    pressure = (Stress_New(1)+Stress_New(2)+Stress_New(3))/3.0D0
    s_tr(1:3) = s_tr(1:3) - pressure
    
    ! Eigendecomposition of anisotropy tensor
    MatM_copy = MatM
    CALL compute_eigendecomposition(MatM_copy, W, Q, Gamma, INFO)

    IF (INFO /= 0) THEN
        WRITE(*,*) 'Error in eigendecomposition, INFO = ', INFO
        RETURN
    END IF
    
    ! Transform trial stress to eigenspace
    s_tilde_tr = MATMUL(TRANSPOSE(Q), s_tr)
    
    ! Initialize trial values
    alpha_tr = STATEV(1)
    s_y_tr = StressYield
    
    ! Compute trial yield function
    f_gamma_tr = 0.5D0 * DOT_PRODUCT(s_tilde_tr, MATMUL(Gamma, s_tilde_tr)) - s_y_tr**2
    

    ! Check yield condition
    IF (f_gamma_tr <= 0.0D0) THEN
        ! Elastic step
        STRESS = Stress_New
        d_ga = 0.0D0
        alpha = alpha_tr
        s_y = StressYield
    ELSE
        ! Plastic step - Newton-Raphson iteration
        tolerance = 5.0D-4
        d_ga = 0.0D0
        alpha = alpha_tr
        converged = .FALSE.
        
        DO iteration = 1, 250
            ! Build K matrix (diagonal)
            K = 0.0D0
            DO K6 = 1, 6
                K(K6,K6) = Gamma(K6,K6) / (1.0D0 + 2.0D0*EG*d_ga*Gamma(K6,K6))**2
            END DO
            
            ! Compute auxiliary variable omega
            omega = SQRT(0.5D0 * DOT_PRODUCT(s_tilde_tr, MATMUL(K, s_tilde_tr)))
            
            IF (omega <= 1.0D-14) THEN
                WRITE(*,*) 'Warning: omega is too small'
                EXIT
            END IF
            
            ! Update alpha
            delta_alpha = 2.0D0 * d_ga * omega
            alpha = alpha_tr + delta_alpha
            
            ! ========================================================
            ! VOCE HARDENING MODEL
            ! σ_y = σ_Y0 + Q∞(1 - exp(-b·α))
            ! ========================================================
            exp_term = EXP(-b_voce * alpha)
            s_y = Sigma_Y0 + Q_inf * (1.0D0 - exp_term)
            
            ! Derivative: ∂σ_y/∂α = Q∞·b·exp(-b·α)
            s_y_D = Q_inf * b_voce * exp_term
            ! ========================================================
            
            ! Compute residual function
            f_gamma_delta = (s_y / omega) - 1.0D0
            
            ! Check convergence
            IF (ABS(f_gamma_delta) <= tolerance) THEN
                StressYield = s_y
                converged = .TRUE.
                EXIT
            END IF
            
            ! Compute derivatives for Newton-Raphson
            dK_dDelta_gamma_matrix = 0.0D0
            DO K7 = 1, 6
                dK_dDelta_gamma_matrix(K7,K7) = -4.0D0*EG*Gamma(K7,K7)**2 &
                    / (1.0D0 + 2.0D0*EG*d_ga*Gamma(K7,K7))**3
            END DO
            
            d_omega_d_d_gamma = (1.0D0 / (4.0D0 * omega)) &
                * DOT_PRODUCT(s_tilde_tr, MATMUL(dK_dDelta_gamma_matrix, s_tilde_tr))
            
            ! Chain rule: ∂σ_y/∂Δγ = (∂σ_y/∂α)·(∂α/∂Δγ)
            d_sy_d_d_gamma = 2.0D0 * s_y_D * (omega + d_ga * d_omega_d_d_gamma)
            
            ! Derivative of residual
            Df_gamma_d_d_gamma = (1.0D0 / omega) * d_sy_d_d_gamma &
                - (s_y / (omega**2)) * d_omega_d_d_gamma
            
            ! Newton-Raphson update
            d_ga = d_ga - f_gamma_delta / Df_gamma_d_d_gamma
            
        END DO
        
        ! Update stress after Newton loop
        Mat_inv = 0.0D0
        DO K8 = 1, 6
            Mat_inv(K8,K8) = 1.0D0 / (1.0D0 + 2.0D0*EG*d_ga*Gamma(K8,K8))
        END DO
        
        dev_corrected = MATMUL(Q, MATMUL(Mat_inv, MATMUL(TRANSPOSE(Q), s_tr)))
        
        ! Assemble full Cauchy stress
        STRESS(1:3) = dev_corrected(1:3) + pressure
        STRESS(4:6) = dev_corrected(4:6)
        
        IF (.NOT. converged) THEN
            WRITE(*,*) 'WARNING: Newton-Raphson did not converge!'
            StressYield = s_y
            STATEV(9) = 1.0D0
        END IF
    END IF
    
    ! Consistent tangent operator
    DDSDDE = MatD
    
    ! Update state variables
    STATEV(1) = alpha                    ! Equivalent plastic strain
    STATEV(2) = s_y                      ! Current yield stress
    strain_p = strain_p + d_ga * MATMUL(MatM, dev_corrected)
    STATEV(3:8) = strain_p(1:6)         ! Plastic strain components
    
    RETURN
    
    CONTAINS
    
    !**********************************************************************
    SUBROUTINE compute_eigendecomposition(A, W, Q, Gamma, INFO)
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: INFO
        REAL(8), DIMENSION(6,6), INTENT(INOUT) :: A
        REAL(8), DIMENSION(6), INTENT(OUT) :: W
        REAL(8), DIMENSION(6,6), INTENT(OUT) :: Q
        REAL(8), DIMENSION(6,6), INTENT(OUT) :: Gamma
        INTEGER :: LWORK, i, j
        REAL(8), DIMENSION(:), ALLOCATABLE :: WORK
        LOGICAL :: is_symmetric
        
        ! Check symmetry
        is_symmetric = .TRUE.
        DO i = 1, 6
            DO j = i + 1, 6
                IF (ABS(A(i,j) - A(j,i)) > 1.0D-12) THEN
                    is_symmetric = .FALSE.
                    EXIT
                END IF
            END DO
            IF (.NOT. is_symmetric) EXIT
        END DO
        
        IF (.NOT. is_symmetric) THEN
            WRITE(*,*) "Error: Matrix not symmetric"
            INFO = -1
            RETURN
        END IF
        
        ! LAPACK eigendecomposition
        LWORK = -1
        ALLOCATE(WORK(1))
        CALL DSYEV('V', 'U', 6, A, 6, W, WORK, LWORK, INFO)
        LWORK = INT(WORK(1))
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
        
        CALL DSYEV('V', 'U', 6, A, 6, W, WORK, LWORK, INFO)
        
        IF (INFO == 0) THEN
            Q = A
            Gamma = 0.0D0
            DO i = 1, 6
                Gamma(i,i) = W(i)
            END DO
        ELSE
            WRITE(*,*) "Eigendecomposition failed, INFO:", INFO
            Gamma = 0.0D0
        END IF
        
        DEALLOCATE(WORK)
        
    END SUBROUTINE compute_eigendecomposition
    
END SUBROUTINE UMAT
!*****************************BOTTOM*****************************