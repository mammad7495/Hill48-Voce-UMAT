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
                JSTEP, KINC)

    IMPLICIT NONE

    ! Abaqus UMAT interface variables
    CHARACTER*80 CMNAME
    INTEGER NDI, NSHR, NTENS, NSTATV, NPROPS
    INTEGER NOEL, NPT, LAYER, KSPT, JSTEP(4), KINC

    REAL(8) STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS)
    REAL(8) SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, DTEMP, PNEWDT, CELENT
    REAL(8) STRAN(NTENS), DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1)
    REAL(8) PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
    REAL(8) DDSDDT(NTENS), DRPLDE(NTENS)


    real(8) :: phase_state(8), trial_stress(6)
    logical :: success

    if (NSTATV < 9) then
        write(*,*) 'Hill48 UMAT requires at least nine state variables.'
        call XIT
        return
    end if
    if (NDI /= 3 .or. NSHR /= 3 .or. NTENS /= 6 .or. NPROPS /= 11) then
        write(*,*) 'Hill48 UMAT requires a 3D stress state and 11 properties.'
        call XIT
        return
    end if

    phase_state(1:6) = STATEV(3:8)
    phase_state(7:8) = STATEV(1:2)
    trial_stress = STRESS

    call hill48_update(trial_stress, DDSDDE, DSTRAN, PROPS, phase_state, success)
    STATEV(9) = 0.0d0
    if (.not. success) then
        STATEV(9) = 1.0d0
        PNEWDT = min(PNEWDT, 0.5d0)
        return
    end if

    STRESS = trial_stress
    STATEV(1:2) = phase_state(7:8)
    STATEV(3:8) = phase_state(1:6)

contains

subroutine hill48_update(stress, ddsdde, dstran, props, statev, success)
implicit none
real(8), intent(inout) :: stress(6), statev(8)
real(8), intent(out) :: ddsdde(6,6)
real(8), intent(in) :: dstran(6), props(11)
logical, intent(out) :: success

!---- Local variables
real(8) :: E, NU, XMU, XLAM
real(8) :: F, G, H, XL, XM, XN, SY0, QINF, BVOCE
real(8) :: AMAT(6,6), QMAT(6,6), GAMMA_EIG0(6,6), GAMMA_EIG(6)
real(8) :: STRESS_TR(6), STRESS_DEV(6), STRESS_TILDE(6), STRESS_DEV_NEW(6)
real(8) :: ALPHA_N, SY_TR, DSYDA, FY_TR, DGAMMA, ALPHA_NEW, SY_NEW
real(8) :: KMAT(6), OMEGA, DOMEGA, DSY_DGAMMA, DFY, RESID, DENOM
real(8) :: HYDRO, TEMP1(6), TEMP2(6), W(6), MatD(6,6), DSTRESS(6), SQRT2
real(8), parameter :: TOL = 1.0D-8
integer :: I, J, K, ITER, INFO
integer, parameter :: MAXITER = 250
logical :: PLASTIC

!======================================================================
!     Read material properties
!======================================================================
E     = props(1)
NU    = props(2)
F     = props(6)
G     = props(7)
H     = props(8)
XL    = props(9)
XM    = props(10)
XN    = props(11)
SY0   = props(3)
QINF  = props(4)
BVOCE = props(5)

XMU  = E / (2.0D0 * (1.0D0 + NU))
XLAM = E * NU / ((1.0D0 + NU) * (1.0D0 - 2.0D0 * NU))
SQRT2 = DSQRT(2.0D0)

!     Elastic stiffness
MatD = 0.0D0
do I = 1, 3
    do J = 1, 3
        MatD(I,J) = XLAM
    end do
    MatD(I,I) = XLAM + 2.0D0*XMU
end do
MatD(4,4) = XMU
MatD(5,5) = XMU
MatD(6,6) = XMU

success = .false.
DDSDDE = MatD
ALPHA_N = statev(7)

!======================================================================
!     Build Hill48 anisotropy tensor (Mandel notation)
!======================================================================
AMAT = 0.0D0
AMAT(1,1) = G + H
AMAT(2,2) = F + H
AMAT(3,3) = F + G
AMAT(1,2) = -H;  AMAT(2,1) = -H
AMAT(1,3) = -G;  AMAT(3,1) = -G
AMAT(2,3) = -F;  AMAT(3,2) = -F
AMAT(4,4) = XN
AMAT(5,5) = XM
AMAT(6,6) = XL

!======================================================================
!     Eigendecomposition: A = Q * GAMMA * Q^T
!======================================================================
call compute_eigen_GRR(AMAT, W, QMAT, GAMMA_EIG0, INFO)
if (INFO /= 0) then
    write(*,*) "Hill48 eigensolve failed, INFO =", INFO
    return
end if
do I = 1, 6
    GAMMA_EIG(I) = GAMMA_EIG0(I,I)
end do

!======================================================================
!     Trial state
!======================================================================
call matvec6(MatD, dstran, DSTRESS)
STRESS_TR = stress + DSTRESS

HYDRO = (STRESS_TR(1) + STRESS_TR(2) + STRESS_TR(3)) / 3.0D0
STRESS_DEV(1:3) = STRESS_TR(1:3) - HYDRO
STRESS_DEV(4:6) = STRESS_TR(4:6) * SQRT2  ! Mandel

call matvec6_T(QMAT, STRESS_DEV, STRESS_TILDE)

call voce_hard(ALPHA_N, SY0, QINF, BVOCE, SY_TR, DSYDA)

FY_TR = -SY_TR**2
do I = 1, 6
    FY_TR = FY_TR + 0.5D0 * GAMMA_EIG(I) * STRESS_TILDE(I)**2
end do

!======================================================================
!     Yield check and plastic correction
!======================================================================
if (FY_TR <= TOL) then
    PLASTIC = .FALSE.
    stress = STRESS_TR
    DGAMMA = 0.0D0
    ALPHA_NEW = ALPHA_N
    SY_NEW = SY_TR
else
    PLASTIC = .TRUE.
    DGAMMA = 0.0D0
    ALPHA_NEW = ALPHA_N

    do ITER = 1, MAXITER
        do I = 1, 6
            if (GAMMA_EIG(I) > 1.0D-12) then
                KMAT(I) = GAMMA_EIG(I) / (1.0D0 + 2.0D0*XMU*DGAMMA*GAMMA_EIG(I))**2
            else
                KMAT(I) = 0.0D0
            end if
        end do

        OMEGA = 0.0D0
        do I = 1, 6
            OMEGA = OMEGA + 0.5D0 * KMAT(I) * STRESS_TILDE(I)**2
        end do
        OMEGA = DSQRT(MAX(OMEGA, 1.0D-30))

        ALPHA_NEW = ALPHA_N + 2.0D0 * DGAMMA * OMEGA
        call voce_hard(ALPHA_NEW, SY0, QINF, BVOCE, SY_NEW, DSYDA)

        RESID = SY_NEW / OMEGA - 1.0D0
        if (ABS(RESID) <= TOL) exit

        DOMEGA = 0.0D0
        do I = 1, 6
            if (GAMMA_EIG(I) > 1.0D-12) then
                TEMP1(I) = -4.0D0*XMU*GAMMA_EIG(I)*KMAT(I) / &
                    (1.0D0 + 2.0D0*XMU*DGAMMA*GAMMA_EIG(I))
                DOMEGA = DOMEGA + 0.25D0/OMEGA * TEMP1(I) * STRESS_TILDE(I)**2
            end if
        end do

        DSY_DGAMMA = 2.0D0 * DSYDA * (OMEGA + DGAMMA * DOMEGA)
        DFY = DSY_DGAMMA/OMEGA - SY_NEW*DOMEGA/OMEGA**2

        if (ABS(DFY) > 1.0D-30) then
            DGAMMA = DGAMMA - RESID / DFY
        end if
        DGAMMA = MAX(DGAMMA, 0.0D0)
    end do

    if (ITER > MAXITER) then
        write(*,*) 'Hill48 return mapping did not converge.'
        return
    end if

    !     Update stress
    do I = 1, 6
        if (GAMMA_EIG(I) > 1.0D-12) then
            TEMP1(I) = STRESS_TILDE(I) / (1.0D0 + 2.0D0*XMU*DGAMMA*GAMMA_EIG(I))
        else
            TEMP1(I) = STRESS_TILDE(I)
        end if
    end do

    call matvec6(QMAT, TEMP1, STRESS_DEV_NEW)
    STRESS_DEV_NEW(4:6) = STRESS_DEV_NEW(4:6) / SQRT2

    stress(1:3) = STRESS_DEV_NEW(1:3) + HYDRO
    stress(4:6) = STRESS_DEV_NEW(4:6)
end if

!======================================================================
!     Update state variables
!     statev(1:6) = plastic strain components  (engineering shear)
!     statev(7)   = alpha (equivalent plastic strain)
!     statev(8)   = SY (current yield stress)
!======================================================================
statev(7) = ALPHA_NEW
statev(8) = SY_NEW

if (PLASTIC) then
    TEMP2(1:3) = STRESS_DEV_NEW(1:3)
    TEMP2(4:6) = STRESS_DEV_NEW(4:6) * SQRT2

    AMAT = 0.0D0
    AMAT(1,1) = G + H;  AMAT(2,2) = F + H;  AMAT(3,3) = F + G
    AMAT(1,2) = -H;  AMAT(2,1) = -H
    AMAT(1,3) = -G;  AMAT(3,1) = -G
    AMAT(2,3) = -F;  AMAT(3,2) = -F
    AMAT(4,4) = XN;  AMAT(5,5) = XM;  AMAT(6,6) = XL

    call matvec6(AMAT, TEMP2, TEMP1)
    ! Mandel flow component 4 is eps_p4 = Dgamma*N*sqrt(2)*sig12;
    ! engineering shear gamma_p = sqrt(2)*eps_p4.
    TEMP1(4:6) = TEMP1(4:6) * SQRT2

    statev(1:6) = statev(1:6) + DGAMMA * TEMP1(1:6)
end if

!======================================================================
!     Tangent stiffness — consistent (algorithmic) tangent for GRR
!
!     C^alg = R · C^e_dev + β · n ⊗ q + (volumetric part of MatD)
!
!     Eigenspace linearization following Versino & Bennett (2018):
!       R = T⁻¹ · Q · D⁻¹ · Q^T · T    (Mandel/Voigt mapping)
!         D_inv(i) = 1/(1 + 2μ Δγ Γ_i)
!       n_i = (T⁻¹ Q D⁻² Γ s̃)(i)       convert Mandel→Voigt: shear÷√2
!       m_i = (T  Q K s̃/(2ω))(i)         convert Voigt←Mandel: shear×√2
!       q   = C^e_dev · m
!       β   = 2μ (DSYDA·2Δγ − 1) / (ω · DFY)
!
!     T scales shear slots by sqrt(2); its inverse scales by 1/sqrt(2).
!     T maps physical stress to Mandel and Mandel strain to engineering strain.
!
!     Scratch usage (arrays re-used after stress update is done):
!       TEMP1(6) = D_inv, TEMP2(6) = D_inv²
!       AMAT(6,6)= R (after eigendecomp; refilled in statev block then reused here)
!       GAMMA_EIG0(6,6) = intermediate scratch (C^e_dev, then Q·D_inv, ...)
!       W(6)     = n_vec (Voigt)
!       DSTRESS(6)= m_vec (Voigt)
!       STRESS_DEV(6) = q_vec
!     The stress / statev update above (lines up to 230) is UNCHANGED.
!======================================================================
if (PLASTIC) then

    !------------------------------------------------------------------
    !  Recompute converged-state KMAT, OMEGA, DOMEGA, DFY (same
    !  formulas as the N-R loop but evaluated once at the final DGAMMA).
    !  Reference: hill48_umat.f:293-320
    !------------------------------------------------------------------
    do I = 1, 6
        if (GAMMA_EIG(I) > 1.0D-12) then
            KMAT(I) = GAMMA_EIG(I) / &
                      (1.0D0 + 2.0D0*XMU*DGAMMA*GAMMA_EIG(I))**2
        else
            KMAT(I) = 0.0D0
        end if
    end do

    OMEGA = 0.0D0
    do I = 1, 6
        OMEGA = OMEGA + 0.5D0 * KMAT(I) * STRESS_TILDE(I)**2
    end do
    OMEGA = DSQRT(MAX(OMEGA, 1.0D-30))

    DOMEGA = 0.0D0
    do I = 1, 6
        if (GAMMA_EIG(I) > 1.0D-12) then
            TEMP1(I) = -4.0D0*XMU*GAMMA_EIG(I)*KMAT(I) / &
                       (1.0D0 + 2.0D0*XMU*DGAMMA*GAMMA_EIG(I))
            DOMEGA = DOMEGA + 0.25D0/OMEGA * TEMP1(I) * STRESS_TILDE(I)**2
        end if
    end do
    DSY_DGAMMA = 2.0D0 * DSYDA * (OMEGA + DGAMMA * DOMEGA)
    DFY = DSY_DGAMMA/OMEGA - SY_NEW*DOMEGA/OMEGA**2

    !------------------------------------------------------------------
    !  Step 1: D_inv(I) = 1/(1+2μΔγΓ_i),  D_inv2(I) = D_inv(I)^2
    !  Stored in TEMP1 and TEMP2 respectively.
    !------------------------------------------------------------------
    do I = 1, 6
        if (GAMMA_EIG(I) > 1.0D-12) then
            DENOM     = 1.0D0 + 2.0D0 * XMU * DGAMMA * GAMMA_EIG(I)
            TEMP1(I)  = 1.0D0 / DENOM
            TEMP2(I)  = TEMP1(I) * TEMP1(I)
        else
            TEMP1(I)  = 1.0D0
            TEMP2(I)  = 1.0D0
        end if
    end do

    !------------------------------------------------------------------
    !  Step 2: Build R = T⁻¹ · (Q D_inv Q^T) · T
    !  Use GAMMA_EIG0 as scratch (Mandel form first, then scaled).
    !  (a) GAMMA_EIG0 = Q · diag(D_inv)  (scale cols of Q)
    !  (b) AMAT = GAMMA_EIG0 · Q^T  (= Q D_inv Q^T in Mandel)
    !  (c) Apply T⁻¹ (÷√2 on rows 4:6) and T (×√2 on cols 4:6)
    !------------------------------------------------------------------
    do I = 1, 6
        do J = 1, 6
            GAMMA_EIG0(I,J) = QMAT(I,J) * TEMP1(J)
        end do
    end do
    AMAT = 0.0D0
    do I = 1, 6
        do J = 1, 6
            do K = 1, 6
                AMAT(I,J) = AMAT(I,J) + GAMMA_EIG0(I,K) * QMAT(J,K)
            end do
        end do
    end do
    !  T on cols 4:6 (×√2), T⁻¹ on rows 4:6 (÷√2) → net: T⁻¹·AMAT·T
    do J = 4, 6
        do I = 1, 6
            AMAT(I,J) = AMAT(I,J) * SQRT2
        end do
    end do
    do I = 4, 6
        do J = 1, 6
            AMAT(I,J) = AMAT(I,J) / SQRT2
        end do
    end do
    !  AMAT now holds R in Voigt

    !------------------------------------------------------------------
    !  Step 3: n = T⁻¹ · Q · (D_inv² · Γ · s̃_tr)
    !  temp = D_inv2(I)*GAMMA_EIG(I)*STRESS_TILDE(I) in eigen-space
    !  then back-transform: Q·temp (Mandel) then ÷√2 on shear
    !  Result stored in W(6).
    !------------------------------------------------------------------
    do I = 1, 6
        GAMMA_EIG0(I,1) = TEMP2(I) * GAMMA_EIG(I) * STRESS_TILDE(I)
    end do
    call matvec6(QMAT, GAMMA_EIG0(1:6,1), W)
    W(4:6) = W(4:6) / SQRT2   ! Mandel→Voigt for shear

    !------------------------------------------------------------------
    !  Step 4: m = (T · Q · K · s̃_tr) / (2ω)
    !  temp = KMAT(I)*STRESS_TILDE(I) in eigen-space
    !  then Q·temp (Mandel), ×√2 on shear, ÷(2ω)
    !  Result stored in DSTRESS(6).
    !------------------------------------------------------------------
    do I = 1, 6
        GAMMA_EIG0(I,1) = KMAT(I) * STRESS_TILDE(I)
    end do
    call matvec6(QMAT, GAMMA_EIG0(1:6,1), DSTRESS)
    DSTRESS(1:3) = DSTRESS(1:3)          / (2.0D0 * OMEGA)
    DSTRESS(4:6) = DSTRESS(4:6) * SQRT2  / (2.0D0 * OMEGA)

    !------------------------------------------------------------------
    !  Step 5: C^e_dev = MatD − vol projection
    !  Use GAMMA_EIG0 as C^e_dev (overwrite all 6 cols).
    !------------------------------------------------------------------
    do I = 1, 6
        do J = 1, 6
            GAMMA_EIG0(I,J) = MatD(I,J)
        end do
    end do
    do J = 1, 6
        HYDRO = (MatD(1,J) + MatD(2,J) + MatD(3,J)) / 3.0D0
        do I = 1, 3
            GAMMA_EIG0(I,J) = GAMMA_EIG0(I,J) - HYDRO
        end do
    end do

    !------------------------------------------------------------------
    !  Step 6: q = C^e_dev · m  (store in STRESS_DEV)
    !------------------------------------------------------------------
    call matvec6(GAMMA_EIG0, DSTRESS, STRESS_DEV)

    !------------------------------------------------------------------
    !  Step 7: β = 2μ (DSYDA·2Δγ − 1) / (ω · DFY)
    !  Stored in DENOM (scalar).
    !------------------------------------------------------------------
    if (ABS(DFY) > 1.0D-30) then
        DENOM = 2.0D0 * XMU * (DSYDA * 2.0D0 * DGAMMA - 1.0D0) / &
                (OMEGA * DFY)
    else
        DENOM = 0.0D0
    end if

    !------------------------------------------------------------------
    !  Step 8: C^alg = R·C^e_dev + β·n⊗q + volumetric
    !  (a) ddsdde = R · C^e_dev   (R=AMAT, C^e_dev=GAMMA_EIG0)
    !  (b) ddsdde += β · n(=W) ⊗ q(=STRESS_DEV)
    !  (c) ddsdde += volumetric col contributions
    !------------------------------------------------------------------
    ddsdde = 0.0D0
    do I = 1, 6
        do J = 1, 6
            do K = 1, 6
                ddsdde(I,J) = ddsdde(I,J) + AMAT(I,K) * GAMMA_EIG0(K,J)
            end do
        end do
    end do

    do I = 1, 6
        do J = 1, 6
            ddsdde(I,J) = ddsdde(I,J) + DENOM * W(I) * STRESS_DEV(J)
        end do
    end do

    do J = 1, 6
        HYDRO = (MatD(1,J) + MatD(2,J) + MatD(3,J)) / 3.0D0
        do I = 1, 3
            ddsdde(I,J) = ddsdde(I,J) + HYDRO
        end do
    end do

else
    ddsdde = MatD
end if

success = .true.
end subroutine hill48_update

subroutine voce_hard(alpha, sy0, qinf, b, sy, dsyda)
implicit none
real(8), intent(in) :: alpha, sy0, qinf, b
real(8), intent(out) :: sy, dsyda
real(8) :: ex
ex = exp(-b * alpha)
sy = sy0 + qinf * (1.0D0 - ex)
dsyda = qinf * b * ex
end subroutine voce_hard

!======================================================================
subroutine compute_eigen_GRR(A, W, Q, Gamma, INFO)
implicit none
real(8), intent(inout) :: A(6,6)
real(8), intent(out) :: W(6), Q(6,6), Gamma(6,6)
integer, intent(out) :: INFO
real(8), allocatable :: WORK(:)
integer :: LWORK, i

LWORK = -1
allocate(WORK(1))
call DSYEV('V', 'U', 6, A, 6, W, WORK, LWORK, INFO)
if (INFO /= 0) return
LWORK = INT(WORK(1))
deallocate(WORK)
allocate(WORK(LWORK))
call DSYEV('V', 'U', 6, A, 6, W, WORK, LWORK, INFO)

if (INFO /= 0) return
Q = A
Gamma = 0.0D0
do i = 1, 6
    Gamma(i,i) = W(i)
end do
deallocate(WORK)
end subroutine compute_eigen_GRR

!======================================================================
subroutine matvec6(A, X, Y)
implicit none
real(8), intent(in) :: A(6,6), X(6)
real(8), intent(out) :: Y(6)
integer :: i
do i = 1, 6
    Y(i) = sum(A(i,:) * X(:))
end do
end subroutine matvec6

!======================================================================
subroutine matvec6_T(A, X, Y)
implicit none
real(8), intent(in) :: A(6,6), X(6)
real(8), intent(out) :: Y(6)
integer :: i
do i = 1, 6
    Y(i) = sum(A(:,i) * X(:))
end do
end subroutine matvec6_T

end subroutine UMAT
