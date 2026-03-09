!======================================================================
!     HillPlasticity: GRR-based Hill48 with Voce Hardening
!     Adapted interface for DNS_Umat.f and ROM_Umat.f
!
!     Material Properties (props array):
!     props(1)  = E      - Young's modulus
!     props(2)  = NU     - Poisson's ratio
!     props(3)  = F      - Hill48 parameter
!     props(4)  = G      - Hill48 parameter
!     props(5)  = H      - Hill48 parameter
!     props(6)  = L      - Hill48 parameter
!     props(7)  = M      - Hill48 parameter
!     props(8)  = N      - Hill48 parameter
!     props(9)  = SY0    - Initial yield stress (Voce)
!     props(10) = QINF   - Saturation stress increment (Voce)
!     props(11) = BVOCE  - Hardening rate parameter (Voce)
!
!     State Variables (statev array):
!     statev(1:6) = plastic strain components
!     statev(7)   = alpha (equivalent plastic strain)
!     statev(8)   = SY (current yield stress)
!
!======================================================================

SUBROUTINE HillPlasticity(stress, stran, ddsdde, dstran, &
                          time, dtime, props, nprops, statev, nstatv)
!
USE NumType
USE lapack95
implicit none
!
integer(kind=ikind), intent(in) :: nstatv, nprops
real(kind=rkind), intent(in) :: dtime
real(kind=rkind), intent(inout) :: stress(6), statev(nstatv)
real(kind=rkind), intent(out) :: ddsdde(6,6)
real(kind=rkind), intent(in) :: stran(6), dstran(6), time(2), props(nprops)
!
!---- Local variables
real(8) :: E, NU, XMU, XLAM
real(8) :: F, G, H, XL, XM, XN, SY0, QINF, BVOCE
real(8) :: AMAT(6,6), QMAT(6,6), GAMMA_EIG0(6,6), GAMMA_EIG(6)
real(8) :: STRESS_TR(6), STRESS_DEV(6), STRESS_TILDE(6), STRESS_DEV_NEW(6)
real(8) :: ALPHA_N, SY_TR, DSYDA, FY_TR, DGAMMA, ALPHA_NEW, SY_NEW
real(8) :: KMAT(6), OMEGA, DOMEGA, DSY_DGAMMA, DFY, RESID, DENOM
real(8) :: HYDRO, TEMP1(6), TEMP2(6), W(6), MatD(6,6), DSTRESS(6), SQRT2
real(8), parameter :: TOL = 1.0D-8
integer :: I, J, ITER, INFO
integer, parameter :: MAXITER = 250
logical :: PLASTIC

!======================================================================
!     Read material properties
!======================================================================
E     = props(1)
NU    = props(2)
F     = props(3)
G     = props(4)
H     = props(5)
XL    = props(6)
XM    = props(7)
XN    = props(8)
SY0   = props(9)
QINF  = props(10)
BVOCE = props(11)

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

    if (ITER >= MAXITER) then
        write(*,*) '***** WARNING: HillPlasticity_GRR Newton iteration did not converge *****'
        write(*,*) '      Max iterations (', MAXITER, ') reached'
        write(*,*) '      Residual: ', ABS(RESID)
        write(*,*) '      Tolerance: ', TOL
        write(*,*) '      DGAMMA: ', DGAMMA
        write(*,*) '      ALPHA_NEW: ', ALPHA_NEW
        write(*,*) '      Time: ', time(2)
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
!     statev(1:6) = plastic strain components
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
    TEMP1(4:6) = TEMP1(4:6) / SQRT2

    statev(1:6) = statev(1:6) + DGAMMA * TEMP1(1:6)
end if

!======================================================================
!     Tangent stiffness (physics-based softening)
!======================================================================
if (PLASTIC) then
    DENOM = 1.0D0 + 2.0D0 * XMU * DGAMMA * MAXVAL(GAMMA_EIG)
    ddsdde = MatD / DENOM
else
    ddsdde = MatD
end if

RETURN
END SUBROUTINE HillPlasticity

!======================================================================
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
LWORK = INT(WORK(1))
deallocate(WORK)
allocate(WORK(LWORK))
call DSYEV('V', 'U', 6, A, 6, W, WORK, LWORK, INFO)

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
