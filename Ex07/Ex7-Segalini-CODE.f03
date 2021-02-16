module quantum_util
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! quantum_util module:
!       contains utility functions for setting a "quantum environment"
! ---------------------------------------------------------------------
! Grid1dim type:
!    a type defined to include all the parameters of a discretized
!    space.
!   INCLUDES:
!       minimum = double precision, lower bound of the interval
!       maximum = double precision, upper bound of the interval
!       N       = integer, size of the lattice (number of points)
!       dx      = double precision, step between two points
!       lattice = double precision, dimension(:), allocatable,
!                 contains the grid points
! ---------------------------------------------------------------------
! GenLattice:
!    function to create a lattice environment, i.e. to
!    discretize an interval. Given the maximum and minimum value of x,
!    @up and @down, the function produces a grid of @N intervals between
!    the two. It returns a type Grid1dim, whose %lattice is a
!    @N+1 dimensional vector including the boundaries of the
!    aforementioned intervals, %dx is the step.
!   INPUT:
!       down  = double precision, lower bound of the lattice
!       up    = double precision, upper bound of the lattice
!       N_int = integer, number of intervasl to divide the interval
!               [down, up] in
!   OUTPUT:
!       grid  = type Grid1dim, structure which includes all the
!               info about the discretized interval
! ---------------------------------------------------------------------
! Potential:
!   function to calculate the potential as a function of the
!   positions @q and of the time @t. The potential is computed
!   accordingly to the problem's requests.
!   INPUT:
!       q  = double precision, position
!       t  = double precision, time
!       TT = double precision, parameter that determines potential
!            shape
!   OUTPUT:
!       V  = double precision, potential computed in point (@q, @t)
! ---------------------------------------------------------------------
! norm:
!   function that computes the discretised norm of a wavefunction
!   @psi, given the width of the discretization @step. normpsi
!   INPUT:
!       psi     = double complex, dimension(:), allocatable,
!                 wavefunction
!       step    = double precision, discretization width
!   OUTPUT:
!       normpsi = double precision, norm of @psi
! ---------------------------------------------------------------------
! ExpectedValue
    !function that computes the expected value @E of some observable
!   @obs wrt wavefunction @psi.
!   INPUT:
!       psi  = double complex, dimension(:), allocatable,
!              wavefunction
!       obs  = double precision, dimension(:), allocatable,
!              observable to which we want to compute the expected
!              value
!       step = double precision, discretization width
!   OUTPUT:
!       E    = double precision, expected value
! ---------------------------------------------------------------------
! printmatrix:
!    subroutine to print matrices in a clearer form
!   INPUT:
!       A = matrix to be printed
!   OUTPUT:
!       print on screen matrix A arranged in a "grid" form (elements
!       arranged in rows and columns)
! ---------------------------------------------------------------------
! str:
!    function that convert the integer @k to a string, in order to
!    properly print/concatenate it to other strings (especially for
!    handling output file names)
!   INPUT:
!       k = integer, to be converted to string
!   OUTPUT:
!       return @k converted into a string
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

    implicit none
    integer aa

    type Grid1dim

        double precision :: minimum
        double precision :: maximum
        integer          :: N
        double precision :: dx

        double precision, dimension(:), allocatable :: lattice

    end type Grid1dim

    contains

    function GenLattice(down, up, N_int)result(grid)

        double precision, intent(IN) :: down, up
        integer, intent(IN)          :: N_int
        type(Grid1dim)               :: grid

        grid%minimum = down
        grid%maximum = up
        grid%N       = N_int+1
        grid%dx      = (up-down)/N_int

        allocate(grid%lattice(grid%N))

        do aa = 1, grid%N
            grid%lattice(aa) = down + (aa-1)*grid%dx
        end do

    end function

    function Potential(q, t, TT)result(V)

        implicit none

        double precision, intent(IN) :: TT
        double precision, intent(IN) :: q
        double precision, intent(IN) :: t
        double precision             :: V

        if (TT == 0d0) then
            V = 0.5d0*(q)**2d0
        else
            V = 0.5d0*(q - t/TT)**2d0
        end if

    end function Potential

    function norm(psi, step)result(normpsi)

        implicit none
        double complex, dimension(:), allocatable :: psi
        double precision                          :: step, normpsi

        normpsi = SUM(ABS(psi)**2 * step)

    end function norm

    function ExpectedValue(psi, obs, step)result(E)

        implicit none

        double complex, dimension(:), allocatable   :: psi
        double precision, dimension(:), allocatable :: obs
        double precision                            :: step, E

        E = 0

        do aa = 1, size(psi)
            E = E + abs(psi(aa))**2 *obs(aa)* step
        end do

    end function ExpectedValue

    subroutine printmatrix(A)

        implicit none
        double precision, dimension(:,:) :: A

        ! cycle to print matrix in the requested form
        do aa = 1, ubound(A, 1)
            print *, A(aa, :)
        end do
    end subroutine

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

end module

program timedepSE

    use quantum_util
    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

    double complex, dimension(:), allocatable   :: psi_now, psi_next, psi_tmp

    type(Grid1dim)                              :: xgrid, tgrid, pgrid

    double complex, dimension(:), allocatable   :: vect_U_V, vect_U_T, A_q, A_p, B_q, B_p
    double precision, dimension(:), allocatable :: q2, p2, ps

    double precision, dimension(3)              :: TT

    double precision                            :: L, hbar, m, omega, pi
    double precision                            :: ex, ex2, ep, ep2

    integer                                     :: ii, jj, n_x, n_t, pp
    integer*8                                   :: dfft_plan, idfft_plan

    ! define main constants
    pi    = 4.d0 * datan(1.d0)
    m     = 1.d0
    omega = 1.d0
    hbar  = 1.d0

    ! grid parameters
    L   = 50.d0
    n_x = 10000
    n_t = 1000

    ! different values of TT
    TT  = (/ 1.d0, 20.d0, 100.d0 /)

    do pp = 1, 3
        ! create lattice
        xgrid = GenLattice(-L, L, n_x)
        tgrid = GenLattice(0.d0, TT(pp), n_t)
        pgrid = GenLattice(-pi/xgrid%dx, pi/xgrid%dx, n_x)

        allocate(ps(pgrid%N))

        ! reordering necessary due to how the FFTW works
        ! positive frequencies in the first half, then negative frequencies in the second half
        do ii = 1, int(pgrid%N/2)
            ps(ii) = pgrid%dx * dble(ii-1)
        end do

        do ii = int(pgrid%N/2) + 1, pgrid%N
            ps(ii) = -1d0 * dble(pgrid%N - ii + 1) * pgrid%dx
        end do

        ! now ps is [0, pgrid%dx, ..., (pgrid%N/2-1)*pgrid%dx, -(pgrid%N/2+1)*pgrid%dx, ..., -pgrid%dx]

        allocate( psi_now(xgrid%N), psi_next(xgrid%N) )
        allocate( vect_U_V(xgrid%N), vect_U_T(xgrid%N), psi_tmp(xgrid%N) )
        allocate( A_q(xgrid%N), A_p(xgrid%N), B_q(xgrid%N), B_p(xgrid%N), q2(xgrid%N), p2(pgrid%N) )

        ! set filenames
        open(11, file = "T"//trim(str(int(TT(pp))))//"E_x.dat")
        open(12, file = "T"//trim(str(int(TT(pp))))//"sigma_x.dat")
        open(13, file = "T"//trim(str(int(TT(pp))))//"E_p.dat")
        open(14, file = "T"//trim(str(int(TT(pp))))//"sigma_p.dat")

        ! DFFTW plans generation
        call dfftw_plan_dft_1d(dfft_plan, xgrid%N, psi_now, psi_tmp, FFTW_FORWARD, FFTW_MEASURE)
        call dfftw_plan_dft_1d(idfft_plan, xgrid%N, psi_now, psi_tmp, FFTW_BACKWARD, FFTW_MEASURE)

        ! set initial state
        do ii = 1, xgrid%N
            psi_now(ii) = pi**(-0.25d0) * EXP(-(xgrid%lattice(ii)**2d0)/2d0)
        end do
        psi_now = psi_now/norm(psi_now, xgrid%dx)

        !perform step-by-step time simulation, with split operator method
        do ii = 1, tgrid%N
            do jj = 1, xgrid%N
                ! kinetic part
                vect_U_T(jj) = EXP(dcmplx(0d0,-0.5d0) * tgrid%dx * ps(jj)**2d0 )
                ! potential part
                vect_U_V(jj) = EXP(dcmplx(0d0,-0.5d0) * tgrid%dx * Potential(xgrid%lattice(jj), tgrid%lattice(ii), TT(pp)) )
            end do

            do jj = 1, xgrid%N
                A_q(jj) = vect_U_V(jj) * psi_now(jj)
            end do

            call dfftw_execute_dft(dfft_plan, A_q, A_p)

            do jj = 1, xgrid%N
                B_p(jj) = vect_U_T(jj) * A_p(jj)/sqrt( dble(xgrid%N) )
            end do

            call dfftw_execute_dft(idfft_plan, B_p, B_q)

            do jj = 1, xgrid%N
                psi_next(jj) = vect_U_V(jj) * B_q(jj)/sqrt( dble(xgrid%N) )
            end do

            psi_next = psi_next/norm(psi_next, xgrid%dx)

            call dfftw_execute_dft(dfft_plan, psi_next, psi_tmp)

            psi_tmp = psi_tmp/norm(psi_tmp, pgrid%dx)

            q2 = xgrid%lattice**2
            p2 = ps**2

            ex  = ExpectedValue(psi_next, xgrid%lattice, xgrid%dx)
            ex2 = ExpectedValue(psi_next, q2, xgrid%dx)
            ep  = ExpectedValue(psi_tmp, ps, pgrid%dx)
            ep2 = ExpectedValue(psi_tmp, p2, pgrid%dx)

            ! write outputs
            ! position expected value
            write(11, *) tgrid%lattice(ii), ex
            ! variance of position
            write(12, *) tgrid%lattice(ii), sqrt(ex2 - ex**2)
            ! momentum expected value
            write(13, *) tgrid%lattice(ii), ep
            ! variance of momentum
            write(14, *) tgrid%lattice(ii), sqrt(ep2 - ep**2)

            psi_now = psi_next

        end do

        deallocate( psi_now, psi_next, psi_tmp, ps)
        deallocate( vect_U_V, vect_U_T, q2, p2)
        deallocate( A_q, A_p, B_q, B_p )
        deallocate( xgrid%lattice, pgrid%lattice, tgrid%lattice )

    end do

end program
