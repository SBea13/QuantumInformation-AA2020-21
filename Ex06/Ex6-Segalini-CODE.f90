module matrix_utilities
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! matrix_utilities module:
!       contains utility functions for matrix handling
! ---------------------------------------------------------------------
! HermRandMat:
!    function to initialize a N-by-N hermitian matrix (a) with random
!    complex numbers with real and imaginary part drawn uniformly
!    in (-1,1)
!   INPUT:
!       N = integer, matrix size
!   OUTPUT:
!       a = complex*16, dimension(N,N), hermitian matrix
! ---------------------------------------------------------------------
! DiagRandMat:
!    function to initialize a N-by-N diagonal matrix (a) with random
!    real numbers drawn uniformly in (-1,1)
!   INPUT:
!       N = integer, matrix size
!   OUTPUT:
!       a = double precision, dimension(N,N), diagonal matrix
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
    integer :: aa, bb ! indexes for loops

    contains

    function HermRandMat(N)result(h)

        integer                     :: N
        complex*16, dimension(N, N) :: h

        ! cycle to fill in hermitian matrix a
        do aa = 1,N
            do bb = 1, aa
                h(aa, bb) = cmplx(RAND(0)*2 - 1, RAND(0)*2 - 1)
                if (aa /= bb) then
                    h(bb, aa) = conjg(h(aa, bb))
                end if
            end do
        end do

    end function HermRandMat

    function DiagRandMat(N, evl)result(d)

    integer                              :: N, INFO
    double precision, dimension(N), intent(OUT)    :: evl
    double precision, dimension(N, N)              :: d

    ! set all elements to 0
    d = 0

    ! cycle to fill diagonal matrix a
    do aa = 1,N
        d(aa, aa) = RAND(0)*2 - 1
        evl(aa) = d(aa, aa)
    end do

    call dlasrt('I', N, evl, INFO)

    end function DiagRandMat

    subroutine printmatrix(A)

        implicit none
        real*4, dimension(:,:) :: A

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

end module matrix_utilities

module quantum_util
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! quantum_util module:
!       contains utility functions for setting a "quantum environment"
! ---------------------------------------------------------------------
! GenLattice:
!    function to create a lattice environment, i.e. to
!    discretize the space. Given the maximum and minimum value of x,
!    the function produces a grid of @N intervals between @xmax and
!    @xmin. It returns a N+1 dimensional vector containing the
!    boundaries of the aforementioned intervals.
!   INPUT:
!       xmax = double precision, right bound of the lattice (maximum value)
!       N    = integer, number of interval to divide the interval
!           [xmin, xmax] in
!   OUTPUT:
!       lat  = double precision, dimension(N+1), lattice
! ---------------------------------------------------------------------
! H_d:
!    subroutine to compute the discretized hamiltonian.
!    It takes as input an array of discretized positions @lattice and
!    a the value @omega for defining the potential shape.
!    It generates the elements of a tridiagonal matrix with diagonal
!    elements equal to 2d0 + (omega*dx*lattice(aa))**2
!    with @dx a support variable representing the spacing of the
!    spacial grid. The subdiagonal elements are set equal to
!    -1. The diagonal elements are stored in @D, the subdiagonal in @E
!    The hamiltonian is not scaled with the grid spacing @dx**2 for
!    avoiding divisions by potentially very small numbers.
!   INPUT:
!       lattice = double precision, dimension(:), discretized space
!                 grid
!       omega   = double precision, value of the potential constant
!       D       = double precision, dimension(@N), diagonal elements
!                 of the hamiltonian
!       E       = double precision, dimension(@N), subdiagonal elements
!                 of the hamiltonian (E(N) is just a support element
!                 needed for further computations)
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

    implicit none
    integer :: aa ! index for loops

    contains

    function GenLattice(xmax, xmin, N) result(grid)

        integer                          :: N
        double precision, intent(IN)     :: xmax, xmin
        double precision                 :: dx
        double precision, dimension(N+1) :: grid

        dx = (xmax - xmin)/N

        do aa = 1, N+1
            grid(aa) = xmin+(aa-1)*dx
        end do

    end function

    subroutine H_d(lattice, omega, D, E)

        double precision, dimension(:), intent(IN)  :: lattice
        double precision, intent(IN)                :: omega
        integer                                     :: N
        double precision                            :: dx
        double precision, dimension(:), intent(OUT) :: D, E

        N = size(lattice)

        dx = lattice(2) - lattice(1)

        do aa= 1, N
            D(aa) = 2d0/(dx**2) + (omega*lattice(aa))**2
            E(aa) = -1d0/(dx**2)
        end do

    end subroutine

end module quantum_util


module debugging_mod
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! debugging_mod module:
!       contains subroutine for checkpoints/debugging
! ---------------------------------------------------------------------
! debugging:
!    subroutine to enable debug mode. Checks if the @condition is true,
!    and produce an output on screen according to the other optional
!    inputs.
!   MANDATORY INPUT:
!       condition = logical, condition to be verified for debugging.
!                   If True, no bug is detected.
!   OPTIONAL INPUTS:
!       msg       = string, it is a text label for expressing the
!                   condition. It can be printed on screen.
!       verbose   = logical, express condition on "long" output,
!                   i.e. if True the full output is printed on screen
!       stopprg   = logical, if True the program is stopped at the
!                   first error detected
!       content   = general variable on which the condition is verified,
!                   can be printed as output if present and if we are
!                   selecting the "full" output
!   OUTPUT:
!       print on screen messages/checkpoints for debugging
!   SUPPORT VARIABLES:
!       full_text = logical, if True "full" output
!                   is displayed
!       stp       = logical, if True program is
!                   stopped whenever an error is
!                   found
!       msg_yes   = logical, if True the argument
!                   @msg is provided, hence it can
!                   be printed
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

    implicit none

    contains

    subroutine debugging(condition, msg, verbose, stopprg, content)

        logical, intent(IN)                 :: condition
        character(*), intent(IN), optional  :: msg
        logical, intent(IN), optional       :: verbose, stopprg
        class(*), intent(IN), optional      :: content

        logical                             :: full_text
        logical                             :: stp
        logical                             :: msg_yes

!       support variables
        full_text = (present(verbose).AND.verbose).OR.(.NOT.present(verbose))
        stp = present(stopprg).AND.stopprg
        msg_yes = present(msg)

        if (msg_yes .AND. full_text) then
!       "full" output mode
            if (condition .EQV. .TRUE.) then
!               check contion
                if (present(content)) then
!                   if variable is present, print it on screen
                    select type(content)
                        type is (integer(1))
                            print*, msg, " => [OK], Variable = ", content
                        type is (integer(2))
                            print*, msg, " => [OK], Variable = ", content
                        type is (integer(4))
                            print*, msg, " => [OK], Variable = ", content
                        type is (integer(8))
                            print*, msg, " => [OK], Variable = ", content
                        type is (real(4))
                            print*, msg, " => [OK], Variable = ", content
                        type is (real(8))
                            print*, msg, " => [OK], Variable = ", content
                        type is (logical)
                            print*, msg, " => [OK], Variable = ", content
                    end select
                else
                    print*, msg, " => [OK]"
                end if
            elseif (stp) then
!               exit program if stop condition is True
                print*, msg, " => [ERROR], Abort execution"
                stop
            else
!               just print ERROR string if an error is detected but stp is False
                print*, msg, " => [ERROR]"
            end if
        elseif (msg_yes) then
!           less verbose output, just with msg print
            if (condition .EQV. .TRUE.) then
                print*, msg, " => [OK]"
            elseif (stp) then
               print*, msg, " => [ERROR], Abort execution"
               stop
            else
                print*, msg, " => [ERROR]"
            end if
       else
!           short output
            if (condition .EQV. .TRUE.) then
                print*, "[OK]"
            elseif (stp) then
                stop
            else
                print*, "[ERROR]"
            end if

        end if

    end subroutine

end module

program QHO

    use quantum_util
    use matrix_utilities

    implicit none

    integer                                       :: ii, jj, nn, kk
    integer, dimension(3)                         :: N
    double precision, dimension(3)                :: dx
    double precision, dimension(:), allocatable   :: lattice, D, E, evl, supp
    double precision                              :: xmax, xmin, omega
    character(30)                                 :: outfile

    ! variables for dstemr LAPACK subroutine
    integer, PARAMETER                            :: LWMAX = 100000
	integer  	                                  :: IL, LWORK, LIWORK, INFO, points
	double precision                              :: WORK(LWMAX), IWORK(LWMAX)
    double precision  	                          :: norm
	double precision, dimension(:,:), allocatable :: evct
	integer, dimension(:), allocatable            :: ISUPPZ
	logical  	                                  :: TRYRAC

    ! define system parameters
    N = (/ 1000, 2000, 5000 /)
    xmax = 25.
    xmin = -25.

    dx = (xmax-xmin)/N

    omega = 1d0

    IL = 1
    kk = 1000
    TRYRAC = .TRUE.
    WORK = 1
    IWORK = 10000000

    ! compute eigenvalues and eigenvectors for different spacings dx
    do ii = 1, size(dx)

        points = N(ii)+1

        allocate(lattice(points))
        allocate(D(points))
        allocate(E(points))
        allocate(evl(points))
        allocate(ISUPPZ(2*kk))
        allocate(evct(points, kk))
        allocate(supp(points))

        ! Initialize support vector ISUPPZ
        do jj = 1, kk
            ISUPPZ(2 * jj - 1) = 1
            ISUPPZ(2 * jj) = points
        end do

        print*, "spacing dx=", dx(ii)
    !   compute lattice
        lattice = GenLattice(xmax, xmin, N(ii))

    !   compute hamiltonian
        print*, "Computing Hamiltonian..."
        call H_d(lattice, omega, D, E)

    !   compute optimal size of workspace
        print*, "Finding optimal workspace..."
        call DSTEMR('V', 'I', points, D, E, 0d0, 0d0, IL, kk, kk, &
                    evl, evct, points, kk, ISUPPZ, TRYRAC, &
                    WORK, -1, IWORK, -1, INFO)
!       print*, "LWORK:", WORK(1)
!       print*, "LIWORK:", IWORK(1)
        LIWORK = int(IWORK(1))
        LWORK = int(WORK(1))

    !   compute eigenvalues and eigenvactors
        print*, "Computing Eigevalues and Eigenvectors..."
        call DSTEMR('V', 'I', points, D, E, 0d0, 0d0, IL, kk, kk, &
                    evl, evct, points, kk, ISUPPZ, TRYRAC, &
                    WORK, LWORK, IWORK, LIWORK, INFO)

    !   save eigenvalues on file
        outfile = trim(str(N(ii)))//"points_evl"//".dat"

        open(42, file=outfile)
        do nn = 1, kk
            write(42, *) nn-1, evl(nn), abs(evl(nn)-2*nn+1)/(2*nn-1)
        end do
        close(42)

    !   normalize the eigenvectors
        do jj = 1, kk
            supp = evct(:, jj) ** 2

            !Trapezoidal rule
            norm = (dx(ii) / 2d0) * (supp(1) + 2d0 * sum(supp(2:points-1)) + supp(points))

            evct(:, jj) = evct(:, jj) / sqrt(norm)
        end do

    !   save eigenvectors on file
        outfile = trim(str(N(ii)))//"points_evc"//".dat"

        open(13, file=outfile)
        do jj = 1, points
            write(13, '(F14.7, a)', advance='no') lattice(jj), achar(9)
            write(13, *) evct(jj, :)
        end do
        close(13)


        deallocate(lattice)
        deallocate(D)
        deallocate(E)
        deallocate(evl)
        deallocate(ISUPPZ)
        deallocate(evct)
        deallocate(supp)

    end do

    open(22, file="true_evl.dat")
    do nn = 1, kk
        write(22, *) nn-1, 2*nn-1
    end do
    close(22)
    print*, "Execution complete."

end program
