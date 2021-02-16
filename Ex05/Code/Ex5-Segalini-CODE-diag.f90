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
!       a = real*8, dimension(N,N), diagonal matrix
! ---------------------------------------------------------------------
! printmatrix:
!    subroutine to print matrices in a clearer form
!   INPUT:
!       A = matrix to be printed
!   OUTPUT:
!       print on screen matrix A arranged in a "grid" form (elements
!       arranged in rows and columns)
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
    real*8, dimension(N), intent(OUT)    :: evl
    real*8, dimension(N, N)              :: d

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

    function spacings(evl, step)result(ss)

        implicit none
        double precision, dimension(:), intent(in)  :: evl
        integer, intent(in)                         :: step
        integer                                     :: N, left, right, INFO
        double precision                            :: avg
        double precision, dimension(:), allocatable :: ss, devl

        N = size(evl)

        allocate(ss(N-1))
        allocate(devl(N-1))

        devl = evl(2:N) - evl(1:(N-1))
        !s = devl/sum(devl)*(N-1)

        do aa = 1, (N-1)
            left = max(1, aa - step)
            right = min(N-1, aa + step)
            avg = sum(devl(left:right))/(right-left+1)
            ss(aa)= devl(aa)/avg
        end do

        call dlasrt('I', N-1, ss, INFO)

!        do aa = 1, (N-1), levels
!            avg = sum(devl(aa:(aa+levels-1))/levels
!            do bb = aa, (levels+aa)
!                s = devl(bb)/avg
!            end do
!        end do
        deallocate(devl)

    end function

    character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
    end function str

end module matrix_utilities

module histograms

    implicit none

    contains

    function make_hist(values, nbins, bottom, top)result(counts)

        implicit none
        real*8, dimension(:), intent(IN)   :: values
        integer, intent(IN)                :: nbins
        real*4, intent(IN)                 :: bottom, top
        integer                            :: N, jj, ii

        real*8                             :: temp, limit, dx
        integer, dimension(:), allocatable :: counts
        logical                            :: flag

        flag = .True.
        N = size(values)

        allocate(counts(nbins))

        dx = (top - bottom)/nbins

        counts = 0
        limit = bottom + dx
        jj = 1

        do ii = 1, nbins
            flag = .True.
            do while((jj <= N) .AND. flag)
                temp = values(jj)
                if (temp > limit) then
                    limit = limit + dx
                    flag=.FALSE.
                else
                    counts(ii) = counts(ii) + 1
                    jj = jj + 1
                end if
            end do
        end do

    end function


end module

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

program eigenvalues

    use matrix_utilities
    use histograms

    implicit none

    integer                                     :: N, nbin
    real*8, dimension(:,:), allocatable         :: D
    integer, dimension(:), allocatable          :: levels
    double precision, dimension(:), allocatable :: evl, s
    integer                                     :: kk, it, zz
    real*4                                      :: dx, xmax, xmin, norm
    character(30)                               :: outfile
    character(4)                                :: name
    integer, dimension(:), allocatable          :: counts, total_hist

    N = 1500
    nbin = 100
    xmax = 5.
    xmin = 0.

    dx = (xmax-xmin)/nbin

    allocate(D(N,N))
    allocate(levels(2))
    allocate(evl(N))
    allocate(total_hist(nbin))

    levels = (/ 1500, 10 /)

    do kk = 1, size(levels)

        total_hist = 0
        print*, "avg on:", levels(kk), "values."

        do it = 1, 20
            print*, "iteration:", it
        !   generate diagonal matrix and compute its sorted eigevalues
            D = DiagRandMat(N, evl)

        !   compute spacings with local/global averages
            s = spacings(evl, levels(kk))
        !   create histogram for single iteration
            counts = make_hist(s, nbin, xmin, xmax)
        !   add to total histogram
            total_hist = total_hist + counts

        end do

       !   normalize
            norm = sum(total_hist)*dx

       !   save on file
            outfile = "diag_avg"//trim(str(levels(kk)))//".dat"

            open(42, file=outfile)
            do zz=1,nbin
                write(42, *) total_hist(zz)/norm
            end do
            close(42)

    end do

        deallocate(D)
        deallocate(levels)
        deallocate(evl)
        deallocate(s)
        deallocate(counts)
        deallocate(total_hist)

    print*, "Execution complete."

end program
