module matrix_operations
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! matrix_operations module:
!       contains utility functions for matrix handling
! ---------------------------------------------------------------------
! RandMat:
!    function to initialize a N-by-M matrix (a) with random float
!    numbers in (0,10)
!   INPUT:
!       N = integer, number of matrix rows
!       M = integer, number of matrix columns
!   OUTPUT:
!       a = real*4, dimension(N,M), N-by-M matrix
! ---------------------------------------------------------------------
! ComputeError:
!    function that compute the difference (element-wise) between two
!    matrices, then adds up the results in order to evaluate a certain
!    kind of error. Includes a check on matrix shapes.
!   INPUT:
!       mat1, mat2 = matrices of which we want to compute the "error"
!   OUTPUT:
!       e = real, defined as the sum of the element-wise differences
!           between @mat1 and @mat2
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
    integer :: aa, bb ! indexes for initializing matrix

    contains

    function RandMat(N, M)result(a)

        integer :: N, M
        real*4, dimension(N,M) :: a

        call random_seed()

        ! cycle to fill in matrix a
        do aa = 1,N
            do bb = 1,M
                a(aa, bb) = RAND(0)*10
            end do
        end do

    end function RandMat

    function ComputeError(mat1, mat2)result(e)

        real*4, dimension(:,:), intent(IN) :: mat1, mat2
        integer                            :: ii, jj
        integer, dimension(2)              :: shape1, shape2
        real*8                             :: e

!       initialize error
        e = 0

!       check matrix sizes
        shape1 = shape(mat1)
        shape2 = shape(mat2)

!       compute error only if matrix size is correct
        if ((shape1(1) == shape2(1)) .AND. (shape1(2) == shape2(2))) then
            do ii = 1, shape1(1)
                do jj = 1, shape1(2)
                    e = e + ABS(mat1(ii, jj) - mat2(ii, jj))
                end do
            end do
        else
            print*, "Incorrect matrix shape!"
        end if

    end function

    subroutine printmatrix(A)

        implicit none
        real*4, dimension(:,:) :: A
        integer :: pp

        ! cycle to print matrix in the requested form
        do pp = 1, ubound(A, 1)
            print *, A(pp, :)
        end do
    end subroutine

end module matrix_operations

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


program mulmat
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! mulmat program:
!       contains three algorithms to compute matrix multiplication
!       and evaluate computation time. It uses the @debugging_mod and
!       @matrix_operations modules.
! ---------------------------------------------------------------------
! PRE-CONDITIONS:
!       - requires the user to insert 4 integers, representing matrix
!         @A and @B sizes
!       - requires the user to select parameters of debug-mode
!       - all matrices are allocated
! ---------------------------------------------------------------------
! After pre-conditions are verified, the program use the @RandMat
! function to randomly initialize matrices @A and @B with numbers
! in (0,10). Then, it performs matrix multiplication:
!   1) with FORTRAN intrinsic function @matmul
!   2) with "manual" method, by cycling on indexes kk-jj-ii
!   2) with "manual" method, by cycling on indexes ii-jj-kk
! For all the three methods the exectution time is computed using
! the @cpu_time FORTRAN built-in function.
!
! Matrices can be printed on screen with the @printmatrix subroutine.
! Of course, this is strongly discouraged if the matrix size is big.
! In this test program, small sizes are considered, hence matrices
! are displayed as outputs.
! ---------------------------------------------------------------------
! POST-CONDITIONS:
!   - check correctness of matrix product by computing the differences
!     between the manually derived results and the FORTRAN intrinsic
!     function (element-wise) and then summing them
!   - deallocation of matrices
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

!   introduce modules
    use matrix_operations
    use debugging_mod

    implicit none
!   dimension of matrices
    integer                             :: N
!   matrices
    real*4, dimension(:,:), allocatable :: A, B, prodF, myprod1, myprod2
!   support variables to compute execution time and errors
    real*8                              :: start, finish, e1, e2
!   indexes for computing matrix product manually
    integer                             :: ii, jj, kk
!   variables for debugging process
    character(1)                        :: debugstr, verbosestr, stopstr
    logical                             :: debug, verbose, stopprg
!   Output variables
    character(99)                       :: inputfile, out_F, out_1, out_2

!   open input file and read matrix dimension
!    print*, "Enter input filename:"
!    read*, inputfile
    inputfile = "N.dat"
    open(1, file=inputfile)

!   loop on all the sizes included in inputfile
    do
    read(1, *, end=2) N

!   set automatically debug parameters, don't ask user
    debug = .TRUE.
    verbose = .FALSE.
    stopprg = .FALSE.

!  UNCOMMENT THIS FOR USER-DEFINED DEBUGGING PROCESS
! ---------------------------------------------------
!    print*, "Check matrix sizes? [y/n]:"
!    read*, debugstr
!    print*, "Print full output? [y/n]:"
!    read*, verbosestr
!    print*, "Stop on error? [y/n]:"
!    read*, stopstr

!    if (debugstr .NE. "y") then
!        debug = .FALSE.
!        print*, "---------------------"
!        print*, "Un-setting debug mode"
!        print*, "---------------------"
!    else
!        debug = .TRUE.
!    end if
!
!    if (verbosestr .NE. "y") then
!        verbose = .FALSE.
!        print*, "---------------------"
!        print*, "Short output selected"
!        print*, "---------------------"
!    else
!        debug = .TRUE.
!    end if
!
!    if (stopstr .NE. "y") then
!        stopprg = .FALSE.
!        print*, "-----------------------"
!        print*, "Continue despite errors"
!        print*, "-----------------------"
!    else
!        stopprg = .TRUE.
!    end if
! -----------------------------------------------

!   debugging process: check all matrix sizes and if the product is legit
!    if (debug) then
!        call debugging(N > 0, "size > 0", verbose, stopprg, N)
!    end if

    ! allocate matrices, PRE-CONDITION
    allocate(A(N, N))
    allocate(B(N, N))
    allocate(prodF(N, N))
    allocate(myprod1(N, N))
    allocate(myprod2(N, N))

!   generate matrices A and B
    A = RandMat(N, N)
    B = RandMat(N, N)

!   commenting these lines hides matrix on-screen output
! ---------------------------------------------------------------------
!    print*, "A matrix:"
!    call printmatrix(A)
!    print*, "B matrix:"
!    call printmatrix(B)
! ---------------------------------------------------------------------

!   support variables for computing execution time
    start = 0
    finish = 0

! ###### Fortran intrinsic function ######

    out_F = "out_f.dat"
    open(unit=10, file=out_F, position='append')

    call cpu_time(start)
    prodF = matmul(A, B)
    call cpu_time(finish)

!   commenting these lines hides matrix on-screen output
! ---------------------------------------------------------------------
!    print*, "A*B Fortran intrinsic function:"
!    call printmatrix(prodF)
! ---------------------------------------------------------------------
    print *, "Time taken:"
    print*, finish-start

    write(10, *) N, finish-start
    close(10)

! ########################################

! ######## kk-jj-ii manual method ########

    out_1 = "out_k.dat"
    open(unit=11, file=out_1, position='append')

    call cpu_time(start)

    do ii = 1, N
        do jj = 1, N
            do kk = 1, N
                myprod1(ii, jj) = myprod1(ii, jj) + A(ii, kk)*B(kk, jj)
            end do
        end do
    end do

    call cpu_time(finish)

!   commenting these lines hides matrix on-screen output
! ---------------------------------------------------------------------
!    print*, "A*B my function (kk-jj-ii):"
!    call printmatrix(myprod1)
! ---------------------------------------------------------------------
    print *, "Time taken:"
    print *, finish-start

    write(11, *) N, finish-start
    close(11)

! ########################################

! ######## ii-jj-kk manual method ########

    out_2 = "out_i.dat"
    open(unit=12, file=out_2, position='append')

    call cpu_time(start)

    ! matrix multiplication with "manual" method:
    ! first cycling on ii, then jj, then kk
    do kk = 1, N
        do jj = 1, N
            do ii = 1, N
                myprod2(ii, jj) = myprod2(ii, jj) + A(ii, kk)*B(kk, jj)
            end do
        end do
    end do

    call cpu_time(finish)

!   commenting these lines hides matrix on-screen output
! ---------------------------------------------------------------------
!    print*, "A*B my function (ii-jj-kk):"
!    call printmatrix(myprod2)
! ---------------------------------------------------------------------
    print *, "Time taken:"
    print *, finish-start

    write(12, *) N, finish-start
    close(12)

!   computing the "error", POST-CONDITION
!    e1 = ComputeError(myprod1, prodF)
!    e2 = ComputeError(myprod2, prodF)
!
!    print*, "Error(kk-jj-ii)= ", e1
!    print*, "Error(ii-jj-kk)= ", e2

!   deallocation of matrices, POST-CONDITION
    deallocate(A)
    deallocate(B)
    deallocate(prodF)
    deallocate(myprod1)
    deallocate(myprod2)

    end do

!   finishing message, after complete execution
    2 print *, "Execution complete."
    close(1)

end program
