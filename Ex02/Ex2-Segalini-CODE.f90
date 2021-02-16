module MATRICES

    implicit none

    type DMATRIX
        ! dimension of matrix
        integer, dimension(2)                   :: N
        ! elements
        complex*16, dimension(:,:), allocatable :: elements
        ! trace and determinant
        complex*16                              :: trace, det
    end type

    interface operator(.Adj.)
        module procedure MatAdjoint
    end interface

    interface operator(.Trace.)
        module procedure MatTrace
    end interface

    contains

    ! function to initialize a DMATRIX
    ! complex random numbers with Re, Im in (0,10)
    ! given the number of rows and columns
    function InitMatrix(rows, cols) result(A)

        integer                                  :: rows, cols
        real*8, dimension(:, :), allocatable     :: xx, yy
        type(DMATRIX)                            :: A

        allocate(A%elements(rows, cols))
        allocate(xx(rows, cols))
        allocate(yy(rows, cols))

        call random_number(xx)
        call random_number(yy)

        A%N = (/ rows, cols /)
        A%elements = cmplx(xx, yy, kind(1d0))
        A%det = 1 ! not implemented in this code
        A%trace = MatTrace(A)

    end function

    function MatTrace(A) result(t)

        type(DMATRIX), intent(IN) :: A
        complex*16                :: t
        integer                   :: ii

        t = (0d0, 0d0)
        ! if the matrix is non-squared matrix, trace is undefined
        ! trace is initialized to 0,
        ! function changes its value only for squared matrices

        if ( A%N(1) .EQ. A%N(2) ) then
            do ii=1,A%N(1)
                t = t + A%elements(ii, ii)
            end do
        end if

    end function

    function MatAdjoint(A) result(Aplus)

        type(DMATRIX), intent(IN) :: A
        type(DMATRIX)             :: Aplus

        Aplus%N(1) = A%N(2)
        Aplus%N(2) = A%N(1)

        ! defined with intrinsic fortran functions
        Aplus%elements = conjg(transpose(A%elements))

        Aplus%det = conjg(A%det)
        Aplus%trace = conjg(A%trace)

    end function

    ! to print matrices nicely
    subroutine printmatrix(A)

        implicit none
        type(DMATRIX) :: A
        integer       :: pp, nn

        do pp = 1, ubound(A%elements, 1)
            print '(*(F0.4, 1x, SP, F0.4, "i", 2x))', A%elements(pp, :)
        end do

    end subroutine

    ! to write on file DMATRICES nicely
    ! it writes all the attributes of a DMATRIX
    subroutine writematrix(A, filename)

        implicit none
        type(DMATRIX)        :: A
        integer              :: pp
        character(len = 99)  :: filename

        open(1, file=filename)

        write(1, *) "Matrix:"

        do pp = 1, ubound(A%elements, 1)
            write(1, '(*(F0.4, 1x, SP, F0.4, "i", 2x))')  A%elements(pp, :)
        end do

        write(1, *) "Trace:"
        write(1, '(F0.4, 1x, SP, F0.4, "i")') A%trace

        write(1, *) "Determinant:"
        write(1, '(F0.4, 1x, SP, F0.4, "i")') A%det

        close(1)

    end subroutine

end module MATRICES


program testDMATRIX

    use MATRICES

    implicit none
    integer             :: nrow, ncol
    type(DMATRIX)       :: A, Aplus
    complex*16          :: trace, traceplus
    character(len = 99) :: filenameA
    character(len = 99) :: filenameAplus

    print*, "Enter the size of matrix A (rows, columns):"
    read*, nrow, ncol

    A = InitMatrix(nrow, ncol)

    print *, "Matrix A:"
    call printmatrix(A)

    Aplus = .Adj.A

    print *, "A+:"
    call printmatrix(Aplus)

    trace = .Trace.A
    traceplus = .Trace.Aplus

    print *, "Trace of A:"
    print '(F0.4, 1x, SP, F0.4, "i")', trace

    print *, "Trace of A+:"
    print '(F0.4, 1x, SP, F0.4, "i")', traceplus

    print*, "Enter filename for A:"
    read*, filenameA

    call writematrix(A, filenameA)

    print*, "Enter filename for A+:"
    read*, filenameAplus

    call writematrix(Aplus, filenameAplus)

end program
