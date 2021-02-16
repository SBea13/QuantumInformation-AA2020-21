module matrix_operations

    implicit none
    integer :: aa, bb

    contains

    ! fill a N by M matrix with random values between 0 and 10
    function RandMat(N, M)result(a)

        integer :: N, M
        real*4, dimension(N,M) :: a

        do aa = 1,N
            do bb = 1,M
                a(aa, bb) = RAND(0)*10
            end do
        end do

    end function RandMat

    ! check if dimensions of two matrices is compatible with multiplication
    function checksize(a, b)result(check)

        real*4, dimension(:, :) :: a
        real*4, dimension(:, :) :: b
        logical :: check

        if (size(a,2) .eq. size(b,1)) then
            check = .TRUE.
        else
            check = .FALSE.
        end if

    end function checksize

    ! to print matrices nicely
    subroutine printmatrix(A)

        implicit none
        real*4, dimension(:,:) :: A
        integer :: pp

        do pp = 1, ubound(A, 1)
            print *, A(pp, :)
        end do
    end subroutine

end module matrix_operations


program mulmat

    use matrix_operations

    implicit none
    integer :: rowa, rowb, cola, colb
    real*4, dimension(:,:), allocatable :: A, B, prodF, myprod1, myprod2
    real*8 :: start, finish
    integer :: ii, jj, kk
    logical :: check

    print*, "Enter the size of matrix A:"
    read*, rowa, cola

    print*, "Enter the size of matrix B:"
    read*, rowb, colb

    ! allocate matrices
    allocate(A(rowa, cola))
    allocate(B(rowb, colb))
    allocate(prodF(rowa, colb))
    allocate(myprod1(rowa, colb))
    allocate(myprod2(rowa, colb))

    A = RandMat(rowa, cola)
    B = RandMat(rowb, colb)

   ! print *, "Matrix A:"
   ! call printmatrix(A)
   ! print *, "Matrix B:"
   ! call printmatrix(B)

    check = checksize(A, B)

    if (check .eqv. .FALSE.) then
        print *, "Wrong sizes, can't multiply these matrices!"
        stop
    end if

    start = 0
    finish = 0

    ! matrix multiplication with Fortran intrinsic function
    call cpu_time(start)
    prodF = matmul(A, B)
    call cpu_time(finish)

    print*, "A*B Fortran intrinsic function:"
   ! call printmatrix(prodF)
   ! print *, "Time taken:"
    print*, finish-start

    call cpu_time(start)

    ! matrix multiplication with "manual" method:
    ! first cycling on kk, then jj, then ii
    do ii = 1, rowa
        do jj = 1, colb
            do kk = 1, rowb
                myprod1(ii, jj) = myprod1(ii, jj) + A(ii, kk)*B(kk, jj)
            end do
        end do
    end do

    call cpu_time(finish)

    print*, "A*B my function (kk-jj-ii):"
    ! call printmatrix(myprod1)
    ! print *, "Time taken:"
    print *, finish-start

    call cpu_time(start)

    ! matrix multiplication with "manual" method:
    ! first cycling on ii, then jj, then kk
    do kk = 1, rowb
        do jj = 1, colb
            do ii = 1, rowa
                myprod2(ii, jj) = myprod2(ii, jj) + A(ii, kk)*B(kk, jj)
            end do
        end do
    end do

    call cpu_time(finish)

    print*, "A*B my function (ii-jj-kk):"
    ! call printmatrix(myprod2)
    ! print *, "Time taken:",
    print *, finish-start

    deallocate(A)
    deallocate(B)
    deallocate(prodF)
    deallocate(myprod1)
    deallocate(myprod2)

end program
