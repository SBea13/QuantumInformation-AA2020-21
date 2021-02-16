module quantum_state
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! quantum_state module:
!       contains
! ---------------------------------------------------------------------
! to_mat:
!   function that converts indexes from tensor form to matricial form
!   INPUT:
!       idx_t = integer, dimension(:), intent(IN), indexes to be
!               converted
!       d     = integer, intent(IN), dimension of Hilbert space
!   OUTPUT:
!       mm    = integer*8, index converted for matricial notation
! ---------------------------------------------------------------------
! to_tensor:
!   function that converts indexes from matricial form to tensor form
!   INPUT:
!       mm    = integer*8, index to be converted from matricial notation
!       N     = integer, intent(IN), number of subsystems
!       d     = integer, intent(IN), dimension of Hilbert space
!   OUTPUT:
!       idx_t = integer, dimension(:), intent(IN), indexes converted
!               to tensor form
! ---------------------------------------------------------------------
! id_tens_mat:
!   function that performs the tensor (Kronecher) product of between
!   an identity matrix of dimension @i_dim x @i_dim and a given
!   @matrix of dimension @N x @N.
!   INPUT:
!       i_dim  = integer, intent(IN), dimension of identity matrix
!       N      = integer, intent(IN), dimension of input @matrix
!       matrix = double complex, dimension(:,:), allocatable, intent(IN),
!                input matrix
!   OUTPUT:
!       txm    = double complex, dimension(i_dim*N, N*i_dim), result of
!                tensor product
! ---------------------------------------------------------------------
! mat_tens_id:
!   function that performs the tensor (Kronecher) product between
!   a given @matrix of dimension @N x @N and an identity matrix of
!   dimension @i_dim x @i_dim.
!   INPUT:
!       i_dim  = integer, intent(IN), dimension of identity matrix
!       N      = integer, intent(IN), dimension of input @matrix
!       matrix = double complex, dimension(:,:), allocatable, intent(IN),
!                input matrix
!   OUTPUT:
!       mxt    = double complex, dimension(i_dim*N, N*i_dim), result of
!                tensor product
! ---------------------------------------------------------------------
! id_tens_mat_tens_id:
!   function that performs the tensor product between an identity matrix
!   of dimension @i_dim1 x @i_dim1, a given @matrix of dimension @N x @N
!   and an identity matrix of dimension @i_dim2 x @i_dim2, in this
!   respective order.
!   INPUT:
!       i_dim  = integer, intent(IN), dimension of first identity matrix
!       N      = integer, intent(IN), dimension of input @matrix
!       matrix = double complex, dimension(:,:), allocatable, intent(IN),
!                input matrix
!       i_dim  = integer, intent(IN), dimension of second identity matrix
!   OUTPUT:
!       prod   = double complex, dimension(i_dim1*i_dim2*N, i_dim1*i_dim2*N),
!                result of tensor product
! ---------------------------------------------------------------------
! findkEigenvalue:
!   subroutine to compute the first @kk eigenvalues of a given @matrix
!   and to store them in the @eigenvalues array, using the LAPACK subroutine
!   DSYEVR.
!   INPUT:
!       matrix      = double complex, dimension(:,:), allocatable,
!                     intent(IN), input matrix
!       eigenvalues = double precision, dimension(:), allocatable,
!                     intent(INOUT), vector containing the eigenvalues in
!                     ascending order
!       kk          = integer, number of eigenvalues to be computed
!   OUTPUT:
!       @eigenvalues entries updated with the first @kk eigenvalues in
!       increasing order
! ---------------------------------------------------------------------
! str:
!    function that convert the integer @k to a string, in order to
!    properly print/concatenate it to other strings (especially for
!    handling output file names). Format has been set to an integer with
!    3 digits.
!   INPUT:
!       k = integer, to be converted to string
!   OUTPUT:
!       return @k converted into a formatted string
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

    implicit none
    integer aa, bb ! indexes for loops

    contains

    function to_mat(idx_t, d)result(mm)
        ! converts indexes from tensor form to matricial form
        integer, dimension(:), intent(IN) :: idx_t
        integer, intent(IN)               :: d
        integer*8                         :: mm

        mm = 1

        do aa = 1, size(idx_t)
            mm = mm + idx_t(aa)*( d**(aa - 1) )
        end do

    end function to_mat

    function to_tensor(mm, N, d) result(idx_t)
        ! converts indexes from matricial form to tensor form
        integer*8             :: mm
        integer, intent(IN)   :: N, d
        integer, dimension(N) :: idx_t

        do aa = 1, N
            idx_t(aa) = modulo( ( (mm - 1)/(d)**(aa-1) ), d) + 1
        end do

    end function to_tensor

    function id_tens_mat(i_dim, N, matrix)result(txm)
        ! tensor product identity (X) matrix
        integer, intent(IN)                                      :: i_dim
        double complex, dimension(:,:), allocatable, intent(IN)  :: matrix
        integer, intent(IN)                                      :: N

        double complex, dimension(i_dim*N, N*i_dim) :: txm

        txm = (0d0, 0d0)

        txm(1:N, 1:N) = matrix

        do aa =2, i_dim
            bb = N*(aa-1) + 1
            txm( bb:(bb+N-1), bb:(bb+N-1) ) = matrix
        end do

    end function id_tens_mat

    function mat_tens_id(i_dim, N, matrix)result(mxt)
        ! tensor product matrix (X) identity
        integer, intent(IN)                         :: i_dim
        double complex, dimension(:,:), intent(IN)  :: matrix
        integer, intent(IN)                         :: N

        double complex, dimension(:,:), allocatable :: i_mat
        double complex, dimension(i_dim*N, N*i_dim) :: mxt

        ! identity matrix
        allocate(i_mat(i_dim, i_dim))

        mxt = (0d0, 0d0)
        i_mat = (0d0, 0d0)

        do aa = 1, i_dim
            i_mat(aa, aa) = (1d0, 0d0)
        end do

        do aa = 1, N
            do bb = 1, N
            mxt( (i_dim*(aa - 1)+1):(i_dim*aa), (i_dim*(bb - 1)+1):(i_dim*bb) ) = matrix(aa, bb)*i_mat
            end do
        end do

        deallocate(i_mat)

    end function mat_tens_id

    function id_tens_mat_tens_id(i_dim1, N, matrix, i_dim2)result(prod)
        ! identity (X) matrix (X) identity
        integer, intent(IN)                                        :: i_dim1, i_dim2, N
        double complex, dimension(:,:), intent(IN)                 :: matrix

        double complex, dimension(i_dim1*i_dim2*N,i_dim1*i_dim2*N) :: prod
        double complex, dimension(:,:), allocatable                :: supp

        ! matrix (X) id2
        if( i_dim2 > 0 ) then
            allocate( supp(N*i_dim2, N*i_dim2) )
            supp = mat_tens_id( i_dim2, N, matrix )
        end if

        ! id1 (X) (matrix (X) id2)
        if( i_dim1 > 0 ) then
            prod = id_tens_mat( i_dim1, N*i_dim2, supp )
            deallocate(supp)
        end if

    end function id_tens_mat_tens_id

    subroutine findkEigenvalue(matrix, eigenvalues, kk)

        double complex, dimension(:,:), allocatable, intent(IN)    :: matrix
        double precision, dimension(:), allocatable, intent(INOUT) :: eigenvalues

        double precision, dimension(:,:), allocatable              :: Z, supp
        integer :: LWORK, LIWORK, INFO, N, M, LWMAX, dims(2), kk
        double precision :: VL, VU, ABSTOL
        double precision, dimension(:), allocatable                :: WORK
        integer, dimension(:), allocatable                         :: ISUPPZ, IWORK

        dims = shape(matrix)
        N    = dims(1)
        LWMAX = 100000
        ABSTOL = 0.
        VU = 0.
        VL = 0.

        allocate(WORK(LWMAX))
        allocate(IWORK(LWMAX))
        allocate(ISUPPZ(2*kk))
        allocate(Z(N, kk))
        allocate(supp(N, N))

        supp = real(matrix)

    !   compute optimal size of workspace
        LWORK = -1
        LIWORK = -1
        !print*, "Computing optimal workspace..."
        call DSYEVR('N', 'I', 'U', N, supp, N, VL, VU, 1, kk, &
                    ABSTOL, M, eigenvalues, Z, N, ISUPPZ, &
                    WORK, LWORK, IWORK, LIWORK, INFO)
        LWORK = int(WORK(1))
        LIWORK = int(IWORK(1))

    !   compute eigenvalues
        !print*, "Computing eigenvalues..."
        call DSYEVR('N', 'I', 'U', N, supp, N, VL, VU, 1, kk, &
                    ABSTOL, M, eigenvalues, Z, N, ISUPPZ, &
                    WORK, LWORK, IWORK, LIWORK, INFO)

        if (INFO .NE. 0) then
            print*, 'Failure!'
            stop
        end if

        deallocate(IWORK, WORK, ISUPPZ, Z)

    end subroutine findkEigenvalue

    character(len=20) function str(k)
    !   Convert an integer to string
        integer, intent(in) :: k
        write (str, '(I3.3)') k
        str = adjustl(str)
    end function str

    subroutine printmatrix(A)

        implicit none
        double precision, dimension(:,:) :: A

        ! cycle to print matrix in the requested form
        do aa = 1, ubound(A, 1)
            print '(*(F0.4, 2x))', A(aa, :)
        end do
    end subroutine

!    subroutine writematrix(A, filename)
!
!        character(*), intent(IN)         :: filename
!        double precision, dimension(:,:) :: A
!
!        open(100, file=filename)
!        ! cycle to print matrix in the requested form
!        do aa = 1, ubound(A, 1)
!            write (100,'(*(F0.4, 2x))') A(aa, :)
!        end do
!
!        close(100)
!
!    end subroutine

end module

program ising

    use quantum_state
    implicit none

    integer   :: N, d, ii, n_evls
    integer*8 :: hh, jj, kk
    ! N = Number of spins
    ! d = Dimension of the each Hilbert space
    ! n_evls = number of eigenvalues to be computed
    ! ii, jj, kk, hh = indexes for loops

    ! external interaction strength
    double precision :: lambda
    ! time support variables
    double precision :: start, finish, time_H, time_diag

    ! Hamiltonian
    double complex, dimension(:,:), allocatable :: H
    ! eigenvalues
    double precision, dimension(:), allocatable :: evls
    ! Pauli matrices
    double complex, dimension(2, 2)             :: sigma_x, sigma_z!, sigma_y
    ! support matrices for computing tensor products
    double complex, dimension(:,:), allocatable :: temp, temp_int
    ! output file
    character(99)                               :: outfile

    !Pauli matrices definition
    sigma_x = 0.d0
    sigma_x(1,2) = 1.d0
    sigma_x(2,1) = 1.d0

    sigma_z = 0.d0
    sigma_z(1,1) = 1.d0
    sigma_z(2,2) = -1.d0

!    sigma_y = 0.d0
!    sigma_y(1,2) = -(0.d0,1.d0)
!    sigma_y(2,1) = (0.d0,1.d0)

    ! fix dimension d
    d = 2
    ! fix number of eigenvalues to compute
    n_evls = 4

    ! time support variables initialization
    start  = 0
    finish = 0

    ! set system parameters
    print*, "Insert number of subsystems N = "
    read*, N

    print*, "Interaction strength lambda = "
    read*, lambda

    ! create output file name according to the parameters set
    outfile = "N"//trim(str(N))//"lambda"//trim(str(nint(100*lambda)))//".dat"
    print*, outfile

    ! ----- Hamiltonian calculation -----
    ! allocation of H and support variables
    allocate(H(d**N , d**N), temp(d**N, d**N), temp_int(d**N, d**N))
    ! allocate vector for the eigenvalues
    allocate(evls(n_evls))

    ! Initialize main variables
    H = (0d0, 0d0)

    temp = (0d0, 0d0)
    temp_int = (0d0, 0d0)

    ! write H matrix representation and compute execution time
    call cpu_time(start)

    ! Interaction with the external field
    do ii=1, N
        temp = id_tens_mat_tens_id(d**(ii-1), d, sigma_z, d**(N-ii))

        H = H + lambda*temp
    end do
    temp = (0d0, 0d0)

    ! Interaction between nearest neighbors
    do ii =1, N-1
        ! on sigma_x^ii
        temp = id_tens_mat_tens_id(d**(ii-1), d, sigma_x, d**(N-))
        ! on sigma_x^(ii+1)
        temp_int = id_tens_mat_tens_id(d**ii, d, sigma_x, d**(N-ii-1))
        ! combine the two interactions
        temp = MATMUL(temp, temp_int)

        H = H + temp
    end do

    deallocate(temp, temp_int)
    call cpu_time(finish)

    !call printmatrix(real(H))

    time_H = finish - start

    print*, "Time to compute H = "
    print*, time_H

    evls = 0.d0

    ! solve eigenproblem for H and compute execution time
    ! find only first kk = 4 eigenvalues
    call cpu_time(start)
    call findkEigenvalue(H, evls, n_evls)
    call cpu_time(finish)

    time_diag = finish - start

    print*, "Time to solve eigenproblem = "
    print*, time_diag

    ! final output & write on file
    open(11, file = outfile)
    !print*, "Eigenvalues:"
    do ii = 1, n_evls
        write(11, *) evls(ii)
        !print*, evls(ii)
    end do

    close(11)

    deallocate(H, evls)

end program
