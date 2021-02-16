module quantum_util
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! quantum_util module:
!       contains functions and subroutines to compute tensor products
!       and to handle matrix outputs
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
! tens_prod:
!   general function to compute Kronecker product @tp between matrices @M1
!   and @M2
!   INPUT:
!       M1 = double complex, dimension(:,:), intent(IN), first input matrix
!       M2 = double complex, dimension(:,:), intent(IN), second input
!            matrix
!   OUTPUT:
!       tp = double complex, dimension(:,:), allocatable, result of
!            tensor product
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
! printmatrix:
!    subroutine to print real matrices in a clearer form
!   INPUT:
!       A = double precision, dimension(:,:), matrix to be printed
!   OUTPUT:
!       print on screen matrix @A arranged in a "grid" form (elements
!       arranged in rows and columns), formatted properly
! ---------------------------------------------------------------------
! writematrix:
!    subroutine to write on file @filename a real matrix @A
!    in a clearer form
!   INPUT:
!       A = double precision, dimension(:,:), matrix to be written
!       filename = character(*), intent(IN), name of output file
!   OUTPUT:
!       write on file @filename a real matrix @A arranged in a "grid"
!       form (elements arranged in rows and columns), formatted properly
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------

    implicit none
    integer aa, bb ! indexes for loops

    contains

    function id_tens_mat(i_dim, N, matrix)result(txm)
        ! tensor product identity (X) matrix
        integer, intent(IN)                                      :: i_dim
        double complex, dimension(:,:), allocatable, intent(IN)  :: matrix
        integer, intent(IN)                                      :: N

        double complex, dimension(:,:), allocatable              :: txm

        allocate(txm(i_dim*N, N*i_dim))
        txm = (0d0, 0d0)

        txm(1:N, 1:N) = matrix

        do aa = 2, i_dim
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
        double complex, dimension(:,:), allocatable :: mxt

        ! identity matrix
        allocate(i_mat(i_dim, i_dim), mxt(i_dim*N, N*i_dim))

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
        integer, intent(IN)                         :: i_dim1, i_dim2, N
        double complex, dimension(:,:), intent(IN)  :: matrix

        double complex, dimension(:,:), allocatable :: prod
        double complex, dimension(:,:), allocatable :: supp

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

    function tens_prod(M1, M2)result(tp)
        ! Kronecker product M1 (X) M2
        double complex, dimension(:,:), intent(IN)  :: M1
        double complex, dimension(:,:), intent(IN)  :: M2

        integer                                     :: dim_1, dim_2
        double complex, dimension(:,:), allocatable :: tp

        dim_1 = size(M1, 1)
        dim_2 = size(M2, 1)

        allocate(tp(dim_1*dim_2, dim_1*dim_2))

        ! Kronecker product
        do aa = 1, dim_1
            do bb = 1, dim_1
                tp(dim_2*(aa-1)+1:dim_2*aa, dim_2*(bb-1)+1:dim_2*bb ) = M1(aa, bb)*M2
            end do
        end do

    end function tens_prod

    character(len=20) function str(k)
    !   Convert an integer to string
        integer, intent(in) :: k
        write (str, '(I3.3)') k
        str = adjustl(str)
    end function str

!    subroutine printmatrix(A)
!
!        implicit none
!        double complex, dimension(:,:) :: A
!
!        ! cycle to print matrix in the requested form
!        do aa = 1, ubound(A, 1)
!            print '(*(F0.2, 1x, SP, F0.2, "i", 2x))', A(aa, :)
!        end do
!    end subroutine

    subroutine printmatrix(A)

        implicit none
        double precision, dimension(:,:) :: A

        ! cycle to print matrix in the requested form
        do aa = 1, ubound(A, 1)
            print '(*(F0.1, 2x))', A(aa, :)
        end do
    end subroutine

    subroutine writematrix(A, filename)

        character(*), intent(IN)         :: filename
        double precision, dimension(:,:) :: A

        open(100, file=filename)
        ! cycle to print matrix in the requested form
        do aa = 1, ubound(A, 1)
            write (100,'(*(F0.4, 2x))') A(aa, :)
        end do

        close(100)

    end subroutine

end module

module ising_util
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! ising_util module:
!       contains functions and subroutines to perform RSRG algorithm for
!       the Ising model
! ---------------------------------------------------------------------
! identity:
!   function that computes the the identity matrix @ide of given
!   dimension @N
!   INPUT:
!       N      = integer, intent(IN), dimension of identity matrix
!   OUTPUT:
!       matrix = double complex, dimension(:,:), allocatable,
!                identity matrix
! ---------------------------------------------------------------------
! Projection:
!   subroutine that given the projector @P and a matrix @M, computes
!   @M_tilde = transpose(P) M P
!   The function checks sizes of input matrices and returns an error
!   message if the operation is not possible, moreover it checks for
!   previous allocation of output matrix @M_tilde, and reallocate it
!   with the correct dimensions.
!   INPUT:
!       P       = double complex, dimension(:,:), intent(IN),
!                 projection operator
!       M       = double complex, dimension(:,:), intent(IN),
!                 matrix to be projected
!   OUTPUT:
!       M_tilde = double complex, allocatable, dimension(:,:), intent(OUT)
! ---------------------------------------------------------------------
! findkEigenvalue:
!   subroutine to compute the first @kk eigenvalues and eigenvectors of a
!   given @matrix and to store them in the @eigenvalues, @Z array, using
!   the LAPACK subroutine DSYEVR.
!   INPUT:
!       matrix      = double complex, dimension(:,:), allocatable,
!                     intent(IN), input matrix
!       eigenvalues = double precision, dimension(:), allocatable,
!                     intent(INOUT), vector containing the eigenvalues in
!                     ascending order
!       kk          = integer, number of eigenvalues to be computed
!       Z           = double precision, dimension(:,:), allocatable,
!                     intent(OUT), first @kk computed eigenvectors
!   OUTPUT:
!       @eigenvalues entries updated with the first @kk eigenvalues in
!       increasing order, @Z containing the first @kk eigenvectors
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
    use quantum_util
    implicit none

    contains

    function identity(N)result(ide)

        integer, intent(IN)                          :: N
        double complex, dimension(:, :), allocatable :: ide

        allocate(ide(N, N))
        ide = (0d0, 0d0)

        do aa = 1, N
            ide(aa, aa) = (1d0, 0d0)
        end do

    end function identity

    subroutine Projection(P, M, M_tilde)

        double complex, dimension(:,:), intent(IN)               :: P, M
        double complex, allocatable, dimension(:,:), intent(OUT) :: M_tilde

        if ( ( (size(P, 1) .eq. size(M, 2)) .and. (size(M, 2) .eq. size(M, 1)) ) )then
            if ( allocated(M_tilde) ) then
                deallocate(M_tilde)
            end if
            allocate( M_tilde( size(P,2), size(P,2) ) )
            M_tilde = matmul( matmul( transpose(P), M), P )
        else
            print*, shape(P)
            print*, shape(M)
            print*, "Incompatible dimensions"
            stop
        end if

    end subroutine Projection

    subroutine findkEigenvalue(matrix, eigenvalues, kk, Z)

        double complex, dimension(:,:), allocatable, intent(IN)    :: matrix
        double precision, dimension(:), allocatable, intent(INOUT) :: eigenvalues
        double precision, dimension(:,:), allocatable, intent(OUT) :: Z

        double precision, dimension(:,:), allocatable              :: supp
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
        call DSYEVR('V', 'I', 'U', N, supp, N, VL, VU, 1, kk, &
                    ABSTOL, M, eigenvalues, Z, N, ISUPPZ, &
                    WORK, LWORK, IWORK, LIWORK, INFO)
        LWORK = min(LWMAX, int(WORK(1)))
        LIWORK = min(LWMAX, int(IWORK(1)))

    !   compute eigenvalues
        !print*, "Computing eigenvalues..."
        call DSYEVR('V', 'I', 'U', N, supp, N, VL, VU, 1, kk, &
                    ABSTOL, M, eigenvalues, Z, N, ISUPPZ, &
                    WORK, LWORK, IWORK, LIWORK, INFO)

        if (INFO .NE. 0) then
            print*, 'Failure!'
            stop
        end if

        deallocate(IWORK, WORK, ISUPPZ)

    end subroutine findkEigenvalue

end module

program RSRG

    use quantum_util
    use ising_util
    implicit none

    integer   :: N, d, ii, jj, mm, n_it
    ! N = Number of spins
    ! d = Dimension of the each Hilbert space
    ! mm = dimension of subspace
    ! ii, jj = indexes for loops
    ! n_it = number for iterations in RSRG algorithm

    double precision :: lambda, step, ground_state, supp_gs, eps, fsize, start, finish
    ! lambda = external interaction strength
    ! step = lambda increment for final plotting
    ! ground_state = E0, first energy level
    ! supp_gs = support variable for computing E0
    ! eps = threshold for RGRG algorithm
    ! fsize = final system size in RSRG algorithm
    ! start, finish = variables for computing execution time

    ! Hamiltonians
    double complex, dimension(:,:), allocatable   :: H, H_tilde
    ! Interaction Hamiltonians and projector
    double complex, dimension(:,:), allocatable   :: P, A, B
    double complex, dimension(:,:), allocatable   :: A_tilde, B_tilde
    ! eigenvalues and eigenvectors
    double precision, dimension(:,:), allocatable :: evct
    double precision, dimension(:), allocatable   :: evls
    ! Pauli matrices
    double complex, dimension(:,:), allocatable   :: sigma_x, sigma_z!, sigma_y
    ! support matrices for initial H
    double complex, dimension(:,:), allocatable   :: temp, temp_int
    ! output file name
    character(99)                                 :: outfile

    !Pauli matrices definition
    allocate(sigma_x(2,2), sigma_z(2,2))

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

    ! set system parameters
    print*, "Insert number of spins N = "
    read*, N

    ! create and open output file
    outfile = "N"//trim(str(N))//".dat"
    open(11, file = outfile, position="append")

!    print*, "Interaction strength lambda = "
!    read*, lambda
    lambda = 0
    step = 0.01

    ! allocation of H and support variables
    allocate(H(d**N , d**N), temp(d**N, d**N), temp_int(d**N, d**N))
    allocate(H_tilde(d**(2*N), d**(2*N)), evls(d**(2*N)))

    ! fix dimension of subspace for RSRG algorithm
    mm = d**N
    ! condition for stopping loop of RSRG
    eps = 1d-12

    ! loop over range of lambda in [0;3]
    do jj = 0, 300

        lambda = dble(jj)*step

        ! ----- Hamiltonian initilalization -----
        print*, "Initialize Hamiltonian..."

        H = (0d0, 0d0)

        temp = (0d0, 0d0)
        temp_int = (0d0, 0d0)

        ! Interaction with the external field
        do ii=1, N
            temp = id_tens_mat_tens_id(d**(ii-1), d, sigma_z, d**(N-ii))

            H = H + lambda*temp
        end do

        temp = (0d0, 0d0)

        ! Interaction between nearest neighbors
        do ii =1, N-1
            ! on sigma_x^ii
            temp = id_tens_mat_tens_id(d**(ii-1), d, sigma_x, d**(N-ii))
            ! on sigma_x^(ii+1)
            temp_int = id_tens_mat_tens_id(d**ii, d, sigma_x, d**(N-ii-1))
            ! combine the two interactions
            temp = MATMUL(temp, temp_int)

            H = H + temp
        end do

        !deallocate(temp, temp_int)

        ! initialize quantities for RSRG
        print*, "Initialize RSRG algorithm..."

        ! time variables
        start = 0
        finish = 0
        ! initialize Hamiltonian of subsystems interaction
        A = id_tens_mat(d**(N-1), d, sigma_x)
        B = mat_tens_id(d**(N-1), d, sigma_x)
        ! starting values
        ground_state = -1
        supp_gs = 0
        n_it = 1
        fsize = dble(N)

        print*, "RSRG algorithm..."
        ! starting time for iteration
        call cpu_time(start)
        ! iteration
        do while(abs(ground_state - supp_gs) > eps)

            print*, "ITERATION:", n_it
            ground_state = supp_gs

            ! build new Hamiltonian and interaction operators

            H_tilde = tens_prod(identity(d**N), H) + tens_prod(H, identity(d**N)) + tens_prod(A, B)
            A_tilde = tens_prod(A, identity(d**N))
            B_tilde = tens_prod(identity(d**N), B)

            ! double system size
            fsize = 2 * fsize

            ! diagonalization and computation of eigenvalues/eigenvectors
            call findkEigenvalue(H_tilde, evls, mm, evct)
            supp_gs = evls(1)/dble(fsize)

            ! build operator P
            P = dcmplx(evct)

            !projection on the subspace
            call Projection(P, H_tilde, H)
            call Projection(P, A_tilde, A)
            call Projection(P, B_tilde, B)

            ! increment iteration number
            n_it = n_it + 1

        end do
        ! update final ground state
        ground_state = supp_gs
        ! finishing time
        call cpu_time(finish)
        ! write on file final output
        write(11, *) lambda, ground_state, fsize, n_it, finish-start

    end do

    close(11)

end program
