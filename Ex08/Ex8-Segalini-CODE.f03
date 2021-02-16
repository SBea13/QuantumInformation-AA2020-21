module quantum_state
! ---------------------------------------------------------------------
! --------------------------- DOCUMENTATION ---------------------------
! ---------------------------------------------------------------------
! quantum_state module:
!       contains
! ---------------------------------------------------------------------
! qstate type:
!    a type defined to include all the parameters for the quantum
!    systems considered.
!   INCLUDES:
!       d    = integer, dimension of a single Hilbert space
!       N    = integer, number of subsystems
!       sep  = logical, .TRUE. if separable
!       memo = integer*8, number of bytes occupied in memory by
!              the coefficients
!       psi  = double complex, dimension(:), allocatable,
!              stores the coefficients
! ---------------------------------------------------------------------
! allocate_pure_sep:
!    function to allocate in memory the vector containing the
!    coefficients of a pure separable system. It properly initialise
!    a qstate type @state. @state%sep is set equal to .TRUE.
!   INPUT:
!       dimH = integer, intent(IN), provided dimension of a single
!              Hilbert space
!       ns   = integer, intent(IN) provided number of syubsystems
!   OUTPUT:
!       state = type qstate, properly allocated and initialised
! ---------------------------------------------------------------------
! allocate_pure:
!    function to allocate in memory the vector containing the
!    coefficients of a general pure system. It properly initialise
!    a qstate type @state. @state%sep is set equal to .FALSE.
!   INPUT:
!       dimH = integer, intent(IN), provided dimension of a single
!              Hilbert space
!       ns   = integer, intent(IN) provided number of syubsystems
!   OUTPUT:
!       state = type qstate, properly allocated and initialised
! ---------------------------------------------------------------------
! getElement:
!   function that computes the element of the tensor representing the
!   total system, given the array @psi_elem, corresponding to the
!   indexes in @idx.
!   INPUT:
!       psi_elem = double complex, dimension(:, :), intent(IN),
!                  wave-function coefficients
!       idx      = integer, dimension(:), intent(IN), list of indexes
!       N        = integer, intent(IN), number of subsystems
!   OUTPUT:
!       element = double complex, corresponding tensor element
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
    integer aa ! index for loops

    type qstate

        integer   :: d
        integer   :: n
        logical   :: sep
        integer*8 :: memo

        double complex, dimension(:), allocatable :: psi

    end type

    contains

    function allocate_pure_sep(dimH, ns)result(state)

        integer, intent(in) :: dimH, ns
        type(qstate)        :: state

        state%d   = dimH
        state%n   = ns
        state%sep =.TRUE.

        allocate(state%psi(ns*dimH))

        state%memo = sizeof(state%psi)

    end function

    function allocate_pure(dimH, ns)result(state)

        integer, intent(in) :: dimH, ns
        type(qstate)        :: state

        state%d   = dimH
        state%n   = ns
        state%sep =.FALSE.

        allocate(state%psi(dimH**ns))

        state%memo = sizeof(state%psi)

    end function

    function getElement(psi_elem, idx, N)result(element)
        ! computes the Cs coefficients of the tensor, given the
        ! array psi_elem
        integer, dimension(:), intent(IN)          :: idx
        integer, intent(IN)                        :: N
        double complex, dimension(:,:), intent(IN) :: psi_elem
        double complex                             :: element

        element = dcmplx(1, 0)

        do aa = 1, N
            element = element * psi_elem(aa, idx(aa))
        end do

    end function getElement

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

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

end module

program ex8

    use quantum_state
    implicit none

    integer   :: N, d, ii, jj, kk
    integer*8 :: ll, hh, sz
    ! N = Number of subsystems
    ! d = Dimension of the each Hilbert space
    ! ii, jj, kk, ll, hh = indexes for loops
    ! sz = size support variable

    ! support variables for systems generation
    double precision :: re_psi, im_psi, norm

    ! element of a system
    double complex, dimension(:,:), allocatable :: psi_elem
    ! density matrix variables
    double complex, dimension(:,:), allocatable :: rho, rho_1, rhosq, rho_1sq
    double precision   :: trace

    ! general system psi
    type(qstate) :: psi_sys

    ! support variable for checking separability condition
    character(1) :: check_sep

    !call RANDOM_SEED()

    print*, "Number of subsystems N="
    read*, N

    print*, "Dimension of the each Hilbert space d="
    read*, d

    print*, "Separable? (y/n)"
    read*, check_sep

    allocate(psi_elem(N, d))

    if (check_sep == "y") then
        psi_sys = allocate_pure_sep(d, N)
       !call exit(0)
    else
        psi_sys = allocate_pure(d, N)
    end if

    sz = sizeof(psi_elem)

    print*, "Size of psi_elem in memory: ", sz, " bytes"
    print*, "Size of psi in memory: ", psi_sys%memo, " bytes"

    ! generation of random states
    do ii = 1, N
        norm = 0
        do jj = 1, d
            call random_number(re_psi)
            call random_number(im_psi)
            psi_elem(ii, jj) = dcmplx(re_psi, im_psi)
            norm = norm + ABS(psi_elem(ii, jj))**2
        end do
        do jj = 1, d
            psi_elem(ii, jj) = psi_elem(ii, jj)/(sqrt(norm))
        end do
    end do

    ! tensor of the whole system
    do ll = 1, d**N
            psi_sys%psi(ll) = getElement(psi_elem, to_tensor(ll, N, d), N)
    end do

    !------- DENSITY MATRIX -------!
    allocate( rho(d**N, d**N), rhosq(d**N, d**N) )
    sz = sizeof(rho)
    print*, "Size of rho in memory: ", sz, " bytes"

    ! compute density matrix
    do ll = 1, d**N
        do hh = 1, d**N
            rho(ll, hh)  = psi_sys%psi(ll)*CONJG(psi_sys%psi(hh))
         end do
    end do

    ! compute density matrix squared
    rhosq = MATMUL(rho, rho)

    trace = 0.d0

    do ll = 1, d**N
        trace = trace + real(rho(ll, ll))
    end do

    print*, "Trace rho = ", trace

    trace = 0.d0

    do ll = 1, d**N
        trace = trace + real(rhosq(ll, ll))
    end do

    print*, "Trace rhosq = ", trace

    !Computing rho_1 and rho_1^2

    allocate( rho_1(d**(N-1), d**(N-1)), rho_1sq(d**(N-1), d**(N-1)) )
    sz = sizeof(rho_1)
    print*, "Size of rho_1 in memory: ", sz, " bytes"

    do ii = 1, d**(N-1)
        do jj = 1, d**(N-1)
            rho_1(ii, jj) = 0
            do kk = 1, d
                rho_1(ii, jj) = rho_1(ii, jj) + rho(to_mat([ii-1, kk-1], d), to_mat([jj-1, kk-1], d))
            end do
        end do
    end do

    rho_1sq = MATMUL(rho_1, rho_1)

    trace = 0.d0

    do ll = 1, d**(N-1)
        trace = trace + real(rho_1(ll, ll))
    end do

    print*, "Trace rho_1 = ", trace

    trace = 0.d0

    do ll = 1, d**(N-1)
        trace = trace + real(rho_1sq(ll, ll))
    end do

    print*, "Trace rho_1sq = ", trace

    !deallocate(psi_sys%psi)
    !deallocate(psi_elem)
    !deallocate(rho, rhosq)
    !deallocate(rho_1, rho_1sq)

end program
