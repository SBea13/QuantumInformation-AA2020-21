program test_limits

    integer*2 :: aa, bb
    integer*4 :: cc, dd
    real*4    :: ee, ff
    real*8    :: gg, hh

    aa = 2000000
    bb = 1

    cc = 2000000
    dd = 1

    ee = 4e0*atan(1e0)*1e32
    ff = sqrt(2e0)*1e21

    gg = 4d0*datan(1d0)*1d32
    hh = sqrt(2d0)*1d21

    print*, "Sum of 2000000 and 1"
    ! here we observe an overflow phenomenon
    print*, "With INTEGER*2 types = ", aa+bb
    ! here the sum works fine
    print*, "With INTEGER*4 types = ", cc+dd

    print*, "Sum of big real numbers"
    ! here we observe an overflow
    print*, "With single precision = ", ee+ff
    ! here it works properly
    print*, "With double precision = ", gg+hh


end program
