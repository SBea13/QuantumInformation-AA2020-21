program test

    implicit none
    double precision, dimension(:), allocatable::vect
    integer *2 ii

    allocate(vect(10))

    do ii = 1,10
        vect(ii) = ii
    end do

print *, vect, "It works!"
deallocate(vect)
end program test
