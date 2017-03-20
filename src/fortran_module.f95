!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Instructions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To compile for use with Python/f2py, please run the following command, which
! requires NumPy and a Fortran compiler:
!
!    f2py -c fortran_module.f95 -m fortran_module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! condense_adjacency_matrix(B, A, V, k, m, n)
!
!   condenses the graph given by the weighted adjacency matrix A into a graph
!   given by the weighted adjacency matrix B by contracting the vertices in V
!

subroutine condense_adjacency_matrix(B, A, V, k, m, n)

    implicit none

    integer, intent(in) :: m, n
    integer, intent(in) :: V(m), k(n+1)
    double precision, intent(in) :: A(m, m)

    double precision, intent(out) :: B(n, n)

    integer :: i, j

    ! m: number of vertices; dimension of A
    ! n: number of SCCs; dimension of B
    ! V: vertices of the components; component(k(i): k(i+1)-1) are the vertices
    !    in component i
    ! k: vertex indices of the components; k(i) is the index of the first
    !    vertex in component i
    ! A: weighted adjacency matrix
    ! B: condensed weighted adjacency matrix

    do i=1,n
        do j=1,n
            if (i /= j) then
                B(i,j) = minimum_slice(V(k(i): k(i+1)-1), V(k(j): k(j+1)-1), k(i+1)-k(i), k(j+1)-k(j))
            else
                B(i,j) = 0.d0
            end if
        end do
    end do

contains

    function minimum_slice(columns, rows, p, q)

        integer, intent(in) :: p, q
        integer, intent(in) :: columns(p), rows(q)

        double precision :: minimum_slice

        double precision :: weight
        integer :: k, l

        minimum_slice = huge(0.d0)

        do l = 1, q
            do k = 1, p
                weight = A(columns(k), rows(l))
                if (weight > 0.d0) then
                    if (weight < minimum_slice) then
                        minimum_slice = weight
                    end if
                end if
            end do
        end do

        if (minimum_slice == huge(0.d0)) then
            minimum_slice = 0.d0
        end if

    end function minimum_slice

end subroutine condense_adjacency_matrix

!
! slice_array(B, A, columns, rows, m, n, p, q)
!
!   finds the slice of a (weighted) adjacency matrix given by the columns and
!   rows
!

subroutine slice_array(B, A, columns, rows, m, n, p, q)

    implicit none

    integer, intent(in) :: m, n, p, q
    integer, intent(in) :: columns(p), rows(q)
    double precision, intent(in) :: A(m, n)

    double precision, intent(out) :: B(p, q)

    integer :: i, j

    do j = 1, q
        do i = 1, p
            B(i, j) = A(columns(i), rows(j))
        end do
    end do

end subroutine slice_array

!
! strongly_connected_components(component, A, n)
!
!   finds the strongly connected components of the graph given by the
!   (weighted) adjacency matrix A
!

subroutine strongly_connected_components(component, A, n)

    implicit none

    integer, intent(in) :: n
    double precision, intent(in) :: A(n, n)

    integer, intent(out) :: component(n)

    integer :: i, j, k, l, u, v, w

    integer :: idx(n), lowlink(n), queue(n), subqueue(n)
    logical :: found(n), updated_queue

    idx = -1
    lowlink = -1
    found = .false.
    queue = 0
    subqueue = 0
    component = 0

    i = 0
    j = 0
    k = 0
    l = 0

    do u = 1, n
        if (found(u) .eqv. .false.) then
            k = k + 1
            queue(k) = u
        end if

        do while (k > 0)
            v = queue(k)
            if (idx(v) == -1) then
                i = i + 1
                idx(v) = i
            end if

            updated_queue = .false.
            do w = 1, n
                if (A(v, w)/=0.d0) then
                    if (idx(w) == -1) then
                        k = k + 1
                        queue(k) = w
                        updated_queue = .true.
                        exit
                    end if
                end if
            end do

            if (updated_queue .eqv. .false.) then
                lowlink(v) = idx(v)
                do w = 1, n
                    if (A(v, w)/=0.d0) then
                        if (found(w) .eqv. .false.) then
                            if (idx(w)>idx(v)) then
                                lowlink(v) = min(lowlink(v), lowlink(w))
                            else
                                lowlink(v) = min(lowlink(v), idx(w))
                            end if
                        end if
                    end if
                end do
                k = k-1

                if (lowlink(v) == idx(v)) then
                    found(v) = .true.
                    j = j + 1
                    component(v) = j
                    do while (l>=1 .and. idx(subqueue(l))>idx(v))
                        w = subqueue(l)
                        l = l - 1
                        found(w) = .true.
                        component(w) = j
                    end do
                else
                    l = l + 1
                    subqueue(l) = v
                end if
            end if
        end do
    end do

end subroutine strongly_connected_components

!
! threshold_matrix(B, A, weight, m, n)
!
!   set entries of a matrix exeeding a threshold equal to zero
!

subroutine threshold_matrix(B, A, weight, m, n)

    implicit none

    integer, intent(in) :: m, n
    double precision, intent(in) :: A(m, n), weight

    double precision, intent(out) :: B(m, n)

    integer :: i, j

    do j = 1, n
        do i = 1, m
            if (A(i,j) <= weight) then
                B(i,j) = A(i,j)
            else
                B(i,j) = 0.d0
            end if
        end do
    end do

end subroutine threshold_matrix

!
! unique_entries(B, k, A, m, n)
!
!   find k unique, sorted entries of a matrix A
!

subroutine unique_entries(B, l, A, m, n)

    implicit none

    integer, intent(in) :: m, n
    double precision, intent(in) :: A(m, n)

    integer, intent(out) :: l
    double precision, intent(out) :: B(m*n)

    integer :: i, j, k

    k = 1
    B(k) = 0.d0
    do j = 1, n
        do i = 1, m
            if (A(i, j) /= 0.d0) then
                k = k + 1
                B(k) = A(i, j)
            end if
        end do
    end do

    call quicksort(B, 1, k, k)

    l = 1
    do i = 2, k
        if (B(i-1) /= B(i)) then
            l = l + 1
            B(l) = B(i)
        end if
    end do

end subroutine unique_entries

!
! quicksort(x, lo, hi, n)
!
!   sort entries of a vector
!

recursive subroutine quicksort(x, lo, hi, n)

    implicit none

    integer, intent(in) :: lo, hi, n
    double precision, intent(inout) :: x(n)

    integer :: i, j, cutoff=8
    double precision :: pivot, tmp

    ! Use insertion sort for small lists.
    if (hi-lo+1>cutoff) then
        pivot = x((lo+hi)/2)
        i = lo
        j = hi
        do
            do while (x(i)<pivot)
                i = i + 1
            end do
            do while (pivot<x(j))
                j = j - 1
            end do

            if (i>=j) then
                exit
            end if

            tmp = x(i)
            x(i) = x(j)
            x(j) = tmp

            i = i + 1
            j = j - 1
        end do

        if (lo<i-1) then
            call quicksort(x, lo, i-1, n)
        end if
        if (j+1<hi) then
            call quicksort(x, j+1, hi, n)
        end if
    else
        call insertion_sort(x(lo:hi), hi-lo+1)
    end if

end subroutine quicksort

!
! insertion_sort(x, n)
!
!   sort entries of a vector
!

recursive subroutine insertion_sort(x, n)

    implicit none

    integer, intent(in) :: n
    double precision, intent(inout) :: x(n)

    integer :: i, j
    double precision :: tmp

    do i = 2, n
        tmp = x(i)
        j = i - 1
        do while (j>=1 .and. x(j)>tmp)
            x(j+1) = x(j)
            j = j - 1
        end do
        x(j+1) = tmp
    end do

end subroutine insertion_sort
