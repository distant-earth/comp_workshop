program main
use matrix_solver

implicit none
real(8), allocatable :: A_matrix(:, :), b_vector(:), res_vector(:)
real(8) :: test
integer :: N, i, A_error, b_error, res_error

! This program solves matrix equation A * u = b, where A is supposed to be 
! a tridiagonal square matrix of a given dimension N.

! Get the dimension of the matrix while reading input file:
N = 0
open(1, file = 'INPUT_matrix')
do
    read(1, *, end=10) test 
    ! Real variable 'test' lets us abandon empty lines.
    N = N + 1
enddo
10 close (1)

! 'INPUT_matrix' contains the matrix A itself.
open(1, file = 'INPUT_matrix')
! Allocate 2-dimensional array of size (N, N) for the matrix A:
allocate (A_matrix(0:N-1, 0:N-1), stat=A_error)
if (A_error.ne.0) stop 'Matrix array: allocation failed!'
! Read matrix A from the file row by row:
do i = 0, N-1
	read(1, *) A_matrix(i, :)
enddo

! 'INPUT_right' contains the right parts of the equations.
open(2, file = 'INPUT_right')
! Allocate 1-dimensional array of size N for the vector b:
allocate (b_vector(0:N-1), stat=b_error)
if (b_error.ne.0) stop 'Right parts array: allocation failed!'
! Read right parts vector from the file:
read(2, *) b_vector

! Allocate 1-dimensional array of size N for the result vector:
allocate (res_vector(0:N-1), stat=res_error)
if (res_error.ne.0) stop 'Result array: allocation failed!'

! Call a subroutine for tridiagonal matrix solving:
call tridiagonal(N, A_matrix, b_vector, res_vector)

end
