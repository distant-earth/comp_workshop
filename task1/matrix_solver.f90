module matrix_solver

contains

subroutine tridiagonal(N, A_matrix, b_vector, res_vector)
implicit none
integer :: N, i
real(8) :: A_matrix(0:N-1, 0:N-1), b_vector(0:N-1), res_vector(0:N-1)
real(8) :: alpha_vector(0:N-2), beta_vector(0:N-1), gamma_vector(1:N-1)
real(8) :: x_vector(0:N-2), y_vector(0:N-2)
real(8) :: det, det_1st, det_2nd
character(80) :: frmt

! Extract the non-zero diagonals:
do i = 0, N-2
	beta_vector(i) = A_matrix(i, i) ! the main diagonal
	alpha_vector(i) = A_matrix(i, i+1) ! the upper subdiagonal
	gamma_vector(i+1) = A_matrix(i+1, i) ! the lower subdiagonal
enddo
beta_vector(N-1) = A_matrix(N-1, N-1) ! one remaining element

! Find the determinant of a tridiagonal matrix using a recurrence relation
! det_N= beta(N-1) * det_N-1 - alpha(N-2) * gamma(N-1) * det_N-2
! (in terms of three diagonals extracted above), where det_N is the determinant 
! of the NxN matrix, det_1 = beta(0), det_0 = 1 and det_-1 = 0:
det_2nd = 1.d0
det_1st = beta_vector(0)
do i = 2, N
	det = beta_vector(i-1) * det_1st - alpha_vector(i-2) * gamma_vector(i-1) * det_2nd
	det_2nd = det_1st
	det_1st = det
enddo

! Prepare the output file:
open(3, file = 'RESULT')

! When the matrix A is degenerate, the system has either no or infinitely many solutions.
! In order to use Rouché–Capelli (Kronecker–Capelli) theorem, one has to find ranks of both
! matrix A and its augmented matrix A|b. This task is rather complicated and may require
! Gaussian elimination to turn the matrices into the row echelon form. It provides another
! stand-alone method to solve the system. Therefore, I decided not to implement it here.

! Check if the determinant of a real-valued matrix A is zero. There, epsilon(X) returns 
! the smallest number E of the same kind as X such that 1 + E > 1:
if (abs(det) <= epsilon(det)) then
	! File output:
	write(3,*) 'Warning: single solution is not available. Method does not work.'
	! Screen output:
	write(*,*) 'Warning: single solution is not available. Method does not work.'
else
	! Use Thomas algorithm to solve the tridiagonal matrix.
	! Compute the coefficients x(i), y(i) of the reccurence relation res(i+1) = x(i) * res(i) + y(i):
	x_vector(N-2) = - gamma_vector(N-1) / beta_vector(N-1)
	y_vector(N-2) = b_vector(N-1) / beta_vector(N-1)
	do i = N-2, 1, -1
		x_vector(i-1) = - gamma_vector(i) / (alpha_vector(i) * x_vector(i) + beta_vector(i))
		y_vector(i-1) = (b_vector(i) - alpha_vector(i) * y_vector(i)) / (alpha_vector(i) * x_vector(i) + beta_vector(i))
	enddo
! Compute the solution res(i):
	res_vector(0) = (b_vector(0) - alpha_vector(0) * y_vector(0)) / (beta_vector(0) + alpha_vector(0) * x_vector(0))
	do i = 0, N-2
		res_vector(i+1) = x_vector(i) * res_vector(i) + y_vector(i)
	enddo
	! Print the solution into the 'RESULT' file using formatted output.
	! Here, 'frmt' is a string with required format specification for a vector of length N.
	write(frmt, '(a,i5,a)') '(',N,'f10.3)'
	write(3, frmt) res_vector	
	! Screen output:
	write(*, frmt) res_vector	
endif

end subroutine tridiagonal

end module matrix_solver
