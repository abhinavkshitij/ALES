program testLU

use omp_lib

implicit none
integer, parameter :: M=5, N=3, P=2  ! Define array size here.

real(8),allocatable,dimension(:,:):: V,T,h_ij      ! Non-linear combination matrix
real(8) :: lambda = 0.1d0         ! lambda, damping factor
real(8) :: random
real(8) :: tic, toc
integer :: i,j,nthread
logical :: printval=.false.
logical :: randomV =.false.


allocate (V(M,N), T(M,P), h_ij(N,P))
V=0.;T=0.;h_ij=0.
! ASSERT (M>N)
if (M.lt.N) then
   print*, "SVD computation needs M >= N "
   stop
end if

!$ call omp_set_num_threads(8)
! INITIALIZE ARRAY:
   do j = 1,N
   do i = 1,M
      V(i,j) =  i**2 + 2*j
   end do
   end do
   T = reshape((/1, 5, 10, 2, 3, 1, 10, 5, 2, 3 /), shape(T))



print*,'V:'
call printmatrix(V,size(V,dim=1))
print*,'T:'
call printmatrix(T,size(T,dim=1))


! SVD DECOMPOSITION:
if (M.gt.8.or.N.gt.8)  printval=.false. !suppress printing large matrices

call LU(V,T,h_ij,lambda)

print*,'V:'
call printmatrix(V,size(V,dim=1))
print*,'T:'
call printmatrix(T,size(T,dim=1))

! PRINT h_ij:
print*,'h_ij:'
call printmatrix(h_ij,size(h_ij,dim=1))





contains



subroutine LU(V, T_ij, h_ij, lambda)
 implicit none
 
 ! dgesv destroys the original matrix. So V is copied into 'a' matrix. This is the LU product matrix.
 ! dgesv returns x vector. First it copies the values from 'T_ij' vector and then computes 'h_ij'.
 ! 'h_ij' is the h_ij vector.
 
 ! ARGUMENTS:
 real(8), dimension(:,:), intent(in) :: V 
 real(8), dimension(:,:), intent(in) :: T_ij
 real(8), dimension(:,:), intent(out) :: h_ij
 real(8), intent(in) :: lambda
 
 ! DGESV ARGUMENTS:
 integer, parameter        :: LDA = N, LDB = N, nrhs = P
 real(8), dimension(:,:),allocatable :: A,VT, eye ! EYE - IDENTITY MATRIX
 real(8), dimension(:,:),allocatable   :: b
 integer, dimension(:), allocatable     :: ipiv
 integer                   :: info
 
 integer :: i
 
 allocate (A(LDA,N), eye(LDA,N), b(LDB,nrhs),ipiv(N),VT(N,LDA))
 A=0.d0;eye=0.d0;b=0.d0

 forall(i = 1:N) eye(i,i) = 1.d0 ! Identity matrix
 ! Use the SAVE attribute or something to avoid repeated construction.
 
 A = matmul(transpose(V),V) +  (lambda * eye)
 print*,'A:'
 call printmatrix(A,size(A,dim=1)) 
 VT = transpose(V)
 b = matmul(VT,T_ij) 
 print*,'b:'
 call printmatrix(b,size(b,dim=1))
 
 call DGESV(N, nrhs, A, LDA, ipiv, b, LDB, info)
 h_ij = b
 
 deallocate (A,eye,b,ipiv,VT)
 return
 end subroutine LU





! SUBROUTINE PRINTMATRIX:
subroutine printmatrix(A,LDA)
implicit none

! ARGUMENTS:
real(8), dimension(:,:) :: A
integer :: LDA
integer :: i,j

if (LDA.gt.8) LDA=8
do i=1,LDA
   print*, A(i,1:size(A,dim=2))
end do

return
end subroutine printmatrix

end program testLU


