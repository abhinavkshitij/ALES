program testSVD

use omp_lib

implicit none
integer, parameter :: M=5, N=3, P=2  ! Define array size here.

real(8),allocatable,dimension(:,:):: V,T,h_ij      ! Non-linear combination matrix
real(8) :: lam = 0.1d0         ! lambda, damping factor
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
if (randomV) then
! CREATE A PSEUDO-RANDOM MATRIX:
call cpu_time(tic)
   call init_random_seed()
!$omp parallel do
   do j = 1,N
   do i = 1,M
      call random_number(random)
      V(i,j) =  random
   end do
   call random_number(random)
   T(j,:) = random              !Copies the first column over to all six columns
   end do
call cpu_time(toc)
print*, 'Elapsed time:', toc-tic,'s'
else

! CREATE A USER-DEFINED MATRIX:

   do j = 1,N
   do i = 1,M
      V(i,j) =  i**2 + 2*j
   end do
   end do
   T = reshape((/1, 5, 10, 2, 3, 1, 10, 5, 2, 3 /), shape(T))
end if



print*,'V:'
call printmatrix(V,size(V,dim=1))
print*,'T:'
call printmatrix(T,size(T,dim=1))

! SVD DECOMPOSITION:
if (M.gt.8.or.N.gt.8)  printval=.false. !suppress printing large matrices
call SVD(V,T,h_ij,printval)


! PRINT h_ij:
   print*,'h_ij:'
if (printval) then
   print*, h_ij
else
      call printmatrix(h_ij,size(h_ij,dim=1))
end if

contains



subroutine SVD(A,T,h_ij,printval)

! ARGUMENTS:
real(8),dimension(:,:),intent(in)  :: A
real(8),dimension(:,:),intent(in)  :: T
real(8),dimension(:,:),intent(out) :: h_ij

! DGESVD ARGUMENTS:                                                                                                
integer, parameter :: LDA = M
integer, parameter :: LDU = M
integer, parameter :: LDVT = N
integer :: LWMAX
integer :: i,j, info, LWORK

real(8),dimension(:,:),allocatable :: U
real(8),dimension(:,:),allocatable :: VT 
real(8),dimension(:,:),allocatable :: D, Vinv
real(8),dimension(:), allocatable :: S, work

logical :: printval
character(1) :: N_char 


LWMAX = M*N
allocate (U(LDU,M), VT(N,N), S(N), D(N,M), Vinv(N,M), work(LWMAX))

if (printval) then
print*,'A'                      ! A MATRIX
call printmatrix(A,size(A,dim=1))
end if


! CALL SVD ROUTINE:
LWORK = -1
call dgesvd('All','All', M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info)
LWORK = min( LWMAX, int(WORK(1)) ) 
!$omp parallel
!$omp critical
call dgesvd('All','All',M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info)
!$omp end critical
!$omp end parallel


! Convergence check:
if (info.gt.0) then
   print*, 'Failed to converge'
   stop
end if


if (printval) then
print*,'A'                      ! A MATRIX
call printmatrix(A,size(A,dim=1))
print*,'U'                      ! U MATRIX
call printmatrix(U,size(U,dim=1))
print*,'S'                      ! S VECTOR
print*, S
print*,'VT'                     ! VT MATRIX
call printmatrix(VT,size(VT,dim=1))
end if


! COMPUTE PSEUDOINVERSE:
D = 0.d0
forall(i=1:N) D(i,i) = S(i) / (S(i)**2 + lam**2)
Vinv = matmul(matmul(transpose(VT),D),transpose(U))

if(printval) then
print*,'D'                     ! D MATRIX - DIAG MATRIX
call printmatrix(D,size(D,dim=1))


print*,'Vinv'                     ! Vinv MATRIX - DIAG MATRIX
call printmatrix(Vinv,size(Vinv,dim=1))
end if

h_ij = matmul(Vinv,T)
return
end subroutine SVD


! SUBROUTINE INIT_RANDOM_SEED: RANDOM NUMBER GENERATION 
subroutine init_random_seed()

integer, dimension(:), allocatable :: seed
integer :: i, n_seed, clock

call RANDOM_seed(size = n_seed) !Generate random seeds                                                                   
allocate(seed(n_seed))
call SYSTEM_clock(COUNT=clock)
seed = clock + 37 * (/ (i-1, i=1,n_seed) /)
call RANDOM_seed(PUT = seed)
deallocate(seed)
return
end subroutine init_random_seed


! SUBROUTINE PRINTMATRIX:
subroutine printmatrix(A,LDA)
implicit none

! ARGUMENTS:
real(8), dimension(:,:) :: A
integer :: LDA
integer :: i,j

if (LDA.gt.8) LDA=8
do i=1,LDA
   print*, A(i,1:4)
end do

return
end subroutine printmatrix

end program testSVD


