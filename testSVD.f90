program testSVD

use omp_lib

implicit none
!integer, parameter :: M=5, N=3, P=2  ! Define array size here.
integer,parameter :: M=17576, N=3403, P=6

real(8),allocatable,dimension(:,:):: V,T,h_ij      ! Non-linear combination matrix
real(8) :: lambda = 0.1d0         ! lambda, damping factor
real(8) :: random
real(8) :: tic, toc
integer :: i,j,nthread
logical :: printval=.false.
logical :: randomV =.true.

character(3) :: solver = 'lud' 


!$ call omp_set_num_threads(2)

! ASSERT (M >= N) BEFORE ALLOCATION:
if (M.lt.N) then
   print*, "Error: M < N ...  Need square or over-determined system"
   stop
else
   allocate (V(M,N), T(M,P), h_ij(N,P))
end if

!
! INITIALIZE ARRAY:
!
if (randomV) then
! PSEUDO-RANDOM MATRIX:
call cpu_time(tic)              !### START TIMING
   call init_random_seed()      
   !$omp parallel do
   do j = 1,N
   do i = 1,M
      call random_number(random)
      V(i,j) =  random
   end do
   end do
   do j=1,P
   do i=1,M
      call random_number(random)
      T(i,j) = random              !Copies the first column over to all six columns
   end do
   end do
   call cpu_time(toc)           !### STOP TIMING
   print*, 'Elapsed time:', toc-tic,'s'
else 
! USER-DEFINED MATRIX:
   do j = 1,N
   do i = 1,M
      V(i,j) =  i**2 + 2*j
   end do
   end do
   T = reshape((/1, 5, 10, 2, 3, 1, 10, 5, 2, 3 /), shape(T))
end if

! PRINT V,T:
print*,'V:'
call printmatrix(V,size(V,dim=1))
print*,'T:'
call printmatrix(T,size(T,dim=1))

if (M.gt.8.or.N.gt.8)  printval=.false. !suppress printing large matrices

!
! SOLVER:
!
print*,'Using solution method:',solver
select case (solver)
case('svd')
   call SVD(V,T,h_ij,printval)
case('lud')
   call LU(V,T,h_ij,lambda)
case default
   print*,'Select svd/lud'
   stop
end select

print*,'V:'
call printmatrix(V,size(V,dim=1))
print*,'T:'
call printmatrix(T,size(T,dim=1))
print*,'h_ij:'
call printmatrix(h_ij,size(h_ij,dim=1))


contains

!****************************************************************
!                              SVD                              !
!****************************************************************

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

real(8) :: tic,toc
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
call cpu_time(tic)
call dgesvd('All','All',M,N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, info)
call cpu_time(toc)
!$omp end critical
!$omp end parallel
print*,'Elapsed time:',toc-tic

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
forall(i=1:N) D(i,i) = S(i) / (S(i)**2 + lambda**2)
call cpu_time(tic)
Vinv = matmul(matmul(transpose(VT),D),transpose(U))
call cpu_time(toc)
print*,"Matmul() Elapsed time:",toc-tic

if(printval) then
print*,'D'                     ! D MATRIX - DIAG MATRIX
call printmatrix(D,size(D,dim=1))


print*,'Vinv'                     ! Vinv MATRIX - DIAG MATRIX
call printmatrix(Vinv,size(Vinv,dim=1))
end if

h_ij = matmul(Vinv,T)
return
end subroutine SVD


!****************************************************************
!                                LU                             !
!****************************************************************
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
  real(8), dimension(:,:),allocatable :: A,VT ! EYE - IDENTITY MATRIX
  real(8), dimension(:,:),allocatable :: b
  integer, dimension(:), allocatable  :: ipiv
  integer                   :: info
  real(8) :: tic, toc
  integer :: i
 
  !DGEMM ARGUMENTS:
  real(8):: alpha=1.d0, beta=0.d0

 allocate (A(N,N), b(N,P), ipiv(N), VT(N,M))

 VT = transpose(V)
 call cpu_time(tic)
 call dgemm('T','N',N,N,M,alpha,V,M,V,M,beta,A,N)
 call cpu_time(toc)
 print*,'Elapsed time', toc-tic

forall(i=1:N) A(i,i) = A(i,i) + lambda 

!A
 print*,'A:'
 call printmatrix(A,size(A,dim=1)) 
!b
! b = matmul(VT,T_ij)
 call dgemm('T','N',N,P,M,alpha,V,M,T_ij,M,beta,b,N)
 print*,'b:'
 call printmatrix(b,size(b,dim=1))
 
call cpu_time(tic)
!$omp parallel 
call DGESV(N, nrhs, A, LDA, ipiv, b, LDB, info)
!$omp end parallel 
call cpu_time(toc)
print*,'Elapsed time:', toc-tic

h_ij = b

deallocate (A,b,ipiv,VT)
return
end subroutine LU


!****************************************************************
!    SUBROUTINE INIT_RANDOM_SEED: RANDOM NUMBER GENERATION      !
!****************************************************************

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


!****************************************************************
!                     SUBROUTINE PRINTMATRIX                    !            
!****************************************************************
subroutine printmatrix(A,LDA)
implicit none

! ARGUMENTS:
real(8), dimension(:,:) :: A
integer :: LDA,width
integer :: i

if (LDA.gt.8) then
   LDA=8
   width=8
end if

if (size(A,dim=2).lt.8) width = size(A,dim=2)

do i=1,LDA
   print*, A(i,1:width)
end do

return
end subroutine printmatrix

end program testSVD


