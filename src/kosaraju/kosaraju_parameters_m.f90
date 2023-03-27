module parameters
implicit none
integer, parameter :: alen=12
type info
	character(len=alen) ::  tag
	integer :: ji
  integer, dimension(2) :: jpar
  integer :: jx,ns,nd,irec,isw
end type info
type(info), dimension(:), allocatable :: v7out
!
type net
	character(len=alen) :: tag
	integer :: ji,npat,nmat
	integer, dimension(2) :: jpar
	integer, dimension(:), pointer :: joff => null()
end type net
type(net), dimension(:), pointer :: graph => null()
type(net), pointer :: z => null()
!
type path
  integer :: node
  type(path), pointer :: next => null()
end type path
type(path), pointer :: my_path
logical, dimension(:), pointer :: visited, tour
logical, dimension(:), allocatable :: in_graph,not_s,not_d
integer, dimension(:), pointer :: stack,store	
integer, dimension(:), allocatable :: ivals,s_inv
integer :: i,ierr,ii,iset,j,jend,kk,lset,maxp,nd,ns,nsub
logical :: iask
character(len=5) :: acha
character(len=10) :: atoken
character(len=50) :: ifile
character(len=120) :: atxt
!
end module parameters
