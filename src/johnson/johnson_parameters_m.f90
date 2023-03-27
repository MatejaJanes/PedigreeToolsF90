module parameters
implicit none
integer, parameter :: alen=12
type inset
  integer :: jset,ji,js,jd
  character(len=alen) :: tag
end type inset
type(inset), dimension(:), allocatable, target :: my_sets
type(inset), pointer :: zin
type sprog 
	integer :: ji
	logical :: xji
end type sprog
type net
  character(len=alen) :: tag
  integer :: ji,noff
  integer, dimension(2) :: parent
  integer, dimension(2) :: ioff
  type(sprog), dimension(:), pointer :: child
end type net
type(net), dimension(:), pointer :: graph
type(net), pointer :: znet
type path
  integer :: node,code
  type(path), pointer :: next=>null()
end type path
type(path), pointer :: my_path
logical, dimension(:), pointer :: blocked,removed
logical :: iask
integer, dimension(:), allocatable :: jinv,moff,cfreq,iperm
integer :: maxn,maxp,nset,ssize,istart,icirc,jlast,i,ierr,ii,jdum,ip,iset 
real :: u
character(len=1) :: acha
character(len=10) :: atoken
character(len=50) :: ifile
character(len=120) :: atxt,adum1,adum2
end module parameters