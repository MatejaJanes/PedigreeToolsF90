module parameters
implicit none
integer, parameter :: alen=12
integer, parameter :: dp=kind(0.0d0)
type info ! for reading in the mating list
  integer :: ji,js,jd,jx
  real(kind=dp) :: f
  character(len=alen) :: tag
  real(kind=dp) :: y,y1,x 
  logical :: actv
end type info
type(info), dimension(:), pointer :: my_ped => null()
type(info), pointer :: z => null()
type grec   ! for information accumulated in binary search trees
    type(grec), pointer :: go_l => null(), go_r => null()
    character(len=alen) :: tag
    integer :: ji
end type grec
type(grec), pointer :: ped_tree => null()
integer, dimension(:), allocatable :: jord
integer :: ierr, i,iact,j,maxp
real(kind=dp), dimension(:), allocatable :: fdum
real(kind=dp) :: gcoa,t1,t2,ff
character(len=alen) :: adum
character(len=10) :: ctime,cdate,atoken
character(len=:), allocatable :: ifile,ofile,afile
character(len=120) :: atxt
logical :: iask,jask,iget 
end module parameters