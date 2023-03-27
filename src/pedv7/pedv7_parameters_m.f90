module parameters
implicit none
integer, parameter :: alen=10
integer, parameter :: dp=kind(0.0d0)
type b_cert ! for reading in the mating list
  character(len=alen) :: ai, as, ad
end type b_cert
type(b_cert), dimension(:), allocatable :: my_rec
type grec   ! for information accumulated in binary search trees
  type(grec), pointer :: go_l => null(), go_r => null()
  character(len=alen) :: tag
  integer :: ti,ts,td,irec
  logical :: serr,derr
end type grec
type(grec), pointer :: sd_tree => null()
type(grec), pointer :: r0_tree => null()
type info   ! for information extracted from trees
  character(len=alen) :: tag
  integer :: ji,js,jd,jx ! jx is usage as sire and dam
  integer :: ns,nd ! number of offspring as sire and dam
	integer :: isw,nyi,irec ! record number, sweep number, and the number allocated in consistent pedigree
  logical :: serr,derr ! flags for discrepancies in sire, and dam across records
end type info
type(info), dimension(:), pointer :: my_info => null(), zinfo => null()
type(info), pointer :: z
type prec   ! for an ordered consistent pedigree
  integer :: ji,js,jd,jx
  character(len=alen) :: tag
end type prec
type(prec), dimension(:), allocatable :: ped_list
integer, dimension(:), allocatable :: jord,sweep,err1,err2,err3,err4,err4x,err5,emax
integer, dimension(:), allocatable, target :: noff
integer, dimension(:), pointer :: zoff => null()
integer :: addi,flag1,flag2,flag3,flag4,flag5
integer :: i,ierr,ii,istep,j,jj,jdam,jlast,jpath,jsir,iskip
integer :: klast,kmax,kpar,kodd,ksweep,maxp,ndam,npar,nr0,nsir
character(len=alen) :: adum
character(len=1) :: acha
character(len=2), dimension(2) :: nkn
character(len=10) :: ctime,cdate,atoken
character(len=50) :: ifile,efile
character(len=120) :: atxt
logical :: iask 
real :: t1,t2
end module parameters