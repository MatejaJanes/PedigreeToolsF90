module parameters
implicit none
integer, parameter :: alen=12
integer, parameter :: dp=kind(0.0d0)
type info ! for reading in the mating list
    integer :: ji,js,jd,jx
    character(len=alen) :: tag
    real(kind=dp) :: f
end type info
type(info), dimension(:), allocatable, target :: my_ped
type(info), pointer :: z => null()
integer :: ierr,maxp,i,np
character(len=5) :: adum
character(len=10) :: atoken
character(len=:), allocatable :: ifile,ofile
character(len=120) :: atxt
logical :: iask
real :: t1,t2
end module parameters