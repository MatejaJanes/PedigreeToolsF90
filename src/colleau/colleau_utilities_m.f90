module utilities
use parameters
implicit none
    contains
!
subroutine report_line(ich,atxt)
integer, dimension(:), intent(in) :: ich
character(len=*), intent(in) :: atxt
integer :: j
	do j=1,size(ich)
	    write(ich(j),'(a)') trim(adjustl(atxt))
    end do
end subroutine report_line
!
subroutine report_time(ich,ctxt)
integer, dimension(:), intent(in) :: ich
character(len=*) :: ctxt
character(len=10) :: cdate, ctime
character(len=140) :: atxt 
	call date_and_time(cdate,ctime)
	write(atxt,'(a)') cdate//ctime//' '//ctxt
	call report_line(ich,atxt)
end subroutine report_time
!
recursive subroutine add_ped_tree(list,adum,idum)
character(len=*), intent(in) :: adum ! the new identity
integer, intent(in) :: idum
type (grec), pointer :: list
	if(.not.associated(list)) then ! new ai
		allocate (list)
		list=grec(null(),null(),trim(adjustl(adum)),idum)
	else if(lgt(adum,list%tag)) then
		call add_ped_tree(list%go_r,adum,idum)
	else if(llt(adum,list%tag)) then
		call add_ped_tree(list%go_l,adum,idum)
	end if
end subroutine add_ped_tree
!
!
recursive subroutine find_in_tree(list,my_info,adum,iget)
type(grec), pointer :: list
type(info), dimension(:), pointer :: my_info
character(len=alen), intent(in) :: adum
logical, intent(inout) :: iget
    if(associated(list)) then
        if(lgt(adum,list%tag)) then
            call find_in_tree(list%go_r,my_info,adum,iget)
	    else if(llt(adum,list%tag)) then
		    call find_in_tree(list%go_l,my_info,adum,iget)
        else
            my_info(list%ji)%actv=.true.
            iget=.true.
        end if
    end if
end subroutine find_in_tree
!
recursive subroutine delete_tree(list)
type(grec), pointer :: list
	if(associated(list%go_l)) call delete_tree(list%go_l)
	if(associated(list%go_r)) call delete_tree(list%go_r)
	nullify(list%go_l); nullify(list%go_r); deallocate(list)
end subroutine delete_tree
!
function shuffle(n) ! Fisher Yates Durstenfeld
integer, intent(in) :: n
integer :: item, i, j
integer, dimension(n) :: shuffle
real :: u
	shuffle = (/ (i,i=1,n) /)
	swapping: do i = n, 2, -1
		call random_number(u)
		j = floor(u*i) + 1
		item = shuffle(j)
		shuffle(j) = shuffle(i)
		shuffle(i) = item
	end do swapping
end function shuffle
!
function fname(atoken,prog,finfo,fext) result(fnout)
character(len=*), intent(in) :: atoken,prog,finfo,fext
character(len=len_trim(atoken)+len_trim(prog)+len_trim(finfo)+len_trim(fext)+3) :: fnout
  fnout=trim(atoken)//'_'//trim(prog)//'_'//trim(finfo)//'.'//trim(fext)
end function fname
!
end module utilities    