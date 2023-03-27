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
	write(atxt,'(a)') cdate//ctime//' '//trim(adjustl(ctxt))
	call report_line([ich],atxt)
end subroutine report_time
! 
function fname(atoken,prog,finfo,fext) result(fnout)
character(len=*), intent(in) :: atoken,prog,finfo,fext
character(len=len_trim(atoken)+len_trim(prog)+len_trim(finfo)+len_trim(fext)+3) :: fnout
  fnout=trim(atoken)//'_'//trim(prog)//'_'//trim(finfo)//'.'//trim(fext)
end function fname
!
end module utilities    