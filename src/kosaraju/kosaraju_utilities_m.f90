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
recursive subroutine  follow_pass_1(graph,visited,stack,jend,jpath)
type(net), dimension(:), intent(in), pointer :: graph
logical, dimension(:), pointer :: visited !(intents and pointers)
integer, dimension(:), pointer :: stack
integer, intent(inout) :: jend
type(path), pointer :: jpath
integer :: i,kk
	visited(jpath%node)=.true.
  do i=1,2 ! there are two parents
    kk=graph(jpath%node)%jpar(i)
    if(kk==0) cycle
    if(.not.visited(kk)) then
      allocate(jpath%next)
      jpath%next%node=kk
      call follow_pass_1(graph,visited,stack,jend,jpath%next)
      deallocate(jpath%next) ! clean up
      nullify(jpath%next)    ! 
    end if
  end do ! having explored both parents
  jend=jend+1 ! now put on end of stack
	stack(jend)=jpath%node
end subroutine follow_pass_1
!
recursive subroutine follow_pass_2(graph,visited,tour,jpath)
type(net), dimension(:), intent(in), pointer :: graph
logical, dimension(:), pointer :: visited
logical, dimension(:), pointer :: tour
type(path), pointer :: jpath ! is current path
integer :: i,kk
	visited(jpath%node)=.true.
  tour(jpath%node)=.true.
  if(associated(graph(jpath%node)%joff)) then
    do i=1,size(graph(jpath%node)%joff)
      kk=graph(jpath%node)%joff(i)	    
      if(.not.visited(kk)) then
        allocate(jpath%next)
        jpath%next%node=kk
        call follow_pass_2(graph,visited,tour,jpath%next)
        deallocate(jpath%next) ! clean up
        nullify(jpath%next)    ! 
      end if
    end do ! having explored all offspring
  end if
end subroutine follow_pass_2
!
subroutine report_set(graph,tour,iset,store)
	type(net), dimension(:), intent(in), pointer :: graph
	logical, dimension(:), intent(in), pointer :: tour
    integer, intent(inout) :: iset
    integer, dimension(:), pointer :: store
	integer :: i,j
	!character(len=20) :: atrack
    iset=iset+1
    j=0
	!atrack=''
	do i=1,size(tour)
        if(.not.tour(i)) cycle
		j=j+1
        !atrack(j:j)=trim(adjustl(graph(i)%tag))
        store(i)=iset
	end do
	!write(atxt,'(a,i4,a,i8,a)') '... set ',iset,' of length ',j,' with self-ancestry'
  !call report_line([6,14],atxt)
end subroutine report_set
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
function fname(atoken,prog,finfo,fext) result(fnout)
character(len=*), intent(in) :: atoken,prog,finfo,fext
character(len=len_trim(atoken)+len_trim(prog)+len_trim(finfo)+len_trim(fext)+3) :: fnout
  fnout=trim(atoken)//'_'//trim(prog)//'_'//trim(finfo)//'.'//trim(fext)
end function fname
!
subroutine subg_info(ich)
integer, intent(in) :: ich
  write(14,'(a)') '... ... subgraph file has no headers, as record lengths vary'
  write(14,'(a)') '... ==> tag'
  write(14,'(a)') '... ==> the graph node number'
  write(14,'(a)') '... ==> the node number of the sire, or 0'
  write(14,'(a)') '... ==> the node number of the dam, or 0'
  write(14,'(a)') '... ==> npat, the number of offspring in graph as sire'
  write(14,'(a)') '... ==> nmat, the number of offspring in graph as dam'
  write(14,'(a)') '... ==> joff, list of nodes of paternal offspring, then maternal offspring'  
end subroutine subg_info 
!
subroutine sets_info(ich)
integer, intent(in) :: ich
  write(14,'(a)') '... ... uses subgraph node numbers'
  write(14,'(a)') '... ==> the graph node number'
  write(14,'(a)') '... ==> the node number of the sire, or 0'
  write(14,'(a)') '... ==> the node number of the dam, or 0'
  write(14,'(a)') '... ==> tag' 
end subroutine sets_info
!
end module utilities