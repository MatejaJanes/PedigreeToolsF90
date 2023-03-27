module utilities
use parameters
implicit none
!
contains
!
subroutine report_line(ich,ctxt)
integer, dimension(:), intent(in) :: ich
character(len=*), intent(in) :: ctxt
integer :: j
	do j=1,size(ich)
   write(ich(j),'(a)') trim(adjustl(ctxt))
  end do
end subroutine report_line
!
subroutine report_time(ich,ctxt)
integer, dimension(:), intent(in) :: ich
character(len=*) :: ctxt
character(len=10) :: cdate, ctime
character(len=(len(ctxt)+21)) :: dtxt 
	call date_and_time(cdate,ctime)
	write(dtxt,'(a)') cdate//ctime//' '//ctxt
	call report_line(ich,dtxt)
end subroutine report_time
!
recursive subroutine add_step(graph,removed,blocked,jpath,istart,icirc,ich,iset) 
type(net), dimension(:), intent(in), pointer :: graph
logical, dimension(:), pointer :: removed
logical, dimension(:), pointer :: blocked
type(path), pointer :: jpath
integer, intent(in) :: istart,ich,iset
integer, intent(inout) :: icirc
integer :: i,jcode,joff,jpar ! joff and jpar are code simplifiers
	blocked(jpath%node)=.true.
  do i=1,2
    if(graph(jpath%node)%parent(i)==0) cycle ! path does not exist
    if(removed(graph(jpath%node)%parent(i))) cycle ! path no longer exists
    if(graph(jpath%node)%parent(i)==istart) then ! first check on next step
      jcode=2*jpath%code+(i-1)
      call circuit(graph,istart,jcode,icirc,cfreq,ich,iset)
      cycle
    end if
		if(blocked(graph(jpath%node)%parent(i))) then ! second check on next step
			jpar=graph(jpath%node)%parent(i) ! parent id
			joff=graph(jpath%node)%ioff(i)   ! position in jpar's family
			graph(jpar)%child(joff)%xji=.true.
			cycle
    end if
    jcode=2*jpath%code+(i-1)
    allocate(jpath%next)
    jpath%next=path(graph(jpath%node)%parent(i),jcode,null())
    call add_step(graph,removed,blocked,jpath%next,istart,icirc,ich,iset)
		!print *, associated(jpath%next)
		!if(associated(jpath%next)) deallocate(jpath%next)
    nullify(jpath%next)    !
  end do
  call unblock(graph,jpath%node,blocked)
  deallocate(jpath) ! itself
end subroutine add_step
!       
recursive subroutine unblock(graph,jnode,blocked)
type(net), dimension(:), intent(in), pointer :: graph
integer, intent(in) :: jnode
logical, dimension(:), pointer :: blocked
integer :: i,j ! j is a code simplifier
do i=1,size(graph(jnode)%child) ! it must have 1+
    if(graph(jnode)%child(i)%xji) then
	    j=graph(jnode)%child(i)%ji
        if(blocked(j)) call unblock(graph,j,blocked)
        graph(jnode)%child(i)%xji=.false.
    end if
end do
blocked(jnode)=.false.
end subroutine unblock
!
subroutine circuit(graph,istart,jcode,icirc,cfreq,ich,iset)
type(net), dimension(:), intent(in), pointer :: graph
integer, intent(in) :: istart,jcode,ich,iset
integer, intent(inout) :: icirc
integer, dimension(:), intent(inout) :: cfreq
integer :: i,j
character(len=20) :: acirc
    icirc=icirc+1
    acirc=decode_path(jcode)
    !write(14,'(3i8,2a)')  icirc,jcode,istart,graph(istart)%tag,acirc
    write(ich,'(*(g0,:,","))')  iset,icirc,istart,jcode,trim(acirc),trim(graph(istart)%tag)
    j=istart
    do i=1,len_trim(acirc)
        if(acirc(i:i)=='p') j=graph(j)%parent(1)
        if(acirc(i:i)=='m') j=graph(j)%parent(2)
        cfreq(j)=cfreq(j)+1
    end do
end subroutine circuit
!    
character(len=10) function decode_path(i) result(a)
integer, intent(in) :: i
integer :: j,k,q
k=bit_size(i)-leadz(i)-1
a=''
do j=k-1,0,-1
    q=k-j
    if(btest(i,j)) then; a(q:q)='m'; else; a(q:q)='p'; end if
end do    
end function decode_path
!
recursive subroutine check_path(jpath)
type(path), pointer :: jpath
    write(14,'(2i8,l2)') jpath%node,jpath%code,associated(jpath%next)
    if(associated(jpath%next)) call check_path(jpath%next)
end subroutine check_path
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
subroutine swap(i,j)
integer, intent(inout) :: i,j
integer :: kdum
kdum=i; i=j; j=kdum
end subroutine swap
!
function fname(atoken,prog,finfo,fext) result(fnout)
character(len=*), intent(in) :: atoken,prog,finfo,fext
character(len=len_trim(atoken)+len_trim(prog)+len_trim(finfo)+len_trim(fext)+3) :: fnout
  fnout=trim(atoken)//'_'//trim(prog)//'_'//trim(finfo)//'.'//trim(fext)
end function fname
!
subroutine subg_info(ich)
integer, intent(in) :: ich
  write(14,'(a)') '... ... renumbered subgraph file has no headers, as record lengths vary'
  write(14,'(a)') '... ==> Kosaraju set number'
  write(14,'(a)') '... ==> tag of individual'
  write(14,'(a)') '... ==> the graph node number'
  write(14,'(a)') '... ==> the node number of the sire, or 0'
  write(14,'(a)') '... ==> the node number of the dam, or 0'
  write(14,'(a)') '... ==> total offspring'
  write(14,'(a)') '... ==> total labelled as sire offspring'
  write(14,'(a)') '... ==> total labelled as dam offspring'  
  write(14,'(a)') '... ==> node numbers of offspring' 
end subroutine subg_info
!
subroutine freq_info(ich)
integer, intent(in) :: ich
  write(14,'(a)') '... ... summary of distinct circuits found, with header line'
  write(14,'(a)') '... ==> iset, Kosaraju set number'
  write(14,'(a)') '... ==> tag, tag of individual'
  write(14,'(a)') '... ==> ctot, the total number of circuits in which tag appears'
  write(14,'(a)') '... ==> cfreq, the fraction of all circuits in which tag appears'
end subroutine freq_info
!
subroutine circ_info(ich)
integer, intent(in) :: ich
  write(14,'(a)') '... ... circuits found, with header line'
  write(14,'(a)') '... ==> iset, Kosaraju set number'
  write(14,'(a)') '... ==> icirc, index of circuit for the set'
  write(14,'(a)') '... ==> istart, the starting node for circuit description' 
  write(14,'(a)') '... ==> jcode, integer code defining circuit'
  write(14,'(a)') '... ==> acirc, geneological pathway for circuit'
  write(14,'(a)') '... ==> tag, tag of starting node'
end subroutine circ_info
!
end module utilities