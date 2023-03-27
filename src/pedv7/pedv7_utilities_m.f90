module utilities
use parameters
implicit none
contains
!    
subroutine compare_sd(recp,newp,errp)
integer, intent(in) :: newp
integer, intent(inout) :: recp
logical, intent(out) :: errp
  if( min(recp,newp)==0 .and. max(recp,newp)>0 ) then
    recp=max(recp,newp)
  else if(.not.(recp==newp)) then
    errp=.true.
    recp=0
  end if
end subroutine compare_sd
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
recursive subroutine add_sd_tree(list,apar,last)
character(len=*), intent(in), target :: apar ! the new identity
integer, intent(inout) :: last ! pedigree length to date
type (grec), pointer :: list
  if(.not.associated(list)) then ! new ai
    last=last+1
    allocate (list)
    list=grec(null(),null(),apar,last,0,0,0,.false.,.false.)
  else if(lgt(apar,list%tag)) then
    call add_sd_tree(list%go_r,apar,last)
  else if(llt(apar,list%tag)) then
    call add_sd_tree(list%go_l,apar,last)
  end if
end subroutine add_sd_tree
!
recursive subroutine update_sd_tree(list,aid,js,jd,addi)
character(len=*), intent(in), target :: aid
integer, intent(in) :: js,jd ! parents
integer, intent(out) :: addi 
type (grec), pointer :: list
  if(aid==list%tag) then
    list%irec=list%irec+1
    if(list%irec==1) then ! first entry as offspring
      list%ts=js
      list%td=jd
    else ! compare information
      call compare_sd(list%ts,js,list%serr)
      call compare_sd(list%td,jd,list%derr)
    end if
		addi=list%ti
	else if(llt(aid,list%tag)) then
		if(associated(list%go_l)) then; call update_sd_tree(list%go_l,aid,js,jd,addi); else; addi=-1; end if
	else if(lgt(aid,list%tag)) then
		if(associated(list%go_r)) then; call update_sd_tree(list%go_r,aid,js,jd,addi); else; addi=-1; end if
	end if
end subroutine update_sd_tree  
!
recursive subroutine add_r0_tree(list,apar,last,js,jd)
character(len=*), intent(in), target :: apar ! the new identity
integer, intent(in) :: js,jd
integer, intent(inout) :: last ! pedigree length to date
type (grec), pointer :: list
	if(.not.associated(list)) then ! new ai
		last=last+1
		allocate (list)
		list=grec(null(),null(),apar,last,js,jd,1,.false.,.false.)
	else if(llt(apar,list%tag)) then
		call add_r0_tree(list%go_l,apar,last,js,jd)
	else if(lgt(apar,list%tag)) then
		call add_r0_tree(list%go_r,apar,last,js,jd)
	else ! found identity so duplicate entry
		list%irec=list%irec+1
		call compare_sd(list%ts,js,list%serr)
		call compare_sd(list%td,jd,list%derr)
	end if
end subroutine add_r0_tree
!
recursive subroutine extract_tree(list,my_info)
type(grec), pointer :: list
type(info), dimension(:), pointer :: my_info
	if(associated(list%go_l)) call extract_tree(list%go_l,my_info)
	if(associated(list%go_r)) call extract_tree(list%go_r,my_info)
	my_info(list%ti)=info(list%tag,list%ti,list%ts,list%td,0,0,0,0,0,list%irec,list%serr,list%derr)
	nullify(list%go_l,list%go_r); deallocate(list)
end subroutine extract_tree
!
recursive subroutine delete_tree_list(list)
type(grec), pointer :: list
	if(associated(list%go_l)) call delete_tree_list(list%go_l)
	if(associated(list%go_r)) call delete_tree_list(list%go_r)
	nullify(list%go_l); nullify(list%go_r); deallocate(list)
end subroutine delete_tree_list
!
function shuffle(n) ! Fisher Yates Durstenfeld
integer, intent(in) :: n
integer :: item, i, j
integer, dimension(n) :: shuffle
real :: u
	shuffle=[(i,i=1,n)]
	swapping: do i=n,2,-1
		call random_number(u)
		j=floor(u*i) + 1
		item=shuffle(j)
		shuffle(j)=shuffle(i)
		shuffle(i)=item
	end do swapping
end function shuffle
!
integer recursive function find_in_tree(list,apar) result(addi)
character(len=*), intent(in) :: apar
type(grec), pointer :: list
	if(apar==list%tag) then
		addi=list%ti
	elseif(llt(apar,list%tag)) then
		if(associated(list%go_l)) then; addi=find_in_tree(list%go_l,apar); else; addi=-1; end if
	elseif(lgt(apar,list%tag)) then
		if(associated(list%go_r)) then; addi=find_in_tree(list%go_r,apar); else; addi=-1; end if
	end if
end function find_in_tree
!
function fname(atoken,prog,finfo,fext) result(fnout)
character(len=*), intent(in) :: atoken,prog,finfo,fext
character(len=len_trim(atoken)+len_trim(prog)+len_trim(finfo)+len_trim(fext)+3) :: fnout
  fnout=trim(atoken)//'_'//trim(prog)//'_'//trim(finfo)//'.'//trim(fext)
end function fname
!
subroutine info_info(ich)
integer, intent(in) :: ich
  write(ich,'(a)') '... ==> tag    : the identifier from input'
  write(ich,'(a)') '... ==> ji     : a unique index number of id, NOT the renumbered id'
  write(ich,'(a)') '... ==> js     : the index number of the sire, 0 unknown'
  write(ich,'(a)') '... ==> jd     : the index number of the dam, 0 unknown'
  write(ich,'(a)') '... ==> jx     : the use as parent: 1 (sire), 2(dam), 3 (both) or 0 (unknown)'
  write(ich,'(a)') '... ==> ns     : the number of offspring as sire'
  write(ich,'(a)') '... ==> nd     : the number of offspring as dam'
  write(ich,'(a)') '... ==> isw    : the sweep number to obtain renumbered id, 0 if unswept'
  write(ich,'(a)') '... ==> nyi    : the renumbered id appearing in pedigree file, if valid pedigree'
  write(ich,'(a)') '... ==> irec   : the number of birth records found in mating list'
  write(ich,'(a)') '... ==> serr   : error flag (T/F) for conflicting sire information in multiple records'
  write(ich,'(a)') '... ==> derr   : error flag (T/F) for conflicting dam information in multiple records'
end subroutine info_info
!
subroutine info_ped(ich)
integer, intent(in) :: ich
  write(ich,'(a)') '... ==> nyi    : renumbered id'
  write(ich,'(a)') '... ==> nys    : renumbered sire, 0 unknown'
  write(ich,'(a)') '... ==> nyd    : renumbered dam, 0 unknown'
  write(ich,'(a)') '... ==> nyx    : use as parent: 1(sire), 2(dam), 3(both), or 0(unused)'
  write(ich,'(a)') '... ==> tag    : identifier from input'
end subroutine info_ped
!
subroutine write_list(ich,wfile,gen,list,append)
type(info), dimension(:), pointer :: gen
type(info), pointer :: z => null()
integer, dimension(:), intent(in) :: list
integer, intent(in) :: ich
character(len=*) :: wfile
character(len=13) :: fcsv='(*(g0,:,","))'
logical, intent(in) :: append
  if(append) then; open(ich,file=wfile,position='append'); 
  else; open(ich,file=wfile,status='replace')
  end if
  select case (ich) ! for header
  case(31); write(ich,'(a)') 'tag' ! birth records added
  case(32); write(ich,'(a)') 'tag,stag,dtag' ! selfing
  case(33); write(ich,'(a)') 'tag,ns,nd' ! multiple role parenting
  case(34); write(ich,'(a)') 'tag,irec,serr,derr' ! incompatible multiple records
  case(35); write(ich,'(a)') 'tag,stag,dtag' ! inferred sexes
  end select
  do i=1,size(list)
    z => gen(list(i))
    select case(ich)
    case(31); write(ich,'(a)') trim(z%tag)
    case(32); write(ich,fcsv) trim(z%tag),trim(gen(z%js)%tag),trim(gen(z%jd)%tag)
    case(33); write(ich,fcsv) trim(z%tag),z%ns,z%nd
    case(34); write(ich,fcsv) trim(z%tag),z%irec,z%serr,z%derr
    case(35); write(ich,fcsv) trim(z%tag),trim(gen(z%js)%tag),trim(gen(z%jd)%tag)
    end select
  end do
  nullify(z)
  close(ich)
end subroutine write_list
end module utilities    