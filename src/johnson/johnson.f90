program johnson ! indebted to https://www.youtube.com/watch?v=johyrWospv0 by Tushar Roy
use parameters
use utilities
use jaw_seeds
implicit none
!
! organize files
!
inquire(file='johnson_in.txt', exist=iask)
if(.not.iask) stop '... input file johnson_in.txt does not exist'
open(10,file='johnson_in.txt') 
read(10,*) atoken; atoken=adjustl(atoken)
open(14,file=fname(atoken,'johns','log','txt'), status='replace')
ifile=fname(atoken,'kosar','sets','csv')
inquire(file=ifile, exist=iask)
if(.not.iask) then
    call report_time([6,14],'... file of strongly connected sets '//trim(ifile)//' does not exist')
    read(*,*); stop
end if
! can proceed
open(11,file=ifile); read(11,*) acha ! skip header line 
open(15,file=fname(atoken,'johns','circ','csv'), status='replace')
open(16,file=fname(atoken,'johns','freq','csv'), status='replace')
open(17,file=fname(atoken,'johns','subg','csv'), status='replace')
write(14,'(2a)') '... input file ...                ',ifile
write(14,'(2a)') '... output for circuits found ... ',fname(atoken,'johns','circ','csv')
call circ_info(14)
write(14,'(2a)') '... output for frequencies ...    ',fname(atoken,'johns','freq','csv')
call freq_info(14)
write(14,'(2a)') '... output for subgraphs ...      ',fname(atoken,'johns','subg','csv')
call subg_info(14)
write(15,'(*(g0,:,","))') 'iset,icirc,istart,jcode,acode,tag' ! header line 
write(16,'(*(g0,:,","))') 'iset,tag,ctot,cfreq' ! header line
!
write(14,'(a)') '... random number seed read from jaw_seed.txt, and replaced by seed for next run' 
call seed_set(14) ! for shuffle
!
maxn=0
list_length: do
  read(11,*,iostat=ierr) acha
  if(ierr<0) exit list_length
  maxn=maxn+1
  !if(mod(maxn,100000)==0) print *, maxn,' counted'
end do list_length
if(maxn==0) stop '... no pedigree info in file! ... stopping!'
rewind(11); read(11,*) acha ! skip header line
call report_time([6,14],'... size of problem established')
allocate(my_sets(maxn))
read_sets: do i=1,maxn  
  read(11,*,iostat=ierr) my_sets
  if(ierr<0) exit read_sets 
end do read_sets  
nset=my_sets(maxn)%jset
write(atxt,'(a,i0,a,i0,a)') '... ',maxn,' nodes found for ',nset,' sets'
call report_time([6,14],atxt)
!
! cycle through sets
!
jlast=0
graph_list: do iset=1,nset	
  write(atxt,'(a,i4)') '... starting set ',iset
  ssize=count(my_sets(:)%jset==iset)
  !if(ssize==1) then
  !  jlast=jlast+ssize
  !  acha='m' ! must be sire or dam
  !  if(my_sets(jlast)%ji==my_sets(jlast)%jd) acha='f'
  !  write(15,'(*(g0,:,","))')  iset,"1","1","1",trim(acha),trim(my_sets(jlast)%tag)
  !  write(atxt,'(a,i8,2a)') '... ... finishing ',iset,' *** a self parent !!! ',my_sets(jlast)%tag
	!	 call report_line([6,14],atxt)
  !  cycle graph_list
  !end if
  allocate(graph(ssize),iperm(ssize))
  graph=net('NA',0,0,(/0,0/),(/0,0/),null())
  iperm=shuffle(ssize)
  read_next: do i=1,ssize ! remove second parent in a selfing as searching for simple circuits 
    ip=jlast+i
    zin => my_sets(ip)    
    znet => graph(iperm(i))
    znet=net(zin%tag,zin%ji,0,(/zin%js,zin%jd/),(/0,0/),null())
    call random_number(u)
    !if(u<0.5) call swap(znet%parent(1),znet%parent(2)) ! the random swapping allows testing
    if(znet%parent(1)==znet%parent(2)) znet%parent(2)=0
  end do read_next
  nullify(zin,znet)
  ! create subgraph for this set ... in strongly connected sets, everyone is an offspring & parent
  !maxp=maxval((/maxval(graph(:)%ji),maxval(graph(:)%parent(1)),maxval(graph(:)%parent(2))/))
  maxp=maxval(graph(:)%ji)
  allocate(jinv(0:maxp), source=0) ! lower bound 0 to simplify missing
  jinv(graph(:)%ji)=(/(i,i=1,size(graph))/)
  graph(:)%ji=(/(i,i=1,size(graph))/)
  renumber: do i=1,size(graph) ! renumber parents and count offspring
    jdum=jinv(graph(i)%parent(1))
    if(jdum>0) graph(jdum)%noff=graph(jdum)%noff+1
    jdum=jinv(graph(i)%parent(2))
    if(jdum>0) graph(jdum)%noff=graph(jdum)%noff+1
    graph(i)%parent(:)=jinv(graph(i)%parent(:)) ! now can renumber
  end do renumber
  do i=1,size(graph) ! noff should be >0 or else it is not strongly connected
    allocate(graph(i)%child(graph(i)%noff), source=sprog(0,.false.))
  end do
  allocate(moff(ssize),source=0)
  children: do i=1,size(graph) ! add offspring information to renumbered graph
    if(graph(i)%parent(1)>0) then
      jdum=graph(i)%parent(1)
      moff(jdum)=moff(jdum)+1
      graph(jdum)%child(moff(jdum))=sprog(i,.false.)
      graph(i)%ioff(1)=moff(jdum)
    end if
    if(graph(i)%parent(2)>0) then
      jdum=graph(i)%parent(2)
      moff(jdum)=moff(jdum)+1
      graph(jdum)%child(moff(jdum))=sprog(i,.false.)
      graph(i)%ioff(2)=moff(jdum)
    end if
  end do children
  deallocate(moff,jinv)
  do i=1,size(graph)
    znet => graph(i)
    write(17,'(*(g0,:,","))') iset,trim(znet%tag),znet%ji,znet%parent(:),znet%noff,znet%ioff(:),znet%child(:)%ji
  end do
  nullify(znet)
  ! carry out search
  icirc=0
  allocate(removed(size(graph)),blocked(size(graph)), source=.false.)
  allocate(cfreq(size(graph)), source=0)
  circuit_search: do i=1,size(graph)
    allocate(my_path)
    my_path=path(i,1,null())
    blocked(:)=.false.
    do ii=1,size(graph)
      graph(i)%child(:)%xji=.false.
    end do
    istart=i
    call add_step(graph,removed,blocked,my_path,istart,icirc,15,iset)
    !call check_path(my_path)
		!print *, associated(my_path)
		!if(associated(my_path)) deallocate(my_path)
    nullify(my_path)
    removed(i)=.true.
  end do circuit_search
  do i=1,size(graph)
    write(16,'(*(g0,:,","))') iset,trim(graph(i)%tag),cfreq(i),real(cfreq(i))/icirc
  end do
  do i=1,size(graph)
    deallocate(graph(i)%child)
  end do
  deallocate(graph,removed,blocked,cfreq,iperm)
	write(atxt,'(a,i8,a,i8,a)') '... ... finishing ',iset,' with ',icirc,' circuits found'
	call report_line([6,14],atxt)
  jlast=jlast+ssize
end do graph_list
call report_time([6,14],'... job finished')
!
end program johnson