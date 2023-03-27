program kosaraju ! indebted to https://www.youtube.com/watch?v=RpgcYiky7uw by Tushar Roy
use parameters
use utilities
implicit none
!
! setting up files
!
inquire(file='kosaraju_in.txt', exist=iask)
if(.not.iask) stop '... input file kosaraju_in.txt does not exist'
open(10,file='kosaraju_in.txt')
read(10,*) atoken; atoken=adjustl(atoken) ! input file from pedv7_out ...
read(10,*) acha; acha=adjustl(acha)
open(unit=14,file=fname(atoken,'kosar','log','txt'),status='replace')
ifile=fname(atoken,acha,'info','csv')
inquire(file=ifile, exist=iask)
if(.not.iask) then
    call report_time([6,14],'... pedigree information file '//trim(ifile)//' does not exist')
    read(*,*); stop
end if
open(unit=11,file=ifile); read(11,*) acha ! skip header line
write(14,'(2a)') '... input file ... ',ifile
!
!size of problem and initialisation
!
maxp=0
list_length: do
  read(11,*,iostat=ierr) acha
  if(ierr<0) exit list_length
  maxp=maxp+1
  if(mod(maxp,100000)==0) print *, maxp,' counted'
end do list_length
write(atxt,'(a,i8,a)') '... ',maxp,' records found'
call report_line([6,14],atxt)
if(maxp==0) stop '... no pedigree info in file! ... stopping!'
rewind(11); read(11,*) acha ! skip header line
call report_time([6,14],'... size of problem established')
!
! create in_graph
!
allocate(v7out(maxp),source=info('NA',0,(/0,0/),0,0,0,0,0))
allocate(in_graph(maxp),source=.true.)
allocate(not_s(maxp),not_d(maxp),source=.false.)
do i=1,size(v7out)
  read(11,*) v7out(i)
end do
if(all(v7out(:)%isw>0)) then
  write(atxt,'(a)') '... all individuals were swept by v7 and no further processing needed'
  call report_line([6,14],atxt); call exit
end if
where(v7out(:)%isw>0) in_graph(:)=.false.
nsub=count(in_graph(:))
write(14,'(a,i8)') '... nodes in maximum graph ',nsub
where(v7out(:)%jpar(1)==0) not_s=.true.
where(v7out(:)%jpar(2)==0) not_d=.true.
exclude_base: do
  where(not_s(:).and.not_d(:)) in_graph(:)=.false.
  if(.not.count(in_graph(:))<nsub) exit exclude_base
  nsub=count(in_graph(:))
  do i=1,maxp
    if(v7out(i)%jpar(1)>0) then
      if(.not.in_graph(v7out(i)%jpar(1))) not_s(i)=.true.
    end if
    if(v7out(i)%jpar(2)>0) then
      if(.not.in_graph(v7out(i)%jpar(2))) not_d(i)=.true.
    end if
  end do
end do exclude_base 
write(atxt,'(a,i8)') '... after removing from top, nodes in subgraph ',nsub 
call report_line([6,14],atxt)
!
! special case of 1 node
!
if(nsub==1) then
  write(atxt,'(a)') '... *** only ONE node left in graph and listed below and stopping'
  call report_line([6,14],atxt)
  do i=1,size(in_graph)
    if(.not.in_graph(i)) cycle
    write(14,'(a,6i8)') trim(v7out(i)%tag),v7out(i)%ji,v7out(i)%jpar(:),v7out(i)%jx,v7out(i)%ns,v7out(i)%nd
  end do
  call exit
end if
!
! prepare connected graph, requires parents and offspring in graph
!
allocate(graph(nsub), source=net('NA',0,0,0,(/0,0/),null()))	
allocate(s_inv(maxp), source=0)
ivals=pack(v7out(:)%ji,in_graph(:)) ! mapping relevant identities from graph to v7out
do i=1,size(ivals)
  s_inv(ivals(i))=i ! an inverse mapping
end do
graph(:)%tag=v7out(ivals(:))%tag ! labels graph with tags
graph(:)%ji=(/(i,i=1,nsub)/)     ! now can renumber graph nodes 
trace_parents: do i=1,nsub 
  if(v7out(ivals(i))%jpar(1)>0) then ! has known sire
    if(in_graph(v7out(ivals(i))%jpar(1))) then ! if that is in the graph
      graph(i)%jpar(1)=s_inv(v7out(ivals(i))%jpar(1)) ! set the node number of sire
      graph(graph(i)%jpar(1))%npat=graph(graph(i)%jpar(1))%npat+1 ! increment sire's offspring count
    end if
  end if
  if(v7out(ivals(i))%jpar(2)>0) then ! has known dam
    if(in_graph(v7out(ivals(i))%jpar(2))) then ! if that is in the graph 
      graph(i)%jpar(2)=s_inv(v7out(ivals(i))%jpar(2)) ! set the node number of dam
      graph(graph(i)%jpar(2))%nmat=graph(graph(i)%jpar(2))%nmat+1 ! increment dam's offspring count
    end if
  end if
end do trace_parents
get_offspring: do i=1,nsub ! allocate and populate joffren
  if((graph(i)%npat+graph(i)%nmat)>0) then ! now know how many joffren for nodes
    ns=graph(i)%npat 
    nd=graph(i)%nmat
    allocate(graph(i)%joff(ns+nd), source=0)
    graph(i)%joff(1:ns)=pack(graph(:)%ji,graph(:)%jpar(1)==i) ! identify the first set of joffren
    graph(i)%joff((ns+1):(ns+nd))=pack(graph(:)%ji,graph(:)%jpar(2)==i) ! ... and second set
  end if
end do get_offspring
open(unit=15,file=fname(atoken,'kosar','subg','csv'), status='replace')
do i=1,nsub
  z => graph(i)
  write(15,'(*(g0,:,","))'), trim(z%tag),z%ji,z%jpar(:),z%npat,z%nmat,z%joff(:)
end do
nullify(z)
deallocate(in_graph,v7out,not_s,not_d,s_inv)
write(14,'(2a)') '... connected graph formed and written out to file ',fname(atoken,'kosar','subs','csv')
call subg_info(14)
call report_time([6,14],'... date and time after completing graph')
!
! pass 1 of Kosaraju
!
allocate(visited(nsub),tour(nsub),source=.false.)
allocate(stack(nsub),store(nsub),source=0)
jend=0
do i=1,size(graph)
	if(.not.visited(i)) then
		allocate(my_path)
		my_path%node=i
		call follow_pass_1(graph,visited,stack,jend,my_path)
		deallocate(my_path)
	end if
end do
write(atxt,'(a,i8,a,i8,a)') '... pass 1 completed with stack ending at ',jend,' and ',count(.not.visited(:)),' unvisited'
call report_line([6,14],atxt)
call report_time([6,14],'... pass 1 of 2 completed')
!
! pass 2 of Kosaraju
!
visited(:)=.false. ! reset visited
iset=0
do i=size(graph),1,-1
	ii=stack(i)
	if(.not.visited(ii)) then ! will be start of a new set 
		allocate(my_path)
		my_path%node=ii
    tour(:)=.false.
		call follow_pass_2(graph,visited,tour,my_path)
    lset=count(tour(:))
    if(lset>1.or.any(graph(ii)%jpar(:)==ii)) call report_set(graph,tour,iset,store)
		deallocate(my_path)
  end if
end do
write(14,'(a,i4,a,i8,a)') '... pass 2 completed with ',iset,' sets and ',count(.not.visited(:)),' unvisited'
call report_time([6,14],'... pass 2 of 2 completed')
!
! report on sets
!
open(unit=16,file=fname(atoken,'kosar','sets','csv'), status='replace')
write(16,'(a)') 'node,snode,dnode,tag'
do i=1,iset
  ivals=pack((/(j,j=1,size(store))/),store(:)==i)
  write(14,'(a,i8,a,i8)') '... self-ancestry set ... ',i,' ... length ... ',size(ivals)
  do j=1,size(ivals)
    z => graph(ivals(j))
    ! must have a joff to be strongly connected, otherwise parent could not reach it
    write(16,'(*(g0,:,","))') i,z%ji,z%jpar(:),trim(z%tag)
  end do
end do
nullify(z)
write(14,'(2a)') '... connected sets written out to file ',fname(atoken,'kosar','sets','csv')
call sets_info(14)
!
! clean up by deallocating
!
do i=1,size(graph)
  if(associated(graph(i)%joff)) deallocate(graph(i)%joff) 
end do
deallocate(graph,visited,tour,stack,store,ivals)
write(14,'(a)') '... ... use Johnson algorithm to find all the circuits within the sets'
open(45,file='johnson_in.txt',status='replace') 
write(45,'(a)') trim(atoken)
write(14,'(2a)') '... ... by using input files johnson_in.txt and ',fname(atoken,'kosar','sets','csv')
call report_time([6,14],'... job finished')
end program kosaraju
