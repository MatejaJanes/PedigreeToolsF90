program pedigree
!
! ... creates an integer valued pedigree file from an unordered mating list
! ... may contain duplicate entries which are checked for inconsistency
! ... may be monoecious
!
use parameters
use utilities
use jaw_seeds
implicit none
!
! setting up, input and output files ... unit 6 is standard output 
!
inquire(file='pedv7_in.txt', exist=iask)
if(.not.iask) stop '... input file pedv7_in.txt does not exist'
open(unit=10,file='pedv7_in.txt')
read(10,*) ifile ; ifile=adjustl(ifile) ! (unordered) pedigree mating list, input file
read(10,*) iskip ! # header lines in input to skip
read(10,*) atoken; atoken=adjustl(atoken) ! tag for files
open(unit=14,file=fname(atoken,'pedv7','log','txt'),status='replace')  ! log file
inquire(file=trim(ifile), exist=iask)
if(.not.iask) then
  call report_time([6,14],'... pedigree input file '//trim(ifile)//' does not exist')
  read(*,*); stop
end if
read(10,*) flag1 ! flag1 for parents without birth records ... 0 (no) or 1 (yes) ... '
read(10,*) flag2 ! flag2 for selfings ... 0 (no) or 1 (yes) ... '
read(10,*) flag3 ! flag3 for mixed use parents ... 0 (no) or 1 (yes) ... '
read(10,*) flag4 ! flag4 for incompatability of multiple records ... 0 (no) or 1 (yes) ... '
read(10,*) flag5 ! flag5 for self parenting ... 0 (no) or 1 (yes) ... '
open(unit=11,file=trim(ifile)); do i=1,iskip; read(11,*) acha; end do
write(14,'(a,t65,a)') '... (unordered) pedigree mating list, input file ',ifile
write(14,'(a,t65,a)') '... output file for valid pedigree (with header line) ',fname(atoken,'pedv7','ped','csv')
call info_ped(14) 
write(14,'(a,t65,a)') '... mating list information file (with header line) ',fname(atoken,'pedv7','info','csv')
call info_info(14)
!
write(14,'(a)') '... random number seed read from jaw_seed.txt, and replaced by seed for next run' 
call seed_set(14) ! for shuffle
nkn=['0 ','NA']  ! trailing blank added to clarify comparisons
write(14,'(4a)') '... unknown characters ',trim(nkn(1)),' and ',nkn(2)
if(flag1/=0) then; atxt='enabled with output '//fname(atoken,'pedv7','err1','csv') 
else; atxt='disabled'; end if
write(14,'(a,t45,a)') '... check additional birth records ',trim(adjustl(atxt))
if(flag2/=0) then; atxt='enabled with output '//fname(atoken,'pedv7','err2','csv')
else; atxt='disabled'; end if
write(14,'(a,t45,a)') '... check selfing records ',trim(adjustl(atxt))
if(flag3/=0) then; atxt='enabled with output '//fname(atoken,'pedv7','err3','csv')
else; atxt='disabled'; end if
write(14,'(a,t45,a)') '... check mixed use parents ',trim(adjustl(atxt))
if(flag4/=0) then; atxt='enabled with output '//fname(atoken,'pedv7','err4','csv')
else; atxt='disabled'; end if
write(14,'(a,t45,a)') '... check multiple record incompatibility ',trim(adjustl(atxt))
if(flag5/=0) then; atxt='enabled with output '//fname(atoken,'pedv7','err5','csv')
else; atxt='disabled'; end if
write(14,'(a,t45,a)') '... check self parenting ',trim(adjustl(atxt))
call report_time([14],'... files present, flags set, and ready to start processing')
!
! size of problem and initialisation ... this just counts the mating list
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
if(maxp==0) stop '... no mating list in file! ... stopping!'
rewind(11); do i=1,iskip; read(11,*) acha; end do
call report_time([14],'... date and time after counting')
!
! read in mating list and create a shuffle for efficient searching
!
allocate(my_rec(maxp),jord(maxp))
do i=1,maxp
  read(11,*,iostat=ierr) my_rec(i) 
  my_rec(i)%ai=adjustl(my_rec(i)%ai)
  my_rec(i)%as=adjustl(my_rec(i)%as)
  my_rec(i)%ad=adjustl(my_rec(i)%ad)
end do
jord(:)=shuffle(size(jord))
!
! read in sires to place at head of the tree, and then dams
!
jlast=0
build_list_sires: do i=1,maxp
  j=jord(i)
  adum=my_rec(j)%as
  if(.not.any(nkn==adum)) call add_sd_tree(sd_tree,adum,jlast)
  if(mod(i,100000)==0) print *, i,' sire records processed to list'
end do build_list_sires
nsir=jlast
write(atxt,'(a,i8,a)') '... list of sires built with ',nsir,' entries'
call report_line([6,14],atxt)
call report_time([14],'... date and time after adding sires ')
!
build_list_dams: do i=1,maxp ! done after sires assuming that sires are less numerous than dams
  j=jord(i)
  adum=my_rec(j)%ad 
  if(.not.any(adum==nkn)) call add_sd_tree(sd_tree,adum,jlast)
  if(mod(i,100000)==0) print *, i,' dam records processed to list'
end do build_list_dams
ndam=jlast-nsir
npar=jlast
write(atxt,'(a,i8,a)') '... list of dams built with  ',ndam,' additional entries'
call report_line([6,14],atxt)
write(atxt,'(a,i8,a)') '... ... with a total of .......  ',npar,' parents'
call report_line([6,14],atxt)
call report_time([14],'... date and time after adding dams')
!
! deal with the individuals
!
build_offspring: do i=1,maxp
  j=jord(i)
  if(any(nkn==my_rec(j)%ai)) then
    write(14,'(a,i8,a)') '... entry ',j,' has an unknown as its identity and is skipped'
    cycle build_offspring
  end if
  if(any(nkn==my_rec(j)%as)) then; jsir=0
  else 
    jsir=find_in_tree(sd_tree,my_rec(j)%as)
    if(jsir<0) write(14,'(a)') '... this should not happen ... in build_offspring with jsir ... stopping'
    if(jsir<0) stop '... this should not happen ... in build_offspring with jsir ... stopping'
  end if
  if(any(nkn==my_rec(j)%ad)) then; jdam=0
  else
    jdam=find_in_tree(sd_tree,my_rec(j)%ad)
    if(jdam<0) write(14,'(a)') '... this should not happen ... in build_offspring with jdam ... stopping'
    if(jdam<0) stop '... this should not happen ... in build_offspring with jdam ... stopping'
  end if
  adum=my_rec(j)%ai
  call update_sd_tree(sd_tree,adum,jsir,jdam,addi) ! update if identity is sire or dam
  if(addi<0) call add_r0_tree(r0_tree,adum,jlast,jsir,jdam)
  if(mod(i,100000)==0) print *, i,' offspring records processed to list'
end do build_offspring
!
! summarise the processing of the mating list
!
nr0 = jlast-npar ! number not processed as a parent
write(atxt,'(a,i8,a)') '... list of r0s built with   ',  nr0,' entries'
call report_line([6,14],atxt)
write(atxt,'(a,i8,a)') '... pedigree list has total  ',jlast,' entries'
call report_line([6,14],atxt)
call report_time([14],'... date and time after adding offspring')
deallocate(my_rec,jord) ! information extracted
if(jlast==0) then
	write(14,'(a)') '... no mating records with recognised identities! ... stopping!'
  stop            '... no mating records with recognised identities! ... stopping!'
end if
!
! now extract info for all individuals encountered in mating list
!
allocate(my_info(jlast))
my_info(:)=info('NA',0,0,0,0,0,0,0,0,0,.false.,.false.)
if(associated(sd_tree)) call extract_tree(sd_tree,my_info) ! it deletes the tree as it leaves
if(associated(r0_tree)) call extract_tree(r0_tree,my_info) ! it deletes the tree as it leaves
nullify(sd_tree,r0_tree)
write(14,'(a,l2)') '... sd_tree deallocated ... ', .not.associated(sd_tree)
write(14,'(a,l2)') '... r0_tree deallocated ... ', .not.associated(r0_tree)
call report_time([14],'... date and time after forming info')
!
! count offspring and add to information
!
mating_stats: do i=1,jlast
  if(my_info(i)%js>0) my_info(my_info(i)%js)%ns=my_info(my_info(i)%js)%ns+1
  if(my_info(i)%jd>0) my_info(my_info(i)%jd)%nd=my_info(my_info(i)%jd)%nd+1
end do mating_stats
zinfo => my_info(1:npar) ! points to parents
where(zinfo%ns>0) zinfo%jx=1 ! sire usage
where(zinfo%nd>0) zinfo%jx=zinfo%jx+2 ! dam usage
call report_time([14],'... offspring counted')
! 
! report on mating list anomalies
! 
write(14,'(a)') '... mating list anomalies report ... '
!
! no birth record
err1=pack(zinfo%ji,zinfo%irec==0)
write(14,'(a,i8,a)') '... *** ',size(err1),' parent(s) have no record as an offspring'
if(flag1/=0.and.size(err1)>0) then
  efile=fname(atoken,'pedv7','err1','csv')
  call write_list(31,efile,my_info,err1,.false.)
end if
deallocate(err1)
! selfings   
err2=pack(my_info(:)%ji,my_info(:)%js>0.and.my_info(:)%js==my_info(:)%jd)
write(14,'(a,i8,a)') '... *** ',size(err2),' offspring from selfing records'
if(flag2/=0.and.size(err2)>0) then
  efile=fname(atoken,'pedv7','err2','csv')
  call write_list(32,efile,my_info,err2,.false.)
end if
deallocate(err2)
! mixed use parenting   
err3=pack(zinfo%ji,zinfo%jx==3)
write(14,'(a,i8,a)') '... *** ',size(err3),' parents used as both sire and dam'
if(flag3/=0.and.size(err3)>0) then
  efile=fname(atoken,'pedv7','err3','csv')
  call write_list(33,efile,my_info,err3,.false.)
end if
deallocate(err3)
! multiple records
err4=pack(my_info(:)%ji,my_info(:)%irec>1)
write(14,'(a,i8,a)') '... *** ',size(err4),' identities with multiple parentage records'
if(size(err4)>0) then
  err4x=pack(err4(:),my_info(err4(:))%serr.or.my_info(err4(:))%derr)
  write(14,'(a,i8,a)') '... ==> ',size(err4x),' identities with incompatible multiple records'
  if(flag4/=0.and.size(err4x)>0) then
    efile=fname(atoken,'pedv7','err4','csv')
    call write_list(34,efile,my_info,err4x,.false.)
  end if
  deallocate(err4x)
end if
deallocate(err4)
! self parenting 
err5=pack(my_info(:)%ji,(my_info(:)%ji==my_info(:)%js).or.(my_info(:)%ji==my_info(:)%jd))
write(14,'(a,i8,a)') '... *** ',size(err5),' self parent(s), if >0 renumbering must FAIL as pedigree invalid'
if(flag5/=0.and.size(err5)>0) then
  efile=fname(atoken,'pedv7','err5','csv')
  call write_list(35,efile,my_info,err5,.false.)
end if
deallocate(err5)
!
! summary
!
allocate(noff(size(my_info))) ! used for sweeping
noff(:)=my_info(:)%ns+my_info(:)%nd
zoff => noff(1:npar) ! ... parallel pointer to zinfo
kpar=count(zoff(:)>0)
kmax=maxval(zoff(:)) 
emax=pack(zinfo(:)%ji,zoff(:)==kmax)
write(atxt,'(a,i8,a,i8)') '... pedigree now has ',kpar,' parent(s) with maximum offspring ',kmax
call report_line([6,14],atxt)
write(atxt,'(a,i8,a)') '... with ',size(emax),' parent(s) with maximum offspring numbers'
call report_line([6,14],atxt)
write(14,'(a,t19,3a)') '...','... e.g. (up to two) identities ... ',(zinfo(emax(i))%tag, i=1,min(2,size(emax)))
write(atxt,'(a,i8,a)') '... with ',npar-kpar,' original parent(s) lost due to conflicting info'
call report_line([6,14],atxt)
nullify(zoff,zinfo)
deallocate(emax)
call report_time([14],'... summaries printed')
!
! now sweep and renumber 
!
klast=jlast
ksweep=0
sweeper: do 
  ksweep=ksweep+1
  sweep=pack(my_info(:)%ji,noff(:)==0)
  if(size(sweep)==0) then
    deallocate(sweep)
    exit sweeper
  end if
  tick_list: do i=size(sweep),1,-1
    ii=sweep(i)
    noff(ii)=noff(ii)-1 ! goes to -1 so it does not get re-swept
    my_info(ii)%nyi=klast
    my_info(ii)%isw=ksweep
    if(my_info(ii)%js>0) noff(my_info(ii)%js)=noff(my_info(ii)%js)-1
    if(my_info(ii)%jd>0) noff(my_info(ii)%jd)=noff(my_info(ii)%jd)-1
    klast=klast-1
  end do tick_list
  write(14,'(a,i4,a,i8,a,i8,a)') '... at sweep ',ksweep,': ',size(sweep),' entries were swept, leaving ',klast,' for sweeping' 
  deallocate(sweep)
  if(klast==0) exit sweeper
end do sweeper
call report_time([14],'... sweeping completed') 
if(klast==0) then ! all is well
	allocate(ped_list(size(my_info)))
	ped_list(:)=prec(0,0,0,0,'NA')
  unscramble: do i=1,size(my_info)
    ped_list(my_info(i)%nyi)%ji=my_info(i)%nyi
    if(my_info(i)%js>0) ped_list(my_info(i)%nyi)%js=my_info(my_info(i)%js)%nyi
    if(my_info(i)%jd>0) ped_list(my_info(i)%nyi)%jd=my_info(my_info(i)%jd)%nyi
    ped_list(my_info(i)%nyi)%jx=my_info(i)%jx
    ped_list(my_info(i)%nyi)%tag=my_info(i)%tag
  end do unscramble
  open(unit=15,file=fname(atoken,'pedv7','ped','csv'), status='replace')
  write(15,'(a)') 'nyi,nys,nyd,nyx,tag'
  do i=1,size(ped_list)
    write(15,'(*(g0,:,","))') ped_list(i)
  end do
  deallocate(ped_list) ! done
  call report_time([14],'... renumbered pedigree written to file')
  write(14,'(a)') '... ... use modified Meuwissen & Luo algorithm to find inbreeding coefficients'
  open(45,file='wmlv5_in.txt',status='replace') 
  write(45,'(a/a)') trim(atoken),'pedv7'
  write(14,'(2a)') '... ... with input files wmlv5_in.txt and ',fname(atoken,'pedv7','ped','csv')
else if(klast>0) then ! all is not well
  write(atxt,'(a)') '... *** circular pedigrees detected'
  call report_line([6,14],atxt)
  write(14,'(a)') '... ... use Kosaraju algorithm to find distinct subsets with self-ancestry'
  open(45,file='kosaraju_in.txt',status='replace') 
  write(45,'(a/a)') trim(atoken),'pedv7'
  write(14,'(2a)') '... ... with input files kosaraju_in.txt and ',fname(atoken,'pedv7','info','csv')
end if
open(unit=16,file=fname(atoken,'pedv7','info','csv'), status='replace') ! list info file
write(16,'(a)') 'tag,ji,js,jd,jx,ns,nd,isw,nyi,irec,serr,derr'
do i=1,size(my_info) ! report the full information
  z => my_info(i)
  write(16,'(*(g0,:,","))') trim(z%tag),z%ji,z%js,z%jd,z%jx,z%ns,z%nd,z%isw,z%nyi,z%irec,z%serr,z%derr
end do
nullify(z)
call report_time([14],'... accumulated information written to file')
deallocate(my_info,noff)
write(14,'(a,l2,l2)') '... my_info, noff deallocated ... ',.not.associated(my_info),.not.allocated(noff)
call report_time([14],'... job is finished')
!read(*,*)
end program pedigree
