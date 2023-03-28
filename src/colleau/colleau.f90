program colleau
use parameters
use utilities
use jaw_seeds
implicit none
!
! setting up files
!
inquire(file='colleau_in.txt', exist=iask)
if(.not.iask) stop '... input file colleau_in.txt does not exist'
open(unit=10,file='colleau_in.txt')
read(10,*) atoken; atoken=adjustl(atoken)
ifile=fname(atoken,'wmlv5','pedf','csv')
afile=fname(atoken,'colle','actv','csv')
open(unit=14,file=fname(atoken,'colle','log','txt'))
inquire(file=ifile, exist=iask); inquire(file=afile, exist=jask);
if(.not.iask.or..not.jask) then
  if(.not.iask) call report_time([6,14],'... pedigree input file '//trim(ifile)//' does not exist')
  if(.not.jask) call report_time([6,14],'... active input file '//trim(afile)//' does not exist')
  stop
end if
! can proceed
open(unit=11,file=ifile); read(11,*) adum ! skip header line
open(unit=12,file=afile)
ofile=fname(atoken,'colle','pedc','csv')
open(unit=15,file=ofile)
call report_time([6,14],'... files present and opened')
write(14,'(2a)') '... pedigree input file  ',ifile
write(14,'(3a)') '... pedigree active file ',afile,' assumed to have NO header line'
write(14,'(3a)') '... pedigree output file ',ofile,' with 2 active set values added before tag'
call report_time([6,14],'... files present')
write(14,'(a)') '... random number seed read from jaw_seed.txt, and replaced by seed for next run' 
call seed_set(14)
!
! size of problem and initialisation
!
maxp=0
list_length: do
  read(11,*,iostat=ierr) adum
  if(ierr<0) exit list_length
  maxp=maxp+1
  if(mod(maxp,100000)==0) print *, maxp,' counted'
end do list_length
print *, maxp,' entries found in pedigree'
write(14,'(a,i8,a)') '... ',maxp,' entries found in pedigree'
if(maxp==0) stop '... no pedigree in file! ... stopping!'
rewind(11); read(11,*) adum ! skip header line 
call report_time([6,14],'... date & time after counting')
!
! read in and create a shuffle for efficient searching
!
allocate(my_ped(maxp),jord(maxp))
my_ped(:)=info(0,0,0,0,0.,'NA',0.,0.,0.,.false.)
do i = 1,maxp
  read(11,*,iostat=ierr) my_ped(i)%ji,my_ped(i)%js,my_ped(i)%jd,my_ped(i)%jx,my_ped(i)%f,my_ped(i)%tag
  my_ped(i)%tag=adjustl(my_ped(i)%tag)
end do
jord(:)=shuffle(size(jord))
build_binary_tree: do i=1,maxp
  j=jord(i)
  call add_ped_tree(ped_tree,my_ped(j)%tag,my_ped(j)%ji)
  if(mod(i,100000)==0) print *, i,' pedigree records processed to tree'
end do build_binary_tree
deallocate(jord)
call report_time([6,14],'... date & time after building pedigree tree')
!
! add information on active list
!
iact=0
read_active: do 
  read(12,'(a)',iostat=ierr) adum
  if(ierr<0) exit read_active
  adum=trim(adjustl(adum))
  iact=iact+1
  iget=.false.
  call find_in_tree(ped_tree,my_ped,adum,iget)
  if(.not.iget) write(14,'(3a)') '... identity ',trim(adum),' in active list not found in pedigree'
  if(mod(iact,1000)==0) print '(i8,a)', iact,' individuals read from active file'
end do read_active
print '(i8,a)', iact,' individuals read from active file, including those not found'
if(iact==0) then
  call report_time([6,14],'... no active list found ... stopping!')
  stop
else 
  write(14,'(a,i8,a,i8,a)') '... ',iact,' active identities read, ',count(my_ped(:)%actv),' found in pedigree'
end if
call delete_tree(ped_tree)
nullify(ped_tree)
write(14,'(a,l2)') '... ped_tree deallocated ... ', .not.associated(ped_tree)
call report_time([6,14],'... date & time after forming info with active list, starting Colleau procedure')
!
! start calculation with y1 = (I-T)'x, y1 is set to 0 at start
!
where(my_ped(:)%actv) my_ped(:)%x=1. ! x added for generality of other subsets being desired 
do i=size(my_ped),1,-1
  my_ped(i)%y1=my_ped(i)%x+my_ped(i)%y1
  if(my_ped(i)%js>0) my_ped(my_ped(i)%js)%y1=my_ped(my_ped(i)%js)%y1+0.5*my_ped(i)%y1
  if(my_ped(i)%jd>0) my_ped(my_ped(i)%jd)%y1=my_ped(my_ped(i)%jd)%y1+0.5*my_ped(i)%y1
end do
!
! finish calculation with (I-T)y = Dy1
!
allocate(fdum(0:maxp))
fdum(1:maxp)=my_ped(:)%f; fdum(0)=-1.
do i=1,size(my_ped)
  ff=0.5-0.25*fdum(my_ped(i)%js)-0.25*fdum(my_ped(i)%jd) ! fdum has -1 if js or jd is 0
  my_ped(i)%y=ff*my_ped(i)%y1
  if(my_ped(i)%js>0) my_ped(i)%y=my_ped(i)%y+0.5*my_ped(my_ped(i)%js)%y
  if(my_ped(i)%jd>0) my_ped(i)%y=my_ped(i)%y+0.5*my_ped(my_ped(i)%jd)%y
end do
my_ped(:)%y1=my_ped(:)%y-1.-my_ped(:)%f
my_ped(:)%y=my_ped(:)%y/real(iact)
if(iact>1) my_ped(:)%y1=my_ped(:)%y1/real(iact-1.)
call report_time([6,14],'... Colleau procedure finished')
!
! write out & clean up
!
do i=1,maxp
  if(.not.my_ped(i)%actv) cycle 
  z => my_ped(i)
  write(15,'(*(g0.6,:,","))') z%ji,z%js,z%jd,z%jx,z%f,z%y,z%y1,trim(z%tag)
end do
nullify(z)
write(14,'(a)') '... relationship parameters for active set written to file '
gcoa=sum(my_ped(:)%y,my_ped(:)%actv)/real(iact)/2.
write(14,'(a,f10.5)') '... group coancestry (F after random mating within active set) ', gcoa
deallocate(my_ped,fdum)
call report_time([6,14],'... job is finished')
end program colleau
