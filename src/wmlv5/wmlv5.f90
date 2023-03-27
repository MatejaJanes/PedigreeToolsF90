program meuwissen_luo
use parameters
use utilities
use jaw_meu_luo
implicit none
!
! setting up files
!
inquire(file='wmlv5_in.txt', exist=iask)
if(.not.iask) stop '... input file wmlv5_in.txt does not exist'
open(unit=10,file='wmlv5_in.txt')
read(10,*) atoken; atoken=adjustl(atoken)
read(10,*) adum; adum=adjustl(adum)
open(unit=14,file=fname(atoken,'wmlv5','log','txt'), status='replace')
ifile=fname(atoken,adum,'ped','csv')
inquire(file=ifile, exist=iask)
if(.not.iask) then
  call report_time([6,14],'... pedigree input file '//trim(ifile)//' does not exist')
  stop
end if
! can proceed
open(unit=11,file=ifile); read(11,*) adum ! skip the header
ofile=fname(atoken,'wmlv5','pedf','csv')
open(unit=15,file=ofile, status='replace')
write(14,'(2a)') '... pedigree input file  ',ifile
write(14,'(3a)') '... pedigree output file ',ofile,' with F added ahead of tag'
call report_time([6,14],'... files present and opened')
! size of problem and initialisation
maxp=0
list_length: do
  read(11,*,iostat=ierr) adum
  if(ierr<0) exit list_length
  maxp=maxp+1
  if(mod(maxp,100000)==0) print *, maxp,' counted'
end do list_length
write(atxt,'(a,i8,a)') '... ',maxp,' entries found in pedigree'
call report_time([6,14],atxt) 
if(maxp==0) stop '... no pedigree in file! ... stopping!'
rewind(11); read(11,*) adum ! skip the header
!
! create space and read in
!
allocate(my_ped(maxp))
my_ped(:)=info(0,0,0,0,'NA',0.) ! F set to 0
do i=1,maxp
  read(11,*,iostat=ierr) my_ped(i)%ji,my_ped(i)%js,my_ped(i)%jd,my_ped(i)%jx,my_ped(i)%tag
  my_ped(i)%tag=adjustl(my_ped(i)%tag)
end do
np=max(maxval(my_ped(:)%js),maxval(my_ped(:)%jd))
write(14,'(a,i10)') '... maximum index of parent ',np 
!
!call report_time([6,14],'... finished reading, testing Meuwissen & Luo')
!call test_meu_luo(6)
!read(*,*)
!
! calculate f
!
call report_time([6,14],'... finished setting up, starting modified Meuwissen & Luo')
call meu_luo(my_ped,maxp,np) ! does not look for consecutive FS
call report_time([6,14],'... finished modified Meuwissen & Luo, starting output')
!
! write out
!
write(15,'(a)') 'ji,js,jd,jx,f,tag'
do i=1,maxp
  z => my_ped(i)
  write(15,'(*(g0.6,:,","))') z%ji,z%js,z%jd,z%jx,z%f,trim(z%tag)
end do
nullify(z)
write(14,'(a)') '... ... use Colleau algorithm for further functions of relationships'
open(45,file='colleau_in.txt',status='replace') 
write(45,'(a)') trim(atoken)
write(14,'(2a)') '... ... with input files colleau_in.txt and ',trim(ofile)
call report_time([6,14],'... job is finished')
deallocate(my_ped,ifile,ofile)
!
end program meuwissen_luo
