module jaw_meu_luo
!
! this uses the standard formula used in the Meuwissen & Luo algorithm
!
use parameters
use utilities
implicit none
type fdat
  logical :: qq             ! ancestor flag
  real(kind=dp) :: dd,rr ! F, var(a) & contribution for ancestor
end type fdat
!
contains
!
subroutine meu_luo(ped,maxp,np)
type(info), dimension(:), pointer, intent(in) :: ped
integer, intent(in) :: maxp,np
type(fdat), dimension(0:np) :: wp
real(kind=dp), dimension(0:maxp) :: ff
integer :: i,ks,kd,kx,ka,kas,kad,na,t1,t2
real(kind=dp) :: di,fi,rk
  wp(:)=fdat(.false.,0._dp,0._dp)
  ff(:)=0._dp; ff(0)=-1._dp
  call system_clock(t1)
  do i=1,maxp ! F one by one
    ks=ped(i)%js; kd=ped(i)%jd
    di=0.50_dp-0.25_dp*(ff(ks)+ff(kd)) ! Mendelian variance needed by descendants
    if(i.le.np) wp(i)%dd=di ! needed for descendants
    kx=min(ks,kd)
    if(kx==0) cycle ! F=0 so move on
    fi=-1._dp ! offset, could loop to 0
    fi=fi+di  ! the individual contribution to F
    wp(ks)%qq=.true.; wp(ks)%rr=wp(ks)%rr+0.5_dp ! flag sire as ancestor and pass r to ks 
    wp(kd)%qq=.true.; wp(kd)%rr=wp(kd)%rr+0.5_dp ! ... and ditto kd, both will be in 1:np 
    na=ks+kd-kx ! the latest ancestor
    do ka=na,1,-1 ! add the ancestors to F calculation 
      if(.not.wp(ka)%qq) cycle
      rk=wp(ka)%rr
      fi=fi+wp(ka)%dd*rk**2.  ! contribution to F from ka
      kas=ped(ka)%js; kad=ped(ka)%jd ! parents of ka
      wp(kas)%qq=.true.; wp(kas)%rr=wp(kas)%rr+0.5_dp*rk ! passing r up pedigree to kas
      wp(kad)%qq=.true.; wp(kad)%rr=wp(kad)%rr+0.5_dp*rk ! ... and kd
      wp(ka)%qq=.false.; wp(ka)%rr=0._dp ! resetting ka elements
    end do
    wp(0)%rr=0._dp ! tidier to stop accumulation over identities, although does not enter into F
    ff(i)=fi ! set F
    if(mod(i,100000)==0) then
      call system_clock(t2)
      write(atxt,'(a,i10,a,i10)') '... ',i,' F calculations completed, increment in system time ',t2-t1
      call report_line([6,14],atxt) 
      t1=t2
    end if
  end do
  ped(:)%f=ff(1:)
end subroutine meu_luo
!
subroutine test_meu_luo(runit)
integer, intent(in) :: runit
type(info), dimension(:), pointer :: jpar
real(kind=dp), dimension(:),allocatable :: ff,dd
integer :: i,isz,np
  isz=13
  allocate(jpar(isz),ff(isz),dd(isz))
  jpar(:)=info(0,0,0,0,'NA',0.) ! F set to 0
  jpar(:)%js=[0,0,1,1,3,2,3,5,6,8,10,0,11] !    sire or dam's ID
  jpar(:)%jd=[0,0,2,2,4,4,5,6,6,9,0,10,12] 
  ff(:)=[0.,0.,0.,0.,0.25,0.25,0.375,0.3125,0.625,0.46875,0.,0.,0.18359375]
  dd(:)=[1.,1.,0.5,0.5,0.5,0.5,0.4375,0.375,0.375,0.265625,0.6328125,0.6328125,0.5]
  np=max(maxval(jpar(:)%js),maxval(jpar(:)%jd))
  call meu_luo(jpar,isz,np)
  write(runit,'(a)') '... ... Meuwissen & Luo test ...'
  do i=1,isz
    write(runit,'(a,3i3,3f12.8)') '... ... ... ',i,jpar(i)%js,jpar(i)%jd,jpar(i)%f,ff(i),dd(i)
  end do
  write(runit,'(a)') '... ... finished test'
  deallocate(jpar)
end subroutine test_meu_luo
!
end module jaw_meu_luo
! 
! ... the algorithm works because all ancestors are encountered in the correct order in the reversed loop
! ... an ancestor cannot be missed (%qq is false): suppose j is a missed ancestor
! ... then j cannot be sire and dam as they are explicitly set (%qq is true) before loop
! ... and j cannot be parents of sire or parents of dam, must be encountered after sire or dam,
! ... as lower index, and %qq is explicitly set to true in loops for the sire and the dam 
! ... therefore ancestor must have an offspring that is not the sire or dam and also has %qq set to false
! ... proof follows by recursion leading to a contradiction when unidentified ancestors exceed list spaces
!