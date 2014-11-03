module domain
use mpi
implicit none
include './resglo1.h'

integer, parameter :: i1=i0-1,i2=i0-2,i3=i0-3,& 
                      j1=j0-1, j2=j0-2,j3=j0-3,nb1=nb0-1
!!! Define domain
real*8,parameter  :: by =1.,ay=0.,bx=1.,ax=0.
!!!intensive used variables
integer :: myid,ierr,istat(mpi_status_size)
integer :: left, right
character*50 :: fname


end module domain


!================================================!
!               PREP for rep                     !
!================================================!
PROGRAM MAIN
use mpi
use domain

implicit none

!!! DEFINE BLOCK PARTITION
integer, dimension(nb0)  :: id,ip,im,js,jf
!!! Coefficent matrix
real*8, dimension(i2,j2) :: al,ab,ac,ar,at
!!! Coefficent matrix
real*8                   :: ccol(nbg0*j2,nb0*j2)
real*8                   ::time
integer                  :: n, nprocs,i,j


call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nprocs,ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)

if (nprocs .ne. ngrid) then
  if(myid == 0) print*,'number of processors are not equal to ngrid'
  call mpi_finalize(ierr)
  stop
endif

left=myid-1
right=myid+1
if (myid .eq. 0) left=ngrid-1
if (myid .eq. ngrid-1) right=0

!!! Prepocessing 
time=mpi_wtime()
call pre(al,ab,ac,ar,at,ccol,id,ip,im,js,jf)
write(*,*) 'myid',myid, 'id2 ',id
time=mpi_wtime()-time
print*,'myid', myid, 'pretime =', time/1000

!!! Output efficient matrixs

if(myid == 0) then
  print*,'EVP PREP has finished successfully.'
endif

call mpi_barrier(mpi_comm_world,ierr)

write(*,*) 'myid',myid, 'id3 ',id
!!! Solver equations with EVP 
call solver(al,ab,ac,ar,at,ccol,id,ip,im,js,jf)


call mpi_finalize(ierr)

end PROGRAM MAIN

subroutine pre(al,ab,ac,ar,at,ccol,id,ip,im,js,jf)
  use mpi  
  use domain
  implicit none

  !INPUT VARIABLES
  real*8, dimension(i2,j2),intent(inout) :: al,ab,ac,ar,at
  real*8, intent(inout) :: ccol(nbg0*j2,nb0*j2)
  integer, dimension(nb0),intent(inout) :: id,ip,im,js,jf
  !OUTPUT VARIABLES

  !LOCAL VARIABLES
  integer   :: ipvt(nbg0*j2)
  integer :: i,j,k,l,jj,ii1,ii2,info
  integer :: nb,nny,nsend
  real*8 :: binv(nbg0*j2,nbg0*j2),cinv(nbg0*j2,nbg0*j2), &
        etmp(j2,j2,2,nb0,2),ebuf(j2,j2,4),ehat(j2,j2,i0,2)
  real*8  :: fw, bw
  real*8 dx,dy

  dy=(by-ay)/dble(i0*nbg0)
  dx=(bx-ax)/dble(j0)


  do j=1,j2
    do i=1,i2
      at(i,j)=1/dy
      ab(i,j)=1/dy
      ar(i,j)=1/dx
      al(i,j)=1/dx
      ac(i,j)=-2*(1/dx+1/dy) 
    end do 
  end do 
!   ---------------------
!          AT
!          |
!   AL --- AC ---  AR
!          |
!          AB
!   ---------------------


  do nb=1,nb0
    id(nb)=(-1)**(nb+1)
    ip(nb)=mod(nb,2)
    im(nb)=mod(nb+1,2)
  enddo



  !determined the step of the marching size

  nny=i1/nb0
  if (nny .le. 2) then
    print*,'nblk=',nb0,'nnx=',nny
    print*,'nblk is so small that evp can''t works. prep stops'
    call mpi_finalize(ierr)
    elseif (nny .gt. 10) then
    print*,'nblk=',nb0,'nnx=',nny
    print*,'nblk is so large that evp can''t works. prep stops'
    call mpi_finalize(ierr)
  endif 
  print*,'marching steps =',nny


  !c       determine start line of each block. forward & backward
  do nb=1,nb0-1
    js(nb)=nny*(nb-ip(nb))+1+ip(nb)
  enddo
  js(nb0)=i2

  do nb=1,nb0
    jf(nb)=(js(nb)+js(nb+id(nb))-id(nb))/2
  enddo

  write(fname,'(a6,i2.2,a6)')'idipim',myid,'.90log'
  open(99,file=fname,form='formatted')
  rewind 99
  write(99,*) id,ip,im,js,jf
  close(99)

  do nb=1,nb0
    do j=1,j2
      ehat(j,j,js(nb),ip(nb)+1)=1.d0
    enddo 
  enddo

  do j=1,j2
    ehat(j,j,i1,2)=1.d0
    ehat(j,j,1,1)=1.d0
  enddo

  do nb=1,nb0
    do l=1,2
      do i=js(nb),jf(nb)-id(nb),id(nb)
        do jj=1,j2
          do j=2,j2-1
            fw=real(ip(nb))*ar(i-1,j)+real(im(nb))*al(i-1,j)
            bw=real(ip(nb))*al(i-1,j)+real(im(nb))*ar(i-1,j)
            ehat(j,jj,i+id(nb),l)=(-ac(i-1,j)*ehat(j,jj,i,l)-at(i-1,j)*ehat(j+1,jj,i,l) &
              -ab(i-1,j)*ehat(j-1,jj,i,l)-bw*ehat(j,jj,i-id(nb),l))/fw
          enddo

          fw=real(ip(nb))*ar(i-1,1)+real(im(nb))*al(i-1,1);
          bw=real(ip(nb))*al(i-1,1)+real(im(nb))*ar(i-1,1);
          ehat(1,jj,i+id(nb),l)=(-ac(i-1,1)*ehat(1,jj,i,l)-at(i-1,1)* &
            ehat(2,jj,i,l)-bw*ehat(1,jj,i-id(nb),l))/fw;
          fw=real(ip(nb))*ar(i-1,j2)+real(im(nb))*al(i-1,j2);
          bw=real(ip(nb))*al(i-1,j2)+real(im(nb))*ar(i-1,j2);
          ehat(j2,jj,i+id(nb),l)=(-ac(i-1,j2)*ehat(j2,jj,i,l)-ab(i-1,j2)* &
            ehat(j2-1,jj,i,l)-bw*ehat(j2,jj,i-id(nb),l))/fw;
        enddo
      enddo


      do jj=1,j2
        do j=2,j2-1 
          etmp(j,jj,1,nb,l)=ehat(j,jj,jf(nb),l)
          fw=real(ip(nb))*ar(jf(nb)-1,j)+real(im(nb))*al(jf(nb)-1,j)
          bw=real(ip(nb))*al(jf(nb)-1,j)+real(im(nb))*ar(jf(nb)-1,j)
          etmp(j,jj,2,nb,l)=(-ac(jf(nb)-1,j)*ehat(j,jj,jf(nb),l) &
            -at(jf(nb)-1,j)*ehat(j+1,jj,jf(nb),l)-ab(jf(nb)-1,j)* &
            ehat(j-1,jj,jf(nb),l)-bw*ehat(j,jj,jf(nb)-id(nb),l))/fw
        enddo
        etmp(1,jj,1,nb,l)=ehat(1,jj,jf(nb),l)
        fw=real(ip(nb))*ar(jf(nb)-1,1)+real(im(nb))*al(jf(nb)-1,1)
        bw=real(ip(nb))*al(jf(nb)-1,1)+real(im(nb))*ar(jf(nb)-1,1)
        etmp(1,jj,2,nb,l)=(-ac(jf(nb)-1,j)*ehat(1,jj,jf(nb),l) &
          -at(jf(nb)-1,j)*ehat(2,jj,jf(nb),l) &
          -bw*ehat(1,jj,jf(nb)-id(nb),l))/fw;
        etmp(j2,jj,1,nb,l)=ehat(j2,jj,jf(nb),l)
        fw=real(ip(nb))*ar(jf(nb)-1,j2)+real(im(nb))*al(jf(nb)-1,j2)
        bw=real(ip(nb))*al(jf(nb)-1,j2)+real(im(nb))*ar(jf(nb)-1,j2)
        etmp(j2,jj,2,nb,l)=(-ac(jf(nb)-1,j2)*ehat(j2,jj,jf(nb),l) &
          -ab(jf(nb)-1,j2)*ehat(j2-1,jj,jf(nb),l) &
          -bw*ehat(j2,jj,jf(nb)-id(nb),l))/fw;  
      enddo

    enddo
  enddo

  nsend=j2*j2
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_sendrecv(etmp(1,1,1,nb0,2),nsend,mpi_real8,right, &
    1,ebuf(1,1,1),nsend,mpi_real8,left,1,mpi_comm_world,istat,ierr)
  call mpi_sendrecv(etmp(1,1,2,nb0,2),nsend,mpi_real8,right, &
    1,ebuf(1,1,2),nsend,mpi_real8,left,1,mpi_comm_world,istat,ierr)
  call mpi_sendrecv(etmp(1,1,1,1,1),nsend,mpi_real8,left,1, &
    ebuf(1,1,3),nsend,mpi_real8,right,1,mpi_comm_world,istat,ierr)
  call mpi_sendrecv(etmp(1,1,2,1,1),nsend,mpi_real8,left,1, &
    ebuf(1,1,4),nsend,mpi_real8,right,1,mpi_comm_world,istat,ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  do ii1=1,j2
    do ii2=1,j2
      etmp(ii1,ii2,1,1,1)=ebuf(ii1,ii2,1)
      etmp(ii1,ii2,2,1,1)=ebuf(ii1,ii2,2)
      etmp(ii1,ii2,1,nb0,2)=ebuf(ii1,ii2,3)
      etmp(ii1,ii2,2,nb0,2)=ebuf(ii1,ii2,4)     
    enddo
  enddo

  l=(myid*nb0)*j2
  k=0
  call arr(nbg0,nb0,j2,etmp(1,1,1,1,2),ccol,l,k,-1.d0)
  call arr(nbg0,nb0,j2,etmp(1,1,2,1,2),ccol,l+j2,k,-1.d0)
  if (myid .eq. 0) then
    l=(nbg0-2)*j2
    !      call arr(nbg0,nb0,j2,etmp(1,1,2,1,1),ccol,l,k,-1.d0)
    !      call arr(nbg0,nb0,j2,etmp(1,1,1,1,1),ccol,l+j2,k,-1.d0)
  else
    call arr(nbg0,nb0,j2,etmp(1,1,2,1,1),ccol,l-2*j2,k,1.d0)
    call arr(nbg0,nb0,j2,etmp(1,1,1,1,1),ccol,l-j2,k,1.d0)
  endif

  if (nb0 .gt. 2) then
    do nb=2,nb0-1
      l=((myid*nb0)+(nb-1)+im(nb))*j2
      k=(nb-1)*j2
      call arr(nbg0,nb0,j2,etmp(1,1,2,nb-ip(nb),ip(nb)+1), &
        ccol,l-2*j2,k,1.d0)
      call arr(nbg0,nb0,j2,etmp(1,1,1,nb-ip(nb),ip(nb)+1), &
        ccol,l-1*j2,k,1.d0)
      call arr(nbg0,nb0,j2,etmp(1,1,1,nb+im(nb),ip(nb)+1), &
        ccol,l,k,-1.d0)
      call arr(nbg0,nb0,j2,etmp(1,1,2,nb+im(nb),ip(nb)+1), &
        ccol,l+j2,k,-1.d0)
    enddo
  endif

  l=((myid+1)*nb0)*j2
  k=(nb0-1)*j2
  call arr(nbg0,nb0,j2,etmp(1,1,2,nb0,1),ccol,l-2*j2,k,1.d0)
  call arr(nbg0,nb0,j2,etmp(1,1,1,nb0,1),ccol,l-j2,k,1.d0)
  if (myid .eq. 3) then
    !c      call arr(nbg0,nb0,j2,etmp(1,1,1,nb0,2),ccol,0,k,1.d0)
    !c      call arr(nbg0,nb0,j2,etmp(1,1,2,nb0,2),ccol,j2,k,1.d0)
  else
    call arr(nbg0,nb0,j2,etmp(1,1,1,nb0,2),ccol,l,k,-1.d0)
    call arr(nbg0,nb0,j2,etmp(1,1,2,nb0,2),ccol,l+j2,k,-1.d0)
  endif
  nsend=nb0*nbg0*j2*j2
  call mpi_barrier(mpi_comm_world,ierr)
  print*,myid
  call mpi_gather(ccol(1,1),nsend,mpi_real8,cinv(1,1),nsend, &
    mpi_real8,0,mpi_comm_world,ierr)
  if (myid .eq. 0) then
    do i=1,nbg0*j2
      binv(i,i)=1.d0
    enddo
    call dgesv(nbg0*j2,nbg0*j2,cinv,nbg0*j2,ipvt,binv,nbg0*j2,info)
    do i=1,nbg0*j2
      do j=1,nbg0*j2
        cinv(j,i)=binv(i,j)
      enddo 
    enddo
  endif

  !c      if (.false.) then

  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(info,1,mpi_integer,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  !c      if (.false.) then
  if (info .ne. 0) then
    print*,'the influnce matrix is singular.'  
  endif
  nsend=nbg0*nb0*j2*j2
  call mpi_scatter(cinv(1,1),nsend,mpi_real8,ccol(1,1),nsend, &
    mpi_real8,0,mpi_comm_world,ierr)
  write(fname,'(a3,i2.2,a6)')'inf',myid,'.90out'
  open(99,file=fname,form='unformatted')
  rewind 99
  write(99) al,ab,ac,ar,at,ccol,id,ip,im,js,jf
  close(99)
  print*,'pre has finished successfully.'

end subroutine pre

subroutine pre_check(al,ab,ac,ar,at,ccol,id,ip,im,js,jf)
  use domain
  implicit none

  !INPUT VARIABLES
  real*8, dimension(i2,j2),intent(in) :: al,ab,ac,ar,at
  real*8, intent(in) :: ccol(nbg0*j2,nb0*j2)
  integer, dimension(nb0),intent(in) :: id,ip,im,js,jf

  !OUTPUT VARIABLES

  !LOCAL VARIABLES
  
  !!! TO DO


end subroutine pre_check

subroutine solver(al,ab,ac,ar,at,ccol,id,ip,im,js,jf)
  use mpi  
  use domain
  implicit none
  !INPUT VARIABLES
  real*8, dimension(i2,j2),intent(in) :: al,ab,ac,ar,at
  real*8, intent(in) :: ccol(nbg0*j2,nb0*j2)
  integer, dimension(nb0),intent(in) :: id,ip,im,js,jf

  !LOCAL VARIABLES
  integer :: i,j
  real*8 :: f(i2,j2),x(i0,j0)
  real*8 :: time



  ! set initial and left-hand value
  do j=2,j1
    do i=2,i1
      x(i,j)=0.d0
      f(i-1,j-1)=1
    enddo 
  enddo

  do j=1,j0
    x(i0,j)=0.d0
    x(1,j)=0.d0
  enddo

  do i=1,i0
    x(i,1)=0.d0
    x(i,j0)=0.d0
  enddo

  !write(fname,'(a3,i3,a6)')'inf',myid,'.90out'
  !open(99,file=fname,form='unformatted')
  !read(99) al,ar,ab,at,ac,ccol,ehat,ip,id,im,js,jf
  !close(99)
  ! test case 

  time=mpi_wtime()
  call rep(al,ab,ac,ar,at,ccol,id,ip,im,js,jf,f,x)
  call rep_check(al,ab,ac,ar,at,f,x)
  
  time=mpi_wtime()-time
  call mpi_barrier(mpi_comm_world,ierr)
  print*,'myid', myid, 'runtime =', time/1000
  


end subroutine solver



! --------------------------REP2 from Global.F--------------------------------------------
subroutine rep(al,ab,ac,ar,at,ccol,id,ip,im,js,jf,f,x) 
  ! ----------------------------------------------------------------------
  use mpi  
  use domain
  implicit none
  !INPUT VARIABLES
  integer, dimension(nb0),intent(in) :: id,ip,im,js,jf
  real*8, dimension(i2,j2),intent(in) :: al,ab,ac,ar,at
  real*8, intent(in) :: ccol(nbg0*j2,nb0*j2)
  real*8, intent(in) :: f(i2,j2)
  !OUTPUT VARIABLES
  real*8, intent(inout):: x(i0,j0)

  !!!LOCAL
  integer:: nny, nsend
  real*8 :: etmp(j2,j2,2,nb0,2),ebuf(j2,j2,4)
  real*8 :: r(nb0*j2),rb(nbg0*j2)
  real*8 :: fw,bw
  integer :: nb,i,j,jl 


  !!! for read
  write(fname,'(a3,i2.2,a6)')'inf',myid,'.90log'
  open(99,file=fname,form='formatted')
  rewind 99
  write(99,*) al,ab,ac,ar,at,ccol,id,ip,im,js,jf
  close(99)


  do nb=1,nb0
    !write(*,*) 'js jf id ',js(nb),jf(nb)-id(nb),id(nb) 
    do i=js(nb),jf(nb)-id(nb),id(nb)
      do j=1,j2
        fw=ip(nb)*ar(i-1,j)+im(nb)*al(i-1,j)
        bw=ip(nb)*al(i-1,j)+im(nb)*ar(i-1,j)
        x(i+id(nb),j+1)=(f(i-1,j)-ac(i-1,j)*x(i,j+1)-at(i-1,j)*x(i,j+2) &
          -ab(i-1,j)*x(i,j)-bw*x(i-id(nb),j+1))/fw
      enddo
    enddo
    jl=(nb+id(nb)-1)*j2
    do j=1,j2
      fw=ip(nb)*ar(jf(nb)-1,j)+im(nb)*al(jf(nb)-1,j)
      bw=ip(nb)*al(jf(nb)-1,j)+im(nb)*ar(jf(nb)-1,j)
      r(j+jl)=(f(jf(nb)-1,j)-ac(jf(nb)-1,j)*x(jf(nb),j+1) &
        -at(jf(nb)-1,j)*x(jf(nb),j+2)-ab(jf(nb)-1,j)*x(jf(nb),j) &
        -bw*x(jf(nb)-id(nb),j+1))/fw*(dble(id(nb)))
    enddo
  enddo

  do nb=1,nb0
    jl=(nb-1)*j2
    do j=1,j2
      r(j+jl)=dble(id(nb))*x(jf(nb),j+1)+r(j+jl)
    enddo 
  enddo


  nsend=j2*nb0
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_allgather(r(1),nsend,mpi_real8, &
    rb(1),nsend,mpi_real8,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)

  do j=1,nb0*j2
    r(j)=0.d0
    do i=1,nbg0*j2
      r(j)=rb(i)*ccol(i,j)+r(j)
    enddo 
  enddo

  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_allgather(r(1),nsend,mpi_real8, &
    rb(1),nsend,mpi_real8,mpi_comm_world,ierr)

  do nb=1,nb0
    do j=1,j2
      x(js(nb),j+1)=x(js(nb),j+1)+rb((myid*nb0+nb-1)*j2+j)
    enddo 
  enddo

  !c for dbc

  if (myid .gt. 0) then
    do j=1,j2
      x(js(1)-1,j+1)=x(js(1)-1,j+1)+rb(((left+1)*nb0-1)*j2+j)
    enddo
  endif      
  if (myid .lt. ngrid-1) then
    do j=1,j2
      x(js(nb0)+1,j+1)=x(js(nb0)+1,j+1)+rb((right*nb0)*j2+j)
    enddo
  endif

  do nb=1,nb0
    do i=js(nb),jf(nb)-id(nb),id(nb)
      do j=1,j2
        fw=ip(nb)*ar(i-1,j)+im(nb)*al(i-1,j)
        bw=ip(nb)*al(i-1,j)+im(nb)*ar(i-1,j)
        x(i+id(nb),j+1)=(f(i-1,j)-ac(i-1,j)*x(i,j+1)-at(i-1,j)*x(i,j+2) &
          -ab(i-1,j)*x(i,j)-bw*x(i-id(nb),j+1))/fw
      enddo
    enddo
  enddo


  write(fname,'(a3,i2.2,a6)')'ans',myid,'.90out'
  open(1233,file=fname)
  write(1233,*),((x(i,j),i=2,i2),j=2,j1)
  close(1233)


end subroutine rep

subroutine rep_check(al,ab,ac,ar,at,f,x) 
  ! ----------------------------------------------------------------------
  use domain
  implicit none
  !!!INPUT 
  real*8, dimension(i2,j2), intent(in) :: al,ab,ac,ar,at,f
  real*8, intent(in) :: x(i0,j0)

  !!!LOCAL
  integer :: i,j,flag
  real*8  :: rtemp
  flag = 0
  do j = 1, j2
    do i = 1, i2
      rtemp=f(i,j)-al(i,j)*x(i,j+1)-ab(i,j)*x(i+1,j)-ac(i,j)* &
        x(i+1,j+1)-ar(i,j)*x(i+2,j+1)-at(i,j)*x(i+1,j+2)
      if (rtemp > 1.0e-6 )  then
        write(*,*) myid, i, j, rtemp
        flag = flag +1
      endif
    enddo
  enddo
  if (flag == 0 )  then
    write(*,*) 'PRE_CHECK PASS on ', myid
  else 
    write(*,*) 'PRE_CHECK FAIL on ', myid, 'flag ', flag
  endif


end subroutine rep_check
!c     submatrix to full matrix
subroutine arr(n1,n2,m,c,a,ix,iy,d)
  implicit none
  !!!INPUT
  integer,intent(in) ::  n1,n2,m,ix,iy
  real*8, intent(in) ::  d
  real*8, intent(in) :: c(m,m)
  real*8, intent(inout) :: a(m*n1,m*n2)
  !!!LOCAL 
  integer i,j

  do i=1,m
    do j=1,m
      a(ix+i,iy+j)=d*c(i,j)
    enddo
  enddo
end subroutine arr

