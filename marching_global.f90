module domain
implicit none
include './resglo.h'

integer,parameter :: i01=i0+1,i0t1=i0t+1,j01=j0+1,j0t1=j0t+1,&
  nbg0=nb0*ngy0,nbg1=nbg0-1,ibir=i2t/ngy0

end module domain

module MPI_GRID
use mpi
implicit none

character, public :: grd*7
integer, public   :: mpi_comm_2d,mpi_comm_lon,mpi_comm_lat,nproc, myid,& 
  mylon, mylat,ierr,mpi_n,mpi_e,mpi_s,mpi_w, istat(MPI_STATUS_SIZE)
public :: mpi_grid_gen

contains
subroutine mpi_grid_gen
  use domain, only : ngx0, ngy0

  integer, dimension(2) ::  ndim
  logical, dimension(2) ::  peri,rdim1,rdim2
  integer :: n

  ndim(1)=ngx0
  ndim(2)=ngy0
  peri(1)=.true.
  peri(2)=.false.
  rdim1(1)=.true.
  rdim1(2)=.false.
  rdim2(1)=.false.
  rdim2(2)=.true.
  call mpi_cart_create(mpi_comm_world,2,ndim,peri,.true.,mpi_comm_2d,ierr)
  call mpi_cart_shift(mpi_comm_2d,0,1,mpi_w,mpi_e,ierr)
  call mpi_cart_shift(mpi_comm_2d,1,1,mpi_s,mpi_n,ierr)
  call mpi_cart_sub(mpi_comm_2d,rdim1,mpi_comm_lon,ierr)
  call mpi_cart_sub(mpi_comm_2d,rdim2,mpi_comm_lat,ierr)
  call mpi_comm_rank(mpi_comm_lon,mylon,ierr)
  call mpi_comm_rank(mpi_comm_lat,mylat,ierr)

  ! cartisean remark
  !  -------------------------------
  ! |   2   |   5   |   8   |   11  |
  ! | (0,2) | (1,2) | (2,2) | (3,2) |
  ! |-------|-------|-------|-------|
  ! |   1   |   4   |   7   |   10  |
  ! | (0,1) | (1,1) | (2,1) | (3,1) |
  ! |-------|-------|-------|-------|
  ! |   0   |   3   |   6   |   9   |
  ! | (0,0) | (1,0) | (2,0) | (3,0) |
  !  -------------------------------

  write(grd,'(i3.3,a1,i3.3)')mylon,'_',mylat
  if(myid.eq.0) write(*,*) '      myid         up           dwon           left    right'
  write(*,*) myid, mpi_n, mpi_s, mpi_w, mpi_e

end subroutine MPI_GRID_GEN

end module MPI_GRID



!================================================!
!               PREP for rep                     !
!================================================!
PROGRAM MAIN
use mpi
use domain
use mpi_grid

implicit none

real*8  :: rinv(ibir,i0,nbg0),rinv1(ibir,i0,nbg1),h(i0,j0)
integer :: ie(nb0)
real*8, dimension(i2,j2) :: al,ab,ac,ar
real*8, dimension(i2,0:j2) :: at
real*8 by,ay,bx,ax,dx,dy
real*8 time
integer :: n, nprocs, i,j


call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nprocs,ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)

if (nprocs .ne. ng0) then
  if(myid == 0) print*,'number of processors are not equal to ng0'
  call mpi_finalize(ierr)
  stop
endif

do n=1,nb0
  ie(n)=min(1+n0*n,j1)
end do
write(*,*) 'ie :', ie(:)
if (ie(nb0).le.ie(nb1).or.j1.ge.ie(nb1)+9) then
  if(myid.eq.0) write(*,252) ie(nb1),j1
  252  format('ie(nb1),j1=',2i5,' not properly defined. ' &
    ,'fix above do 250 loop and restart.') 
endif

call mpi_grid_gen

ie(nb0)=j1

by=1
bx=1
ax=0
ay=0      

dy=(by-ay)/dble(j0t)
dx=(bx-ax)/dble(i0t)

open(199,file='./log'//grd,form='formatted')

time=mpi_wtime()

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

call pre(al,ab,ac,ar,at,rinv,rinv1,h,ie)


open(99,file='./data/evp'//grd,form='unformatted')
rewind 99
write(99) al,ar,ab,at,ac,rinv,rinv1,ie
close(99)

if(myid == 0) then
  print*,'EVP PREP has finished successfully.'
endif
call solver

close(199)

call mpi_finalize(ierr)

end PROGRAM MAIN

subroutine pre(ax,ay,bb,cx,cy,rinv,rinv1,h,ie)
  use mpi  
  use domain
  use mpi_grid
  implicit none

  !INPUT VARIABLES
  integer, dimension(nb0),intent(in)        :: ie
  real*8, dimension(i2,j2),intent(in)      :: ax, ay, bb, cx
  real*8, dimension(i2,0:j2),intent(inout) :: cy

  !OUTPUT VARIABLES
  real*8, dimension(ibir, i0, nbg0),intent(inout) :: rinv
  real*8, dimension(ibir, i0, nbg1),intent(inout) :: rinv1
  real*8, dimension(i0,j0),intent(inout) :: h

  !LOCAL VARIABLES
  integer :: i,j,k,jl,jh,jg,jhm,ii,ig,jhp,it,jt,jss,jff,info,inmod
  integer :: n,nsd,ng,nb,nbg,idxy, nsend,nx,ny,ngt,ngxt,nt
  real*8,dimension(i2,j2)   :: axt,ayt,bbt,cxt,cyt
  real*8,dimension(i2t,i0t) :: rtemp
  real*8,dimension(i2t,i2t) :: rtp
  real*8,dimension(i0t*i0t) :: rp
  real*8,dimension(i2)      :: ct
  real*8,dimension(i2t)     :: ipvt
  real*8                    :: ctemp


  ct(:) = 0.d0
  do ng = 1, ngy0 !loop 100 : proceed by mpi_comm_lat
    jl =1
    if(mylat == ng-1) then
      do j = 1,j2
        do i = 1,i2
          axt(i,j)=ax(i,j)
          ayt(i,j)=ay(i,j)
          bbt(i,j)=bb(i,j)
          cxt(i,j)=cx(i,j)
          cyt(i,j)=cy(i,j)        
        end do
      end do
    endif


    if (mylat == ng-2) then !!! ??? boundary along y 
      do i =1, i2
        ct(i) = cy(i,j2)
      end do
    end if

    nsd=i2*j2
    call mpi_bcast(axt,nsd,mpi_real8,ng-1,mpi_comm_lat,ierr)
    call mpi_bcast(ayt,nsd,mpi_real8,ng-1,mpi_comm_lat,ierr)
    call mpi_bcast(bbt,nsd,mpi_real8,ng-1,mpi_comm_lat,ierr)
    call mpi_bcast(cxt,nsd,mpi_real8,ng-1,mpi_comm_lat,ierr)
    call mpi_bcast(cyt,nsd,mpi_real8,ng-1,mpi_comm_lat,ierr)
    if (ng /= 1)  &
      call mpi_bcast(ct,i2,mpi_real8,ng-2,mpi_comm_lat,ierr)

    if(mylat == ng-1) then
      do i = 1, i2
        cy(i,0) = ct(i)
      end do
    end if

    do nb = 1, nb0  !!! proceed by block on ng
      nbg = (ng-1)*nb0+nb !!! global block index
      if(myid == 0 )  write(*,111) nbg,nbg0
      111  format ('processing block #',i3,' out of ',i3,' total evp solver')
      jh=ie(nb)  !!! last line of #nb block, initial line of netx block
      jhp=jh+1   !!! up boundary of #nb block second(guess) line of #(nb+1) block
      jhm=jh-2
      jg=jl+1    !!! guess line

      !!!! different comm_lat processors handle different initial value case for block ng
      !!!! each comm_lat has i2t/ngy0 cases
      do ii =1, ibir
        ig=mylat*ibir+ii
        idxy=(ig-1)/i2

        do j=jl,jhp
          do i=1,i0
            h(i,j)=0.d0
          end do
        end do
        !!! set the only nonzero value on each comm_lat group
        if (idxy .eq. mylon) then
          inmod=mod(ig-1,i2)+2
          h(inmod,jg)=1.d0
        endif

        !!! exchange boundary
        call mpi_sendrecv(h(2,jg),1,mpi_real8,mpi_w,1,h(i0,jg),1, &
          mpi_real8,mpi_e,1,mpi_comm_world,istat,ierr)
        call mpi_sendrecv(h(i1,jg),1,mpi_real8,mpi_e,1,h(1,jg),1, &
          mpi_real8,mpi_w,1,mpi_comm_world,istat,ierr)

        !!! get coefficient value of the bottom boundary line
        !!! if nb not equal to 1, it comes from itself
        !!! otherwise, it comes from south processor
        if (nbg .ne. 1) then
          if (idxy .eq. mylon) then
            inmod = mod(ig-1,i2)+1
            if (nb .ne. 1) then
              ctemp=cyt(inmod,jg-2)
            else
              ctemp=ct(inmod)
            endif
          endif

          !!! form bottom boundary value according to the F-error of this
          !!!   boundary
          !!! RINV1 :
          call mpi_bcast(ctemp,1,mpi_real8,idxy,mpi_comm_lon,ierr)
          do n=1,i0
            h(n,jl)=rinv1(ii,n,nbg-1)*ctemp
          end do
        endif

        !!! marching from the line above guess line
        do j=jl,jhm
          do  i=1,i2
            h(i+1,j+2)=-(axt(i,j)*h(i,j+1)+ayt(i,j)*h(i+1,j)+ &
              bbt(i,j)*h(i+1,j+1)+cxt(i,j)*h(i+2,j+1))/cyt(i,j)
          end do 
          !!! exchange boundary everytime. Optimaztion needed
          call mpi_sendrecv(h(2,j+2),1,mpi_real8,mpi_w,1,h(i0,j+2),1, &
            mpi_real8,mpi_e,1,mpi_comm_world,istat,ierr)
          call mpi_sendrecv(h(i1,j+2),1,mpi_real8,mpi_e,1,h(1,j+2),1, &
            mpi_real8,mpi_w,1,mpi_comm_world,istat,ierr)
        end do


        !!! F error
        j=jh-1
        do i=1,i2
          rinv(ii,i+1,nbg)=axt(i,j)*h(i,j+1)+ayt(i,j)*h(i+1,j)+bbt(i,j)* &
            h(i+1,j+1)+cxt(i,j)*h(i+2,j+1)
        end do

        !!! save E error of last line
        if (nbg.ne.nbg0) then
          j=ie(nb)
          do n=1,i0
            rinv(ii,n,nbg0)=h(n,j)
          end do
          if (idxy .eq. mylon) then
            inmod = mod(ig-1,i2) +2
            h(inmod,jg) = 0.d0
          endif
        endif

      enddo


      !!! gather rinv for each block
      nsend=i2*ibir

      call mpi_gather(rinv(1,2,nbg),nsend,mpi_real8,rp(1),nsend, &
        mpi_real8,0,mpi_comm_world,ierr)


      !!! reshape rinv
      if (myid .eq. 0) then
        do  nx=1,ngx0
          ngxt=(nx-1)*ngy0
          do  ny=1,ngy0
            ngt=((ny-1)+ngxt)*ibir*i2
            do  j=1,i2
              nt=ngt+(j-1)*ibir
              jt=j+(nx-1)*i2
              do  i=1,ibir
                it=i+(ny-1)*ibir
                rtemp(it,jt)=rp(nt+i)
              end do
            end do
          end do
        end do
        do i=1,i2t
          do j=1,i2t
            rtp(j,i)=0.d0
          end do
        end do

        do i=1,i2t
          rtp(i,i)=1.d0
        end do
        call dgesv(i2t,i2t,rtemp,i2t,ipvt,rtp,i2t,info)
        do i=1,i2t
          do j=1,i2t
            rtemp(j,i)=rtp(i,j)
          end do
        enddo
      endif

      call mpi_bcast(info,1,mpi_integer,0,mpi_comm_world,ierr)
      if (info .ne. 0) then
        if(myid == 0) print*,'preporcessor failed because rinv is singular at nbg =',nbg
        call mpi_finalize(ierr)
        stop
      endif
      nsend=i2t*i2t
      call mpi_bcast(rtp,nsend,mpi_real8,0,mpi_comm_world,ierr)

      it=mylon*i2
      jt=mylat*ibir
      jss=1
      jff=i0
      if (mylon .eq.0) then
        do  i=1,ibir
          rinv(i,1,nbg)=rtp(jt+i,i2t)
        enddo
        jss=2
      endif

      if (mylon .eq. ngx1) then
        do i=1,ibir
          rinv(i,i0,nbg)=rtp(jt+i,1)
        enddo
        jff=i1
      endif

      do  j=jss,jff
        do  i=1,ibir
          rinv(i,j,nbg)=rtp(jt+i,it+j-1)
        enddo
      enddo
      if (nbg.eq.nbg0) return


      nsend=i0*ibir
      call mpi_allgather(rinv(1,1,nbg0),nsend,mpi_real8,rp,nsend, &
        mpi_real8,mpi_comm_lat,ierr)
      do  ny=1,ngy0
        ngt=(ny-1)*ibir*i0
        do  j=1,i0
          nt=ngt+(j-1)*ibir
          do  i=1,ibir
            it=i+(ny-1)*ibir
            rtemp(it,j)=rp(nt+i)
          enddo
        enddo
      enddo

      !!! RINV1 : from F-error of last line to its Error
      it=ibir*mylat
      do  i=1,ibir
        do  j=1,i0
          rinv1(i,j,nbg)=0.d0
          do  k=1,i2t
            rinv1(i,j,nbg)=rinv1(i,j,nbg)-rtp(it+i,k)*rtemp(k,j)
          enddo
        enddo
      enddo
      jl=jh

      call mpi_barrier(mpi_comm_world,ierr)

    end do 
  end do ! loop 100
end subroutine pre

subroutine pre_check(ax,ay,bb,cx,cy,rinv,rinv1,ie)
  use domain
  use mpi_grid
  implicit none

  !INPUT VARIABLES
  integer, dimension(nb0),intent(in)        :: ie
  real*8, dimension(i2,j2),intent(in)      :: ax, ay, bb, cx
  real*8, dimension(i2,0:j2),intent(inout) :: cy

  !OUTPUT VARIABLES
  real*8, dimension(ibir, i0, nbg0),intent(inout) :: rinv
  real*8, dimension(ibir, i0, nbg1),intent(inout) :: rinv1

  !LOCAL VARIABLES
  
  !!! TO DO


end subroutine pre_check

subroutine solver
  use mpi  
  use domain
  use mpi_grid
  implicit none
  real*8 :: rinv(ibir,i0,nbg0),rinv1(ibir,i0,nbg1)
  real*8 :: dum0(i0,nb0), dum1(ibir),dum2(i0),dumg(i2t)
  real*8 :: f(i2,j2),x(i0,j0)
  real*8 :: al(i2,j2),ab(i2,j2),ac(i2,j2),ar(i2,j2),at(i2,0:j2)
  integer:: ie(nb0)

  integer :: i,j
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

  open(99,file='./data/evp'//grd,form='unformatted')
  if(myid.eq.0) write(*,*)'evp solver'
  read(99) al,ar,ab,at,ac,rinv,rinv1,ie
  close(99)
  ! test case 

  time=mpi_wtime()
  call rep(al,ab,ac,ar,at,rinv,rinv1,f,x,ie)
  call rep_check(al,ab,ac,ar,at,f,x)
  
  time=mpi_wtime()-time
  call mpi_barrier(mpi_comm_world,ierr)
  print*,'myid', myid, 'runtime =', time/1000
  
  open(99,file='./glo'//grd//'.out',form='formatted')
  write(99,*) x(2:i1,2:j1)
  close(99)


end subroutine solver



! --------------------------REP2 from Global.F--------------------------------------------
subroutine rep(ax,ay,bb,cx,cy,rinv,rinv1,f,x,ie) 
  ! ----------------------------------------------------------------------
  use mpi  
  use domain
  use mpi_grid
  implicit none
  !!!INPUT 
  real*8, dimension(i2,j2), intent(in) :: ax,ay,bb,cx,f
  real*8, intent(in) :: rinv(ibir,i0,nbg0),rinv1(ibir,i0,nbg1)
  integer,intent(in) :: ie(nb0)
  real*8, intent(inout) :: cy(i2,0:j2),x(i0,j0)

  !!!LOCAL
  real*8  :: dum0(i0,nb0),dum1(ibir),dum2(i0),dumg(i2t)
  integer :: i,j,jsp,jfp
  integer :: n,m,ng,nb,nbs,nbg
  real*8  :: dd


  do  ng=0,ngy0-1
    jsp=1
    do nb=1,nb0
      nbg=ng*nb0+nb
      if (mylat .eq. ng) then
        jfp=ie(nb)-2
        do  j=jsp,jfp
          do  i=1,i2
            x(i+1,j+2)=(f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* &
              x(i+1,j+1)-cx(i,j)*x(i+2,j+1))/cy(i,j)
          enddo
          call mpi_sendrecv(x(2,j+2),1,mpi_real8,mpi_w,1,x(i0,j+2),1, &
            mpi_real8,mpi_e,1,mpi_comm_2d,istat,ierr)
          call mpi_sendrecv(x(i1,j+2),1,mpi_real8,mpi_e,1,x(1,j+2),1,&
            mpi_real8,mpi_w,1,mpi_comm_2d,istat,ierr)
        enddo
        if (nbg.ne.nbg0) then
          j=ie(nb)-1
          do i=1,i2
            dum2(i)=f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* &
              x(i+1,j+1)-cx(i,j)*x(i+2,j+1)-cy(i,j)*x(i+1,j+2)
          enddo
          j=ie(nb)
          call mpi_allgather(dum2,i2,mpi_real8,dumg,i2,mpi_real8, &
            mpi_comm_lon,ierr)
          do n=1,i0
            dum2(n)=x(n,j)
            dum0(n,nb)=x(n,j)
          end do
          dd=1.d0
        endif
      else
        if (nbg .ne. nbg0) then
          dd=0.d0
        endif
      endif
      if (nbg.ne.nbg0) then
        call mpi_scatter(dumg,ibir,mpi_real8,dum1,ibir,mpi_real8,ng, &
          mpi_comm_lat,ierr)
        call dgemv('t',ibir,i0,-1.d0,rinv1(1,1,nbg),ibir,dum1,1,dd,dum2,1)
        if (nb .ne. nb0) then
          call mpi_reduce(dum2,x(1,ie(nb)),i0,mpi_real8,mpi_sum,ng, &
            mpi_comm_lat,ierr)
          jsp=ie(nb)
        endif
      endif

    enddo 
    if (nbg.ne.nbg0) then
      call mpi_reduce(dum2,x(1,1),i0,mpi_real8,mpi_sum,ng+1, &
        mpi_comm_lat,ierr)
    endif
  enddo

  do ng=ngy0-1,0,-1
    do  nbs=1,nb0
      nb=nb0-nbs+1
      nbg=nb0*ng+nb

      jsp=1
      if (nb.ne.1) jsp=ie(nb-1)
      if (mylat .eq. ng) then
        if (nbg.ne.nbg0) then
          j=ie(nb)
          do  n=1,i0
            x(n,j)=dum0(n,nb)
          enddo
        endif
        j=ie(nb)-1
        do  i=1,i2
          dum2(i)=f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* &
            x(i+1,j+1)-cx(i,j)*x(i+2,j+1)-cy(i,j)*x(i+1,j+2)
        enddo
        call mpi_allgather(dum2,i2,mpi_real8,dumg,i2,mpi_real8,&
          mpi_comm_lon,ierr)
      endif

      call mpi_scatter(dumg,ibir,mpi_real8,dum1,ibir,mpi_real8,ng, &
        mpi_comm_lat,ierr)
      call dgemv('t',ibir,i0,1.d0,rinv(1,1,nbg),ibir,dum1,1,0.d0,dum2,1)
      call mpi_reduce(dum2,dum0(1,nb),i0,mpi_real8,mpi_sum,ng, &
        mpi_comm_lat,ierr)

      if (mylat .eq. ng) then
        do n=1,i0
          x(n,jsp+1)=x(n,jsp+1)+dum0(n,nb)
        enddo
      endif
      if (nbg.ne.1) then
        if (mylat .eq. ng) then
          do m=1,i2
            dum2(m)=dum0(m+1,nb)*cy(m,jsp-1)
          enddo
          call mpi_allgather(dum2,i2,mpi_real8,dumg,i2,mpi_real8,&
            mpi_comm_lon,ierr)
        endif
        call mpi_scatter(dumg,ibir,mpi_real8,dum1,ibir,mpi_real8,ng,&
          mpi_comm_lat,ierr)
        call dgemv('t',ibir,i0,1.d0,rinv1(1,1,nbg-1),ibir,dum1,&
          1,0.d0,dum2,1)
        call mpi_reduce(dum2,dum0(1,nb),i0,mpi_real8,mpi_sum,ng,&
          mpi_comm_lat,ierr)
        if (mylat .eq. ng) then
          do n=1,i0
            dum0(n,nb)=x(n,jsp)+dum0(n,nb)
          enddo
        endif
      endif
    enddo
    if (nbg.ne.1) then
      if (mylat .eq. ng) &
        call mpi_send(x(1,2),i0,mpi_real8,mpi_s,1,mpi_comm_2d,ierr)
      if (mylat .eq. ng-1) &
        call mpi_recv(x(1,j0),i0,mpi_real8,mpi_n,1,mpi_comm_2d,istat,ierr)
    endif
  enddo
  if (mylat .eq. 0) then
    do i=1,i0
      dum0(i,1)=0.d0
    enddo
  endif

  jsp=1
  do nb=1,nb0
    jfp=ie(nb)-2
    if (nb.eq.nb0) jfp=ie(nb)-2
    do i=1,i0
      x(i,jsp)=dum0(i,nb)
    enddo
    do j=jsp,jfp
      do i=1,i2
        x(i+1,j+2)=(f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* &
          x(i+1,j+1)-cx(i,j)*x(i+2,j+1))/cy(i,j)
      enddo
      call mpi_sendrecv(x(2,j+2),1,mpi_real8,mpi_w,1,x(i0,j+2),1, &
        mpi_real8,mpi_e,1,mpi_comm_2d,istat,ierr)
      call mpi_sendrecv(x(i1,j+2),1,mpi_real8,mpi_e,1,x(1,j+2),1, &
        mpi_real8,mpi_w,1,mpi_comm_2d,istat,ierr)
    enddo
    jsp=ie(nb)
  enddo 

end subroutine rep

subroutine rep_check(ax,ay,bb,cx,cy,f,x) 
  ! ----------------------------------------------------------------------
  use domain
  use mpi_grid
  implicit none
  !!!INPUT 
  real*8, dimension(i2,j2), intent(in) :: ax,ay,bb,cx,f
  real*8, dimension(i2,0:j2), intent(in) :: cy
  real*8, intent(in) :: x(i0,j0)

  !!!LOCAL
  integer :: i,j,flag
  real*8  :: rtemp
  flag = 0
  do j = 1, j2
    do i = 1, i2
      rtemp=f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* &
        x(i+1,j+1)-cx(i,j)*x(i+2,j+1)-cy(i,j)*x(i+1,j+2)
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
