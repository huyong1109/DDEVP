module MPI_GRID
use mpi
implicit none

character, public :: grd*7
integer, public   :: mpi_comm_2d,mpi_comm_lon,mpi_comm_lat,nproc,myid, & 
  mylon, mylat,ierr,mpi_n,mpi_e,mpi_s,mpi_w, istat(MPI_STATUS_SIZE)

public :: mpi_grid_gen

contains
subroutine mpi_grid_gen

  include './resglo.h90'
  integer, dimension(2) ::  ndim
  logical, dimension(2) ::  peri,rdim1,rdim2

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
  call mpi_barrier(mpi_comm_world,ierr)
  if (myid == 0) then
    write(*,*) '      myid         up           dwon           left    right'
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  write(*,*) myid, mpi_n, mpi_s, mpi_w, mpi_e
  call mpi_barrier(mpi_comm_world,ierr)

end subroutine MPI_GRID_GEN

end module MPI_GRID

module domain
  implicit none
  include './resglo.h90'
  
  integer,parameter :: i01=i0+1,i0t1=i0t+1,j01=j0+1,j0t1=j0t+1,&
    nbg0=nb0*ngy0,nbg1=nbg0-1,ibir=i2t/ngy0

end module domain


!================================================!
!               PREP for rep                     !
!================================================!
PROGRAM MAIN
  use mpi_grid
  use mpi
  use domain
  
  implicit none
  
  
  real*8            :: rinv(ibir,i0,nbg0),rinv1(ibir,i0,nbg1),h(i0,j0),ie(nb0)
  real*8, dimension(i2,j2) :: al,ab,ac,ar
  real*8, dimension(i2,0:j2) :: at
  real*8 by,ay,bx,ax,dx,dy
  real*8 time
  integer :: n, nprocs, i,j


  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,nprocs,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  
  if (nprocs .ne. ng0) then
    print*,'number of processors are not equal to ng0'
    call mpi_finalize(ierr)
    stop
  endif
  
  do n=1,nb0
    ie(n)=min(1+n0*n,j1)
  end do
  write(*,*) nb0, ie
  write(*,251) ie
  251  format('ie'/(20i4))
  if (ie(nb0).le.ie(nb1).or.j1.ge.ie(nb1)+9) then
    write(*,252) ie(nb1),j1
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
  
  
  time=mpi_wtime()
  
  do j=1,j2
    do i=1,i2
      at(i,j)=1/dy/dy 
      ab(i,j)=1/dy/dy
      ar(i,j)=1/dx/dx
      al(i,j)=1/dx/dx
      ac(i,j)=-2*(1/dx/dx+1/dy/dy) 
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
  
  
  
  print*,'EVP has finished successfully.'
  call mpi_finalize(ierr)

end PROGRAM MAIN

subroutine pre(ax,ay,bb,cx,cy,rinv,rinv1,h,ie)
  use mpi  
  use domain
  use mpi_grid
  implicit none

  !INPUT VARIABLES
  real*8, dimension(i2,j2),intent(in)      :: ax, ay, bb, cx
  real*8, dimension(i2,0:j2),intent(inout) :: cy
  real*8, dimension(nb0),intent(in)        :: ie

  !OUTPUT VARIABLES
  real*8, dimension(ibir, i0, nbg0),intent(inout) :: rinv
  real*8, dimension(ibir, i0, nbg1),intent(inout) :: rinv1
  real*8, dimension(i0,j0),intent(inout) :: h

  !LOCAL VARIABLES
  integer :: i,j,k,jl,jh,jg,jhm,ii,ig,jhp,it,jt,jss,jff,info
  integer :: n,nsd,ng,nb,nbg,idxy, nsend,nx,ny,ngt,ngxt,nt
  real*8,dimension(i2,j2)   :: axt,ayt,bbt,cxt,cyt
  real*8,dimension(i2t,i0t) :: rtemp
  real*8,dimension(i2t,i2t) :: rtp
  real*8,dimension(i0t*i0t) :: rp
  real*8,dimension(i2)      :: ct
  real*8,dimension(i2t)     :: ipvt
  real*8                    :: ctemp

  do ng = 1, ngy0 !loop 100
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
    if (mylat == ng-2) then
      do i =1, i2
        ct(i) = cy(i,j2)
      end do
    end if
    
    nsd=i2*j2
    call mpi_bcast(axt,nsd,mpi_real,ng-1,mpi_comm_lat,ierr)
    call mpi_bcast(ayt,nsd,mpi_real,ng-1,mpi_comm_lat,ierr)
    call mpi_bcast(bbt,nsd,mpi_real,ng-1,mpi_comm_lat,ierr)
    call mpi_bcast(cxt,nsd,mpi_real,ng-1,mpi_comm_lat,ierr)
    call mpi_bcast(cyt,nsd,mpi_real,ng-1,mpi_comm_lat,ierr)
    if (ng /= 1)  &
      call mpi_bcast(ct,i2,mpi_real,ng-2,mpi_comm_lat,ierr)
    if(mylat == ng-1) then
      do i = 1, i2
        cy(i,0) = ct(i)
      end do
    end if

    do nb = 1, nb0
      nbg = (ng-1)*nb0+nb
      write(*,*) nb, nb, nb0
      write(*,111) nbg,nbg0
      111  format ('processing block #',i3,' out of ',i3,' total evp solver')
      jh=ie(nb)
      jhp=jh+1
      jhm=jh-2
      jg=jl+1

      do ii =1, ibir
        ig=mylat*ibir+ii
        idxy=(ig-1)/i2
      end do
      do j=jl,jhp
        do i=1,i0
          h(i,j)=0.d0
        end do
      end do

      if (idxy .eq. mylon) then
        if (mod(ig,i2) .eq. 0) then
          h(i1,jg)=1.d0
        else
          h(mod(ig,i2)+1,jg)=1.d0
        endif
      endif

      write(*,*)  "before sendrecv" 

      call mpi_sendrecv(h(2,jg),1,mpi_real8,mpi_w,1,h(i0,jg),1, &
        mpi_real8,mpi_e,1,mpi_comm_world,istat,ierr)
      call mpi_sendrecv(h(i1,jg),1,mpi_real8,mpi_e,1,h(1,jg),1, &
        mpi_real8,mpi_w,1,mpi_comm_world,istat,ierr)
      write(*,*)  "after sendrecv" 
      if (nbg .ne. 1) then
        if (idxy .eq. mylon) then
          if (nb .ne. 1) then
            if (mod(ig,i2) .eq. 0) then
              ctemp=cyt(i2,jg-2)
            else
              ctemp=cyt(mod(ig,i2),jg-2)
            endif
          else
            if (mod(ig,i2) .eq. 0) then
              ctemp=ct(i2)
            else
              ctemp=ct(mod(ig,i2))
            endif
          endif
        endif
        write(*,*)  "before bcast" 

        call mpi_bcast(ctemp,1,mpi_real,idxy,mpi_comm_lon,ierr)
        do n=1,i0
          h(n,jl)=rinv1(ii,n,nbg-1)*ctemp
        end do
      endif
      do j=jl,jhm
        do  i=1,i2
          h(i+1,j+2)=-(axt(i,j)*h(i,j+1)+ayt(i,j)*h(i+1,j)+ &
            bbt(i,j)*h(i+1,j+1)+cxt(i,j)*h(i+2,j+1))/cyt(i,j)
        end do 
        call mpi_sendrecv(h(2,j+2),1,mpi_real8,mpi_w,1,h(i0,j+2),1, &
          mpi_real8,mpi_e,1,mpi_comm_world,istat,ierr)
        call mpi_sendrecv(h(i1,j+2),1,mpi_real8,mpi_e,1,h(1,j+2),1, &
          mpi_real8,mpi_w,1,mpi_comm_world,istat,ierr)
      end do

      j=jh-1
      do i=1,i2
        rinv(ii,i+1,nbg)=axt(i,j)*h(i,j+1)+ayt(i,j)*h(i+1,j)+bbt(i,j)* &
          h(i+1,j+1)+cxt(i,j)*h(i+2,j+1)
      end do

      if (nbg.ne.nbg0) then
        j=ie(nb)
        do n=1,i0
          rinv(ii,n,nbg0)=h(n,j)
        end do
        if (idxy .eq. mylon) then
          if (mod(ig,i2) .eq. 0) then
            h(i1,jg)=0.d0
          else
            h(mod(ig,i2)+1,jg)=0.d0
          endif
        endif
      endif

      nsend=i2*ibir

      write(*,*)  "before gather" ,myid

      call mpi_gather(rinv(1,2,nbg),nsend,mpi_real8,rp(1),nsend, &
        mpi_real8,0,mpi_comm_world,ierr)

      write(*,*)  "after gather", myid

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
      print*,'info =',info
      if (info .ne. 0) then
        print*,'preporcessor failed because rinv is singular at nbg =',nbg
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





