module MPI_GRID
  use mpi
  implicit none

  character, public :: grd*7
  integer, public   :: mpi_comm_2d,mpi_comm_lon,mpi_comm_lat,nproc,myid, & 
              mylon, mylat,ierr,mpi_n,mpi_e,mpi_s,mpi_w, istat(MPI_STATUS_SIZE)

  public :: MPI_GRID_GEN

contains
  subroutine MPI_GRID_GEN

    include './resglo.h90'
    integer, dimension(2) ::  ndim,mycrd
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
 end subroutine MPI_GRID_GEN
end module MPI_GRID



!================================================!
!               PREP for rep                     !
!================================================!
PROGRAM MAIN
  use MPI_GRID
  use mpi

  implicit none
  include './resglo.h90'
  
  
  integer,parameter :: i01=i0+1,i0t1=i0t+1                        &
                  ,j01=j0+1,j0t1=j0t+1                            &
                  ,nbg0=nb0*ngy0,nbg1=nbg0-1,ibir=i2t/ngy0
  real*8            :: rinv(ibir,i0,nbg0),rinv1(ibir,i0,nbg1),h(i0,j0),ie(nb0)
  real*8, dimension(i2,j2) :: al,ab,ac,ar,at
  integer :: n, nprocs
   
  
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(mpi_comm_world,nprocs,ierr)
  call MPI_COMM_RANK(mpi_comm_world,myid,ierr)

  IF (NPROCS .ne. NG0) THEN
        PRINT*,'NUMBER of PROCESSORs are not equal to NG0'
        CALL MPI_FINALIZE(IERR)
        stop
  ENDIF

  DO N=1,NB0
    IE(N)=MIN(1+N0*N,J1)
  end do
  WRITE(*,*) NB0, IE
  WRITE(*,251) IE
  251  FORMAT('IE'/(20I4))
  IF (IE(NB0).le.IE(NB1).or.J1.ge.IE(NB1)+9) then
  WRITE(*,252) IE(NB1),J1
  252  FORMAT('IE(NB1),J1=',2I5,' not properly defined. ' &
           ,'Fix above DO 250 loop and restart.') 
  endif
  CALL MPI_GRID_GEN


end PROGRAM MAIN

