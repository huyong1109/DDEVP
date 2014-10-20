C ----------------------------------------------------------------------
      PROGRAM PREP
C ----------------------------------------------------------------------
      INCLUDE 'mpif.h'
C NGRID : total proc; NB0: blocks in each proc; NBG0: total
C blocks
      PARAMETER(NGRID=4,I0=19,J0=66,NB0=4,NBG0=NB0*NGRID)
      PARAMETER(IJ0=I0+J0,IJ1=IJ0-1,I1=I0-1,I2=I0-2,I3=I0-3,J1=J0-1,
     1 J2=J0-2,J3=J0-3,NB1=NB0-1)
      REAL*8 RINV,RINV1,DUM0,DUM1,DUM2,X,H
      REAL*8 EHAT,CCOL,BINV,CINV,ETMP,EBUF
      REAL AL,AB,AC,AR,AT,DX,DY
      REAL*8 BY,AY,BX,AX
      INTEGER ID,IP,IM,JS,JF,IDM,NNY
      INTEGER   LEFT,RIGHT,MYID,IERR,NSEND,ISTAT,IREQ
      INTEGER   INFO,IPVT
      DIMENSION IE(NB0),ID(NB0),IP(NB0),IM(NB0),JS(NB0),JF(NB0),
     1          IDM(NB0,2),ISTAT(MPI_STATUS_SIZE),IPVT(NBG0*J2)
      DIMENSION EHAT(J2,J2,I0,2),CCOL(NBG0*J2,NB0*J2),
     1          BINV(NBG0*J2,NBG0*J2),CINV(NBG0*J2,NBG0*J2),
     2          ETMP(J2,J2,2,NB0,2),EBUF(J2,J2,4)
      DIMENSION ALSBUF(J2),ALRBUF(J2)
      DIMENSION AL(I2,J2),AB(I2,J2),AC(I2,J2),AR(I2,J2),AT(I2,J2)

      CHARACTER*50 fname

      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(mpi_comm_world,nprocs,ierr)
      call MPI_COMM_RANK(mpi_comm_world,myid,ierr)


      BY=1
      BX=1
      AX=0
      AY=0      

      dy=(by-ay)/dble(I0*NBG0)
      dx=(bx-ax)/dble(J0)
C matrix A*x = b coefficient
      DO 100 J=1,J2
         DO 100 I=1,I2
      AT(I,J)=1/DY              !Top 
      AB(I,J)=1/DY              !Bottom
      AR(I,J)=1/DX              !Right
      AL(I,J)=1/DX              !Left
 100  AC=-2*(1/DX+1/DY)         !Center



      DO 310 NB=1,NB0
      ID(NB)=(-1)**(NB+1)
      IP(NB)=mod(NB,2)
 310  IM(NB)=mod(NB+1,2)
      DO 313 NB=1,NB0
      IDM(NB,1)=(MYID*NB0+NB-IP(NB)-1)*J2
 313  IDM(NB,2)=(MYID*NB0+NB-IP(NB))*J2
      IF (MYID .EQ. 0) IDM(1,1)=(NBG0-1)*J2
      IF (MYID .EQ. NGRID-1) IDM(NB0,1)=0
      LEFT=MYID-1
      RIGHT=MYID+1
      IF (MYID .EQ. 0) LEFT=NGRID-1
      IF (MYID .EQ. NGRID-1) RIGHT=0


C Determined the step of the marching size

      NNY=I1/NB0
      if (NNY .LE. 2) then
      print*,'NBLK=',NB0,'NNX=',NNY
      print*,'NBLK is so small that EVP can''t works. PREP STOPS'
      CALL MPI_FINALIZE(IERR)
      elseif (NNY .GT. 10) then
      print*,'NBLK=',NB0,'NNX=',NNY
      print*,'NBLK is so large that EVP can''t works. PREP STOPS'
      CALL MPI_FINALIZE(IERR)
      endif 
      print*,'Marching steps =',NNY

C       Determine start line of each block. forward & backward
      DO 315 NB=1,NB0-1
 315  JS(NB)=NNY*(NB-IP(NB))+1+IP(NB)
      JS(NB0)=I2

      DO 316 NB=1,NB0
 316  JF(NB)=(JS(NB)+JS(NB+ID(NB))-ID(NB))/2
      DO 320 NB=1,NB0
      DO 320 J=1,J2
 320  EHAT(J,J,JS(NB),IP(NB)+1)=1.d0
      DO 325 J=1,J2
      EHAT(J,J,I1,2)=1.d0
 325  EHAT(J,J,1,1)=1.d0

      DO 350 NB=1,NB0
      DO 350 L=1,2
      DO 335 I=JS(NB),JF(NB)-ID(NB),ID(NB)
      DO 335 JJ=1,J2
      DO 330 J=2,J2-1
      FW=real(IP(NB))*AR(I-1,J)+real(IM(NB))*AL(I-1,J)
      BW=real(IP(NB))*AL(I-1,J)+real(IM(NB))*AR(I-1,J)
 330  EHAT(J,JJ,I+ID(NB),L)=(-AC(I-1,J)*EHAT(J,JJ,I,L)-AT(I-1,J)*EHAT(J+
     1 1,JJ,I,L)-AB(I-1,J)*EHAT(J-1,JJ,I,L)-BW*EHAT(J,JJ,I-ID(NB),L))/FW

      FW=real(IP(NB))*AR(I-1,1)+real(IM(NB))*AL(I-1,1);
      BW=real(IP(NB))*AL(I-1,1)+real(IM(NB))*AR(I-1,1);
      EHAT(1,JJ,I+ID(NB),L)=(-AC(I-1,1)*EHAT(1,JJ,I,L)-AT(I-1,1)*
     1 EHAT(2,JJ,I,L)-BW*EHAT(1,JJ,I-ID(NB),L))/FW;
      FW=real(IP(NB))*AR(I-1,J2)+real(IM(NB))*AL(I-1,J2);
      BW=real(IP(NB))*AL(I-1,J2)+real(IM(NB))*AR(I-1,J2);
 335  EHAT(J2,JJ,I+ID(NB),L)=(-AC(I-1,J2)*EHAT(J2,JJ,I,L)-AB(I-1,J2)*
     1 EHAT(J2-1,JJ,I,L)-BW*EHAT(J2,JJ,I-ID(NB),L))/FW;
 

      DO 345 JJ=1,J2
      DO 340 J=2,J2-1 
      ETMP(J,JJ,1,NB,L)=EHAT(J,JJ,JF(NB),L)
      FW=real(IP(NB))*AR(JF(NB)-1,J)+real(IM(NB))*AL(JF(NB)-1,J)
      BW=real(IP(NB))*AL(JF(NB)-1,J)+real(IM(NB))*AR(JF(NB)-1,J)
 340  ETMP(J,JJ,2,NB,L)=(-AC(JF(NB)-1,J)*EHAT(J,JJ,JF(NB),L)
     1 -AT(JF(NB)-1,J)*EHAT(J+1,JJ,JF(NB),L)-AB(JF(NB)-1,J)*
     2 EHAT(J-1,JJ,JF(NB),L)-BW*EHAT(J,JJ,JF(NB)-ID(NB),L))/FW
      ETMP(1,JJ,1,NB,L)=EHAT(1,JJ,JF(NB),L)
      FW=real(IP(NB))*AR(JF(NB)-1,1)+real(IM(NB))*AL(JF(NB)-1,1)
      BW=real(IP(NB))*AL(JF(NB)-1,1)+real(IM(NB))*AR(JF(NB)-1,1)
      ETMP(1,JJ,2,NB,L)=(-AC(JF(NB)-1,J)*EHAT(1,JJ,JF(NB),L)
     1 -AT(JF(NB)-1,J)*EHAT(2,JJ,JF(NB),L)
     2 -BW*EHAT(1,JJ,JF(NB)-ID(NB),L))/FW;
      ETMP(J2,JJ,1,NB,L)=EHAT(J2,JJ,JF(NB),L)
      FW=real(IP(NB))*AR(JF(NB)-1,J2)+real(IM(NB))*AL(JF(NB)-1,J2)
      BW=real(IP(NB))*AL(JF(NB)-1,J2)+real(IM(NB))*AR(JF(NB)-1,J2)
 345  ETMP(J2,JJ,2,NB,L)=(-AC(JF(NB)-1,J2)*EHAT(J2,JJ,JF(NB),L)
     1 -AB(JF(NB)-1,J2)*EHAT(J2-1,JJ,JF(NB),L)
     2 -BW*EHAT(J2,JJ,JF(NB)-ID(NB),L))/FW;  

 350  CONTINUE

      NSEND=J2*J2
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_sendrecv(ETMP(1,1,1,NB0,2),NSEND,MPI_REAL8,right,
     11,ebuf(1,1,1),NSEND,MPI_REAL8,left,1,MPI_COMM_WORLD,ISTAT,IERR)
      call mpi_sendrecv(ETMP(1,1,2,NB0,2),NSEND,mpi_real8,right,
     11,ebuf(1,1,2),NSEND,mpi_real8,left,1,MPI_COMM_WORLD,ISTAT,IERR)
      call mpi_sendrecv(ETMP(1,1,1,1,1),NSEND,mpi_real8,left,1,
     1EBUF(1,1,3),NSEND,mpi_real8,right,1,MPI_COMM_WORLD,ISTAT,IERR)
      call mpi_sendrecv(ETMP(1,1,2,1,1),NSEND,mpi_real8,left,1,
     1EBUF(1,1,4),NSEND,mpi_real8,right,1,MPI_COMM_WORLD,ISTAT,IERR)
      call mpi_barrier(mpi_comm_world,ierr)
      do 400 II1=1,J2
      do 400 II2=1,J2
      ETMP(II1,II2,1,1,1)=EBUF(II1,II2,1)
      ETMP(II1,II2,2,1,1)=EBUF(II1,II2,2)
      ETMP(II1,II2,1,NB0,2)=EBUF(II1,II2,3)
 400  ETMP(II1,II2,2,NB0,2)=EBUF(II1,II2,4)     

      L=(myid*NB0)*J2
      k=0
      call ARR(NBG0,NB0,J2,ETMP(1,1,1,1,2),Ccol,L,k,-1.d0)
      call ARR(NBG0,NB0,J2,ETMP(1,1,2,1,2),Ccol,L+J2,k,-1.d0)
      if (MYID .eq. 0) then
      L=(NBG0-2)*J2
C      call ARR(NBG0,NB0,J2,ETMP(1,1,2,1,1),Ccol,L,k,-1.d0)
C      call ARR(NBG0,NB0,J2,ETMP(1,1,1,1,1),Ccol,L+J2,k,-1.d0)
      else
      call ARR(NBG0,NB0,J2,ETMP(1,1,2,1,1),Ccol,L-2*J2,k,1.d0)
      call ARR(NBG0,NB0,J2,ETMP(1,1,1,1,1),Ccol,L-J2,k,1.d0)
      endif

      if (NB0 .gt. 2) then
      do NB=2,NB0-1
      L=((MYID*NB0)+(NB-1)+IM(NB))*J2
      k=(NB-1)*J2
      call ARR(NBG0,NB0,J2,ETMP(1,1,2,NB-IP(NB),IP(NB)+1),
     1 Ccol,L-2*J2,k,1.d0)
      call ARR(NBG0,NB0,J2,ETMP(1,1,1,NB-IP(NB),IP(NB)+1),
     1 Ccol,L-1*J2,k,1.d0)
      call ARR(NBG0,NB0,J2,ETMP(1,1,1,NB+IM(NB),IP(NB)+1),
     1 Ccol,L,k,-1.d0)
      call ARR(NBG0,NB0,J2,ETMP(1,1,2,NB+IM(NB),IP(NB)+1),
     1 Ccol,L+J2,k,-1.d0)
      enddo
      endif

      L=((myid+1)*NB0)*J2
      k=(NB0-1)*J2
      call ARR(NBG0,NB0,J2,ETMP(1,1,2,NB0,1),Ccol,L-2*J2,k,1.d0)
      call ARR(NBG0,NB0,J2,ETMP(1,1,1,NB0,1),Ccol,L-J2,k,1.d0)
      if (MYID .eq. 3) then
C      call ARR(NBG0,NB0,J2,ETMP(1,1,1,NB0,2),Ccol,0,k,1.d0)
C      call ARR(NBG0,NB0,J2,ETMP(1,1,2,NB0,2),Ccol,J2,k,1.d0)
      else
      call ARR(NBG0,NB0,J2,ETMP(1,1,1,NB0,2),Ccol,L,k,-1.d0)
      call ARR(NBG0,NB0,J2,ETMP(1,1,2,NB0,2),Ccol,L+J2,k,-1.d0)
      endif
      NSEND=NB0*NBG0*J2*J2
      call mpi_barrier(mpi_comm_world,ierr)
      print*,myid
      call mpi_gather(CCOL(1,1),NSEND,MPI_REAL8,CINV(1,1),NSEND,
     1                MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      if (myid .eq. 0) then
      do 499 I=1,NBG0*J2
 499  BINV(I,I)=1.d0
      call dgesv(NBG0*J2,NBG0*J2,CINV,NBG0*J2,IPVT,BINV,NBG0*J2,INFO)
      do 500 I=1,NBG0*J2
      do 500 J=1,NBG0*J2
 500  CINV(J,I)=BINV(I,J)
      endif

C      if (.false.) then

      call mpi_barrier(mpi_comm_world,ierr)
      call MPI_BCAST(INFO,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      call mpi_barrier(mpi_comm_world,ierr)
C      if (.false.) then
      if (INFO .ne. 0) then
      print*,'The influnce matrix is singular.'  
      endif
      NSEND=NBG0*NB0*J2*J2
      call mpi_scatter(CINV(1,1),NSEND,mpi_REAL8,Ccol(1,1),NSEND,
     1                 mpi_REAL8,0,mpi_comm_world,ierr)
      write(fname,'(A3,I2.2,A4)')'Inf',MYID,'.out'
      OPEN(99,file=fname,form='unformatted')
      WRITE(99) AL,AR,AB,AT,AC,Ccol,EHAT,IDM,IP,ID,IM,JS,JF
      print*,'EVP has finished successfully.'
      call mpi_finalize(ierr)
      END

C     submatrix to full matrix
      SUBROUTINE ARR(n1,n2,m,C,A,ix,iy,d)
      implicit double precision (a-h,o-z)
      integer n1,n2,m,ix,iy,id
      dimension C(m,m),A(m*n1,m*n2)
      integer i,j
      do i=1,m
         do j=1,m
         A(ix+i,iy+j)=d*C(i,j)
         enddo
      enddo
      end

