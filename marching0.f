      implicit DOUBLE PRECISION (a-h,o-z)
      include 'mpif.h'
      PARAMETER(NGRID=4,I0=8,J0=10,NB0=2,NBG0=NB0*NGRID)
      PARAMETER(IJ0=I0+J0,IJ1=IJ0-1,I1=I0-1,I2=I0-2,I3=I0-3,J1=J0-1,
     1 J2=J0-2,J3=J0-3,NB1=NB0-1)
      REAL*8 RINV,RINV1,DUM0,DUM1,DUM2,X,H
      REAL*8 CCOL,BINV,CINV,ETMP,EBUF
      REAL*8 AL,AB,AC,AR,AT,DX,DY,F
      REAL*8 BY,AY,BX,AX
      REAL*8 R,RB,RBUF,XSB,XRB,XSB1,XRB1
      INTEGER ID,IP,IM,JS,JF,IDM,NNY
      INTEGER   LEFT,RIGHT,MYID,IERR,NSEND,ISTAT,IREQ
      INTEGER   INFO,IPVT
      CHARACTER*50 fname
      DIMENSION IE(NB0),ID(NB0),IP(NB0),IM(NB0),JS(NB0),JF(NB0),
     1          IDM(NB0,2),ISTAT(MPI_STATUS_SIZE),IPVT(NBG0*J2)
      DIMENSION CCOL(NBG0*J2,NB0*J2),
     1          BINV(NBG0*J2,NBG0*J2),CINV(NBG0*J2,NBG0*J2),
     2          ETMP(J2,J2,2,NB0,2),EBUF(J2,J2,4)
      DIMENSION ALSBUF(J2),ALRBUF(J2)
      DIMENSION AL(I2,J2),AB(I2,J2),AC(I2,J2),AR(I2,J2),AT(I2,J2)
      DIMENSION F(I2,J2),X(I0,J0)
      DIMENSION XSB1(J2),XRB1(J2),R(NB0*J2),RB(NBG0*J2),RBUF(J2,2)


      CALL MPI_INIT(ierr)
      call MPI_COMM_SIZE(mpi_comm_world,nprocs,ierr)
      call MPI_COMM_RANK(mpi_comm_world,myid,ierr)


      write(fname,'(A3,I2.2,A4)')'Inf',MYID,'.out'
      OPEN(99,file=fname,form='unformatted')
      READ(99) al,ab,ac,ar,at,ccol,id,ip,im,js,jf

C      if (.false.) then
      DO 100 J=2,J1
          DO 100 I=2,I1
      X(I,J)=0.d0
 100  F(I-1,J-1)=1 !!! out of range for F 

      DO 101 J=1,J0
      X(I0,J)=0.d0
 101  X(1,J)=0.d0

      DO 103 I=1,I0
      X(I,1)=0.d0
 103  X(I,J0)=0.d0

      LEFT=MYID-1
      RIGHT=MYID+1
      IF (MYID .EQ. 0) LEFT=NGRID-1
      IF (MYID .EQ. NGRID-1) RIGHT=0
C      if (.false.) then

      DO 150 NB=1,NB0
      DO 105 I=JS(NB),JF(NB)-ID(NB),ID(NB)
      DO 105 J=1,J2
      FW=IP(NB)*AR(I-1,J)+IM(NB)*AL(I-1,J)
      BW=IP(NB)*AL(I-1,J)+IM(NB)*AR(I-1,J)
 105  X(I+ID(NB),J+1)=(F(I-1,J)-AC(I-1,J)*X(I,J+1)-AT(I-1,J)*X(I,J+2)
     1 -AB(I-1,J)*X(I,J)-BW*X(I-ID(NB),J+1))/FW
      JL=(NB+ID(NB)-1)*J2
      DO 107 J=1,J2
      FW=IP(NB)*AR(JF(NB)-1,J)+IM(NB)*AL(JF(NB)-1,J)
      BW=IP(NB)*AL(JF(NB)-1,J)+IM(NB)*AR(JF(NB)-1,J)
 107  R(J+JL)=(F(JF(NB)-1,J)-AC(JF(NB)-1,J)*X(JF(NB),J+1)
     1 -AT(JF(NB)-1,J)*X(JF(NB),J+2)-AB(JF(NB)-1,J)*X(JF(NB),J)
     2 -BW*X(JF(NB)-ID(NB),J+1))/FW*(dble(ID(NB)))
 150  CONTINUE

      DO 160 NB=1,NB0
      JL=(NB-1)*J2
      DO 160 J=1,J2
 160  R(J+JL)=dble(ID(NB))*X(JF(NB),J+1)+R(J+JL)


      NSEND=J2*NB0
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_allgather(R(1),NSEND,mpi_double_precision,
     &     Rb(1),NSEND,mpi_double_precision,mpi_comm_world,ierr)
      call mpi_barrier(mpi_comm_world,ierr)

      do 180 J=1,NB0*J2
      R(J)=0.d0
      do 180 I=1,NBG0*J2
 180  R(J)=Rb(I)*CCOL(I,J)+R(J)

      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_allgather(R(1),NSEND,mpi_double_precision,
     &     Rb(1),NSEND,mpi_double_precision,mpi_comm_world,ierr)

      DO 200 NB=1,NB0
      DO 200 J=1,J2
 200  X(JS(NB),J+1)=X(JS(NB),J+1)+Rb((MYID*NB0+NB-1)*J2+J)

C For DBC

      IF (MYID .gt. 0) then
      DO 210 J=1,J2
 210  X(JS(1)-1,J+1)=X(JS(1)-1,J+1)+Rb(((LEFT+1)*NB0-1)*J2+J)
      ENDIF      
      IF (MYID .lt. NGRID-1) then
      DO 220 J=1,J2
 220  X(JS(NB0)+1,J+1)=X(JS(NB0)+1,J+1)+Rb((RIGHT*NB0)*J2+J)
      ENDIF
CC
      DO 250 NB=1,NB0
      DO 205 I=JS(NB),JF(NB)-ID(NB),ID(NB)
      DO 205 J=1,J2
      FW=IP(NB)*AR(I-1,J)+IM(NB)*AL(I-1,J)
      BW=IP(NB)*AL(I-1,J)+IM(NB)*AR(I-1,J)
 205  X(I+ID(NB),J+1)=(F(I-1,J)-AC(I-1,J)*X(I,J+1)-AT(I-1,J)*X(I,J+2)
     1 -AB(I-1,J)*X(I,J)-BW*X(I-ID(NB),J+1))/FW
 250  CONTINUE

      write(fname,'(A3,I2.2,A4)')'ANS',MYID,'.out'
      OPEN(1233,file=fname)
      WRITE(1233,*),((X(I,J),I=2,I2),J=2,J1)
      call MPI_FINALIZE(ierr)

      END

