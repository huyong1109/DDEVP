      implicit DOUBLE PRECISION (a-h,o-z)
      INCLUDE 'resglo.h'
      include 'mpif.h'
      PARAMETER(I01=I0+1,J01=J0+1,JG0=J2*NGY0+2,JG1=JG0-1,JG2=JG0-2,
     1 IG0=I2*NGX0+2,IG1=IG0-1,IG2=IG0-2,
     1 NBG0=NB0*NGY0,NBG1=NBG0-1,IBIR=I2T/NGY0)
      REAL*8 RINV,RINV1,DUM0,DUM1,DUM2,DUMG,X
      REAL*8 AL,AR,AB,AT,AC,F
      COMMON/SEVP/RINV(IBIR,I0,NBG0),RINV1(IBIR,I0,NBG1),DUM0(I0,NB0),
     1 DUM1(IBIR),DUM2(I0),DUMG(I2T),F(I2,J2),X(I0,J0)
     + ,AL(I2,J2),AB(I2,J2),AC(I2,J2),AR(I2,J2),AT(I2,0:J2)
     + ,S(I2,J2),IE(NB0)
      CHARACTER GRD*7


      CALL MPI_INIT(ierr)
      call MPI_COMM_SIZE(mpi_comm_world,nprocs,ierr)
      call MPI_COMM_RANK(mpi_comm_world,myid,ierr)

! set initial and left-hand value
      DO 100 J=2,J1
          DO 100 I=2,I1
      X(I,J)=0.d0
 100  F(I,J)=1

      DO 101 J=1,J0
      X(I0,J)=0.d0
 101  X(1,J)=0.d0

      DO 103 I=1,I0
      X(I,1)=0.d0
 103  X(I,J0)=0.d0

      time=MPI_Wtime()
      OPEN(99,file='./EVP'//GRD,form='unformatted')
      IF(MYID.EQ.0) WRITE(*,*)'EVP solver'
      READ(99) AL,AR,AB,AT,AC,RINV,RINV1,IE
      CLOSE(99)

      CALL REP2(AL,AB,AC,AR,AT,RINV,RINV1,DUM0,DUM1,DUM2,DUMG,S,X,IE)

      time=MPI_Wtime()-time
      print*,myid,time/1000,time1/1000

      OPEN(99,file='./EVPANS'//GRD,form='unformatted')
      WRITE(99,*),((X(I,J),I=2,I2),J=2,J1)

C      write(fname,'(A4,I2.2,A4)')'TIME',MYID,'.out'
C      OPEN(1534,file=fname)
C      WRITE(1534,*),time



      call MPI_FINALIZE(ierr)

      END



C ----------------------------------------------------------------------
      SUBROUTINE REP2(AX,AY,BB,CX,CY,RINV,RINV1,DUM0,DUM1,DUM2,DUMG,
     1               F,X,IE)
C ----------------------------------------------------------------------
CIBM  REAL*8 RINV,RINV1,DUM0,DUM1,DUM2,X,H,ROW
      INCLUDE 'mpif.h'
      INCLUDE 'resglo.h'
      PARAMETER(NBG0=NB0*NGY0,NBG1=NBG0-1,IBIR=I2T/NGY0)
      REAL*8 X,RINV,RINV1,DUM1,DUM2,DUM0,DUMG,DD
      DIMENSION AX(I2,J2),AY(I2,J2),BB(I2,J2),CX(I2,J2),CY(I2,0:J2),
     1 F(I2,J2),RINV(IBIR,I0,NBG0),RINV1(IBIR,I0,NBG1),IE(NB0),
     2 DUM0(I0,NB0),DUM1(IBIR),DUM2(I0),X(I0,J0),DUMG(I2T)
      COMMON/MPI/M_CART,M_CLON,M_CLAT,MYID,MYLON,MYLAT,M_N,M_E,
     1 M_S,M_W,M_VLON,M_VLAT,M_V8LON,M_V8LAT,JS,JF,IERR,
     1 NPX(14),NPY(14),ISTAT(MPI_STATUS_SIZE)

      CHARACTER GRD*7
      COMMON/GRD/GRD


      DO 160 NG=0,NGY0-1
      JSP=1
      DO 150 NB=1,NB0
      NBG=NG*NB0+NB
      IF (MYLAT .EQ. NG) THEN
      JFP=IE(NB)-2
      DO 106 J=JSP,JFP
      DO 105 I=1,I2
 105  X(I+1,J+2)=(F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*
     1 X(I+1,J+1)-CX(I,J)*X(I+2,J+1))/CY(I,J)
      CALL MPI_SENDRECV(X(2,J+2),1,MPI_REAL8,M_W,1,X(I0,J+2),1,
     1                  MPI_REAL8,M_E,1,M_CART,ISTAT,IERR)
      CALL MPI_SENDRECV(X(I1,J+2),1,MPI_REAL8,M_E,1,X(1,J+2),1,
     1                  MPI_REAL8,M_W,1,M_CART,ISTAT,IERR)
 106  CONTINUE
      IF (NBG.EQ.NBG0) GO TO 160
      J=IE(NB)-1
      DO 115 I=1,I2
 115  DUM2(I)=F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*
     1 X(I+1,J+1)-CX(I,J)*X(I+2,J+1)-CY(I,J)*X(I+1,J+2)
      J=IE(NB)
      CALL MPI_ALLGATHER(DUM2,I2,MPI_REAL8,DUMG,I2,MPI_REAL8,
     1 M_CLON,IERR)
      DO 116 N=1,I0
      DUM2(N)=X(N,J)
 116  DUM0(N,NB)=X(N,J)
      DD=1.d0
      ELSE
      IF (NBG .EQ. NBG0) GO TO 160
      DD=0.d0
      ENDIF
      CALL MPI_SCATTER(DUMG,IBIR,MPI_REAL8,DUM1,IBIR,MPI_REAL8,NG,
     1                 M_CLAT,IERR)
      CALL DGEMV('T',IBIR,I0,-1.d0,RINV1(1,1,NBG),IBIR,DUM1,1,DD,DUM2,1)
      IF (NB .EQ. NB0) GO TO 150
      CALL MPI_REDUCE(DUM2,X(1,IE(NB)),I0,MPI_REAL8,MPI_SUM,NG,
     1     M_CLAT,IERR)
      JSP=IE(NB)

 150  CONTINUE
      CALL MPI_REDUCE(DUM2,X(1,1),I0,MPI_REAL8,MPI_SUM,NG+1,
     1     M_CLAT,IERR)
 160  CONTINUE

CW WENIEN
CW      OPEN(160,file='OPT/REP160'//GRD,form='unformatted')
CW      WRITE(160) X,F
CW      CLOSE(160)


      DO 300 NG=NGY0-1,0,-1
      DO 250 NBS=1,NB0
      NB=NB0-NBS+1
      NBG=NB0*NG+NB

      JSP=1
      IF (NB.NE.1) JSP=IE(NB-1)
      IF (MYLAT .EQ. NG) THEN
      IF (NBG.EQ.NBG0) GO TO 201
      J=IE(NB)
      DO 200 N=1,I0
 200  X(N,J)=DUM0(N,NB)
 201  J=IE(NB)-1
      DO 210 I=1,I2
 210  DUM2(I)=F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*
     1 X(I+1,J+1)-CX(I,J)*X(I+2,J+1)-CY(I,J)*X(I+1,J+2)
      CALL MPI_ALLGATHER(DUM2,I2,MPI_REAL8,DUMG,I2,MPI_REAL8,
     1 M_CLON,IERR)
      ENDIF

      CALL MPI_SCATTER(DUMG,IBIR,MPI_REAL8,DUM1,IBIR,MPI_REAL8,NG,
     1 M_CLAT,IERR)
      CALL DGEMV('T',IBIR,I0,1.d0,RINV(1,1,NBG),IBIR,DUM1,1,0.d0,DUM2,1)
      CALL MPI_REDUCE(DUM2,DUM0(1,NB),I0,MPI_REAL8,MPI_SUM,NG,
     1     M_CLAT,IERR)

      IF (MYLAT .EQ. NG) THEN
      DO 220 N=1,I0
 220  X(N,JSP+1)=X(N,JSP+1)+DUM0(N,NB)
      ENDIF
      IF (NBG.EQ.1) GO TO 300
      IF (MYLAT .EQ. NG) THEN
      DO 230 M=1,I2
 230  DUM2(M)=DUM0(M+1,NB)*CY(M,JSP-1)
      CALL MPI_ALLGATHER(DUM2,I2,MPI_REAL8,DUMG,I2,MPI_REAL8,
     1     M_CLON,IERR)
      ENDIF
      CALL MPI_SCATTER(DUMG,IBIR,MPI_REAL8,DUM1,IBIR,MPI_REAL8,NG,
     1 M_CLAT,IERR)
      CALL DGEMV('T',IBIR,I0,1.d0,RINV1(1,1,NBG-1),IBIR,DUM1,
     1 1,0.d0,DUM2,1)
      CALL MPI_REDUCE(DUM2,DUM0(1,NB),I0,MPI_REAL8,MPI_SUM,NG,
     1 M_CLAT,IERR)
      IF (MYLAT .EQ. NG) THEN
      DO 245 N=1,I0
 245  DUM0(N,NB)=X(N,JSP)+DUM0(N,NB)
      ENDIF
 250  CONTINUE
      IF (MYLAT .EQ. NG)
     1 CALL MPI_SEND(X(1,2),I0,MPI_REAL8,M_S,1,M_CART,IERR)
      IF (MYLAT .EQ. NG-1)
     1 CALL MPI_RECV(X(1,J0),I0,MPI_REAL8,M_N,1,M_CART,
     2 ISTAT,IERR)
 300  CONTINUE
      IF (MYLAT .EQ. 0) THEN
      DO 350 I=1,I0
 350  DUM0(I,1)=0.d0
      ENDIF

      JSP=1
      DO 400 NB=1,NB0
      JFP=IE(NB)-2
      IF (NB.EQ.NB0) JFP=IE(NB)-2
      DO 410 I=1,I0
 410  X(I,JSP)=DUM0(I,NB)
      DO 420 J=JSP,JFP
      DO 430 I=1,I2
 430  X(I+1,J+2)=(F(I,J)-AX(I,J)*X(I,J+1)-AY(I,J)*X(I+1,J)-BB(I,J)*
     1 X(I+1,J+1)-CX(I,J)*X(I+2,J+1))/CY(I,J)
      CALL MPI_SENDRECV(X(2,J+2),1,MPI_REAL8,M_W,1,X(I0,J+2),1,
     1                  MPI_REAL8,M_E,1,M_CART,ISTAT,IERR)
      CALL MPI_SENDRECV(X(I1,J+2),1,MPI_REAL8,M_E,1,X(1,J+2),1,
     1                  MPI_REAL8,M_W,1,M_CART,ISTAT,IERR)
 420  CONTINUE
 400  JSP=IE(NB)

CW WENIEN
CW      OPEN(160,file='OPT/REP'//GRD,form='unformatted')
CW      WRITE(160) X,F
CW      CLOSE(160)




      END

