C ----------------------------------------------------------------------
      PROGRAM PREP
C ----------------------------------------------------------------------
      INCLUDE './resglo.h'
      INCLUDE 'mpif.h'
      PARAMETER(I01=I0+1,I0T1=I0T+1
     1,J01=J0+1,J0T1=J0T+1
     1,NBG0=NB0*NGY0,NBG1=NBG0-1,IBIR=I2T/NGY0)
      REAL*8 RINV,RINV1,H,IE
      REAL AL,AB,AC,AR,AT,DX,DY,FW,BW
      REAL*8 BY,AY,BX,AX
      DIMENSION RINV(IBIR,I0,NBG0),RINV1(IBIR,I0,NBG1),H(I0,J0),IE(NB0)
      DIMENSION AL(I2,J2),AB(I2,J2),AC(I2,J2),AR(I2,J2),AT(I2,J2)
      CHARACTER GRD*7
      COMMON/MPI/MPI_COMM_2D,MPI_COMM_LON,MPI_COMM_LAT,NPROC,MYID,MYLON,
     1 MYLAT,IERR,MPI_N,MPI_E,MPI_S,MPI_W,ISTAT(MPI_STATUS_SIZE),GRD
  



      CHARACTER*50 fname

      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(mpi_comm_world,nprocs,ierr)
      call MPI_COMM_RANK(mpi_comm_world,myid,ierr)

      IF (NPROCS .ne. NG0) THEN
      PRINT*,'NUMBER of PROCESSORs are not equal to NG0'
      CALL MPI_FINALIZE(IERR)
      stop
      ENDIF

      DO 250 N=1,NB0
 250  IE(N)=MIN(1+N0*N,J1)
      WRITE(*,*) NB0, IE
      WRITE(*,251) IE
 251  FORMAT('IE'/(20I4))
      IF (IE(NB0).GT.IE(NB1).AND.J1.LT.IE(NB1)+9) GO TO 260
      WRITE(*,252) IE(NB1),J1
 252  FORMAT('IE(NB1),J1=',2I5,' not properly defined. '
     1 'Fix above DO 250 loop and restart.')

      CALL MPI_GRID_GEN
 260  IE(NB0)=J1

      BY=1
      BX=1
      AX=0
      AY=0      

      dy=(by-ay)/dble(J0T)
      dx=(bx-ax)/dble(I0T)


      time=MPI_WTIME()

      DO 100 J=1,J2
         DO 100 I=1,I2
      AT(I,J)=1/DY/DY 
      AB(I,J)=1/DY/DY
      AR(I,J)=1/DX/DX
      AL(I,J)=1/DX/DX
 100  AC(I,J)=-2*(1/DX/DX+1/DY/DY) 
C   ---------------------
C          AT
C          |
C   AL --- AC ---  AR
C          |
C          AB
C   ---------------------


C      if (.false.) then

CTS      write(fname,'(A7,I2.2)')'EDATA_1',MYID+1
C      print*,fname
CTS      OPEN(100,file=fname,form='unformatted')
CTS      READ(100) AL,AR,AB,AT,AC


      CALL PRE(AL,AB,AC,AR,AT,RINV,RINV1,H,IE)


C Determined the step of the marching size


      time=MPI_WTIME()-time
      print*,MYID,time

      write(fname,'(A3,I2.2,A4)')'Inf',MYID,'.out'
      OPEN(99,file='./EVP'//GRD,form='unformatted')
      WRITE(99) AL,AR,AB,AT,AC,RINV,RINV1,IE

      write(fname,'(A5,I2.2,A4)')'./ptime',MYID,'.out'
      OPEN(1534,file=fname)
      WRITE(1534,*),time


C      endif

      print*,'EVP has finished successfully.'
      call mpi_finalize(ierr)
      END

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


      SUBROUTINE PRE(AX,AY,BB,CX,CY,RINV,RINV1,H,IE)
C ----------------------------------------------------------------------
      INCLUDE './resglo.h'
      INCLUDE 'mpif.h'
      PARAMETER(I01=I0+1,I0T1=I0T+1
     1,J01=J0+1,J0T1=J0T+1
     1,NBG0=NB0*NGY0,NBG1=NBG0-1,IBIR=I2T/NGY0)
      CHARACTER GRD*7
      REAL*8 RINV,RINV1,H,RUN,RTP,RTEMP,RP
      DIMENSION AX(I2,J2),AY(I2,J2),BB(I2,J2),CX(I2,J2),CY(I2,0:J2),
     1 RINV(IBIR,I0,NBG0),RINV1(IBIR,I0,NBG1),H(I0,J0),IE(NB0),
     2 RP(I0T*I0T),RTEMP(I2T,I0T),RTP(I2T,I2T),IPVT(I2T)
      DIMENSION AXT(I2,J2),AYT(I2,J2),BBT(I2,J2),CXT(I2,J2),CYT(I2,J2),
     1 CT(I2)

      COMMON/MPI/MPI_COMM_2D,MPI_COMM_LON,MPI_COMM_LAT,NPROC,MYID,MYLON,
     1 MYLAT,IERR,MPI_N,MPI_E,MPI_S,MPI_W,ISTAT(MPI_STATUS_SIZE),GRD

      DO 100 NG=1,NGY0
      JL=1
      IF (MYLAT .EQ. NG-1) THEN
      DO 99 J=1,J2
      DO 99 I=1,I2
      AXT(I,J)=AX(I,J)
      AYT(I,J)=AY(I,J)
      BBT(I,J)=BB(I,J)
      CXT(I,J)=CX(I,J)
 99   CYT(I,J)=CY(I,J)
      ENDIF
      IF (MYLAT .EQ. NG-2) THEN
      DO 98 I=1,I2
 98   CT(I)=CY(I,J2)
      ENDIF
      NSD=I2*J2
      CALL MPI_BCAST(AXT,NSD,MPI_REAL,NG-1,MPI_COMM_LAT,IERR)
      CALL MPI_BCAST(AYT,NSD,MPI_REAL,NG-1,MPI_COMM_LAT,IERR)
      CALL MPI_BCAST(BBT,NSD,MPI_REAL,NG-1,MPI_COMM_LAT,IERR)
      CALL MPI_BCAST(CXT,NSD,MPI_REAL,NG-1,MPI_COMM_LAT,IERR)
      CALL MPI_BCAST(CYT,NSD,MPI_REAL,NG-1,MPI_COMM_LAT,IERR)
      IF (NG .NE. 1)
     1 CALL MPI_BCAST(CT,I2,MPI_REAL,NG-2,MPI_COMM_LAT,IERR)

      IF (MYLAT .EQ. NG-1) THEN
      DO 97 I=1,I2
 97   CY(I,0)=CT(I)
      ENDIF
      DO 100 NB=1,NB0
      NBG=(NG-1)*NB0+NB
      WRITE(*,*) NB, NG, NB0
      WRITE(*,111) NBG,NBG0
 111  FORMAT ('Processing block #',I3,' out of ',I3,' total evp solver')
      JH=IE(NB)
      JHP=JH+1
      JHM=JH-2
      JG=JL+1
      DO 250 II=1,IBIR
      IG=MYLAT*IBIR+II
      IDXY=(IG-1)/I2
      DO 210 J=JL,JHP
      DO 210 I=1,I0
 210  H(I,J)=0.d0
      IF (IDXY .EQ. MYLON) THEN
      IF (MOD(IG,I2) .EQ. 0) THEN
      H(I1,JG)=1.d0
      ELSE
      H(MOD(IG,I2)+1,JG)=1.d0
      ENDIF
      ENDIF

      write(*,*)  "before sendrecv" 

      CALL MPI_SENDRECV(H(2,JG),1,MPI_REAL8,MPI_W,1,H(I0,JG),1,
     1                  MPI_REAL8,MPI_E,1,MPI_COMM_WORLD,ISTAT,IERR)
      CALL MPI_SENDRECV(H(I1,JG),1,MPI_REAL8,MPI_E,1,H(1,JG),1,
     1                  MPI_REAL8,MPI_W,1,MPI_COMM_WORLD,ISTAT,IERR)
      write(*,*)  "after sendrecv" 

      IF (NBG .EQ. 1) GO TO 220
      IF (IDXY .EQ. MYLON) THEN
      IF (NB .NE. 1) THEN
       IF (MOD(IG,I2) .EQ. 0) THEN
       CTEMP=CYT(I2,JG-2)
       ELSE
       CTEMP=CYT(MOD(IG,I2),JG-2)
       ENDIF
      ELSE
       IF (MOD(IG,I2) .EQ. 0) THEN
       CTEMP=CT(I2)
       ELSE
       CTEMP=CT(MOD(IG,I2))
       ENDIF
      ENDIF
      ENDIF

        
      write(*,*)  "before bcast" 
      CALL MPI_BCAST(CTEMP,1,MPI_REAL,IDXY,MPI_COMM_LON,IERR)
      DO 218 N=1,I0
 218  H(N,JL)=RINV1(II,N,NBG-1)*CTEMP
 220  CONTINUE
      DO 226 J=JL,JHM
      DO 225 I=1,I2
 225  H(I+1,J+2)=-(AXT(I,J)*H(I,J+1)+AYT(I,J)*H(I+1,J)+
     1 BBT(I,J)*H(I+1,J+1)+CXT(I,J)*H(I+2,J+1))/CYT(I,J)
      CALL MPI_SENDRECV(H(2,J+2),1,MPI_REAL8,MPI_W,1,H(I0,J+2),1,
     1                  MPI_REAL8,MPI_E,1,MPI_COMM_WORLD,ISTAT,IERR)
 226  CALL MPI_SENDRECV(H(I1,J+2),1,MPI_REAL8,MPI_E,1,H(1,J+2),1,
     1                  MPI_REAL8,MPI_W,1,MPI_COMM_WORLD,ISTAT,IERR)
      J=JH-1
      DO 230 I=1,I2
 230  RINV(II,I+1,NBG)=AXT(I,J)*H(I,J+1)+AYT(I,J)*H(I+1,J)+BBT(I,J)*
     1 H(I+1,J+1)+CXT(I,J)*H(I+2,J+1)

      IF (NBG.EQ.NBG0) GO TO 250
      J=IE(NB)
      DO 240 N=1,I0
 240  RINV(II,N,NBG0)=H(N,J)
      IF (IDXY .EQ. MYLON) THEN
      IF (MOD(IG,I2) .EQ. 0) THEN
      H(I1,JG)=0.d0
      ELSE
      H(MOD(IG,I2)+1,JG)=0.d0
      ENDIF
      ENDIF
 250  CONTINUE

      NSEND=I2*IBIR

      write(*,*)  "before gather" ,MYID

      CALL MPI_GATHER(RINV(1,2,NBG),NSEND,MPI_REAL8,RP(1),NSEND,
     1                MPI_REAL8,0,MPI_COMM_WORLD,IERR)

      write(*,*)  "after gather", MYID

      IF (MYID .EQ. 0) then
      DO 497 NX=1,NGX0
      NGXT=(NX-1)*NGY0
      DO 497 NY=1,NGY0
      NGT=((NY-1)+NGXT)*IBIR*I2
      DO 497 J=1,I2
      NT=NGT+(J-1)*IBIR
      JT=J+(NX-1)*I2
      DO 497 I=1,IBIR
      IT=I+(NY-1)*IBIR
 497  RTEMP(IT,JT)=RP(NT+I)

      DO 498 I=1,I2T
      DO 498 J=1,I2T
 498  RTP(J,I)=0.d0
      DO 499 I=1,I2T
 499  RTP(I,I)=1.d0
      call dgesv(I2T,I2T,RTEMP,I2T,IPVT,RTP,I2T,INFO)
      do 500 I=1,I2T
      do 500 J=1,I2T
 500  RTEMP(J,I)=RTP(I,J)
      ENDIF

      CALL MPI_BCAST(INFO,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      PRINT*,'INFO =',INFO
      IF (INFO .NE. 0) THEN
      PRINT*,'Preporcessor failed because RINV is singular at NBG =',NBG
      CALL MPI_FINALIZE(IERR)
      STOP3
      ENDIF
      NSEND=I2T*I2T
      CALL MPI_BCAST(RTP,NSEND,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      IT=MYLON*I2
      JT=MYLAT*IBIR

      JSS=1
      JFF=I0
      IF (MYLON .EQ.0) THEN
      DO 501 I=1,IBIR
 501  RINV(I,1,NBG)=RTP(JT+I,I2T)
      JSS=2
      ENDIF
      IF (MYLON .EQ. NGX1) THEN
      DO 502 I=1,IBIR
 502  RINV(I,I0,NBG)=RTP(JT+I,1)
      JFF=I1
      ENDIF
      DO 503 J=JSS,JFF
      DO 503 I=1,IBIR
 503  RINV(I,J,NBG)=RTP(JT+I,IT+J-1)

      IF (NBG.EQ.NBG0) RETURN
      NSEND=I0*IBIR
      CALL MPI_ALLGATHER(RINV(1,1,NBG0),NSEND,MPI_REAL8,RP,NSEND,
     1                   MPI_REAL8,MPI_COMM_LAT,IERR)
      DO 504 NY=1,NGY0
      NGT=(NY-1)*IBIR*I0
      DO 504 J=1,I0
      NT=NGT+(J-1)*IBIR
      DO 504 I=1,IBIR
      IT=I+(NY-1)*IBIR
 504  RTEMP(IT,J)=RP(NT+I)

      IT=IBIR*MYLAT
      DO 300 I=1,IBIR
      DO 300 J=1,I0
      RINV1(I,J,NBG)=0.d0
      DO 300 K=1,I2T
 300  RINV1(I,J,NBG)=RINV1(I,J,NBG)-RTP(IT+I,K)*RTEMP(K,J)
      JL=JH

 100  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      END

      SUBROUTINE MPI_GRID_GEN
      INCLUDE 'mpif.h'

      CHARACTER GRD*7
      COMMON/MPI/MPI_COMM_2D,MPI_COMM_LON,MPI_COMM_LAT,NPROC,MYID,MYLON,
     1 MYLAT,IERR,MPI_N,MPI_E,MPI_S,MPI_W,ISTAT(MPI_STATUS_SIZE),GRD

      DIMENSION NDIM(2),MYCRD(2)
      LOGICAL PERI(2),RDIM1(2),RDIM2(2)
      NDIM(1)=NGX0
      NDIM(2)=NGY0
      PERI(1)=.TRUE.
      PERI(2)=.FALSE.
      RDIM1(1)=.TRUE.
      RDIM1(2)=.FALSE.
      RDIM2(1)=.FALSE.
      RDIM2(2)=.TRUE.
      CALL MPI_CART_CREATE(MPI_COMM_WORLD,2,NDIM,PERI,
     1     .TRUE.,MPI_COMM_2D,IERR)
      CALL MPI_CART_SHIFT(MPI_COMM_2D,0,1,MPI_W,MPI_E,IERR)
      CALL MPI_CART_SHIFT(MPI_COMM_2D,1,1,MPI_S,MPI_N,IERR)
      CALL MPI_CART_SUB(MPI_COMM_2D,RDIM1,MPI_COMM_LON,IERR)
      CALL MPI_CART_SUB(MPI_COMM_2D,RDIM2,MPI_COMM_LAT,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_LON,MYLON,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_LAT,MYLAT,IERR)

C CARTISEAN REMARK
C  -------------------------------
C |   2   |   5   |   8   |   11  |
C | (0,2) | (1,2) | (2,2) | (3,2) |
C |-------|-------|-------|-------|
C |   1   |   4   |   7   |   10  |
C | (0,1) | (1,1) | (2,1) | (3,1) |
C |-------|-------|-------|-------|
C |   0   |   3   |   6   |   9   |
C | (0,0) | (1,0) | (2,0) | (3,0) |
C  -------------------------------
      write(GRD,'(I3.3,A1,I3.3)')MYLON,'_',MYLAT

      END
C ----------------------------------------------------------------------
