! *******************
! * elliptic solver *
! *******************
! ----------------------------------------------------------------------
subroutine rep(ax,ay,bb,cx,cy,rinv,rinv1,f,x,ie, i2, j2, nblk)

    !real*8 rinv,rinv1,dum0,dum1,dum2,x,h
    !dimension ax(i2,*),ay(i2,*),bb(i2,*),cx(i2,*),cy(i2,*), rinv(i2,i2,*),rinv1(i2,i2,*),h(i0,*),ie(*),dum0(i2,*),dum1(*), dum2(*),f(i2,*),x(i0,*)
    implicit none
    integer, intent(in) :: i2,j2, nblk
    integer, dimension(0:nblk), intent(in) :: ie
    real(8), dimension(i2, j2), intent(in) :: ax, ay, bb, cx, cy 
    real(8), dimension(i2, i2, nblk), intent(in) :: rinv, rinv1
    real(8), dimension(i2, j2), intent(in) :: f
    real(8), dimension(i2+2, j2+2), intent(inout) :: x

    ! local 
    integer :: i0, js, nb, jf, i, j, m, nbs, n
    real(8), dimension(i2+2, j2+2) :: h 
    real(8), dimension(i2, nblk) :: dum0
    real(8), dimension(i2) :: dum1, dum2
   
    i0 = i2+2
    do nb=1,nblk
      js=ie(nb-1)
      jf=ie(nb)-2
      do j=js,jf
        do i=1,i2
          x(i+1,j+2)=(f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* x(i+1,j+1)-cx(i,j)*x(i+2,j+1))/cy(i,j)
        enddo
      enddo
      if (nb.ne.nblk) then
        j=ie(nb)-1
        do i=1,i2
          dum1(i)=f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* x(i+1,j+1)-cx(i,j)*x(i+2,j+1)-cy(i,j)*x(i+1,j+2)
        enddo
        j=ie(nb)
        do n=1,i2
          dum2(n)=0.
          do m=1,i2
            dum2(n)=dum2(n)+dum1(m)*rinv1(m,n,nb)
            write(*,*) dum1(m), rinv1(m,n,nb) , m ,n ,nb
          enddo
          dum0(n,nb)=x(n+1,j)
          x(n+1,j)=x(n+1,j)-dum2(n)
        enddo
      endif
      !js=ie(nb)
    enddo

    do nbs=1,nblk
      nb=nblk-nbs+1
      js=ie(nb-1)
      jf=ie(nb)-2
      if (nb.ne.nblk) then
        j=ie(nb)
        do n=1,i2
          x(n+1,j)=dum0(n,nb)
        enddo
      endif
      n=ie(nb)
      do j=js,n
        do i=1,i0
          h(i,j)=0.
        enddo
      enddo
      j=ie(nb)-1
      do i=1,i2
        dum1(i)=f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* x(i+1,j+1)-cx(i,j)*x(i+2,j+1)-cy(i,j)*x(i+1,j+2)
      enddo
      do n=1,i2
        dum2(n)=0.
        do m=1,i2
          dum2(n)=dum2(n)+dum1(m)*rinv(m,n,nb)
        enddo
        h(n+1,js+1)=dum2(n)
        x(n+1,js+1)=x(n+1,js+1)+dum2(n)
      enddo
      if (nb.ne.1) then
        do m=1,i2
          dum1(m)=h(m+1,js+1)*cy(m,js-1)
        enddo
        j=ie(nb-1)
        do n=1,i2
          dum2(n)=0.
          do m=1,i2
            dum2(n)=dum2(n)+dum1(m)*rinv1(m,n,nb-1)
          enddo
          h(n+1,j)=dum2(n)
        enddo
      endif
      do j=js,jf
        do i=1,i2
          h(i+1,j+2)=(-ax(i,j)*h(i,j+1)-ay(i,j)*h(i+1,j)-bb(i,j)* h(i+1,j+1)-cx(i,j)*h(i+2,j+1))/cy(i,j)
          x(i+1,j+2)=x(i+1,j+2)+h(i+1,j+2)
        end do
      end do
    enddo




end subroutine rep
! *******************
! * elliptic solver *
! *	DDEVP	    *
! *******************
! ----------------------------------------------------------------------
subroutine ddevp(ax,ay,bb,cx,cy,rinv,rinv1,f,x,ie, i2, j2, nblk)

    !real*8 rinv,rinv1,dum0,dum1,dum2,x,h
    !dimension ax(i2,*),ay(i2,*),bb(i2,*),cx(i2,*),cy(i2,*), rinv(i2,i2,*),rinv1(i2,i2,*),h(i0,*),ie(*),dum0(i2,*),dum1(*), dum2(*),f(i2,*),x(i0,*)
    implicit none
    integer, intent(in) :: i2,j2, nblk
    integer, dimension(0:nblk), intent(in) :: ie
    real(8), dimension(i2, j2), intent(in) :: ax, ay, bb, cx, cy 
    real(8), dimension(i2, i2, nblk), intent(in) :: rinv, rinv1
    real(8), dimension(i2, j2), intent(in) :: f
    real(8), dimension(i2+2, j2+2), intent(inout) :: x

    ! local 
    integer :: i0, js, nb, jf, i, j, m, nbs, bkl,  n, it, maxit
    real(8), dimension(i2+2, j2+2) :: h
    real(8), dimension(i2, j2/nblk, nblk) :: axl, ayl, bbl, cxl, cyl
    real(8), dimension(i2+2, j2/nblk+2, nblk) :: hh, xt
    real(8), dimension(i2, j2/nblk) :: fl
    real(8), dimension(i2, nblk) :: dum0
    real(8), dimension(i2) :: dum1, dum2
    real(8), intent	    :: eps
    integer, dimension(0:1) :: iel
    character(LEN=80)	:: fmt
    real(8), dimension(nblk) :: ratio1, ratio2 ! efficient of left and right contribution

    ratio1(1) = 1.0D0 
    ratio2(nblk) = 1.0D0  
    do i = 1, nblk-1
	ratio1(i+1) = 0.5
	ratio1(i) = 0.5
    end do 


    maxit = 20
    i0 = i2+2
    bkl = j2/nblk
    iel(0) = 1
    iel(1) = bkl+1
    fmt = "(80f8.3)"


    !write(*,*) "iel", iel(0), iel(1)

    do it=1, maxit

	do nb=1,nblk
    	    js = ie(nb-1)
    	    jf = ie(nb)
	    write(*,*) "js - jf", js , jf
    	    axl(:,:,nb) = ax(:, js:jf-1)
    	    ayl(:,:,nb) = ay(:, js:jf-1)
    	    bbl(:,:,nb) = bb(:, js:jf-1)
    	    cxl(:,:,nb) = cx(:, js:jf-1)
    	    cyl(:,:,nb) = cy(:, js:jf-1)
	    xt(:,:,nb)  = x(:, js:jf+1)
    	end do 
	do nb=1,nblk

	    js = ie(nb-1)
	    jf = ie(nb)
	
	    fl(:,:) = f(:, js:jf-1)
	    
    	    !write(*,*) "coefficient axl "
    	    !write(*,fmt) axl
    	    !write(*,*) "coefficient ayl "
    	    !write(*,fmt) ayl
    	    !write(*,*) "coefficient bbl "
    	    !write(*,fmt) bbl
    	    !write(*,*) "coefficient cxl "
    	    !write(*,fmt) cxl
    	    !write(*,*) "coefficient cyl "
    	    !write(*,fmt) cyl
    	    !write(*,*) "rinvl"
    	    !write(*,fmt) rinv(:,:,nb), rinv1(:,:,nb)
    	    !write(*,*) "fffl"
    	    !write(*,fmt) fl
    	    !write(*,*) "i,j "
    	    !write(*,*) i2, bkl+2
    	    !write(*,*) "iel "
    	    !write(*,*) iel(:)

	    call rep(axl(:,:,nb),ayl(:,:,nb),bbl(:,:,nb),cxl(:,:,nb),cyl(:,:,nb),rinv(:,:,nb),rinv1(:,:,nb),fl,xt(:,:, nb),iel,i2, bkl, 1)
    	    !write(*,*) "xt "
    	    !write(*,fmt) xt
	    eps = 1.0D-5
	    CheckRep(axl(:,:,nb),ayl(:,:,nb),bbl(:,:,nb),cxl(:,:,nb),cyl(:,:,nb),fl,xt(:,:,nb),i2,j2,eps)
	end do 

	do nb=1,nblk
    	    js = ie(nb-1)
    	    jf = ie(nb)
	    write(*,*) 'bkl-1', jf-js, bkl-1
	    x(:, js+1:jf) = xt(:,2:bkl+1,nb)
    	    write(*,*) "iter =  ", it, " x, nb", nb
    	    write(*,'(10f8.3)') xt(:, :,nb)
	end do 

	!do j = 2, j2+1
    	!    x(1,j) = x(2,j)
    	!    x(i2+2,j) = x(i2+1,j)
    	!end do 
    	    write(*,*) "iter =  ", it, " x"
    	    write(*,'(10f8.3)') x(:, :)
    end do

    !do it=1, maxit
    !    do nb=1,nblk
    !	    js=ie(nb-1)
    !	    jf=ie(nb)-2
    !	    do j=js,jf
    !		do i=1,i2
    !		    x(i+1,j+2)=(f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* x(i+1,j+1)-cx(i,j)*x(i+2,j+1))/cy(i,j)
    !		enddo
    !	    enddo

    !	    if (nb.ne.nblk) then
    !	    j=ie(nb)-1
    !	    do i=1,i2
    !		dum1(i)=f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* x(i+1,j+1)-cx(i,j)*x(i+2,j+1)-cy(i,j)*x(i+1,j+2)
    !	    enddo
    !        write(*,*) it, nb, 'dum1', sum(abs(dum1(:)))
    !        !write(*,'(8f8.3)') rinv1(:,:, nb)
    !        j=ie(nb)
    !	    do n=1,i2
    !	        dum2(n)=0.
    !	        do m=1,i2
    !    	  dum2(n)=dum2(n)+dum1(m)*rinv1(m,n,nb)
    !	        enddo
    !	        !dum0(n,nb)=x(n+1,j)
    !	        x(n+1,j)=x(n+1,j)-dum2(n)
    !	    end do
    !        write(*,*) it, nb, 'dum2', sum(abs(dum2(:)))
    !        !write(*,*) dum2

    !        endif 

    !	    do n=1,i2
    !	       dum2(n)=0.
    !		   do m=1,i2
    !		       dum2(n)=dum2(n)+dum1(m)*rinv(m,n,nb)
    !		   enddo
    !		   hh(n+1,2, nb)=dum2(n)
    !		   xt(n+1,js+1)=x(n+1,js+1)+dum2(n)
    !	    enddo

    !	    !do m=1,i2
    !	    !    dum1(m)=hh(m+1,2, nb)*cy(m,js-1)
    !	    !enddo

    !        !write(*,*) it, nb, 'dum11', sum(abs(dum1(:)))
    !	    !j=ie(nb-1)
    !	    !do n=1,i2
    !	    !    dum2(n)=0.
    !	    !    do m=1,i2
    !        !      dum2(n)=dum2(n)+dum1(m)*rinv1(m,n,nb-1)
    !	    !    enddo
    !	    !    hh(n+1,1, nb)=dum2(n)
    !	    !enddo
    !        !
    !        !write(*,*) it, nb, 'dum2', sum(abs(dum2(:)))
    !	    !do j=js,jf
    !	    !    do i=1,i2
    !        !        m = j-(nb-1)*j2/nblk
    !        !        !write(*,*) j, m
    !	    !        hh(i+1,m+2, nb)=(-ax(i,j)*hh(i,m+1, nb)-ay(i,j)*hh(i+1,m, nb)-bb(i,j)* hh(i+1,m+1, nb)-cx(i,j)*hh(i+2,m+1, nb))/cy(i,j)
    !	    !        x(i+1,j+2)=x(i+1,j+2)+hh(i+1,m+2, nb)
    !	    !    end do
    !	    !end do
    !	enddo
    !enddo
    
end subroutine ddevp

  !*******************************************************************
  ! elliptic solver preprocessor *
  !*******************************************************************
  subroutine pre(ax,ay,bb,cx,cy,rinv,rinv1,i2,j2, ie, nblk)
    ! ----------------------------------------------------------------------
    ! m2 is the working bir strip width (i2 is evenly divisible by m2)
    implicit none 
    integer, intent(in)                             :: i2,j2, nblk
    integer, dimension(0:nblk), intent(in)            :: ie
    real(8), dimension(i2, j2), intent(in)          :: ax,ay,bb,cx,cy

    real(8), dimension(i2, i2, nblk), intent(inout) :: rinv,rinv1

    logical	:: CheckIfUnitMat

    ! local variables 
    integer     :: jl, nb, jh,  ig, n,ng, i0, i, j, k
    character(LEN=80)	:: fmt
    real(8), dimension(i2+2, j2+2)     :: h
    
    !! test for rep  
    real(8), dimension(i2, i2)     :: htmp

    fmt = "(3f11.4)"
    !write(*,*) "coefficient ax "
    !write(*,fmt) ax
    !write(*,*) "coefficient ay "
    !write(*,fmt) ay
    !write(*,*) "coefficient bb "
    !write(*,fmt) bb
    !write(*,*) "coefficient cx "
    !write(*,fmt) cx
    !write(*,*) "coefficient cy "
    !write(*,fmt) cy
    !write(*,*) "rinv"
    !write(*,fmt) rinv, rinv1
    !write(*,*) "i,j "
    !write(*,*) i2, j2
    !write(*,*) "ie "
    !write(*,*) ie(:)
    i0=i2+2
    jl=1
    do nb= 1, nblk
        jh=ie(nb)
        do ng=1,i2
	    ig = ng +1
            do j = jl, jh+1
            do i = 1,i0
                h(i,j)=0.
            end do 
            end do 

            h(ig,jl+1)=1.0d0

            if (nb.ne.1) then
                do n=1,i2
                    h(n+1,jl)=rinv1(ng,n,nb-1)*cy(ig-1,jl-1)
                end do
            endif 

            do j=jl,jh-2
		do i=1,i2
		    h(i+1,j+2)=-(ax(i,j)*h(i,j+1)+ay(i,j)*h(i+1,j)+bb(i,j)*h(i+1,j+1)+   &
				cx(i,j)*h(i+2,j+1))/cy(i,j)
		    !write(*,*) i, j 
		    !write(*,*) ax(i,j), ay(i,j),bb(i,j), cx(i,j),cy(i,j)
		end do 
		    !write(*,*) 'h on sweep', ng,  j+2
		    !write(*,'(8f8.4)')   h(2:i2+1, j+2)
            end do 

            j=jh-1
            do  i=1,i2
		htmp(ng, i) = h(i+1, j+1)*cy(i,j)
                rinv(ng,i,nb)=ax(i,j)*h(i,j+1)+ay(i,j)*h(i+1,j)+bb(i,j)*             &
				h(i+1,j+1)+cx(i,j)*h(i+2,j+1)
            end do 

            if (nb.ne.nblk) then 
                j=ie(nb)
                do n=1,i2
                    rinv(ng,n,nblk)=h(n+1,j)
                end do
            endif
            h(ig,jl+1)=0.
        end do
	!write(*,*) 'rinv'
	!write(*,'(4f8.3)') rinv(:,:,nb)
	!write(*,*) 'matinv, i2 = ', i2, 'nb =', nb
        call matinv(rinv(:,:,nb),i2)

	!write(*,*) 'htmp'
	!write(*,'(4f8.3)') htmp(:,:)
	!if(CheckIfUnitMat(rinv(:,:,nb), htmp, i2, 1.0D-8)) then 
    	!    print *, "HTMP inverse succeed"
    	!else
    	!    print *, "HTMP inverse failed"
    	!endif 
    
        write(*,111) nb,nblk
        111  format ('processing block #',i3,' out of',i3,' total evp solver blocks')
        
        if (nb.eq.nblk) return
        do i=1,i2
          do j=1,i2
            rinv1(i,j,nb)=0.
            do k=1,i2
              rinv1(i,j,nb)=rinv1(i,j,nb)-rinv(i,k,nb)*rinv(k,j,nblk)
            end do
          end do
        end do
        jl=jh
	!write(*,*) 'rinv1'
	!write(*,'(8f8.3)') rinv1(:,:,nb)

    end do 
  
  end subroutine pre
  !*******************************************************************
  
  subroutine matinv(b,n)
    implicit none 
    logical		       :: CheckIfUnitMat
    integer, intent(in)                       :: n
    real(8), dimension(n,n),intent(inout)     :: b
    
    ! local variables 
    integer                    :: i, j, n1, ip1, ia, ib, im1, astat
    real(8), dimension(n)      :: b1, b2
    real(8), dimension(n,n)    :: bb

   ! write(*,*) 'stat alocate', n
   ! allocate(b1(3), stat =astat)
   ! if(astat /= 0) print *, 'b1 allocate failed' ; stop
    
    bb(:,:) = b(:,:)

    n1=n-1
    do i=1,n1
      b1(1)=1./b(i,i)
!	 write(*,*) 'origin --b1'
!	 write(*,'(4f8.3)') b1
      b(i,i)=1.0
      do j=1,n
        b(i,j)=b(i,j)*b1(1)
      end do
      ip1=i+1
      do ia=ip1,n
        b1(ia)=b(ia,i)
      end do
      do ia=ip1,n
        b(ia,i)=0.
      end do
      do j=1,n
        b2(j)=b(i,j)
      end do
      do ia=ip1,n
        do j=1,n
          b(ia,j)=b(ia,j)-b1(ia)*b2(j)
        end do
      end do
!    write(*,*) 'origin ----bb', i
!    write(*,'(4f8.3)') b
    end do

    b1(1)=1./b(n,n)
    b(n,n)=1.

!    write(*,*) 'origin ----b1'
!    write(*,'(4f8.3)') b1

    do j=1,n
      b(n,j)=b(n,j)*b1(1)
    end do

    do i=2,n
      do ib=1,i
        b1(ib)=b(ib,i)
      end do
      im1=i-1
      do ib=1,im1
        b(ib,i)=0.
      end do
      do j=1,n
        b2(j)=b(i,j)
      end do
      im1=i-1
      do ib=1,im1
        do j=1,n
          b(ib,j)=b(ib,j)-b1(ib)*b2(j)
        end do
      end do
    end do


!    write(*,*) 'inverse'
!    write(*,'(4f8.3)') b
!    write(*,*) 'origin'
!    write(*,'(4f8.3)') bb
    if(CheckIfUnitMat(bb, b, n, 1.0D-5)) then 
	print *, "INV inverse succeed"
    else
	print *, "INV inverse failed"
    endif 
  
  end subroutine matinv

  function CheckIfUnitMat(b, invb, n, eps)

    !implicit none 
    integer, intent(in)	    :: n
    real(8), intent(in)	    :: eps
    real(8), dimension(n,n), intent(in)	    :: b, invb
    logical 	    :: CheckIfUnitMat

    ! local 
    integer	:: i, j, k, errflag
    real(8), dimension(n,n)	    :: cb

    print *, "------Test inverse matrix------"
    print *, eps

    cb(:,:)  = 0.0d0
    errflag = 0
    do i = 1, n
	do j = 1, n

	    do k = 1, n
		cb(i,j) =  cb(i,j)  +  b(i, k)*invb(k,j)
	    end do 
	    if((abs(cb(i,j)) > eps .and. i /=j ) .or. (abs(cb(i,j)-1.0) > eps .and. i ==j )) then
		print *, i, j, cb(i,j)
		errflag = errflag +1
	    endif 
	end do 
    end do 

    if (errflag == 0 ) then 
	print *, "Inverse matrix test passed at tol = ", eps
	CheckIfUnitMat = .true.
    else
	print *, "Inverse matrix test failed at tol = ", eps
	CheckIfUnitMat = .false.
    endif 

    return

  end function CheckIfUnitMat

  function CheckRep(ax,ay,bb,cx,cy,f,x,i2,j2,eps)

    !implicit none 
    integer, intent(in) :: i2,j2
    real(8), dimension(i2, j2), intent(in) :: ax, ay, bb, cx, cy 
    real(8), dimension(i2, j2), intent(in) :: f
    real(8), dimension(i2+2, j2+2), intent(in) :: x
    real(8), intent(in)	    :: eps

    logical 	    :: CheckRep

    ! local 
    integer	:: i, j, k, errflag
    real(8), dimension(i2,j2)	    :: cf

    print *, "------Test rep solution------"

    cf(:,:)  = 0.0d0
    errflag = 0
    do i = 1, i2
	do j = 1, j2
	    cf(i,j) = f(i,j)-ax(i,j)*x(i,j+1)-ay(i,j)*x(i+1,j)-bb(i,j)* x(i+1,j+1)-cx(i,j)*x(i+2,j+1)-cy(i,j)*x(i+1,j+2)
	    if(abs(cf(i,j)) > eps ) then
		print *, i, j, cf(i,j)
		errflag = errflag +1
	    endif 
	end do 
    end do 

    if (errflag == 0 ) then 
	print *, "REP solver test passed at tol = ", eps
	CheckRep = .true.
    else
	print *, "REP solver test failed at tol = ", eps
	CheckRep = .false.
    endif 

    return

  end function CheckRep
  subroutine check_preinv(inv, tmph, n)

    implicit none 
    integer, intent(in)	    :: n
    real(8), dimension(n,n), intent(in)	    ::  inv, tmph

    ! local 
    integer	:: i, j, k, errflag
    real(8), dimension(n,n)	    :: cb

    !print *, "------Test inverse matrix------"

    !cb(:,:)  = 0.0d0
    !errflag = 0
    !do i = 1, n
    !    do j = 1, n

    !        do k = 1, n
    !    	cb(i,j) =  cb(i,j)  +  b(i, k)*invb(k,j)
    !        end do 
    !        if((abs(cb(i,j)) > 1.0e-8 .and. i /=j ) .or. (abs(cb(i,j)-1.0) > 1.0e-8 .and. i ==j )) then
    !    	print *, i, j, cb(i,j)
    !    	errflag = errflag +1
    !        endif 
    !    end do 
    !end do 

    !if (errflag == 0 ) then 
    !    print *, "Inverse matrix test passed!"
    !else
    !    print *, "Inverse matrix test failed!"
    !endif 
    
  end subroutine check_preinv
!********************************************************************
  
