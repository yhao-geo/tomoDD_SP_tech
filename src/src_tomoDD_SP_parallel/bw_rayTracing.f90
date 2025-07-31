! =============================================================================
  subroutine pbr( isp, aas, bbs, hs, aar, bbr, hr, w, ni, tt, &
                  xn,yn,zn,nx,ny,nz,npcom,rak,dak,vpak,vsak,nlayer,xlayer,np15,msg, &
                  maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)    
! =============================================================================
! 3D pseudo-bending for the continuous, spherical earth based on kazuki koketsu's
! algorithm. Ray is tracing from receiver to source. Input coordinates are in 
! degrees and km. Latitude from -90 to 90 and longitude from -180 to 180, depth 
! is radius from the earth's center.
! ni:    the # of segments in the ray path
! w :    stores the coordinates of the segments.
! -----------------------------------------------------------------------------
!  parameters for calculation							
!	xfac   = enhancement factor (see um and thurber, 1987)		
!	nloop  = number of bending iterations
!	n1, n2 = min & max of ray segments
!	mins   = min. length of segment (km)
! -----------------------------------------------------------------------------
	Implicit None
	
	real, intent(in)							::	bld, xl, yl, zl
	integer, intent(in)							::	maxnx,maxny,maxnz,nx, ny, nz, ixkms, iykms, izkms, nlayer, npcom, np15, msg
	real, dimension(maxnx), intent(in)				::	xn
	real, dimension(maxny), intent(in)				::	yn
	real, dimension(maxnz), intent(in)				::	zn
	integer, dimension(ixkms), intent(in)		::	ixloc
	integer, dimension(iykms), intent(in)		::	iyloc
	integer, dimension(izkms), intent(in)		::	izloc
	real, dimension(maxnx,maxny,maxnz*2), intent(in)		::	vel
	real(kind=8), dimension(npcom), intent(in)	::  rak, dak, vpak, vsak
	real(kind=8), dimension(0:nlayer),intent(in)::	xlayer
	
	integer, intent(in)			::	isp
	integer, intent(out)		::	ni
	real, intent(out)			::	tt
	real(kind=8), intent(inout)	::	aas, bbs, hs, aar, bbr, hr
	real(kind=8), intent(out)	::	w(3,msg+1)			! coordinates of segments	
		
	integer				n1, n2, true, nloop, loops, i, j, k, l, kk, idstn
	real(kind=8)		shiftlo	! Hansc
	real(kind=8)		xfac, flim, mins
	real(kind=8)		aas2, as, aar2, ar, ad, bre, bso, dlo, bs, br, rs, rr 
	real(kind=8)		x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3, dx, dy, dz
	real(kind=8)		to, tp, tk
	real(kind=8)		r1, r2, r3, a1, a2, a3, b1, b2, b3
	real(kind=8)		acosa, sina, cosa
	real(kind=8)		dn, ddn, da, db, dr, dseg, ddseg
	real(kind=8)		adv, bdv, rdv 
	real(kind=8)		upz, dwz, v1, v2, v3, vr, vr1, vr2, vb, vb1, vb2, va, va1, va2
	real(kind=8)		pr, pa, pb, vrd, rvr, rva, rvb, rvs, cc, rcur, rdr, rda, rdb, rrp, ap, bp

	real(kind=8), dimension(msg+1)	::	a, b, r
		
	real(kind=8), parameter	::	RNULL = 0.0e10
	real(kind=8), parameter	::	ro = 6371.00	
	real(kind=8), parameter	::	d2r = asin(1.)/ 90.
	real(kind=8), parameter	::	r2d = 90./asin(1.)
	real(kind=8), parameter	::	eps = 0.0
	real(kind=8), external	::	vel2, rtim
    integer         rayError1
! -----------------------------------------------------------------------------
  ! initialization									
	ni = msg+1 ;	n1 = 2 ;		n2 = msg ;	
	xfac = 1.5 ;	flim = 1.e-5 ;	mins = 10. 
	nloop  = 200 ;	dseg = 0.0	
	
  ! Check coordinates 
	if( aas<-90. .or. aas>90. )	then
		write( 6,'(">>> STOP pbr: Latitude of source is out of range ")' )	
		stop
	endif
	if( aar<-90 .or. aar>90.)then
		write( 6,'(">>> STOP pbr: Latitude of station is out of range ")' )
		stop
	endif

    if( bbs<-180 .or. bbs>180.)  then
        write(*,*)'Longitude of source is out of range in pbr.'
        stop
    endif

    if( bbr.lt.-180 .or. bbr.gt.180.)  then
        write(*,*)'Longitude of station is out of range in pbr.'
        stop
    endif
!    write(*,*) aas, bbs, hs, aar, bbr, hr !! lei
!    write(*,*) 'bld,nx,ny,nz,vel(2,2,2):',bld,nx,ny,nz,vel(2,2,2) !! lei
!      do k=1,nz*2
!         do j=1,ny
!            write(*, '(100f7.3)')(vel(i,j,k),i=1,nx)  !! lei
!         enddo
!      enddo
	! Rotate coordinates in order to have longitude and latitude range from 0 to 180. 
	! This program does not work with angles greater than 180. Pass from latitude to colatitude
	!as = (90.00-aas) * d2r !colatitude
	!ar = (90.00-aar) * d2r !colatitude

	! consider the ellipticity correction( coLatitude )
	! -----------------------------------------------------
	!aas2=atan(.99664719*tan(aas*d2r))
	aas2 = atan( 1.00000000*tan(aas*d2r) ) ;		as = asin(1.) - aas2
	aar2 = atan( 1.00000000*tan(aar*d2r) ) ;		ar = asin(1.) - aar2

	! coordinate of longitude
	! -----------------------------------------------------
    if( bbr<0.0 )   then
        bre = 360. + bbr
    else
        bre = bbr
    endif
    
    if( bbs<0.0 )   then
        bso = 360. + bbs
    else
        bso = bbs
    endif

    dlo = abs( bso - bre )

	if( dlo < 180. )then
		shiftlo = 0.0e10
		if( bso < bre )	then
			shiftlo = bso-(180.-dlo)/2.
			bbs = (180.-dlo)/2.
			bbr = bbs+dlo
		else
			shiftlo = bre-(180.-dlo)/2.
			bbr = (180.-dlo)/2.
			bbs = bbr+dlo
		endif
	else
		dlo = 360.0000 - dlo
		shiftlo = 0.0e10
		if( bso < bre )	then
			shiftlo = bso-(dlo+(180.-dlo)/2.)
			bbs = (180.-dlo)/2.+dlo 
			bbr = bbs-dlo
		else    
			shiftlo=bre-(dlo+(180.-dlo)/2.)
			bbr=(180.-dlo)/2.+dlo
			bbs=bbr-dlo
		endif
	endif

	bs = bbs * d2r ;	br = bbr * d2r ;	ad = (as + ar) / 2.								 

	rs = ro - hs ;		rr = ro - hr
									
	! initial straight ray
	! -----------------------------------------------------
	ni = n1
	x1 = rr*sin(ar)*cos(br) ;		y1 = rr*sin(ar)*sin(br) ;		z1 = rr*cos(ar)								 
	x2 = rs*sin(as)*cos(bs) ;		y2 = rs*sin(as)*sin(bs) ;		z2 = rs*cos(as)	 
	dx = (x2-x1) / ni ;				dy = (y2-y1) / ni ;				dz = (z2-z1) / ni

	do j=1,ni+1
		
		x = x1 + dx*(j-1) ;		y = y1 + dy*(j-1) ;		z = z1 + dz*(j-1)								 
		r(j) = sqrt(x**2 + y**2 + z**2) ;		
		acosa= z/r(j)
		
		if( acosa < -1. )	acosa = -1.		
		if( acosa > 1 )		acosa = 1.
		
		a(j) = acos(acosa) ;	acosa = x/r(j)/sin(a(j))
		
		if( acosa < -1. )	acosa = -1.
		if( acosa > 1 )		acosa = 1.
		
		b(j) = acos(acosa)
		if( y < 0.00000 )	b(j) = 360.00000*d2r-b(j)
	
	end do

	to = rtim( isp,shiftlo,msg,ni+1,r,a,b,xn,yn,zn,nx,ny,nz,npcom,rak,dak,vpak,vsak,nlayer,xlayer,np15, &
                 maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
	tp = to

	do i=1,ni+1
		w(1,i) = r(i) ;		w(2,i) = a(i) ;		w(3,i) = b(i)
	end do
	
  ! number of points loop
	loops = 0 ;		true = 1   

	do while( true==1 )
		loops = 0										 								 
		do l=1,nloop				 
			loops = loops + 1			
			do kk=2,ni
			
				if( mod(kk,2)==0 ) then
					k = kk/2 + 1
				else
					k = ni+1 - (kk-1)/2
				endif

				r1 = r(k-1) ;					a1 = a(k-1) ;					b1 = b(k-1)
				x1 = r1*sin(a1)*cos(b1) ;		y1 = r1*sin(a1)*sin(b1) ;		z1 = r1*cos(a1)
				
				r3 = r(k+1) ;					a3 = a(k+1) ;					b3 = b(k+1)
				x3 = r3*sin(a3)*cos(b3) ;		y3 = r3*sin(a3)*sin(b3) ;		z3 = r3*cos(a3)
				dx = x3 - x1 ;					dy = y3 - y1 ;					dz = z3 - z1
				
				x2 = x1 + dx/2 ;				y2 = y1 + dy/2 ;				z2 = z1 + dz/2
				r2 = sqrt(x2**2 + y2**2 + z2**2)				
				acosa = z2 / r2
				
				if( acosa < -1. )	acosa = -1.
				if( acosa > 1 )		acosa = 1.
				a2 = acos(acosa) ;	sina = sin(a2) ;	cosa = cos(a2)								
				
                acosa=x2/r2/sina ;	
                if(acosa < -1.)		acosa=-1.
				if(acosa > 1)		acosa=1.
				
				b2 = acos(acosa) ;	if( y<0.00000 )		b2 = 360.00000*d2r-b2
				
				dn = dx**2 + dy**2 + dz**2 ;			ddn = sqrt(dn)
				
				dr = (r3-r1) / ddn ;	da = (a3-a1) / ddn ;	db = (b3-b1) / ddn

			  ! Begin find the gradients and velocities
			  ! ---------------------------------------------------------------
			  ! first find the length of segment
				dseg  = sqrt( (dx/2)**2 + (dy/2)**2 + (dz/2)**2 )
				ddseg = dseg/2.

			  ! Now ddseg will be a distance to find dV along the coordinates 
!                write(*,*) 'r1,a1,b1:',r1,a1,b1  !! lei	
               	v1 = vel2( isp,shiftlo,r1,a1,b1,1,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				v2 = vel2( isp,shiftlo,r2,a2,b2,2,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				v3 = vel2( isp,shiftlo,r3,a3,b3,3,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				
			  ! Begin to determine coordinates of points surroundibg point a2,b2,r2 at the distance ddseg
				upz = r2 + ddseg ;		dwz = r2 - ddseg
				if( upz > ro )	then
					upz = ro
					dwz = upz-dseg
				endif
				if( dwz<=eps )	then
					dwz = 0.00000001
					upz = ro
				endif

			  ! The following if-endif is just for P & S, thus comment out for SKS & PKP !!!
			  ! This gives the lowermost mantle Vp in the outer core
				if( (dwz-3479.50)<=eps )	then
					dwz = 3479.500
					upz = dwz+dseg
				endif

			  ! find dV along depth, longitude and latitude
			  ! ---------------------------------------------------------------
				vr1 = vel2( isp,shiftlo,upz,a2,b2,4,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				vr2 = vel2( isp,shiftlo,dwz,a2,b2,5,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				vr  = (vr1-vr2) / dseg
		
			  ! find dV along longitude
				call km2deg( a2,b2,r2,ddseg,RNULL,adV,bdV,rdV ) 
				vb2 = vel2( isp,shiftlo,rdV,adV,bdV,6,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				call km2deg( a2,b2,r2,-1.*ddseg,RNULL,adV,bdV,rdV )
				vb1 = vel2( isp,shiftlo,rdV,adV,bdV,7,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				vb  = -1.*(vb1-vb2)/dseg

			  ! find dV along latitude
				call  km2deg( a2,b2,r2,RNULL,ddseg,adV,bdV,rdV )
				va2 = vel2( isp,shiftlo,rdV,adV,bdV,8,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				call  km2deg( a2,b2,r2,RNULL,-1.*ddseg,adV,bdV,rdV )
				va1 = vel2( isp,shiftlo,rdV,adV,bdV,9,xn, yn, zn, nx, ny, nz,&
                 			npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 			maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
				va  = -1.*(va1-va2)/dseg
			 
			  ! spherical velocity gradient
				!va = va / r2
				!vb = vb / r2 / sina
				!(tangential vector) = (slowness vector) / s
				!write(*,*)'perturbations !!!'
				pr = dr ;		pa = r2 * da ;		pb = r2 * sina * db
				
				vrd = pr*vr + pa*va + pb*vb
				rvr = vr - vrd*pr ;		rva = va - vrd*pa ;		rvb = vb - vrd*pb

				rvs = sqrt(rvr*rvr + rva*rva + rvb*rvb)				 

				if( rvs==eps ) then
					r(k) = r2 ;				a(k) = a2 ;				b(k) = b2								 
				else										
					rvr = rvr / rvs ;		rva = rva / rvs ;		rvb = rvb / rvs								 

					cc  = (1./v1+1./v3)/2.						
!                    write(*,*) 'cc,v1,v3:',cc,v1,v3  !! lei	
					rcur = vr*rvr + va*rva + vb*rvb				 
					if( rcur<=eps ) 	rcur=abs(rcur)	
!                    write(*,*) 'cc,rcur:',cc,rcur  !! lei				
					rcur = (cc*v2+1.) / (4.*cc*rcur)				
					rcur = -rcur + sqrt(rcur**2+dn/(8.*cc*v2))     

					rdr = rvr * rcur ;		rda = rva * rcur ;		rdb = rvb * rcur

					rrp  = r2 + rdr ;		ap  = a2 + rda/r2;		bp  = b2 + rdb/(r2*sina)

					r(k) = ( rrp-r(k) )*xfac + r(k)
					a(k) = ( ap-a(k)  )*xfac + a(k)
					b(k) = ( bp-b(k)  )*xfac + b(k)

				endif
			end do								

			idstn = ni
			do j=1,ni+1
				w(1,j) = r(j) ;		w(2,j) = a(j) ;		w(3,j) = b(j)
	 		end do
			ni = idstn
			
		  ! here trace each iteration
			tk = rtim( isp,shiftlo,msg,ni+1,r,a,b,xn,yn,zn,nx,ny,nz,npcom,rak,dak,vpak,vsak,nlayer,xlayer,np15, &
                 maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)

			if( abs(to-tk) <= to*flim )  exit		
			
			to = tk	
											 
		end do									  

	  ! skip increasing of segment number if minimum length of segment is exceed 
	  ! or maximum number of segments was reached
		if( dseg<mins .or. ni>=n2 ) go to 66666
	  
	  ! double the number of points.
		to = tk ;		ni = ni * 2

		do i=1,ni/2+1
			r(i*2-1) = w(1,i) ;		a(i*2-1) = w(2,i) ;		b(i*2-1) = w(3,i)
		enddo
		
		do k=2,ni,2
			r1 = r(k-1) ;					a1 = a(k-1) ;					b1 = b(k-1)						 
			x1 = r1*sin(a1)*cos(b1) ;		y1 = r1*sin(a1)*sin(b1) ;		z1 = r1*cos(a1)						 
			r3 = r(k+1) ;					a3 = a(k+1) ;					b3 = b(k+1)						 
			x3 = r3*sin(a3)*cos(b3) ;		y3 = r3*sin(a3)*sin(b3) ;		z3 = r3*cos(a3)						 
			dx = x3 - x1 ;					dy = y3 - y1 ;					dz = z3 - z1							
			
			x2 = x1 + dx/2 ;				y2 = y1 + dy/2 ;				z2 = z1 + dz/2						
			r2 = sqrt(x2**2 + y2**2 + z2**2)			
			
			acosa = z2/r2
			if( acosa < -1. )		acosa = -1.
			if( acosa > 1 )			acosa = 1.
			a2 = acos(acosa) ;		sina = sin(a2) ;		acosa=x2/r2/sina
			if( acosa < -1. )		acosa=-1.
			if( acosa > 1 )			acosa=1.
			
			b2 = acos(acosa)
			if( y < 0.00000 )		b2 = 360.00000*d2r-b2
			r(k) = r2 ;				a(k) = a2 ;				b(k) = b2
		enddo

		tk = rtim( isp,shiftlo,msg,ni+1,r,a,b,xn,yn,zn,nx,ny,nz,npcom,rak,dak,vpak,vsak,nlayer,xlayer,np15, &
                 maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)

		66666 continue

		if( abs(to-tk) <= to*flim )  exit
		
		to = tk
		 
	enddo

	tt = tk 	! the final output of the traveltime
	idstn = ni
	do i=1,ni+1
		w(1,i) = r(i) ;		w(2,i) = a(i) ;		w(3,i) = b(i)
	enddo
	ni = idstn
  
  ! Return coordinates to the origin
	idstn = ni
	do k=1,ni+1
		w(2,k) = w(2,k)*r2d ;		w(3,k) = w(3,k)*r2d+shiftlo
        if( w(3,k)<0. )             w(3,k) = 360. + w(3,k)
	enddo 
	ni = idstn

	return
  end subroutine pbr	
	

! =============================================================================
  Function rtim( isp, shiftlo, msg, m, r, a, b, &
                 xn,yn,zn,nx,ny,nz,npcom,rak,dak,vpak,vsak,nlayer,xlayer,np15, &
                 maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)	
! =============================================================================
! function to calculate the travel time along a ray path
! -----------------------------------------------------------------------------
	Implicit None

	real, intent(in)							::	bld, xl, yl, zl
	integer, intent(in)							::	maxnx,maxny,maxnz,nx,ny,nz,ixkms,iykms,izkms,nlayer,npcom,np15
	real, dimension(maxnx), intent(in)				::	xn
	real, dimension(maxny), intent(in)				::	yn
	real, dimension(maxnz), intent(in)				::	zn
	integer, dimension(ixkms), intent(in)		::	ixloc
	integer, dimension(iykms), intent(in)		::	iyloc
	integer, dimension(izkms), intent(in)		::	izloc
	real, dimension(maxnx,maxny,maxnz*2), intent(in)		::	vel
	real(kind=8), dimension(npcom), intent(in)	::  rak, dak, vpak, vsak
	real(kind=8), dimension(0:nlayer),intent(in)::	xlayer
    integer rayError1 
	
	integer			::	isp, m, msg
	real(kind=8)	::	shiftlo
	real(kind=8), dimension(msg+1)	::	r, a, b
	
	integer			::	j
	real(kind=8)	::	rtim, rv1, rv2, sm, dl
	real(kind=8)	::	x1, y1, z1, x2, y2, z2
	
	real(kind=8), external	::	vel2
	
	
	if( m > (msg+1) )	write(*,"(/,'--> Pbr: Please increase msg in mod_bodywave.f90.')") 
	
	rtim = 0.		
	
	rv1 = 1./vel2( isp,shiftlo,r(1),a(1),b(1),-1,xn, yn, zn, nx, ny, nz,&
                 	npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 	maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
	
	do j=1,m-1
	
	    x1 = r(j)*sin(a(j))*cos(b(j)) ;	    x2 = r(j+1)*sin(a(j+1))*cos(b(j+1)) 
	    y1 = r(j)*sin(a(j))*sin(b(j)) ;	    y2 = r(j+1)*sin(a(j+1))*sin(b(j+1))
	    z1 = r(j)*cos(a(j)) ;				z2 = r(j+1)*cos(a(j+1))
	    
	    dl = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
	    
	    rv2 = 1./vel2(isp,shiftlo,r(j+1),a(j+1),b(j+1),-2,xn, yn, zn, nx, ny, nz,&
                 		npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 		maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
	    
	    sm = (rv1 + rv2) / 2.	   
	    
	    rtim = rtim + sqrt(dl)*sm
	    
	    rv1 = rv2
	
	end do
									    
	return										
  end function rtim


! =============================================================================
  Function vel2( isp, shiftlo, r, pa, ra, k, xn, yn, zn, nx, ny, nz,&
                 npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                 maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
! =============================================================================
! return the velocity value of given point(r,pa,ra)
! -----------------------------------------------------------------------------
    Implicit None
	
	real, intent(in)							::	bld, xl, yl, zl
	integer, intent(in)							::	maxnx,maxny,maxnz,nx, ny, nz, nlayer, npcom, np15, ixkms, iykms, izkms
	real, dimension(maxnx), intent(in)				::	xn
	real, dimension(maxny), intent(in)				::	yn
	real, dimension(maxnz), intent(in)				::	zn
	real, dimension(maxnx,maxny,maxnz*2), intent(in)		::	vel
	real(kind=8), dimension(npcom), intent(in)	::  rak, dak, vpak, vsak
	real(kind=8), dimension(0:nlayer),intent(in)::	xlayer
	integer, dimension(ixkms), intent(in)		::	ixloc
	integer, dimension(iykms), intent(in)		::	iyloc
	integer, dimension(izkms), intent(in)		::	izloc
    integer rayError1
	
	integer			::	isp, k, nbl
	real(kind=8)	::	r, pa, ra, vel2, shiftlo
	
	real(kind=8)	::	clat, clon, V
	real(kind=8), parameter	::	r2d = 90./asin(1.)									  
	
! -----------------------------------------------------------------------------
	
	! convert to degree and rotate coordinates 
	clat = pa*r2d
	clon = ra*r2d + shiftlo
    if( clon<0 )    clon = 360.+clon  !! lei
!    write(*,*) 'clat,clon:',clat,clon !! lei
	call cellvel( isp, clat, clon, r, V, nbl, xn, yn, zn, nx, ny, nz,&
	              npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                  maxnx,maxny,maxnz,vel, bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)

	vel2 = V
    !write(*,*) 'vel2:',vel2 !! lei
	return
  end function vel2
  

! =============================================================================
  subroutine cellvel( isp, la, lo, era, vl, nbl, xn, yn, zn, nx, ny, nz, &
                      npcom, rak, dak, vpak, vsak, nlayer, xlayer, np15, &
                      maxnx,maxny,maxnz,vel, bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
! =============================================================================
! This subroutine determine the velocity in the block using co-Lat, co-Lon and R.
! Lat and Lon in degrees, era is the radius, dp is the depth
! Note, the coordinates are colongitude and colatitude
! ------------------------------------------------------------------------------
! The inconsistance of coordinate is really a pain in the ass, it took me two days 
! to finally find where is the little bug. -- added by Hongjian @ MIT 2014/08/21.
! -----------------------------------------------------------------------------	
	Implicit None
	
	real, intent(in)							::	bld, xl, yl, zl
	integer, intent(in)							::	maxnx,maxny,maxnz,nx,ny,nz,nlayer,npcom,np15,ixkms,iykms,izkms
	real, dimension(maxnx), intent(in)				::	xn
	real, dimension(maxny), intent(in)				::	yn
	real, dimension(maxnz), intent(in)				::	zn
	real, dimension(maxnx,maxny,maxnz*2), intent(in)		::	vel
	integer, dimension(ixkms), intent(in)		::	ixloc
	integer, dimension(iykms), intent(in)		::	iyloc
	integer, dimension(izkms), intent(in)		::	izloc
	real(kind=8), dimension(npcom), intent(in)	::  rak, dak, vpak, vsak
	real(kind=8), dimension(0:nlayer),intent(in)::	xlayer
    integer rayError1
	
	integer				::	isp, nbl
	real(kind=8)		::	la, lo, era, dp, vl
	
	integer				::	i, nlr, nla, nlo, ilr, num, ip, jp, kp, kpg
	real(kind=8)		::	d2r, r2d, Re, atemp
	real	            ::	tlat, tlon, tdep, x, y, z, v
	real, dimension(8)	::	wv

	integer, parameter		::	nxi=72, nyi=36, nzi=16
	real(kind=8), parameter	::	dx = 5., dy = 5.

  ! --------------------------------------------------
	d2r = asin(1.)/ 90. ;		r2d = 90./asin(1.)

	! Firstly determine number of layer
	Re = 6371.00 ;				dp = Re - era 	! Depth of the point
    
	nlr = 0 ;		nla = 0 ;		nlo = 0

	if( dp<0.00 )	dp = 0.0

    do while(lo>=180) ;     lo = lo-360. ;      enddo    !!lei
    do while(lo<-180) ;   lo = lo+360. ;      enddo      !!lei
    if (abs(lo).gt.180) lo = lo-lo/abs(lo)*360           !!lei
  ! --------------------------------------------------------
  ! change co-latitude back to latitude, it should be from -90 to 90
	la = 90-la 	
	
  ! Now need to consider the ellipticity correction
	!atemp = 1.0/.99664719;
	atemp = 1.0 ;	la = atan( atemp*tan(la*d2r) )*r2d

  ! calculate the layer number for the AK135 model
	do i=2,np15
		if( dp>=dak(i-1) .and. dp<dak(i) )	then
			ilr = i-1
                        !write(*,*) 'ilr,np15,dp,dak(i-1),dak(i):',ilr,np15,dp,dak(i-1),dak(i)  ! lei   
			exit
		endif
           !write(*,*) 'i,np15,dak(i),vpak(i),vsak(i):',i,np15,dak(i),vpak(i),vsak(i)    
	enddo 
    !write(*,*) 'la,lo,dp:',la,lo,dp !! lei
  ! Where is the point --> inside/outside of the region
  ! ------------------------------------------------------
	if( lo>xn(1)+0.1 .and. lo<xn(nx)-0.1 .and. la>yn(1)+0.1 .and. la<yn(ny)-0.1 ) 	then	
	! inside of the study region			
		
		if( dp<zn(nz)-1 )	then
			! Total inside: decide the index number of the box: use Prothero's intmap here
			x = lo ;	y = la ;	z = dp

			call vel3( isp, x, y, z, v, wv, nx, ny, nz, xn, yn, zn, maxnx,maxny,maxnz,vel, &
                   bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,ip,jp,kp,kpg,rayError1)	! Hansc		
			vl = v

!			call intmap( isp,x,y,z,ip,jp,kp,kpg,nz,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc )	! Hansc
			
!			if( isp==1 )	kp = kp+nz		! ???	It has been already updated in vel3, Hansc

			nbl = ip + (jp-1)*nx + (kp-1)*nx*ny
		else
			!write(*,*)'Velocity determination is deeper than regional nodes, go to the global scale'
	   		!write(*,*)'outside the small region!'
			ip = 0 ;	jp = 0 ;	kp = 0
			
			! just divide the Earth into many blocks
			do i=1,nzi
				if( dp>=xlayer(i-1) .and. dp<xlayer(i) )then
					nlr = i
					exit
				endif
			enddo

			do i=1,nxi
				if( lo>=real(dx*real(i-1)) .and. lo<real(dx*real(i)) )	then
					nlo = i
					exit
				endif
			enddo

			! and then Row (lat) number: note here la is latitude, not co-latitude
			do i=1,nyi
				if( la>=dy*real(i-1) .and. la<dy*real(i) )	then
					nla = i
					exit
				endif
			enddo

			num = (nla-1)*nxi + nlo + (nlr-1)*nxi*nyi + nx*ny*nz
	
			if( isp==0 ) 	vl = vpak(ilr)
			if( isp==1 ) 	vl = vsak(ilr)
                        !write(*,*) 'ilr,vl:',ilr,vl  !lei
			nbl = num
			
		end if

	else 	
		! off the region of study
		print *, lo, la, nx, ny, xn(1), xn(nx), yn(1), yn(ny)
	    write(*,'(">>> STOP: Ray is out of the region, change MOD")')
	    stop    
	endif

	return
  end subroutine cellvel
  
  
! =============================================================================
  subroutine vel3( isp, x, y, z, v, wv, nx, ny, nz, xn, yn, zn, maxnx,maxny,maxnz,vel, &
                   bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,ip,jp,kp,kpg,rayError1)
! =============================================================================
! This routine is Cliff Thurber's. Return the velocity value of given coor(x,y,z)
! -----------------------------------------------------------------------------	
	Implicit None
	
	integer, intent(in)						::	maxnx,maxny,maxnz,nx, ny, nz, ixkms, iykms, izkms
	real, intent(in)						::	bld, xl, yl, zl
	real, dimension(maxnx), intent(in)			::	xn
	real, dimension(maxny), intent(in)			::	yn
	real, dimension(maxnz), intent(in)			::	zn
	real, dimension(maxnx,maxny,maxnz*2), intent(in)	::	vel
	integer, dimension(ixkms), intent(in)	::	ixloc
	integer, dimension(iykms), intent(in)	::	iyloc
	integer, dimension(izkms), intent(in)	::	izloc
    integer rayError1
	
	
	integer, intent(in)				::	isp
	real            				::	x, y, z
	real, intent(out)				::	v
	real, dimension(8), intent(out)	::	wv	
	
	integer				::	ip, jp, kp, kpg, ip1, jp1, kp1
	real				::	xf, yf, zf, xf1, yf1, zf1
  ! ------------------------------------------------------------------------------	

    if(abs(z-zn(nz)).lt.1e-06) then
       z=zn(nz)-0.1
    endif

    if(x.ge.180) x=x-360

	call intmap( isp,x,y,z,ip,jp,kp,nz,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc )
    if (ip.lt.1 .or. jp.lt.1 .or. kp.lt.1) then
    rayError1 = 1
    v = 1
    write(*,*) 'vel3',rayError1
    else	
	ip1 = ip+1 ;	jp1 = jp+1 ;	kp1 = kp+1
	wv = 0.0	
  ! -----------------------------------------------------------------------------	
  ! Note: there seems a bug here, xf yf are not at the same level as zf, so this 
  ! 	    kind of interpolation is problematic.
  ! modified by Hongjian @ mit 2014/08/19
  ! -----------------------------------------------------------------------------	
	!yf=(y-yn(jp))/(yn(jp1)-yn(jp))
	!zf=(z-zn(kp))/(zn(kp1)-zn(kp))
	xf = ( x-xn(ip) ) / ( xn(ip1)-xn(ip) )
	yf = ( y-yn(jp) ) / ( yn(jp1)-yn(jp) )
	zf = ( z-zn(kp) ) / ( zn(kp1)-zn(kp) )
	xf1= 1.0 - xf ;		yf1= 1.0 - yf ;		zf1= 1.0 - zf

	wv(1) = xf1*yf1*zf1 ;		wv(2) = xf*yf1*zf1
	wv(3) = xf1*yf*zf1  ;		wv(4) = xf*yf*zf1
	wv(5) = xf1*yf1*zf  ;		wv(6) = xf*yf1*zf
	wv(7) = xf1*yf*zf   ;		wv(8) = xf*yf*zf

  ! calculate velocity
  ! S-velocity is stored after P-velocity (or V*Q if iuseq=1)
	kpg = kp 		
	if( isp==1 ) 	kp = kp+nz
	kp1 = kp+1
	
	v = wv(1)*vel(ip,jp,kp)  + wv(2)*vel(ip1,jp,kp)		&
      + wv(3)*vel(ip,jp1,kp) + wv(4)*vel(ip1,jp1,kp)	&
      + wv(5)*vel(ip,jp,kp1) + wv(6)*vel(ip1,jp,kp1)	&
      + wv(7)*vel(ip,jp1,kp1)+ wv(8)*vel(ip1,jp1,kp1)
       if (v.le.0) then  !! lei
       print *,'ERROR: neg vel in vel3',v
       print *,ip,jp,kp
       print *,wv(1),wv(2),wv(3),wv(4),wv(5),wv(6),wv(7),wv(8)
       print *,vel(ip,jp,kp),vel(ip1,jp,kp),vel(ip,jp1,kp),vel(ip1,jp1,kp)
       print *,vel(ip,jp,kp1),vel(ip1,jp,kp1),vel(ip,jp1,kp1),vel(ip1,jp1,kp1)
!        stop
       endif
    endif

	return
  end subroutine vel3


! =============================================================================
  subroutine intmap( isp, x, y, z, ip, jp, kp, &
                     nz,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
! =============================================================================
! Modified by W. Prothero so a single call can get the indices.
! For given coordinate(x,y,z), return the index of vel(ip,jp,kp)
! If an array element=0, the position is off the map.
! -----------------------------------------------------------------------------
	Implicit None
	
	real, intent(in)						::	bld, xl, yl, zl
	integer, intent(in)						::	nz, ixkms, iykms, izkms
	integer, dimension(ixkms), intent(in)	::	ixloc
	integer, dimension(iykms), intent(in)	::	iyloc
	integer, dimension(izkms), intent(in)	::	izloc
    integer rayError1
	
	integer, intent(in)		::	isp
	real, intent(in)		::	x, y, z
	integer, intent(out)	::	ip, jp, kp	
	real					::	bldk

	ip = int( (x+xl)/bld ) ;
	jp = int( (y+yl)/bld ) ;
	!bldk = 1.0
	bldk = bld
	kp = int( (z+zl)/bldk) ;

    if (ip.lt.1 .or. jp.lt.1 .or. kp.lt.1) then
        ! this ray point is out of the boundary. --Hao Guo.
        write(*,*) 'intmap',ip,jp,kp,x,y,z
        rayError1 = 1
        ip = 0
        jp = 0
        kp = 0
        write(*,*) 'intmap',rayError1
    else
        ip=ixloc(ip)
        jp=iyloc(jp)
        kp=izloc(kp)
    endif
 
!	kpg = kp
!	if( isp==1 )	kp = kp + nz

	return
  end subroutine intmap


! =============================================================================
  subroutine km2deg( ala,alo,adp,dx,dy,bla,blo,bdp )
! =============================================================================
! This subroutine calculate position of new point in polar coordinates basing on 
! the coordinates of main point in radians ( la is colatitude) and dx and dy in km.
! -----------------------------------------------------------------------------
	Implicit None
	
	real(kind=8)	::	ala, alo, adp, dx, dy, bla, blo, bdp
	
	real(kind=8)	::	d2r, dps
	
	d2r = asin(1.)/ 90.
	dps = adp*sin(ala)
	blo = alo + atan2(dx,dps)
	bla = ala + atan2(dy,adp)
	
	if( bla>(180.*d2r) )	then
		bla = 360.*d2r-bla
		blo = blo+180.*d2r
	endif
	if( bla<0.)	then
		bla = abs(bla)
		blo = blo+180.*d2r
	endif
	
	if( blo<0. )			blo = 360.*d2r+blo
	if( blo>(360.*d2r) )	blo = blo-(360.*d2r)
	
	bdp = sqrt(adp**2+dx**2+dy**2)
	
  end subroutine km2deg


! =============================================================================
  subroutine ttmder( isp, rp, nrp, dth1, dtm, nx, ny, nz, xn, yn, zn, maxnx,maxny,maxnz,vel, bld,xl,yl,zl,&
                   	 ixkms,iykms,izkms,ixloc,iyloc,izloc,maxpar,nparvi,wlat, wlon, theta, &
                     nx1, ny1, nz1, nx2, nxy2, stepl, nfix, ndexfx, imerge, jequal,rayError1)
! =============================================================================
! Calculate hypocenter & slowness partial derivatives
! -----------------------------------------------------------------------------
	Implicit None
	integer, intent(in) 					::	isp, nrp
	real, dimension(3,100000), intent(in)	::	rp
	real, dimension(4), intent(out)			::	dth1
	
	integer, intent(in)						::	maxnx,maxny,maxnz,nx,ny,nz,ixkms,iykms,izkms,maxpar,nparvi
	integer, intent(in)						::	nx1, ny1, nz1, nx2, nxy2
	real, intent(in)						::	bld, xl, yl, zl, wlat, wlon, theta, stepl
	real, dimension(maxnx), intent(in)			::	xn
	real, dimension(maxny), intent(in)			::	yn
	real, dimension(maxnz), intent(in)			::	zn
	real,dimension(maxnx,maxny,maxnz*2),intent(in)	::	vel
	integer, dimension(ixkms), intent(in)	::	ixloc
	integer, dimension(iykms), intent(in)	::	iyloc
	integer, dimension(izkms), intent(in)	::	izloc
	integer, dimension(maxpar),intent(in)	::	nfix, ndexfx, imerge, jequal
	real, dimension(maxpar), intent(inout)	::	dtm
    integer rayError1

	integer					::	ip, jp, kp, kpg	! Hansc
	integer					::	nrp1, i, i1, is, izg, kk, kk1, jj, jj1, ii, ii1
	integer					::	ijk, ini, ne, no, j, inp, nseg
	real                	::	xc, yc, zc, tlat, tlon, tdep, us, ds, uds, half, rx, ry, rz, sl, ssl
	real					::	dxs, dys, dzs, xp, yp, zp, v, dt, fnsegi
	real        			::	dx, dy, dz, tt
	
	integer,dimension(8)	::	in
	real, dimension(8)		::	wv	! Hansc
	
  ! ---------------------------------------------------------------------------
  
  ! Calculate travel time derivatives with respect to hypocentral parameters 
  ! --> use coords of first two points on raypath determine slowness at source.
	xc = rp(1,1) ;		yc = rp(2,1) ;		zc = rp(3,1)

  	! change it into lat-lon-dep to find velocity
	call car2sph_ft( xc, yc, zc, wlat, wlon, tlat, tlon, tdep, theta )
      if(tlon .ge. 180) tlon = tlon-360!!!added by lei
      !write(*,*)tlat,tlon,tdep,xl,yl,zl
      if( tlon.le.(xn(1)+0.1) .or. tlon.ge.(xn(nx)-0.1) .or.&
         tlat.le.(yn(1)+0.1) .or. tlat.ge.(yn(ny)-0.1) .or.&
         tdep.le.(zn(1)+1.0) .or. tdep.ge.(zn(nz)-1.0) ) return


	call vel3( isp, tlon, tlat, tdep, v, wv, nx, ny, nz, xn, yn, zn, maxnx,maxny,maxnz,vel, &
                   bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,ip,jp,kp,kpg,rayError1)  ! Hansc             
	us = 1.0/v

  	! determine cartesian derivatives from direction cosines
	dx = rp(1,2)-rp(1,1) ;		dy = rp(2,2)-rp(2,1) ;		dz = rp(3,2)-rp(3,1)
	
	ds = sqrt(dx*dx+dy*dy+dz*dz) ;		                    uds= -us/ds
	
  	! origin time derivative, hypocentral derivatives
	dth1(1)=1.0 ;		dth1(2) = uds*dx ;		dth1(3) = uds*dy ;		dth1(4) = uds*dz
	
	! ---------------------------------------------------------------------------
  	! skip calculating par-deriv if all velocity nodes are fixed (nparvi=0)
	if( nparvi==0 ) 	return
	! ---------------------------------------------------------------------------


  ! Calculate travel time derivatives with respect to slowness
  	! travel time and velocity partial derivatives
	tt = 0.0 ;              half = 0.5
	
  	! initialize the velocity derivative
	dtm = 0.0
	
	! loop over segments comprising the ray path
	nrp1 = nrp-1
	!pl = 0.0	
	do i=1,nrp1

		i1 = i+1	
		
		rx = rp(1,i) ;		ry = rp(2,i) ;		rz = rp(3,i)
		dx = rp(1,i1)-rx ;	dy = rp(2,i1)-ry ;	dz = rp(3,i1)-rz
		sl = sqrt(dx*dx+dy*dy+dz*dz)	! segment length
		!pl = pl+sl

		! decide on number of subsegments and compute length
		nseg = nint(sl/stepl)+1 ;		fnsegi = 1.0/float(nseg) ;		ssl = sl*fnsegi
		dxs = dx*fnsegi ;	            dys = dy*fnsegi ;	            dzs = dz*fnsegi
		xp = rx-half*dxs;	            yp = ry-half*dys;	            zp = rz-half*dzs

		! loop over subsegments
		do is=1,nseg
			
			xp = xp+dxs ;		yp = yp+dys ;	zp = zp+dzs
			
			! now xp, yp, zp are X-Y-Z, change them into lat-lon-dep 
			call car2sph_ft( xp, yp, zp, wlat, wlon, tlat, tlon, tdep, theta )
			if(tlon .ge. 180) tlon = tlon-360!!!added by lei
		  ! ---------------------------------------------------------------------------
		  ! now it is possible that the ray may be outside the effective domain
		  ! first determine if tlat, tlon and tdep are in the domain.
			if( tlon<=(xn(1)+0.1) .or. tlon>=(xn(nx)-0.1 ) .or.		&
     			tlat<=(yn(1)+0.1) .or. tlat>=(yn(ny)-0.1 ) .or.		&
     			tdep<=(zn(1)+1.0) .or. tdep>=(zn(nz)-1.0 ) ) 	return
     	  ! ---------------------------------------------------------------------------
			
!		  	call intmap( isp,tlon,tlat,tdep,ip,jp,kp,kpg,nz,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc )	! Hansc

			call vel3( isp, tlon, tlat, tdep, v, wv, nx, ny, nz, xn, yn, zn, maxnx,maxny,maxnz,vel, &
                   bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,ip,jp,kp,kpg,rayError1)	! Hansc
		
			dt = ssl/v ;	tt = tt + dt

			! The next section is a change from 'block' to 'linear' partial derivatives.
			! Nodes with non-zero weight
            if(ip .lt. 1 .or. jp .lt. 1 .or. kp .lt. 1) then
              in(1)=1
              in(2)=1
              in(3)=1
              in(4)=1
              in(5)=1
              in(6)=1
              in(7)=1
              in(8)=1
              rayError1=1
              write(*,*) 'ttmder regional', rayError1
            else
			in(1) = ip-1 + nx2*(jp-2) + nxy2*(kp-2)-nxy2*(2*isp) ;	in(2) = in(1)+1
			in(3) = in(1)+nx2 ;			in(5) = in(1)+nxy2 ;		in(4) = in(3)+1
			in(6) = in(5)+1 ;			in(7) = in(5)+nx2 ;			in(8) = in(7)+1
            endif

			! Assign zero weight to boundary nodes. ( these nodes are not included in inversion, 
			! but are in the velocity array, thus we want to avoid writing to negative or incorrect 
			! elements of the partial derivative matrix )
			if( ip==1 ) then
				wv(1)=0.0 ;			wv(3)=0.0 ;		wv(5)=0.0 ;		wv(7)=0.0
			else
				if( ip==nx1 ) then
					wv(2)=0.0 ;		wv(4)=0.0 ;		wv(6)=0.0 ;		wv(8)=0.0
				end if
			endif

			if( jp==1 ) then
				wv(1)=0.0 ;			wv(2)=0.0 ;		wv(5)=0.0 ;		wv(6)=0.0
			else
				if( jp==ny1 ) then
					wv(3)=0.0 ;		wv(4)=0.0 ;		wv(7)=0.0 ;		wv(8)=0.0
				endif
			endif

			if( (kpg==1) .or. (kpg==(nz1+1)) ) then
				do izg=1,4 ;		wv(izg)=0.0 ;		end do
			else
				if( (kpg==nz1) .or. (kpg==(2*nz1)) ) then
					do izg=5,8 ;	wv(izg)=0.0 ;		end do
				endif			
			endif

			! Accumulate model partial derivatives
			do kk=1,2
				kk1 = kk-1
				do jj=1,2
					jj1 = jj-1
					do ii=1,2
					
						ii1 = ii-1 ;	ijk = ii+2*jj1+4*kk1
						
						! skip boundary nodes
						if( wv(ijk)<0.05 ) 	cycle	
						
						! ----------------------------------------------------------------------------
						! skip fixed nodes
						if( nfix(in(ijk))==1 )	cycle
						
						ini = ndexfx( in(ijk) )

						if( imerge(in(ijk))==1 ) 	ini = ndexfx( jequal(in(ijk)) )

						! check for writing to an array element that is outside of the inversion solution array
						if( (ini<1) .or. (ini>nparvi) ) then
							write(16,1606) 	ini, ijk
							write(16,1603) 	ne,no,xp,yp,zp,v,ip,jp,kp,kpg
							write(16,1607) 	( j,in(j),j,wv(j),j=1,8 )
							write(16,'(">>> STOP (to avoid writing outside of defined DTM array")')
							stop
						end if

						inp = ini
						
						! ----------------------------------------------------------------------------
						! DEP: Now include weight factor for hit(DWS)
						!if((ne>(neqs+nbls)).and.(nit==0)) then
							!hit(in(ijk))=hit(in(ijk))+wv(ijk)*wtsht
						!else
							!hit(in(ijk))=hit(in(ijk))+wv(ijk)*wtcomb(no,ne)
						!endif
						!hit(in(ijk))=hit(in(ijk))+wv(ijk)*TimeWeight
						! ----------------------------------------------------------------------------
						
						! set up the model derivative matrix
						dtm(inp) = dtm(inp) + wv(ijk)*ssl
						
					end do
				end do
			end do
	 	end do
	end do

  ! format for output
1606    format('>>> Error in TTMDER, accessing', ' gridpoint outside of velocity inversion',    &
               ' gridpoints, ini=',i5,', ijk=',i5,/,22x,'Probably boundary gridpoints are',     &
               ' too close (TTMDER tries to write DTM',' elements with wv >= 0.05)')

1603    format(' ne=',i5,', no=',i5,', xp=',f8.2,', yp=',f8.2,', zp=',f8.2,', v=',f8.3,/,       &
                 21x,'ip=',i6,',   jp=',i6,',   kp=',i6,',   kpg=',i6) 
 
1607    format(' in(',i1,')=',i6,' wv(',i1,')=',e15.5)
 
  end subroutine ttmder


! ================================================================
!  subroutine bldmap
! ================================================================
!	use bw_fwrd, only : bld, xl, yl, zl, xn, yn, zn, nx, ny, nz, &
!                        ixkms, iykms, izkms, ixloc, iyloc, izloc
!	Implicit None

!    include 'RaySPDR.inc'
  ! local variables
!	integer	::	ixmax, iymax, izmax, ix, iy, iz, ix1, iy1, iz1, i
!	real	::	bldk, xnow, ynow, znow

!	bldk= 1.0 	! decide the depth step for ray tracing
	!bldk = bld
!	xl 	= bld - xn(1) ;		ixmax = (xn(nx)+xl)/bld
!	yl	= bld - yn(1) ;		iymax = (yn(ny)+yl)/bld
!	zl	= bldk - zn(1) ;	izmax = (zn(nz)+zl)/bldk

  ! Check for array size overflow
!	ixkms = ixmax + 100 ;	iykms = iymax + 100 ;	izkms = izmax + 100
!	allocate( ixloc(ixkms), iyloc(iykms), izloc(izkms) )
	
!	ix = 1
!	do i=1,ixmax
!		ix1 = ix+1
!		xnow= float(i)*bld-xl
!		if( anint(xnow*100) >= anint(xn(ix1)*100) ) 	ix = ix1
!		ixloc(i) = ix
!	end do	
!	do	i=ixmax,ixkms ;		ixloc(i) = 0 ;		end do		! Fill remainder of array with zeroes.

!	iy=1
!	do i=1,iymax
!		iy1 = iy+1
!		ynow=float(i)*bld-yl
!		if( anint(ynow*100) >= anint(yn(iy1)*100) ) iy = iy1
!		iyloc(i)=iy
!	 end do
!	do i=iymax,iykms ;		iyloc(i) = 0 ;		end do		! Fill remainder of array with zeroes.

!	iz=1
!	do i=1,izmax
!		iz1=iz+1
!		znow=float(i)*bldk-zl
!		if( anint(znow*100) >= anint(zn(iz1)*100) ) iz=iz1
!		izloc(i)=iz
!	end do
!	do i=izmax,izkms ;		izloc(i) = 0 ;		end do		! Fill remainder of array with zeroes.
	
!  end subroutine bldmap

! =============================================================================
!  subroutine find_id( neve, nobs, evesta, staID, evID, k )
! =============================================================================
! find the column index of staID that received the evID-th.
! -----------------------------------------------------------------------------
!    Implicit None
!	integer, intent(in)                         ::  neve, nobs, staID, evID
!	integer, intent(out)                        ::  k
!	integer, dimension(neve,nobs+1), intent(in) ::  evesta

!	integer     ::  i

!	k = 0
!	do i=1,evesta(evID,1)
!		if( staID == evesta(evID,i+1) ) then
!			k = i
!			exit
!		endif
!	enddo

!	return
!  end subroutine

