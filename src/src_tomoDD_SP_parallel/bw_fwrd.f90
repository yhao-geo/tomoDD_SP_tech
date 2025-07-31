! =============================================================================
       subroutine bw_forward(MAXEVE,MAXSTA,MAXOBS,MAXNODE,eve_sta,wlat,wlon,theta,stepl,&
                pgood,tmp_ttp,tmp_xp,tmp_yp,tmp_zp,tmp_vp,tmp_vp_index,&
                sgood,tmp_tts,tmp_xs,tmp_ys,tmp_zs,tmp_vs,tmp_vs_index,&
                nsta,sta_lat,sta_lon,sta_elv,sta_lab,&
                nsrc,src_lat,src_lon,src_dep,src_type,&
                xn,yn,zn,nx,ny,nz,maxnx,maxny,maxnz,vel,&
                bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,&
                nx1, ny1, nz1, nx2, nxy2, nfix, ndexfx, imerge, jequal,&
                npari,nparvi,maxpar,mxpari,msg,np15,npcom,nlayer,MAXSEG1,&
                rak,vpak,vsak,dak,xlayer,&
                iuses,rayErrorp,rayErrors,DISTratio,tmp_sp_close)
! =============================================================================
! Calculate hypocenters and slowness partial derivatives for body wave.
! ---------------------------------------------------------------------------
	use omp_lib
!	use tomoFDD
	Implicit None
!    include 'RaySPDR.inc'
    integer     eve_sta(MAXEVE,MAXOBS+1)
	real		wlat,wlon,theta,stepl
    integer     pgood(MAXOBS,MAXEVE),sgood(MAXOBS,MAXEVE)
    integer     nsta,nsrc
	real		tmp_ttp(MAXOBS,MAXEVE)
	real		tmp_tts(MAXOBS,MAXEVE)
	real		tmp_xp(MAXOBS,MAXEVE)
	real		tmp_yp(MAXOBS,MAXEVE)
	real		tmp_zp(MAXOBS,MAXEVE)
	real		tmp_xs(MAXOBS,MAXEVE)
	real		tmp_ys(MAXOBS,MAXEVE)
	real		tmp_zs(MAXOBS,MAXEVE) 
	character	sta_lab(MAXSTA)*7
	real		sta_lat(MAXSTA)
	real		sta_lon(MAXSTA)
    real        sta_elv(MAXSTA) 
	doubleprecision	src_lat(MAXEVE)
	doubleprecision	src_lon(MAXEVE) 
	real		src_dep(MAXEVE)
    integer     src_type(MAXEVE) 
    real*8 rak(npcom),vpak(npcom),vsak(npcom),dak(npcom)
    real*8 xlayer(0:nlayer)
    integer         rayErrorp(MAXOBS,MAXEVE)    ! first arrival  !!!lei
    integer         rayErrors(MAXOBS,MAXEVE)    ! first arrival
    real            DISTratio
    integer         tmp_sp_close(MAXOBS,MAXEVE)


    integer     npari,nparvi,MAXEVE,MAXSTA,MAXOBS,MAXNODE 
!    real,allocatable ::            tmp_vp(:,:,:)
!    integer,allocatable ::         tmp_vp_index(:,:,:)
!    real,allocatable ::            tmp_vs(:,:,:)
!    integer,allocatable ::         tmp_vs_index(:,:,:) 
    real            tmp_vp(MAXOBS,MAXEVE,MAXNODE)
    integer         tmp_vp_index(MAXOBS,MAXEVE,MAXNODE+1)
    real            tmp_vs(MAXOBS,MAXEVE,MAXNODE)
    integer         tmp_vs_index(MAXOBS,MAXEVE,MAXNODE+1) 
    integer nx, ny, nz,maxnx,maxny,maxnz
    real bld,xn(maxnx),yn(maxny),zn(maxnz),vel(maxnx,maxny,maxnz*2)
    real xl,yl,zl
    integer ixloc(ixkms), iyloc(iykms), izloc(izkms),ixkms, iykms, izkms
    integer nx1, ny1, nz1,nx2,nxy2
    integer maxpar,mxpari,msg,np15, npcom, nlayer,nrp,MAXSEG1
    integer nfix(maxpar),ndexfx(maxpar),imerge(maxpar),jequal(maxpar)
    real rp(3,MAXSEG1)
	real,dimension(mxpari)	::	dtm	
    integer iuses,isp

	
	real			::	tlat, tlon, tdep, xc, yc, zc, dtemp
	real(kind=8)	::	aar, bbr, hr, aar0, bbr0, hr0, aas, bbs, hs, aas0, bbs0, hs0
	integer			::	iunit, trimlen, i, j, k, staID, evID, nseg2, i1, m, n				

	real, dimension(4)							::	dthP, dthS
	real(kind=8), dimension(3,msg+1)			::	rp_geo
	real(kind=8), parameter						::	r2d=90./asin(1.)
    real(kind=8), parameter                 	::  dpi=asin(1.)/90.
	real(kind=8), parameter						::	r0=6371.0
    integer         rayError1
    integer         nray1, nray2,k1,k2
    real            ray_length, haus_dist, seg_length
    real            ray1(3,msg+1), ray2(3,msg+1)
  ! ---------------------------------------------------------------------------	
  ! Calculate travel-time grid for station i, i.e. treating station i as a "source" first

!    	allocate(tmp_vp(MAXOBS,MAXEVE,MAXNODE))
!        if (.not.allocated(tmp_vp)) print *, "Memory allocation ERROR for tmp_vp"
!        allocate(tmp_vp_index(MAXOBS,MAXEVE,MAXNODE+1))
!        if (.not.allocated(tmp_vp_index)) print *, "Memory allocation ERROR for tmp_vp_index"

!        allocate(tmp_vs(MAXOBS,MAXEVE,MAXNODE))
!        allocate(tmp_vs_index(MAXOBS,MAXEVE,MAXNODE+1))
!        if (.not.allocated(tmp_vs)) print *, "Memory allocation ERROR for tmp_vs"
!        if (.not.allocated(tmp_vs_index)) print *, "Memory allocation ERROR for tmp_vs_index"

  ! ---------------------------------------------------------------------------
	   !$omp parallel &
	   !$omp default(private) &
	   !$omp shared( MAXEVE,MAXOBS,eve_sta,npari,nparvi,wlat,wlon,theta,stepl,iuses ) &
	   !$omp shared( pgood, tmp_ttp,tmp_xp,tmp_yp,tmp_zp,tmp_vp,tmp_vp_index ) &
	   !$omp shared( sgood, tmp_tts,tmp_xs,tmp_ys,tmp_zs,tmp_vs,tmp_vs_index ) &
	   !$omp shared( nsta,sta_lat,sta_lon,sta_elv,sta_lab,nsrc,src_lat,src_lon,src_dep,src_type ) &
	   !$omp shared( xn,yn,zn,nx,ny,nz,rak,dak,vpak,vsak,xlayer ) &
	   !$omp shared( maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc ) &
	   !$omp shared( maxpar,mxpari,nx1, ny1, nz1, nx2, nxy2, nfix, ndexfx, imerge, jequal,msg,np15,npcom,nlayer,MAXSEG1 ) &
       !$omp shared( rayErrorp,rayErrors,DISTratio,tmp_sp_close )
	   !$omp do collapse(2)

	   do i=1,nsta
!	   write(*,*)i, sta_lab(i),sta_lat(i),sta_lon(i)
!       Calculate travel-time grid for station i, i.e. treating 
!       station i as a "source" first

           ! limit the range of staion's location ! Bug fixed by Hao Guo

 
	    do j=1, nsrc
!           write(*,*)i, sta_lab(i),sta_lat(i),sta_lon(i),j,src_lat(j),src_lon(j),src_dep(j)
           

!z---- change the program here, do not calculate the partial derivatives
!z---- for every station and event pair. Instead only calculates the 
!z---- partial derivatives between actual event and station pair. 
              
!z---- check if this event is recorded on this station
           if (abs(sta_lat(i)).gt.90) sta_lat(i) = sta_lat(i)-sta_lat(i)/abs(sta_lat(i))*180
           if (abs(sta_lon(i)).gt.180) sta_lon(i) = sta_lon(i)-sta_lon(i)/abs(sta_lon(i))*360

	      tlat=sta_lat(i)
	      tlon=sta_lon(i)
	      tdep=sta_elv(i)
	  
           aar=tlat ! these must be double-precisions
           bbr=tlon
           hr =tdep

           aar0=tlat ! these must be double-precisions
           bbr0=tlon
           hr0 =tdep

              staID=i
              evID=j
              call find_id(MAXEVE,MAXOBS,eve_sta,staID,evID,k)
              if(k.eq.0) goto 1199   ! do not calculate this station-event pair
!              write(*,*)src_cusp(j),src_lat(j),src_lon(j),src_dep(j)
!       convert source coordinates into km (HZ)
	     
              ! limit the range of event's location ! Bug fixed by Hao Guo
              if (abs(src_lat(j)).gt.90) src_lat(j) = src_lat(j)-src_lat(j)/abs(src_lat(j))*180
              if (abs(src_lon(j)).gt.180) src_lon(j) = src_lon(j)-src_lon(j)/abs(src_lon(j))*360 

	      tlat=src_lat(j)
	      tlon=src_lon(j)
	      tdep=src_dep(j)

              aas=tlat ! double-precisions
              bbs=tlon
              hs =tdep

              aas0=tlat ! double-precisions
              bbs0=tlon
              hs0 =tdep
	   
              aar=aar0
              bbr=bbr0
              hr =hr0
              !write(*,*)sta_lab(i),src_cusp(j),aar,bbr,hr,aas,bbs,hs

!---- Call the Spherical Pseudo-bending method
!---  ni is the number of segments in the ray path and w stores the
!---  coordinates of the segments. It starts from receive to the source
              isp=0
             rayError1 = 0
             call pbr( 0, aar, bbr, hr, aas, bbs, hs, rp_geo, nseg2, tmp_ttp(k,j), &
                 xn,yn,zn,nx,ny,nz,npcom,rak,dak,vpak,vsak,nlayer,xlayer,np15,msg, &
                 maxnx,maxny,maxnz,vel,bld,xl,yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
              rayErrorp(k,j) = rayError1  !!!lei
!                write(*,*) rayError1,rayError(k,j)
              if (rayErrorp(k,j) .eq. 1) write(*,*) 'Station:',i,sta_lat(i),sta_lon(i),'Event:',j,src_lat(j),src_lon(j)
              if (rayErrorp(k,j) .eq. 1) goto 1198 ! this ray is off-boundary. -- Hao Guo 

              do i1=1, nseg2+1
                 tlat  = 90-rp_geo(2,i1)
                 !dtemp = 1.0/.99664719;
                 dtemp = 1.0;
                 tlat  = atan(dtemp*tan(tlat*dpi))*r2d
                 tlon  = rp_geo(3,i1)
                 tdep  = r0-rp_geo(1,i1)
                 call  sph2car_ft(tlat,tlon,tdep,wlat,wlon,xc,yc,zc,theta)
                 rp(1,i1)=xc
                 rp(2,i1)=yc
                 rp(3,i1)=zc
                 !write(*,*)i1,tlat,tlon,tdep,xc,yc,zc
              enddo
              nrp = nseg2+1
              pgood(k,j) = 1

                 nray1=nrp      ! # of ray path segments for P wave
!     save ray path information for P waves and the find the ray path length
                 ray_length=0
                 seg_length=0
                
                 do k2=1,nrp
                    ray1(1,k2)=rp(1,k2)
                    ray1(2,k2)=rp(2,k2)
                    ray1(3,k2)=rp(3,k2)
                    k1=k2+1
                    if(k1.gt.nrp) k1=nrp
                    seg_length=(rp(1,k2)-rp(1,k1))*(rp(1,k2)-rp(1,k1))
                    seg_length=seg_length+(rp(2,k2)-rp(2,k1))*(rp(2,k2)-rp(2,k1))
                    seg_length=seg_length+(rp(3,k2)-rp(3,k1))*(rp(3,k2)-rp(3,k1))
                    seg_length=sqrt(seg_length)
                    ray_length=ray_length+seg_length
                 enddo


!       determine the hypocenter derivatives for P-wave
	      call ttmder( 0, rp, nrp, dthP, dtm, nx, ny, nz, xn, yn, zn, maxnx,maxny,maxnz,vel, bld,xl,yl,zl, &
              ixkms,iykms,izkms,ixloc,iyloc,izloc,maxpar,nparvi,wlat, wlon, theta, &
     	       nx1, ny1, nz1, nx2, nxy2, stepl, nfix, ndexfx, imerge, jequal,rayError1)
          rayErrorp(k,j) = rayError1   !!!lei
          if (rayErrorp(k,j) .eq. 1) goto 1198 

	      if(src_type(evID).ne.0) then !shot or Blast data
!		 write(*,*)'Shot data....'
		 tmp_xp(k,j)=0
		 tmp_yp(k,j)=0
		 tmp_zp(k,j)=0
	      else ! earthquake data
		 tmp_xp(k,j)=-dthP(2)
		 tmp_yp(k,j)=-dthP(3)
		 tmp_zp(k,j)=-dthP(4)
	      endif
              
              !write(*,*)tmp_ttp(k,j),tt,tmp_xp(k,j),tmp_yp(k,j),tmp_zp(k,j)
	      n=0
	      do m=1,npari  
		  if (abs(dtm(m)).ge.0.05) then
		    n=n+1
		    tmp_vp_index(k,j,n+1)=m
		    tmp_vp(k,j,n)=-dtm(m)	
		  endif
	      enddo
	      tmp_vp_index(k,j,1)=n ! No. of nonzero model derivatives
!          write(*,*) 'k,j,tmp_vp_index(k,j,1):',k,j,tmp_vp_index(k,j,1),n
1198          continue 	      
!       calculate the S-wave travel time by 3D ray tracing
	      
	      if (iuses.eq.2) then
		
                 aas=aas0 ! double-precisions
                 bbs=bbs0
                 hs =hs0

                 aar=aar0
                 bbr=bbr0
                 hr =hr0
!                write(*,*)aar,bbr,hr,aas,bbs,hs

!---- Call the Spherical Pseudo-bending method
!---  ni is the number of segments in the ray path and w stores the
!---  coordinates of the segments. It starts from receive to the source
                 isp=1
                 rayError1 = 0
                 call pbr( 1,aar, bbr, hr, aas, bbs, hs, rp_geo, nseg2, &
                         tmp_tts(k,j),xn,yn,zn,nx,ny,nz,npcom,rak,dak, &
     			          vpak,vsak,nlayer,xlayer,np15,msg,maxnx,maxny,maxnz,vel,bld,xl, &
                  		  yl,zl,ixkms,iykms,izkms,ixloc,iyloc,izloc,rayError1)
                 rayErrors(k,j) = rayError1  !!!lei
!                write(*,*) rayError1,rayError(k,j)
                 if (rayErrors(k,j) .eq. 1) write(*,*) 'Station:',i,sta_lat(i),sta_lon(i),'Event:',j,src_lat(j),src_lon(j)
                 if (rayErrors(k,j) .eq. 1) goto 1199 ! this ray is off-boundary. -- Hao Guo

                 do i1=1, nseg2+1
                   tlat  = 90-rp_geo(2,i1)
                   !dtemp = 1.0/.99664719;
                   dtemp = 1.0
                   tlat  = atan(dtemp*tan(tlat*dpi))*r2d
                   tlon  = rp_geo(3,i1)
                   tdep  = r0-rp_geo(1,i1)
                   call sph2car_ft(tlat,tlon,tdep,wlat,wlon,xc,yc,zc,theta)
                   rp(1,i1)=xc
                   rp(2,i1)=yc
                   rp(3,i1)=zc
                 enddo
                 nrp = nseg2+1
                 sgood(k,j) = 1

                  nray2=nrp   ! # of ray path segments for S wave
!     save ray path information for S waves
                  do k1=1,nrp
                     ray2(1,k1)=rp(1,k1)
                     ray2(2,k1)=rp(2,k1)
                     ray2(3,k1)=rp(3,k1)
                  enddo

		 call ttmder( 1, rp, nrp, dthS, dtm, nx, ny, nz, xn, yn, zn, maxnx,maxny,maxnz,vel, bld,xl,yl,zl,&
                    ixkms,iykms,izkms,ixloc,iyloc,izloc,maxpar,nparvi,wlat, wlon, theta, &
     				nx1, ny1, nz1, nx2, nxy2, stepl, nfix, ndexfx, imerge, jequal,rayError1)
                 rayErrors(k,j) = rayError1   !!!lei
                 if (rayErrors(k,j) .eq. 1) goto 1199    !--wangfan skip the wrong derivatives 

		 if(src_type(evID).ne.0) then !shot or blast data
		    tmp_xs(k,j)=0
		    tmp_ys(k,j)=0
		    tmp_zs(k,j)=0 
		 else ! Earthquake data
		    tmp_xs(k,j)=-dthS(2)
		    tmp_ys(k,j)=-dthS(3)
		    tmp_zs(k,j)=-dthS(4) 
		 endif
	
		 n=0
		 do m=1,npari
		    if (abs(dtm(m)).ge.0.05) then
		       n=n+1
		       tmp_vs_index(k,j,n+1)=m
		       tmp_vs(k,j,n)=-dtm(m)
		    endif
		 enddo
		 tmp_vs_index(k,j,1)=n ! save the number of non model derivatives
		 
	      endif
!---  find the haudorff distance between P and S wave ray paths
                 call hausdorff(ray1,nray1,ray2,nray2,haus_dist)
!---  determine if the distance between P and S wave rays are close enough;
!---  If hausdorff distance is less than the 5% of the total ray path.
                 if(haus_dist.lt.DISTratio*ray_length) then
                    tmp_sp_close(k,j)=1
                 else
                    !write(*,*)" S and P ray paths are not close"
                    tmp_sp_close(k,j)=0
                 endif                
     
1199          continue
	   enddo
	   enddo
	   !$omp enddo
	   !$omp end parallel
  end subroutine bw_forward

