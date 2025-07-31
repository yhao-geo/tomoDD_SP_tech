c Get input parameters

      subroutine getinpSPDR_SP(MAXEV,log,fn_inp,
     & fn_cc,fn_ct,fn_sta,fn_eve, fn_abs, fn_vel,
     & fn_cc_SP,fn_ct_SP,fn_abs_SP,
     & fn_loc,fn_reloc,fn_initial,fn_initial_SP,fn_res,fn_res_SP,fn_stares,
     & fn_srcpar,fn_vp, fn_vs,fn_vpvs, idata,iphase,
     & minobs_cc,minobs_ct,
     & amaxres_cross,amaxres_net,amaxdcc,amaxdct,
     & noisef_dt,maxdist,
     & awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp, awt_ctd,
     & istart,maxiter,isolv,niter,aiter,
     & iclust,ncusp,icusp,wlat, wlon, 
     & stepl, RayTracing,
     & rota,ajoint, CC_format, thres_vp,thres_vpvs,
     & wt_vp1, wt_vp2, wt_vp3,wt_vpvs1,wt_vpvs2,wt_vpvs3,
     & wt_vs1, wt_vs2, wt_vs3,wt_vpvsconst,DISTratio,
     & air_dep,!lat0,lat1,lon0,lon1,dep0,dep1,
     &  minvp,maxvp,minvs,maxvs,minvpvs,maxvpvs,
     &  maxdvp,maxdvs,maxdvpvs)

	include 'RaySPDR.inc'

c	Parameters:
	integer		maxev		! Array dimension
	integer		log		! Log-file identifier
	character	fn_abs*80       ! Absolute data file
	character	fn_abs_SP*80       ! Absolute data file
	character	fn_inp*80	! File of control info.
	character	fn_cc*80	! File of cross-corr. times
	character	fn_cc_SP*80	! File of cross-corr. times
	character	fn_ct*80	! File of catalog times
	character	fn_ct_SP*80	! File of catalog times
	character	fn_sta*80	! Station file
	character	fn_eve*80	! Event file
	character	fn_loc*80	! Output file of original locs.
	character	fn_reloc*80	! Output file of final locs.
	character       fn_initial*80   ! Output residual file of initial
	character       fn_initial_SP*80   ! Output residual file of initial
	character	fn_res*80	! Output residual file
	character	fn_res_SP*80	! Output residual file
	character	fn_stares*80	! Output station file
	character	fn_srcpar*80	! Output source-parameter file
	character	fn_vel*80	! Velocity output file
	character	fn_vp*80	! Final Vp model
	character	fn_vs*80	! Final Vs model
	character	fn_vpvs*80	! Final Vp/Vs model

	integer		idata		! 0: Synthetics
					! 1: Cross-correlation
					! 2: Catalog
					! 3: Both
	integer		iphase		! 1: P; 2: S; 3: Both
	integer		minobs_cc	! Min. obs./pair for ccor. data
	integer		minobs_ct	! Min. obs./pair for cat. data
	real		amaxres_cross(20)! [1..niter] Ccor. res. thresh.
	real		amaxres_net(20)	! [1..niter] Cat. res. thresh.
	real		amaxdcc(20)	! [1..niter] Ccor. link-dist. limit
	real		amaxdct(20)	! [1..niter] Cat. link-dist. limit
	real		noisef_dt	! Synthetic noise
	real		maxdist		! Max. cluster-station distance
	real		awt_ccp(20)	! [1..niter] Wts. for ccor. P
	real		awt_ccs(20)	! [1..niter] Wts. for ccor. S
	real		awt_ctp(20)	! [1..niter] Wts. for cat. P
	real		awt_cts(20)	! [1..niter] Wts. for cat. S
	real            awt_ctd(20)     ! [1..niter] wts. for abs. vs. diff. cat.
	real		adamp(20)	! [1..niter] Damping (lsqr only)
	integer         ajoint(20)      ! Joint inversion or not
	integer		istart		! 1: From single source
					! 2: From network sources
	integer		maxiter
	integer		isolv		! 1: SVD; 2: LSQR
	integer		niter		! No. of iteration sets
	integer		aiter(0:20)	! [1..niter] Iterations/set
	integer		iclust		! Cluster to relocate (0: all).
	integer		ncusp		! No. of event keys in icusp[]
	integer		icusp(maxev)	! [1..ncusp] Events to relocate
	real            stepl           ! step size for calculating model derivatives
	real            wlat           
	real		wlon          
        integer         RayTracing      ! 1 for Hole's, 2 for Lomax's
        real            rota            ! rotation angle
        real            thres_vp(20)
        real            thres_vpvs(20)
        integer         CC_format
        real            wt_vp1,wt_vp2,wt_vp3
        real            wt_vs1,wt_vs2,wt_vs3
        real            wt_vpvs1,wt_vpvs2,wt_vpvs3
        real            air_dep
!        real            lat0,lat1,lon0,lon1,dep0,dep1

	real            wt_vpvsconst
	real            DISTratio
        real            minvp,maxvp,minvs,maxvs,minvpvs,maxvpvs
        real            maxdvp,maxdvs,maxdvpvs

c	Local variables:
	integer		fu_inp
	integer		i
	integer		ii
	integer		l
	character	line*150
	integer		trimlen

	character rcsid*150
	data rcsid /"$Header:getinpFDD.f $"/
	save rcsid

c--- newest format: 083000 with iteration step dependent weighting
c-- open input file:
      call freeunit (fu_inp)
      open (fu_inp,status='unknown',file=fn_inp,err=998)
      ncusp= 0
      niter= 0  ! number of iteration blocks
      l = 1
      ii= 1

c-- Loop to read each parameter lines, skipping comments
 210  read (fu_inp,'(a)',end=220) line
      if (line(1:1).eq.'*' .or. line(2:2).eq.'*') goto 210
      if (l.eq.1) read (line,'(a)',err=999) fn_cc
      if (l.eq.2) read (line,'(a)',err=999) fn_cc_SP
      if (l.eq.3) read (line,'(a)',err=999) fn_ct
      if (l.eq.4) read (line,'(a)',err=999) fn_ct_SP
      if (l.eq.5) read (line,'(a)',err=999) fn_eve
      if (l.eq.6) read (line,'(a)',err=999) fn_sta
      if (l.eq.7) read (line,'(a)',err=999) fn_loc
      if (l.eq.8) read (line,'(a)',err=999) fn_reloc
      if (l.eq.9) read (line,'(a)',err=999) fn_stares
      if (l.eq.10) read (line,'(a)',err=999) fn_initial
      if (l.eq.11) read (line,'(a)',err=999) fn_initial_SP
      if (l.eq.12) read (line,'(a)',err=999) fn_res
      if (l.eq.13) read (line,'(a)',err=999) fn_res_SP
      if (l.eq.14) read (line,'(a)',err=999) fn_srcpar
      if (l.eq.15) read (line,'(a)',err=999) fn_vel
      if (l.eq.16) read (line,'(a)',err=999) fn_vp
      if (l.eq.17) read (line,'(a)',err=999) fn_vs
      if (l.eq.18) read (line,'(a)',err=999) fn_vpvs
      if (l.eq.19) read (line,'(a)',err=999) fn_abs
      if (l.eq.20) read (line,'(a)',err=999) fn_abs_SP
      if (l.eq.21) read (line,*,err=999) idata, iphase,maxdist
      if (l.eq.22) read (line,*,err=999) minobs_cc,minobs_ct, air_dep
      if (l.eq.23) then 
          read (line,*,err=999) istart, isolv, niter, RayTracing,
     &          wt_vpvsconst, DISTratio
      end if
      if (l.eq.24) then
          read (line,*,err=999) iuses, iuseq, invdel, stepl!,lat0,lat1,
!     &          lon0,lon1,dep0,dep1
      end if
      if (l.eq.25) then 
	 read (line,*,err=999) wlat, wlon, rota, CC_format
      end if
      if (l.eq.26) then
         read (line,*,err=999) minvp,maxvp,minvs,maxvs,minvpvs,maxvpvs,
     &                         maxdvp,maxdvs,maxdvpvs
      end if
      if (l.eq.27) then
         read (line,*,err=999) wt_vp1, wt_vp2, wt_vp3,wt_vs1,wt_vs2,
     &        wt_vs3, wt_vpvs1, wt_vpvs2, wt_vpvs3
      endif
c--Read iteration instructions
      if (l.ge.28 .and. l.le.27+niter) then
         i=l-27
         read (line,*,err=999) aiter(i),
     &  awt_ccp(i), awt_ccs(i), amaxres_cross(i), amaxdcc(i),
     &  awt_ctp(i), awt_cts(i), amaxres_net(i), amaxdct(i), awt_ctd(i),
     &  adamp(i),ajoint(i), thres_vp(i), thres_vpvs(i)
      endif
c--Read specific clusters/events to relocate
      if (l.eq.28+niter) read (line,*,err=999) iclust
      if (l.ge.29+niter) then
         read (line,*,err=999,end=230) (icusp(i),i=ii,ii+7)
230      ii= i
      endif
      l= l+1
      goto 210
220   close (fu_inp)
      ncusp= ii-1

c- rearrange aiter:
      do i=2,niter
        aiter(i)= aiter(i-1)+aiter(i)
      enddo

c- check files
      call exist (fn_eve)
      call exist (fn_sta)
      if ((idata.eq.1 .or.idata.eq.3).and.trimlen(fn_cc).gt.1)
     & call exist(fn_cc)
      if ((idata.eq.2 .or.idata.eq.3).and.trimlen(fn_ct).gt.1)
     & call exist (fn_ct)

      maxiter= aiter(niter)
c synthetic noise:
      noisef_dt= 0.002

c write log output: of newest format
600   if (trimlen(fn_loc).lt.2) fn_loc= 'tomoSPD_SP.loc'
      if (trimlen(fn_reloc).lt.2) fn_reloc= 'tomoSPD_SP.reloc'
      write (6,'("INPUT FILES:",/,
     &"cross dtime data: ",a,/,"catalog dtime data: ",a,/,
     &"cross dtime S-P data: ",a,/,"catalog dtime S-P data: ",a,/,
     &"events: ",a,/,"stations: ",a,/,"OUTPUT FILES:",/,
     &"initial locations: ",a,/,"relocated events: ",a,/,
     &"event pair residuals: ",a,/,"station residuals: ",a,/,
     &"absolute data: ",a,/,"catalog absolute S-P data: ",a,/,
     &"velocity output: ",a,/,
     &"source parameters: ",a,/,"initial residuals:",a)')
     &fn_cc(1:trimlen(fn_cc)),
     &fn_ct(1:trimlen(fn_ct)),
     &fn_cc_SP(1:trimlen(fn_cc_SP)),
     &fn_ct_SP(1:trimlen(fn_ct_SP)),
     &fn_eve(1:trimlen(fn_eve)),
     &fn_sta(1:trimlen(fn_sta)),fn_loc(1:trimlen(fn_loc)),
     &fn_reloc(1:trimlen(fn_reloc)),fn_res(1:trimlen(fn_res)),
     &fn_stares(1:trimlen(fn_stares)),
     &fn_abs(1:trimlen(fn_abs)),fn_abs_SP(1:trimlen(fn_abs_SP)), 
     &fn_vel(1:trimlen(fn_vel)),
     &fn_srcpar(1:trimlen(fn_srcpar)),fn_initial(1:trimlen(fn_initial))

      write (log,'("Input parameters: (from ",a,")",/,
     &"  cross dtime file: ",a,/,"  catalog dtime file: ",a,/,
     &"  catalog absolute file: ",a,/,"  cross S-P dtime file: ",a,/,
     &"  catalog S-P dtime file: ",a,/,"  catalog absolute S-P file: ",a,/,
     &"  station file: ",a,/,"  event file: ",a,/,
     &"  initial locations: ",a,/,
     &"  absolute data file: ",a,/,
     &"  velocity output file: ",a,/,
     &"  relocated events: ",a)')
     &fn_inp(1:trimlen(fn_inp)),
     &fn_cc(1:trimlen(fn_cc)),
     &fn_ct(1:trimlen(fn_ct)),
     &fn_abs(1:trimlen(fn_abs)),
     &fn_cc_SP(1:trimlen(fn_cc_SP)),
     &fn_ct_SP(1:trimlen(fn_ct_SP)),
     &fn_abs_SP(1:trimlen(fn_abs_SP)),
     &fn_sta(1:trimlen(fn_sta)),
     &fn_eve(1:trimlen(fn_eve)),fn_loc(1:trimlen(fn_loc)),
     &fn_abs(1:trimlen(fn_abs)), fn_vel(1:trimlen(fn_vel)),
     &fn_reloc(1:trimlen(fn_reloc))

      write(log,*)'Air depth is:', air_dep
      write (log,'(
     &"  event pair file: ",a,/,"  station residual file: ",a,/,
     &"  source parameter file: ",a,/,
     &"  IDATA= ",i2,2X,"IPHASE= ",i2,2x,"MAXDIST= ",f5.0,/,
     &"  MINOBS_CC= ",i3,2x,"MINOBS_CT= ",i3,/,"  ISTART= ",i1,2x,
     &"ISOLV= ",i1,2x)')
     &fn_res(1:trimlen(fn_res)),
     &fn_stares(1:trimlen(fn_stares)),fn_srcpar(1:trimlen(fn_srcpar)),
     &idata,iphase,maxdist,minobs_cc,minobs_ct,istart,isolv

      aiter(0)=0
      write (log, '("  ITER ",i2,"-",i2,
     &": DAMP= "f7.1,/,"    WT_CCP= ",f7.4,2X,"WT_CCS= ",f7.4,2x,
     &"MAXR_CC= ",f7.4,2X,"MAXD_CC= ",f7.2,2X,/,
     &"    WT_CTP= ",f7.4,2x,"WT_CTS= ",f7.4,2x,"MAXR_CT= ",f7.4,2x,
     &"MAXD_CT= ",f7.2, 2x, "WT_CTD= ",f7.2, 2x, "JOINT= ", i3,2x,
     &"THRE_VP= ",f7.4,2x, "THRE_VPVS= ",f7.4)')
     &(aiter(i-1)+1,aiter(i),adamp(i),awt_ccp(i),awt_ccs(i),
     & amaxres_cross(i),
     & amaxdcc(i),awt_ctp(i),awt_cts(i), amaxres_net(i), amaxdct(i),
     & awt_ctd(i),ajoint(i),thres_vp(i),thres_vpvs(i),
     & i=1,niter)
 
!      write(*,*)'lat0,lat1,lon0,lon1,dep0,dep1',lat0,lat1,lon0,lon1,dep0,dep1

c--Repeat number of clusters, events to relocate
      if (iclust.eq.0) then
        write (*,*) 'Relocate all clusters'
        write (log,*) 'Relocate all clusters'
      else
        write (*,*) 'Relocate cluster number ',iclust
        write (log,*) 'Relocate cluster number ',iclust
      end if

      if (ncusp.eq.0) then
        write (*,*) 'Relocate all events'
        write (log,*) 'Relocate all events'
      else
        write (*,*) 'Relocate ',ncusp,' events'
        write (log,*) 'Relocate ',ncusp,' events'
      end if
      return

c--Input error handling
998   write(*,*)'>>> ERROR OPENING CONTROL PARAMETER FILE'
      goto 1000

999   write (*,*)'>>> ERROR READING CONTROL PARAMETERS IN LINE ',l
      write (*,*) line
1000  stop 'Program run aborted.'
      end  ! of subroutine getinp
