	subroutine dtres_FDD_SP(log, ndt, stdim, nsrc, nsta, 
     &	dt_dt, 
     &	dt_ista, dt_ic1, dt_ic2,
     &	src_cusp, src_t, tmp_ttp, tmp_tts,
     &	dt_res, eve_sta)

	use tomoFDD
	implicit none

c	Parameters:
	integer		log		! Log-file identifier
	integer		ndt		! No. of data
	integer		stdim		! Column dimenson of arrays tmp_tt[ps]
	integer		nsrc		! No. of sources
	integer         nsta
	real		dt_dt(MAXDATA)	! [1..ndt] Observed time differences
	integer		dt_ista(MAXDATA)! [1..ndt] Station indices
	integer		dt_ic1(MAXDATA)	! [1..ndt] Event indices
	integer		dt_ic2(MAXDATA)	! [1..ndt] Event indices
	integer		src_cusp(MAXEVE)! [1..nsrc] Event keys
	real		src_t(MAXEVE)	! [1..nsrc] Event times
	real		tmp_ttp(MAXOBS,MAXEVE)! [1.., 1..MAXOBS]
	real		tmp_tts(MAXOBS,MAXEVE)! [1.., 1..MAXOBS]
	real		dt_cal(MAXDATA)	! [1..ndt] Theoretical time differences
	real		dt_res(MAXDATA)	! [1..ndt] Time-difference residuals

	integer         eve_sta(MAXEVE,MAXOBS+1)
	integer         evID1, evID2, j1, j2

c	Local variables:
	integer		i,j
	integer		k1, k2
	real		ttP1, ttP2
	real            ttS1, ttS2

      write(log,'("~ getting residual vector...")')

      if (nsrc.eq.1) then
c        Single source
         do i=1,ndt
            dt_res(i) = dt_dt(i)
         enddo
      else
c        Mulitple sources
	   ttP1 = 0.0
	   ttP2 = 0.0
	   ttS1 = 0.0
	   ttS2 = 0.0
         do i=1,ndt

cz--- find the sequence number for each event
            evID1 = dt_ic1(i)
            evID2 = dt_ic2(i)

cz--- call the subroutine find_id to look for the sequence number
cz--- corresponding to each station-event pair
            call find_id(MAXEVE,MAXOBS,eve_sta, dt_ista(i), evID1, k1)
            call find_id(MAXEVE,MAXOBS,eve_sta, dt_ista(i), evID2, k2)

c              P phase
	      ttP1 = tmp_ttp(k1,dt_ic1(i)) -
     &               src_t(dt_ic1(i))/1000
	      ttP2 = tmp_ttp(k2,dt_ic2(i)) -
     &               src_t(dt_ic2(i))/1000
              if (ttP1.eq.0 .or. ttP2.eq.0) then
		 !write(*,*)' P Phase'
		 !write(*,'("FATAL ERROR (theor tt).")') 
c                stop
              endif
c       
c       S phase
	      ttS1 = tmp_tts(k1,dt_ic1(i)) -
     &               src_t(dt_ic1(i))/1000
	      ttS2 = tmp_tts(k2,dt_ic2(i)) -
     &               src_t(dt_ic2(i))/1000
	      if (ttS1.eq.0 .or. ttS2.eq.0) then
		 !write(*,*)' S Phase'
		 !write(*,'("FATAL ERROR (theor tt).")')
c                stop 	       
	      endif

	      if (dt_ic1(i) .ne. dt_ic2(i) ) then ! difference time
		 dt_cal(i) = ttS1 - ttS2-(ttP1-ttP2)
		 dt_res(i) = dt_dt(i) - dt_cal(i)
	      else
		 dt_cal(i) =ttS1-ttP1
		 dt_res(i) = dt_dt(i) - dt_cal(i) ! absolute catalog time
	      endif 
	   enddo
	endif

	end			!of subroutine dtres
