c-- This subroutine is used to calculate the Hausdorff distance between two polygonal
c-- curves. Hausdorff distance is defined as the maximum distance of the minimum distance
c-- between two sets. Since the curves for P and S waves are similar to each other, sharing
c-- the same origin and end points, I will not search all the points, rather than choosing only
c-- a limited number of points.

      subroutine hausdorff(ray1, nray1, ray2, nray2, dist)
      
      real                  ray1(3,130)
      real                  ray2(3,130)
      integer               nray1, nray2
      real                  dist
      integer               i, j, k1, k2
      real                  dij,shortest
      
      dist=0
      do i=1,nray1
         k1=max(1,i-10) ! only search the neighboring 20 points
         k2=min(nray2,i+10)
         x1=ray1(1,i)
         x2=ray1(2,i)
         x3=ray1(3,i)
         
         shortest=999999.0
         do j=k1,k2
            y1=ray2(1,j)
            y2=ray2(2,j)
            y3=ray2(3,j)
            dij=sqrt( (x1-y1)*(x1-y1)+(x2-y2)*(x2-y2)+(x3-y3)*(x3-y3) )
            if(dij.lt.shortest) shortest=dij
         enddo

         if(dist.lt.shortest) dist=shortest

      enddo
      return
      end
            
