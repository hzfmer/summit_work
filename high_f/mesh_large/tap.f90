       program tap
       integer(kind=8), parameter    :: nz_ssh = 500, nz_tap = 375
       integer    :: k

       do k=1,1000
       if (k.le.nz_tap) then
	  tapval = 1.0
       elseif (k.gt.nz_ssh) then
	  tapval = 0.
       else
	  tapval = 1. - real(k - nz_tap) / real(nz_ssh - nz_tap)
       endif
       write(*,*) k, tapval

       enddo
       end
