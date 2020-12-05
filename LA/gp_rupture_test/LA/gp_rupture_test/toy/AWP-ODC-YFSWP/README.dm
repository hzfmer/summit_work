Notes for the implementation of DM mesh in AWP
Daniel Roten, Feb 2017 (droten@mail.sdsu.edu)

-For stress and velocity updates, as well as plasticity, same kernels are used for different mesh sizes
-Arrays with stresses and velocities are now allocated as 2-D arrays, where the first dimension controls the grid number.
 E.g.,
 float *u1[0] contains x-velocity for finest grid,
 float *u1[1] contains x-velocity for coarse grid.
 The number of grid sizes can be larger than 2, but also just 1 (uniform mesh).

-GPU device variables, d_DT, d_DH, ..., and d_yline_[1,2] and d_slice_[1,2] are now also arrays.
 d_DH[0] contains grid spacing for fine grid, d_dH[1] for coarse(r) grid, d_dH[2] for coarsest grid, for example
 same with other arrays

-all kernels now accept an additional argument, d_i, which contains the grid index.  This controls which of the
 device variable is used for computing FD updates.  The only exceptions are the free surface B.C., fstr and fvel, as they 
 are always called on grid index 0.
 
-arrays with material properties and velocities/stresses, are still passed to kernels as 1D arrays.
 So we call dvelcx_H using u1[0], v1[0],  for the fine grid, for example.
 Each kernel is invoked once per grid size (and per subdomain partition).


the velocities and yield factors from different grids (fine, coarse, extra coarse, etc ...) are exchanged
using one calls to SendRecvMsgX / PostRecvMsgY for each grid size.
To distinguish between messages pertaining to different grid sizes, the message tags are computed such that
the tag for each grid size is unique.  E.g., for PostSendMsg_Y:

tag = MPIRANKY+rank*gridnum

Therefore, the grid number is passed as extra (last) argument to the following routines:

PostSendMsg_Y, PostRecvMsg_Y, PostSendMsg_X, PostRecvMsg_X

-Updated command-line arguments:
 --NGRIDS or -G: specifies the number of grid sizes (1 = equidistant)
 --NZ or -Z specifies vertical extent of different grids, from top to bottom, separated by comma.  E.g.,
    -NZ 32,36,56 using three different grid sizes.  Note that sizes include the overlap zone.     

Todo: revisit interpolation of normal stresses, u1, v1 and xy from coarse grid at position kh=align+1
 (interpolation in vertical direction may be needed too)
