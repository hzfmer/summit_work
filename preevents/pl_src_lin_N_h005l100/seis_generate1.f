CCC **************************************************************************** 
CCC  Update -- MARCH 23 2018
CCC  Compiler: pgi (module swap PrgEnv-gnu PrgEnv-pgi)
CC             `pgf90 seis_generate8.f -o seis_generate8` 
******************************************************************************* 
CCC ************************************************************************
CCC  MOST RECENT UPDATE -- MARCH 25 2014
CCC  Added the capability of performing averaging across the
CCC  free-surface directly.  Use the old version of the code
CCC  if you want to extract the normal seismograms, depending on
CCC  performance this might be extended to MPI for larger models. --
CCC  William Savran
CCC
C******************************************************************************
CCC  Last update - Aug 30, 2012
CCC    A new parameter (SAVE_STEP) is read to save the seismogram data gathered
CCC    when that number of files are read. - Efecan.
CCC
CCC ****************************************************************************
CC   Kwangyoon updated:
CC   ./seis_generate numsta stationlist path x_prefix y_predix z_prefix tmax dt num_x num_y ntiskp write_step save_step
CC   (example: ./seis_generate 14 w2wstations.txt ../output_sfc/ SX96PS SY96PS SZ96PS 500.0 0.0058 4000, 2000, 20, 100, 2)
CC   num_x should be the number of x elements in output file (nx/NSKPX)
CC   num_y should be the number of y elements in output file (ny/NSKPZ)
CCC ****************************************************************************

C.0 - Variable Declaration & Allocation

      integer :: n_args
      character (len=120) :: output_path, stationlist
      character (len=30) :: numsta_chars, xpref_chars, ypref_chars
      character (len=30) :: zpref_chars, tmax_chars
      character (len=30) :: dt_chars, nx_chars, ny_chars, ntiskp_chars
      character (len=30) :: writestep_chars
      character (len=30) :: savestep_chars

      integer :: i,j,k,m,p,err,idx,idx2,idx3,lastSavedIdx
      integer, parameter :: filen=160
c     ADDED BY WILLIAM SAVRAN -- CAN BE INCLUDED IN INPUT FILE, HOWEVER
C     THIS IS FOR AVERAGING ACROSS THE FREE-SURFACE
      integer :: nx
      integer :: ny
      integer :: ntiskp
      integer :: write_step
      integer, parameter :: floatsize=4    ! change from 4 to 1 for portability
      integer :: numsta
      integer :: save_step
      integer :: current_save_block, lfinish

      integer, dimension (:), allocatable :: xs
      integer, dimension (:), allocatable :: ys

      integer :: fileunit
      real :: totalsteps_f
      integer :: totalsteps
      real :: tstart, tfinish

      character (len=filen) :: fx,fy,fz,fseis
      real, dimension (:), allocatable :: xdata
      real, dimension (:), allocatable :: ydata
      real, dimension (:), allocatable :: zdata

      real, dimension (:), allocatable :: xseis
      real, dimension (:), allocatable :: yseis
      real, dimension (:), allocatable :: zseis

      real :: xtemp1, xtemp2, xtemp3, xtemp4
      real :: ytemp1, ytemp2, ytemp3, ytemp4
      real :: ztemp1, ztemp2, ztemp3, ztemp4
      real :: avgx, avgy, avgz
C.0   Get argument first
      call cpu_time(tstart)
      n_args = iargc()
      if (n_args .ne. 13) then
        print *, 'The program needs 11 arguments, in the order of: '
        print *, '\t numsta stationlist output_path, x_prefix, y_prefix,
     +        z_prefix, tmax, dt, nx, ny, ntiskp, write_step, save_step'
      end if

      call getarg(1, numsta_chars)
      call getarg(2, stationlist)
      call getarg(3, output_path)
      call getarg(4, xpref_chars)
      call getarg(5, ypref_chars)
      call getarg(6, zpref_chars)
      call getarg(7, tmax_chars)
      call getarg(8, dt_chars)
      call getarg(9, nx_chars)
      call getarg(10, ny_chars)
      call getarg(11, ntiskp_chars)
      call getarg(12, writestep_chars)
      call getarg(13, savestep_chars)

      read (numsta_chars,'(I10)')  numsta
      read (tmax_chars,'(F10.3)')  TMAX
      read (dt_chars,'(F10.3)')  DT
      read (nx_chars,'(I10)')  nx
      read (ny_chars,'(I10)')  ny
      read (ntiskp_chars,'(I10)')  ntiskp
      read (writestep_chars,'(I10)')  write_step
      read (savestep_chars,'(I10)')  save_step

      print *, numsta, trim(stationlist)
      print *, trim(output_path), ' ; ', trim(xpref_chars), ' ; ',
     +         trim(ypref_chars), ' ; ', trim(zpref_chars)
      print *, TMAX, DT, nx, ny, ntiskp, write_step
      print *, 'save seismogram after reading ', save_step, ' files.'

      fileunit=ntiskp*write_step
      totalsteps_f=TMAX/DT
      totalsteps=(ifix(totalsteps_f)/fileunit)*fileunit
      if (totalsteps.ne.ifix(totalsteps_f)) then
        totalsteps=totalsteps+fileunit
      endif
      write(*,*) 'fileunit ', fileunit
      write(*,*) 'totalsteps ', int(totalsteps_f),
     +            int(totalsteps_f/ntiskp)

      write(*,*) 'totalsteps (rounded)', totalsteps

      allocate(xdata(nx*ny))
      allocate(ydata(nx*ny))
      allocate(zdata(nx*ny))
      allocate(xs(numsta))
      allocate(ys(numsta))

      allocate(xseis(numsta*write_step*save_step))
      allocate(yseis(numsta*write_step*save_step))
      allocate(zseis(numsta*write_step*save_step))

      open(20,file=stationlist)
      do k=1,numsta
        read(20,*) xs(k),ys(k)
      enddo

C.1 - Initialize variables

      idx=0
      lastSavedIdx=0
      current_save_block=1
      lfinish=0

C.2 - Calculate vmag for each 10th timestep & Calcualte final vcumu for all timesteps

c  use 45600 here for consistent handling, create soft-link
      do k=fileunit,totalsteps,fileunit
C        lfinish=0
        if (k .eq. totalsteps) then
            lfinish=1
        endif
        write(*,*) '* ', k

        write(fx, '(a,a,i7.7)'), trim(output_path), trim(xpref_chars),k
        write(fy, '(a,a,i7.7)'), trim(output_path), trim(ypref_chars),k
        write(fz, '(a,a,i7.7)'), trim(output_path), trim(zpref_chars),k

        write(*,*) 'fx', fx
        write(*,*) 'fy', fy
        write(*,*) 'fz', fz

        open(9,file=fx,form='unformatted',status='old',
     +  access='direct',recl=floatsize*nx*ny,iostat=err)
        open(10,file=fy,form='unformatted',status='old',
     +  access='direct',recl=floatsize*nx*ny,iostat=err)
        open(11,file=fz,form='unformatted',status='old',
     +  access='direct',recl=floatsize*nx*ny,iostat=err)

        do p=1,write_step
          write(*,*) ' -- ', p
          idx=idx+1
          if (idx .gt. int(totalsteps/ntiskp)) then
            print *, 'All data is read.'
            close(9)
            close(10)
            close(11)
            lfinish=1
            goto 333
          end if
c         read z=0 plane and z=1 plane as seperate slices
          read(9,rec=p) xdata
          read(10,rec=p) ydata
c          read(11,rec=nz*p-1) zdata_t
          read(11,rec=p) zdata
          do m=1,numsta
            idx2 = (m-1)*(save_step*write_step)
            idx3 = (ys(m)-1)*nx+xs(m)
c           add other indexes to correspond to the averaging parameters
c           here, we are going to extract the four required parameters
c           and then average them.
            xtemp1 = xdata(idx3)
c           don't extract values on the border's edge, otherwise it will crash
            ytemp1 = ydata(idx3)
            ytemp2 = ydata(idx3-nx)
            ytemp3 = ydata(idx3-1)
            ytemp4 = ydata(idx3-nx-1)
c           we only use the free-surface node for the z-velocity
            ztemp1 = zdata(idx3)
            ztemp2 = zdata(idx3-1)
c           average across the free-surface
            avgx = xtemp1
            avgy = 0.25*(ytemp1+ytemp2+ytemp3+ytemp4)
            avgz = 0.5*(ztemp1+ztemp2)

            if (p==write_step .and. m==1) print *,xdata(idx3)
c            xseis(idx2+idx-lastSavedIdx) = xdata(idx3)
c            yseis(idx2+idx-lastSavedIdx) = ydata(idx3)
c            zseis(idx2+idx-lastSavedIdx) = zdata(idx3)
              xseis(idx2+idx-lastSavedIdx) = avgx
              yseis(idx2+idx-lastSavedIdx) = avgy
              zseis(idx2+idx-lastSavedIdx) = avgz
          enddo

        end do   ! end of "do p=1,write_step" loop
        write(*,*) 'Now idx = ', idx
        close(9)
        close(10)
        close(11)
        if (k .ge. current_save_block*save_step*fileunit) then
c         write seismograms and deallocate arrays
          goto 222

444       continue

          allocate(xdata(nx*ny))
          allocate(ydata(nx*ny))
          allocate(zdata(nx*ny))

          allocate(xseis(numsta*write_step*save_step))
          allocate(yseis(numsta*write_step*save_step))
          allocate(zseis(numsta*write_step*save_step))

          current_save_block = current_save_block + 1
          lastSavedIdx = idx
        end if
      end do    ! end of "do k" loop


C.3 -  Write seismograms

222   continue

      print *, "WRITING SEISMOGRAM"

      print *, "seis_x"
      write(fseis, '(A6,I7.7)')
     +       'seis_x',current_save_block*save_step*fileunit
      open(19,file=fseis,form='unformatted',
     + access='direct',recl=floatsize*numsta*write_step*save_step)
      write(19,rec=1) xseis
      close(19)

      print *, "seis_y"
      write(fseis, '(A6,I7.7)')
     +        'seis_y',current_save_block*save_step*fileunit
      open(19,file=fseis,form='unformatted',
     +  access='direct',recl=floatsize*numsta*write_step*save_step)
      write(19,rec=1) yseis
      close(19)

      print *, "seis_z"
      write(fseis, '(A6,I7.7)')
     +         'seis_z',current_save_block*save_step*fileunit
      open(19,file=fseis,form='unformatted',
     + access='direct',recl=floatsize*numsta*write_step*save_step)
      write(19,rec=1) zseis
      close(19)

C      lfinish=1

C.4 -  Deallocate arrays

333   continue

      deallocate(xdata)
      deallocate(ydata)
      deallocate(zdata)

      deallocate(xseis)
      deallocate(yseis)
      deallocate(zseis)

C     If not finished go back to loop
      if (lfinish.eq.0) goto 444

      deallocate(xs)
      deallocate(ys)

      call cpu_time(tfinish)
      print *, 'Total elapsed time=', tfinish-tstart
666   stop
      end
