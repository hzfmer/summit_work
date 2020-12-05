! author: william savran
! date: 03.02.2015
! this parallelization asks for one processor per depth layer. this is a 1d decomposition where 
! NX = nx, NY = ny, NZ = 1
! for NVE = 3

program add
  implicit none
  include 'mpif.h'
  
  ! Variable declarations
  ! ***********************************
  real, parameter               :: std = 0.05  
  integer(kind=8), parameter    :: nx = 1000, ny = 1000
  integer(kind=8), parameter    :: nrec_len_hom = nx*ny*5, nrec_len_fr = nx*ny, nrec_len_het = nx*ny*5
  integer(MPI_OFFSET_KIND)      :: offset_hom, offset_fr, offset_het
  integer                       :: rank, nthreads
  integer                       :: err
  integer                       :: MCW, fh1, fh2, fh3
  real, dimension(:), allocatable :: buf_fr
  real, dimension(:), allocatable :: buf_hom, buf_het
  real                          :: temp_vp, temp_vs, temp_rho, a, ss, sp
  real                          :: temp_qp, temp_qs, temp_xi, temp_yi, temp_zi
  integer(kind=8)               :: i, j, k, ind
  character(len=50)             :: filename_fr, filename_hom, filename_het
! for 4km depth
!  real, parameter               :: vpmax = 5701, vsmax = 3300., vsmin=140
  real, parameter               :: vpmax = 6000, vsmax = 3400., vsmin=200
  
  filename_fr = 'ssh_4m_h005l100.out'
  filename_hom = 'mesh_smth_gvda_4m.bin'
  filename_het = 'mesh_smth_gvda_4m_het.bin'
    
  call MPI_INIT(err)
  call MPI_COMM_DUP(MPI_COMM_WORLD,MCW,err)
  call MPI_COMM_RANK(MCW,rank,err)
  call MPI_COMM_SIZE(MCW,nthreads,err)

  allocate (buf_fr(nx*ny), buf_hom(nx*ny*5), buf_het(nx*ny*5))

  ! compute offset based on rank
  offset_fr = rank*nrec_len_fr*4
  offset_hom  = rank*nrec_len_hom*4
  offset_het = rank*nrec_len_het*4
  k = rank+1
  
  ! open fractal model
  call MPI_FILE_OPEN(MCW, filename_fr, MPI_MODE_RDONLY, MPI_INFO_NULL, fh1, err)
  call error_check(err, "opening fractal model")
  call MPI_FILE_READ_AT_ALL(fh1, offset_fr, buf_fr, nrec_len_fr, MPI_REAL, MPI_STATUS_IGNORE, err)
  call error_check(err, "reading fractal model")
  call MPI_FILE_CLOSE(fh1, err)
  
  ! open homogeneous mesh
  call MPI_FILE_OPEN(MCW, filename_hom, MPI_MODE_RDONLY, MPI_INFO_NULL, fh2, err)
  call error_check(err, "opening homogeneous model")
  write(*,*) 'offset_hom=', offset_hom
  call MPI_FILE_READ_AT_ALL(fh2, offset_hom, buf_hom, nrec_len_hom, MPI_REAL, MPI_STATUS_IGNORE, err)
  call error_check(err, "reading homogeneous model")
  call MPI_FILE_CLOSE(fh2, err)
  
  ! loop through slice
  do j=1,ny
    do i=1,nx
      ! ordered vp(1), vs(1), rho(1)
      ind = (j-1)*nx+i
      temp_vp = buf_hom(5*(ind-1)+1)
      temp_vs = buf_hom(5*(ind-1)+2)
      temp_rho = buf_hom(5*(ind-1)+3)
      !No SSHs are added to the water layer
      if (temp_vs > 0.) then
         ! convert to slowness
         sp = 1./temp_vp
         ss = 1./temp_vs
         ! add heterogeneities
         sp = sp*(1-std*buf_fr(ind))
         ss = ss*(1-std*buf_fr(ind))
         ! convert back
         buf_het(5*(ind-1)+1)=1./sp
         if (buf_het(5*(ind-1)+1) .gt. vpmax) then
            buf_het(5*(ind-1)+1) = vpmax
         endif
         buf_het(5*(ind-1)+2)=1./ss
         if (buf_het(5*(ind-1)+2) .gt. vsmax) then
            buf_het(5*(ind-1)+2) = vsmax
         endif
         if (buf_het(5*(ind-1)+2) .lt. vsmin) then
            buf_het(5*(ind-1)+2) = vsmin
            if (buf_het(5*(ind-1)+1) .lt. vsmin*sqrt(2.)) then
               buf_het(5*(ind-1)+1) = vsmin*sqrt(2.)
            endif
         endif
         buf_het(5*(ind-1)+3) = temp_rho*(1+std*buf_fr(ind))
      else
         buf_het(5*(ind-1)+1)=temp_vp
         buf_het(5*(ind-1)+2)=temp_vs
         buf_het(5*(ind-1)+3)=temp_rho
      endif 
      buf_het(5*(ind-1)+4) = buf_hom(5*(ind-1)+4)
      buf_het(5*(ind-1)+5) = buf_hom(5*(ind-1)+5)
  ! Get number of threads
    enddo
  enddo
  
  ! open complete file
  call MPI_FILE_OPEN(MCW, filename_het, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh3, err)
  call error_check(err, "opening heterogeneous model")
  call MPI_FILE_WRITE_AT_ALL(fh3, offset_het, buf_het, nrec_len_het, MPI_REAL, MPI_STATUS_IGNORE, err)
  call error_check(err, "writing heterogeneous model")
  call MPI_FILE_CLOSE(fh3, err)
  
  ! open combined model
  call MPI_FINALIZE(err)
end program

subroutine error_check(err, message)
  integer :: err, ierr, errlen
  character (len=*) :: message
  character (len=500) :: errmsg

  if (err .ne. 0) then
     call MPI_Error_string(err, errmsg, errlen, ierr)  
     write(*,*) 'Error', err, 'in ', message, ':'
     write(*,*) errmsg
  endif
  return
end

