!Author: Junyi Cheng Oct.20 2018
!The mapping bewteen GEM and XGC
!Three major subroutines:
!1.mapping_init: setup for mapping and generate the coefficient for linear interpolation
!2.mapping_GEM_XGC: density mapping from GEM to XGC, 
!including 3 procedures: a) generate GEM 3d data b) mapping
!3.mapping_XGC_GEM: field mapping from XGC to GEM
!incuding 2 proceducres: a) read XGC data b) mapping and read the field data for GEM
!Noting that: nphi is the toroidal number in XGC, but phi is field value in GEM
module mapping
  use mpi
  implicit none
  
  type grid_type
    integer :: nnode ! number of nodes in a grid system
    
    !for nodes in .node file
    integer, allocatable :: rgn(:) ! region value for each node point
    real, allocatable :: x(:,:) ! R-Z position for each node point
    real, allocatable :: theta_geo(:) ! the real theta for each node point
    real, allocatable :: theta_flx(:) ! the flux theta for each node point
    
    !start and end indices for flux-surfaces in .ele file
    integer :: npsi_surf ! number of surface
    integer, allocatable :: surf_len(:) ! length of each surface
    integer :: surf_maxlen ! maximal length of the flux-surfaces
    integer, allocatable :: surf_idx(:,:) ! vertex indices for each flux-surface
  end type grid_type

  type mapping_GEM_XGC_linear_coef_type
    !the dimension is 1: density node number, 2: the \phi number

    !the related location for linear interpolation between y and \phi
    integer, dimension(:,:), allocatable :: jy1, jy2

    !the related weight for linear inerpolation between y and \phi
    real, dimension(:,:), allocatable :: wy1, wy2    
  end type mapping_GEM_XGC_linear_coef_type

#ifndef __MAP_PARALLEL
  type mapping_XGC_GEM_linear_coef_type
    !the dimension is 1: imx, 2: jmx 3: \theta number in XGC grid 

    !the related phi (zeta), theta location for linear interpolation between (\theta,\phi) and (y,\theta/z)
    integer, dimension(:,:,:), allocatable :: z1, z2, jl1, jl2, jr1, jr2
    !the related weight for phi (zeta), theta
    real, dimension(:,:,:), allocatable :: wz1, wz2, wl1, wl2, wr1, wr2
  end type mapping_XGC_GEM_linear_coef_type
#else
  type mapping_XGC_GEM_linear_coef_type
    !the dimension is 1: imx, 2: jmx 3: \theta number in XGC grid 

    !the related phi (zeta), theta location for linear interpolation between (\theta,\phi) and (y,\theta/z)
    integer, dimension(:,:,:), allocatable :: z1, z2, jl1, jl2, jr1, jr2
    !the related weight for phi (zeta), theta
    real, dimension(:,:,:), allocatable :: wz1, wz2, wl1, wl2, wr1, wr2
  end type mapping_XGC_GEM_linear_coef_type
#endif
 
  type(grid_type) :: grid
  type(mapping_GEM_XGC_linear_coef_type) :: mapping_GEM_XGC_linear_coef
  type(mapping_XGC_GEM_linear_coef_type) :: mapping_XGC_GEM_linear_coef

  !3d density/field for GEM
  real, dimension(:,:,:), allocatable :: density_3d!,field_3d

  !3d density/field for XGC
  real, dimension(:,:), allocatable :: pot_density_3d!,pot_field_3d

  !MPI decomposition for mapping
  integer :: x_id,gx_id,gxg_id,cce_surface_start,cce_surface_end,x_comm,gx_group,gx_comm,group_world 

  !input parameters
  integer :: cce_first_surface,cce_last_surface,cce_all_field_node_number,cce_all_surface_node_number,nwedge,cce_all_node_number,nphi,coef_opt
  integer :: procsperwriter=1
  character(256) :: eq_filename,node_filename,surf_filename
 
  logical, save :: test_mapping_logical=.false.
contains

!mapping setup and generate the coefficient for linear interpolation
subroutine mapping_init
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_field_node_number,cce_all_surface_node_number,&
  !  cce_all_node_number,nphi,eq_filename,node_filename,surf_filename,coef_opt
  use gem_com, only:MyId,imx,jmx,kmx
  use gem_equil, only:ntheta,dth,thflx
  use EZspline_obj
  use EZspline

  integer :: ierr
  !eq_filename
  character (len=80) :: eq_header
  integer :: eq_mr,eq_mz,eq_mpsi
  real :: eq_min_r,eq_max_r,eq_min_z,eq_max_z,eq_axis_r,eq_axis_z,eq_axis_b
  !nodefilename
  integer :: nodefile=10,n_n,dum
  real :: r_length,z_length,th_tmp,surf_len_local,theta_grid(0:ntheta)
  integer, dimension(:), allocatable :: surf_idx_local
  real, dimension(:), allocatable :: theta_flux
  !surffilename
  integer :: surffile=11,i,ii,iii,j,isize,dum2(1)
  !spline
  real, external :: interpol_spl_1d
  type(EZspline1_r8) :: spl_thflx
  
  namelist /mapping/ cce_first_surface,cce_last_surface,nphi,nwedge, coef_opt,eq_filename,node_filename,surf_filename,procsperwriter
  
  open(unit=20,file='mapping.in',status='old',action='read')
  read(20,nml=mapping)
  close(unit=20)

  open(unit=111,file='mapping.out',status='replace',action='write')

  !read the XGC grid information

  !the .eqd first, only need the eq_axis_r, eq_axis_z for the calculation of angle
  if(MyId==0)then
    write(*,*) 'read eqd file from XGC'
    open(9, file=eq_filename, status='old')
    read(9,300) eq_header
    read(9,200) eq_mr, eq_mz, eq_mpsi
    read(9,100) eq_min_r, eq_max_r, eq_min_z, eq_max_z
    read(9,100) eq_axis_r, eq_axis_z, eq_axis_b ! eq_axis_r, eq_axis_z, the axis R/Z 
    close(9)
    write(*,*) 'end of eqd reader'
  endif
  !broadcast the eq_axis_r, eq_axis_z
  !call mpi_bcast(eq_axis_r, 1, MPI_REAL8, 0, mpi_comm_world, ierr)
  !call mpi_bcast(eq_axis_z, 1, MPI_REAL8, 0, mpi_comm_world, ierr)

  !read XGC each node R/Z data in the plane
  if(MyId==0)then
    write(111,*) 'read .node file from XGC'
    call flush(111)
    open(nodefile, file=node_filename, status='old', form='formatted')
    read(nodefile,*) grid%nnode, dum,dum,dum
  endif
  !broadcast the node number
  call mpi_bcast(grid%nnode, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
  n_n=grid%nnode
  allocate(grid%x(2,n_n),grid%rgn(n_n))
  !read R Z grid information
  if(MyId==0)then
    do i=1,n_n
      read(nodefile,*) dum, grid%x(1,i), grid%x(2,i), grid%rgn(i)
    enddo
    close(nodefile)
    eq_axis_r=grid%x(1,1)
    eq_axis_z=grid%x(2,1)
    write(111,*) 'end of node reader'
    call flush(111)
  endif
  !broadcast the eq_axis_r, eq_axis_z
  call mpi_bcast(eq_axis_r, 1, MPI_REAL8, 0, mpi_comm_world, ierr)
  call mpi_bcast(eq_axis_z, 1, MPI_REAL8, 0, mpi_comm_world, ierr) 
  !broadcast R Z grid data
  call mpi_bcast(grid%x, n_n*2, MPI_REAL8, 0, mpi_comm_world, ierr)
  !broadcast boundary index
  call mpi_bcast(grid%rgn, n_n, MPI_INTEGER, 0, mpi_comm_world, ierr)

  !read the grid information along the flux surface
  if(MyId==0)then
    write(111,*) 'read .aif file from XGC'
    call flush(111)
    open(surffile, file=surf_filename, status='old', form='formatted')
    !read number of surfaces
    read(surffile,*) grid%npsi_surf
    !read length of each surface
    allocate(grid%surf_len(grid%npsi_surf))
    read(surffile,*) grid%surf_len
    !determine maximal length of the flux-surface
    dum2=maxval(grid%surf_len)
    grid%surf_maxlen=dum2(1)
    !read the vertex indices of each surface
    allocate(grid%surf_idx(grid%surf_maxlen,grid%npsi_surf))
    grid%surf_idx=0
    do i=1,grid%npsi_surf
      !do j=1,grid%surf_len(i)
      read(surffile,*) grid%surf_idx(1:grid%surf_len(i),i)
      !enddo
      !write(*,*)'surface number',i,'grid for each surface',grid%surf_idx(1:grid%surf_len(i),i)
    enddo
    close(surffile)
    write(111,*) 'end of .aif reader'
    call flush(111)
  endif
  !broadcast the array size for memory allocation
  call mpi_bcast(grid%npsi_surf, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
  call mpi_bcast(grid%surf_maxlen, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)then
     write(111,*)'before allocate in others'
     call flush(111)
  endif
  if(MyId .ne. 0)then
    allocate(grid%surf_len(grid%npsi_surf))
    allocate(grid%surf_idx(grid%surf_maxlen,grid%npsi_surf))
    grid%surf_idx=0
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)then
     write(111,*)'after allocate in others'
     call flush(111)
  endif
  !broadbcast flux-surface information
  isize=grid%npsi_surf
  call mpi_bcast(grid%surf_len, isize, MPI_INTEGER, 0, mpi_comm_world, ierr)
  isize=grid%surf_maxlen*grid%npsi_surf
  call mpi_bcast(grid%surf_idx, isize, MPI_INTEGER, 0, mpi_comm_world, ierr)

  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)then
     write(111,*)'after broadcast'
     call flush(111)
  endif
 
  !caculate cce_all_field_node_number,cce_all_surface_node_number,cce_all_node_number
  cce_all_surface_node_number=grid%surf_idx(grid%surf_len(cce_last_surface),cce_last_surface)-grid%surf_idx(1,cce_first_surface)+1
  cce_all_node_number=grid%nnode
  !theta in each node
  allocate(grid%theta_geo(n_n),grid%theta_flx(n_n))
  allocate(surf_idx_local(grid%surf_maxlen),theta_flux(grid%surf_maxlen))
  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)then
     write(111,*)'after allocate grid, surf_idx'
     call flush(111)
  endif
  if(MyId==0)then
    do i=1,n_n
      r_length=grid%x(1,i)-eq_axis_r
      z_length=grid%x(2,i)-eq_axis_z
      grid%theta_geo(i)=atan2(z_length,r_length)
    enddo
    theta_grid(ntheta/2)=0.
    do i=ntheta/2+1,ntheta
      th_tmp=real(i-ntheta/2)*dth
      theta_grid(i)=th_tmp
      theta_grid(ntheta-i)=-th_tmp
    enddo
    do ii=cce_first_surface,cce_last_surface
      surf_idx_local=grid%surf_idx(:,ii)
      surf_len_local=grid%surf_len(ii)
      iii=ii-cce_first_surface
      call init_1d_interpolation(spl_thflx,theta_grid,thflx(iii,:),ntheta+1,.false.)
      do i=1,int(surf_len_local)
        th_tmp=grid%theta_geo(surf_idx_local(i))
        theta_flux(i)=interpol_spl_1d(th_tmp,0,spl_thflx)
        grid%theta_flx(surf_idx_local(i))=theta_flux(i)
      enddo
      call finalize_1d_interpolation(spl_thflx) 
    enddo
  endif
  deallocate(surf_idx_local,theta_flux)
  call mpi_bcast(grid%theta_geo, n_n, MPI_REAL8, 0, mpi_comm_world, ierr)
  call mpi_bcast(grid%theta_flx, n_n, MPI_REAL8, 0, mpi_comm_world, ierr)
  !init the 3d data
  !3d density/field for GEM
  allocate(density_3d(0:imx,0:jmx,0:kmx))
  !allocate(field_3d(0:imx,0:jmx,0:kmx))
  !3d density/field for XGC
  allocate(pot_density_3d(cce_all_node_number,nphi))
  !allocate(pot_field_3d(cce_all_node_number,nphi))
  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)then
     write(111,*)'after allocate pot','imx,jmx,kmx',imx,jmx,kmx,'cce_all_surface_node_number,nphi',cce_all_surface_node_number,nphi,'cce_all_node_number,nphi',cce_all_node_number,nphi
     call flush(111)
  endif
  density_3d=0.0
  !field_3d=0.0
  pot_density_3d=0.0
  !pot_field_3d=0.0  
  
#ifdef __MAP_PARALLEL
  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)then
    write(111,*)'before mpi docomposition',imx,x_id,gx_id,cce_surface_start,cce_surface_end,x_comm,gx_comm
    call flush(111)
  endif
  !mapping mpi decomposition
  call mapping_mpi_decomposition(imx,x_id,gx_id,cce_surface_start,cce_surface_end,x_comm,gx_comm)
  !if(gx_id==0)then
    if(coef_opt==0)then
      !GEM to XGC
      call mpi_barrier(mpi_comm_world,ierr)
      if(myid==0)then
        write(111,*)'coef begins',MPI_WTIME()
        call flush(111)
      endif
      call setup_GEM_XGC_coef_new(mapping_GEM_XGC_linear_coef,grid)
      call mpi_barrier(mpi_comm_world,ierr)
      if(myid==0)then
        write(111,*)'GEM to XGC coef',MPI_WTIME()
        call flush(111)
      endif
      !XGC to GEM
      call setup_XGC_GEM_coef_new(mapping_XGC_GEM_linear_coef,grid)
      call mpi_barrier(mpi_comm_world,ierr)
      if(myid==0)then
        write(111,*)'XGC to GEM coef',MPI_WTIME()
        call flush(111)
      endif
       
      !store coef
      !call store_coef_new(mapping_GEM_XGC_linear_coef,mapping_XGC_GEM_linear_coef,grid)
      !if(myid==0)write(111,*)'store coef',MPI_WTIME()
    else
      !restore coef
      !call restore_coef(mapping_GEM_XGC_linear_coef,mapping_XGC_GEM_linear_coef,grid)
      write(*,*)'restore coef',MPI_WTIME()
    endif
  !endif
#else
  !call mapping_mpi_decomposition(imx,x_id,gx_id,cce_surface_start,cce_surface_end,x_comm,gx_comm) 
  !generate the coefficient for the linear interpolation
  if(MyId==0)then
    if(coef_opt==0)then
      !GEM to XGC
      write(111,*)'coef begins'
      call setup_GEM_XGC_coef(mapping_GEM_XGC_linear_coef,grid)
      write(111,*)'GEM to XGC coef'
      !XGC to GEM
      call setup_XGC_GEM_coef(mapping_XGC_GEM_linear_coef,grid)
      write(111,*)'XGC to GEM coef'
      !store coef
      call store_coef(mapping_GEM_XGC_linear_coef,mapping_XGC_GEM_linear_coef,grid)
      write(111,*)'store coef'
    else
      !restore coef
      call restore_coef(mapping_GEM_XGC_linear_coef,mapping_XGC_GEM_linear_coef,grid)
      write(111,*)'restore coef'
    endif
  endif
#endif

  !remind that the index for coupling is different with the index for simulation
  !if(gx_id==0)call mpi_barrier(mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)write(111,*)'after mapping_init'
  !close(111) 
100 format(4(e19.13,1x))
200 format(8I8)
300 format(a80)
end subroutine mapping_init

#ifdef __MAP_PARALLEL
!mapping mpi decomposition in x direction for parallel mapping
subroutine mapping_mpi_decomposition(imx,x_id,gx_id,cce_surface_start,cce_surface_end,x_comm,gx_comm)
  use mpi
  use gem_com, only: myid,numprocs,kmx,ntube,lx,ly,lz
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface
#ifdef EQUIL_OUTPUT
  use gem_equil, only:nr,ntheta,thflx,qhat,yfn,radius,hght,thfnz,zfnth,sf,bfld,dbdth,grcgt,jacob
#endif
#ifdef __ADIOS2__TEST
#include "adios2_macro.h"
  use adios2_comm_module
#endif
  implicit none

  integer, intent(in) :: imx
  integer, intent(out) :: x_id,gx_id,cce_surface_start,cce_surface_end,x_comm,gx_comm

  integer :: i,j,ierr
#ifdef EQUIL_OUTPUT
  integer :: err
  character(512) :: cce_filename
  integer*8 :: buf_id,buf_size,total_size

#ifdef __ADIOS2__TEST
  type(adios2_engine), save :: engine
  type(adios2_io), save :: io
  type(adios2_variable) ::var
  integer, save :: init=0
  integer(8), dimension(2) :: gdims_2d, goffset_2d, ldims_2d
  integer(8), dimension(1) :: gdims_1d, goffset_1d, ldims_1d
#endif

#endif


  x_id=0
  gx_id=0
  cce_surface_start=0
  cce_surface_end=0
  x_comm=0
  gx_comm=0
  
  !total MPI number is larger/equal than x grid number
  if(numprocs/(imx+1) .ge. 1)then
    !the number for each mpi bin
    i=int(numprocs/(imx+1))
    !x_id is the id cooresponding to imx
    x_id=mod(int(myid/i),imx+1)
    !gx_id is the color for the bin
    gx_id=mod(myid,i)
    !the myid is larger than (imx+1)*i, the total resource we want to use, the gx_id=1
    if(myid .ge. (imx+1)*i)gx_id=1
    cce_surface_start=x_id+cce_first_surface
    cce_surface_end=cce_surface_start
    if(cce_surface_end .gt. cce_last_surface)then
      write(*,*)'index in x direction exceeds the boundary,','i=',i,'myid=',myid,'cce_surface_start=',cce_surface_start,'cce_surface_end=',cce_surface_end,'x_id=',x_id,'gx_id=',gx_id
      stop 
    endif
    !write(*,*)'i=',i,'myid=',myid,' cce_surface_start=',cce_surface_start,' cce_surface_end=',cce_surface_end,'x_id=',x_id,'gx_id=',gx_id
    !if(gx_id==0)write(*,*)'gx_id==0,x_id=',x_id,'myid=',myid
  else
  !total MPI number is smaller than x grid number
  !since the imx and numprocs are not exact division, the last myid will have less points number in x direction
    i=int((imx+1)/numprocs)+1
    if(myid*i .le. imx)then
      x_id=myid
      gx_id=0    
      cce_surface_start=x_id*i+cce_first_surface
      if((myid+1)*i .le. imx)then
        cce_surface_end=cce_surface_start+i-1
      else
        cce_surface_end=cce_last_surface 
      endif
      if(cce_surface_end .gt. cce_last_surface)then
        write(*,*)'index in x direction exceeds the boundary'
        stop
      endif
    else
      x_id=myid
      gx_id=1
      cce_surface_start=cce_first_surface
      cce_surface_end=cce_first_surface
    endif
  endif
  !x_comm is the comm with the same gx_id
  !if(myid==0)write(*,*)'rank_num=',rank_num
  call mpi_comm_split(mpi_comm_world,gx_id,x_id,x_comm,ierr)
  call mpi_comm_split(mpi_comm_world,x_id,gx_id,gx_comm,ierr)
  !write(*,*)'myid=',myid,'gx_id=',gx_id,'x_id=',x_id,'mpi_comm_world=',mpi_comm_world,'x_comm=',x_comm,'gx_comm=',gx_comm
  call mpi_barrier(mpi_comm_world,ierr)
  !if(myid==0)write(*,*)'rank=',rank 
  !deallocate(rank)

#ifdef EQUIL_OUTPUT
 if(myid==0)then
    !write(*,*)'sf=',sf
    cce_filename='equil.bp'
#ifdef __ADIOS2__TEST
    
    if(init==0)then
      call adios2_declare_io(io, adios2obj, 'equil', err)

      gdims_2d(1)=nr+1
      gdims_2d(2)=ntheta+1
      goffset_2d(1)=0
      goffset_2d(2)=0
      ldims_2d(1)=nr+1
      ldims_2d(2)=ntheta+1

      call adios2_define_variable(var, io, "lx", adios2_type_dp, err)

      call adios2_define_variable(var, io, "ly", adios2_type_dp, err)

      call adios2_define_variable(var, io, "lz", adios2_type_dp, err)

      call adios2_define_variable(var, io, 'thflx', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      call adios2_define_variable(var, io, 'qhat', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)
      
      call adios2_define_variable(var, io, 'yfn', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      gdims_1d(1)=nr+1
      goffset_1d(1)=0
      ldims_1d(1)=nr+1

      call adios2_define_variable(var, io, 'sf', adios2_type_dp, 1, &
           gdims_1d, &
           goffset_1d, &
           ldims_1d, &
           .true., err)
   
      call adios2_define_variable(var, io, 'radius', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      call adios2_define_variable(var, io, 'hght', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)
      
      call adios2_define_variable(var, io, 'bfld', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      call adios2_define_variable(var, io, 'dbdth', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      call adios2_define_variable(var, io, 'grcgt', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      call adios2_define_variable(var, io, 'jacob', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      gdims_1d(1)=ntheta+1
      goffset_1d(1)=0
      ldims_1d(1)=ntheta+1

      call adios2_define_variable(var, io, 'zfnth', adios2_type_dp, 1, &
           gdims_1d, &
           goffset_1d, &
           ldims_1d, &
           .true., err)

      call adios2_define_variable(var, io, 'thfnz', adios2_type_dp, 1, &
           gdims_1d, &
           goffset_1d, &
           ldims_1d, &
           .true., err)
      
      gdims_1d(1)=grid%npsi_surf
      goffset_1d(1)=0
      ldims_1d(1)=grid%npsi_surf

      
      call adios2_define_variable(var, io, 'surf_len', adios2_type_integer4, 1, &
           gdims_1d, &
           goffset_1d, &
           ldims_1d, &
           .true., err)

      gdims_2d(1)=grid%surf_maxlen
      gdims_2d(2)=grid%npsi_surf
      goffset_2d(1)=0
      goffset_2d(2)=0
      ldims_2d(1)=grid%surf_maxlen
      ldims_2d(2)=grid%npsi_surf

      call adios2_define_variable(var, io, 'surf_idx', adios2_type_integer4, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      gdims_2d(1)=2
      gdims_2d(2)=grid%nnode
      goffset_2d(1)=0
      goffset_2d(2)=0
      ldims_2d(1)=2
      ldims_2d(2)=grid%nnode

      call adios2_define_variable(var, io, 'rz', adios2_type_dp, 2, &
           gdims_2d, &
           goffset_2d, &
           ldims_2d, &
           .true., err)

      gdims_1d(1)=grid%nnode
      goffset_1d(1)=0
      ldims_1d(1)=grid%nnode

      call adios2_define_variable(var, io, 'theta_geo', adios2_type_dp, 1, &
           gdims_1d, &
           goffset_1d, &
           ldims_1d, &
           .true., err)

      call adios2_define_variable(var, io, 'theta_flx', adios2_type_dp, 1, &
           gdims_1d, &
           goffset_1d, &
           ldims_1d, &
           .true., err)

      init=1
    endif

    call adios2_open(engine, io, cce_filename, adios2_mode_write, mpi_comm_self, err)
    call adios2_begin_step(engine, adios2_step_mode_append, err)
    call adios2_put(engine, 'lx', lx, err)
    call adios2_put(engine, 'ly', ly, err)
    call adios2_put(engine, 'lz', lz, err)
    call adios2_put(engine, 'thflx', thflx, err)
    call adios2_put(engine, 'qhat', qhat, err)
    call adios2_put(engine, 'yfn', yfn, err)
    call adios2_put(engine, 'sf', sf, err)
    call adios2_put(engine, 'radius', radius, err)
    call adios2_put(engine, 'hght', hght, err)
    call adios2_put(engine, 'zfnth', zfnth, err)
    call adios2_put(engine, 'thfnz', thfnz, err)
    call adios2_put(engine, 'surf_len', grid%surf_len, err)
    call adios2_put(engine, 'surf_idx', grid%surf_idx, err)
    call adios2_put(engine, 'rz', grid%x, err)
    call adios2_put(engine, 'theta_geo', grid%theta_geo, err)
    call adios2_put(engine, 'theta_flx', grid%theta_flx, err)
    call adios2_put(engine, 'bfld', bfld, err)
    call adios2_put(engine, 'dbdth', dbdth, err)
    call adios2_put(engine, 'grcgt', grcgt, err)
    call adios2_put(engine, 'jacob', jacob, err)
    call adios2_end_step(engine, err)
    call adios2_close(engine, err)
    
#endif
 endif
#endif 
end subroutine mapping_mpi_decomposition
#endif

#ifndef __MAP_PARALLEL
!mapping from GEM to XGC
subroutine mapping_GEM_XGC(density_gem,density_xgc,den_tmp,grid,coef)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi
  use gem_com, only:MyId,imx,jmx,kmx

  integer :: ierr
  real, intent(inout) :: density_gem(0:imx,0:jmx,0:kmx)
  real :: density_xgc_tmp(cce_all_surface_node_number,nphi)
  real, intent(out) :: density_xgc(cce_all_node_number,nphi)
  real, intent(inout) :: den_tmp(0:imx,0:jmx,0:1)
  type(grid_type), intent(in) :: grid
  type(mapping_GEM_XGC_linear_coef_type), intent(in) :: coef
  
  integer: node_start,node_end
  density_gem=0.0
  density_xgc=0.0
  density_xgc_tmp=0.0

  !generate the 3d density
  if(MyId==0)write(*,*)'before 3d GEM',MPI_WTIME()
  call generate_3d_density_GEM(density_gem,den_tmp)

  if(MyId==0)write(*,*)'after 3d GEM, begin mapping GEM XGC',MPI_WTIME()
  if(MyId==0)then
    !mapping
    call mapping_GEM_XGC_core(density_gem,density_xgc_tmp,grid,coef)
  endif
  if(Myid==0)write(*,*)'after mapping GEM XGC',MPI_WTIME()
  node_start=grid%surf_idx(1,cce_first_surface)
  node_end=grid%surf_idx(grid%surf_len(cce_last_surface),cce_last_surface)
  density_xgc(node_start:node_end)=density_xgc_tmp(:)
  call mpi_barrier(mpi_comm_world,ierr)

end subroutine mapping_GEM_XGC
#else
!mapping from GEM to XGC
subroutine mapping_GEM_XGC_new(den_tmp,grid,coef,nor)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi
  use gem_com, only:MyId,imx,jmx,kmx

  integer :: ierr
  !real, intent(inout) :: density_gem(0:imx,0:jmx,0:kmx)
  !real :: density_xgc_tmp(cce_all_surface_node_number,nphi)
  !real, intent(out) :: density_xgc(cce_all_node_number,nphi)
  real, intent(inout) :: den_tmp(0:imx,0:jmx,0:1)
  type(grid_type), intent(in) :: grid
  type(mapping_GEM_XGC_linear_coef_type), intent(in) :: coef
  real, intent(in) :: nor

  integer :: node_start,node_end
  density_3d=0.0
  pot_density_3d=0.0
  !density_xgc_tmp=0.0

  !generate the 3d density
  call mpi_barrier(mpi_comm_world,ierr)
  if(MyId==0)then
     write(111,*)'before 3d GEM',MPI_WTIME()
     call flush(111)
  endif
  call generate_3d_density_GEM(den_tmp)

  call mpi_barrier(mpi_comm_world,ierr)
  if(MyId==0)then
    write(111,*)'after 3d GEM, begin mapping GEM XGC',MPI_WTIME()
    call flush(111)
  endif
  !mapping
  node_start=grid%surf_idx(1,cce_first_surface)
  node_end=grid%surf_idx(grid%surf_len(cce_last_surface),cce_last_surface)
  !call mapping_GEM_XGC_core_new(density_3d,pot_density_3d(node_start:node_end,:),grid,coef)
  call mapping_GEM_XGC_core_new(node_start,node_end,grid,coef,nor)
  call mpi_barrier(mpi_comm_world,ierr)
  if(Myid==0)then
    write(111,*)'after mapping GEM XGC',MPI_WTIME()
    call flush(111)
  endif
  !node_start=grid%surf_idx(1,cce_first_surface)
  !node_end=grid%surf_idx(grid%surf_len(cce_last_surface),cce_last_surface)
  !density_xgc(node_start:node_end,:)=density_xgc_tmp(1:cce_all_surface_node_number,:) 
  call mpi_barrier(mpi_comm_world,ierr)

end subroutine mapping_GEM_XGC_new

#endif

!generate the 3d density data in GEM
subroutine generate_3d_density_GEM(den_tmp)
  use gem_com, only:imx,jmx,kmx,grid_comm,tube_comm,gclr,tclr,glst,tlst

  integer :: ierr,icount
  real, intent(inout) :: den_tmp(0:imx,0:jmx,0:1)
  !real, intent(out) :: density_gem(0:imx,0:jmx,0:kmx)
  real :: density_gem_tmp(0:imx,0:jmx,0:kmx)

  !den_tmp=den(1,:,:,:)/q(1)
  !density_gem=0.0
  !store the density each z domain  
  density_3d(:,:,gclr)=den_tmp(:,:,0)
  if(gclr==0)then
    density_3d(:,:,kmx)=den_tmp(:,:,0)
  endif

  !sum of the 3d density  
  icount=(imx+1)*(jmx+1)*(kmx+1)
  call mpi_allreduce(MPI_IN_PLACE,density_3d,icount,MPI_REAL8,mpi_sum,mpi_comm_world,ierr)

  call mpi_barrier(mpi_comm_world,ierr)

end subroutine generate_3d_density_GEM

#ifndef __MAP_PARALLEL
!the major part for mapping from GEM to XGC
subroutine mapping_GEM_XGC_core(density_gem,density_xgc,grid,coef)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi
  use gem_com, only:imx,jmx,kmx,pi,pi2,lz,myid
  use gem_equil, only:ntheta,delz,thfnz,q0,dth,thflx
  use EZspline_obj
  use EZspline

  type(grid_type), intent(in) :: grid
  type(mapping_GEM_XGC_linear_coef_type), intent(in) :: coef
  type(EZspline1_r8) :: spl_z,spl_thflx
  real, intent(in) :: density_gem(0:imx,0:jmx,0:kmx)
  real, intent(out) :: density_xgc(cce_all_surface_node_number,nphi)
  integer :: i,ii,j,iii,j1,j2,k,node_start
  integer :: surf_len_local,surf_idx_local(grid%surf_maxlen)
  real :: w1,w2,wz0,wz1,th_tmp,theta_z_gem(0:kmx),density_gem_z(0:kmx), &
    theta_local(grid%surf_maxlen),theta_flux(grid%surf_maxlen),density_gem_z_xgc(0:jmx,grid%surf_maxlen), &
    theta_grid(0:ntheta),z_gem(0:kmx),theta_z_gem_inv(0:kmx),density_gem_z_inv(0:kmx),theta_z_gem_flx(0:kmx)
  real, external :: interpol_spl_1d

  !the start node number
  node_start=grid%surf_idx(1,cce_first_surface)

  if(myid==0)write(*,*)'point 1, GEM_XGC,',MPI_WTIME()
  !the theta in z direction
  do i=0,kmx
     z_gem(i)=real(i)*lz/real(kmx)
     k=floor(z_gem(i)/delz)
     wz0=real(k+1)-z_gem(i)/delz
     wz1=1.-wz0
     theta_z_gem(i)=wz0*thfnz(k)+wz1*thfnz(modulo(k+1,ntheta))
  enddo
  if(myid==0)write(*,*)'point 2, GEM_XGC,',MPI_WTIME()
  theta_grid(ntheta/2)=0.
  do i=ntheta/2+1,ntheta
     th_tmp=real(i-ntheta/2)*dth
     theta_grid(i)=th_tmp
     theta_grid(ntheta-i)=-th_tmp
  enddo
  if(myid==0)write(*,*)'point 3, GEM_XGC,',MPI_WTIME()
  if(q0<0)then
    do i=0,kmx
       theta_z_gem_inv(i)=theta_z_gem(kmx-i)
    enddo
    theta_z_gem=theta_z_gem_inv
  endif
  if(myid==0)write(*,*)'point 4, GEM_XGC,',MPI_WTIME()
  !generate each surface
  do ii=cce_first_surface,cce_last_surface
 if(myid==0 .and. ii==150)write(*,*)'point 4.0, GEM_XGC',MPI_WTIME()
    surf_idx_local=grid%surf_idx(:,ii)
    surf_len_local=grid%surf_len(ii)
 if(myid==0 .and. ii==150)write(*,*)'point 4.1, GEM_XGC',MPI_WTIME() 
    theta_local=0.0
    do j=1,surf_len_local
      theta_flux(j)=grid%theta_flx(surf_idx_local(j))
    enddo
 if(myid==0 .and. ii==150)write(*,*)'point 4.2, GEM_XGC',MPI_WTIME() 
    !i is the idex in x direction
    i=ii-cce_first_surface
    
    call init_1d_interpolation(spl_thflx,theta_grid,thflx(i,:),ntheta+1,.false.)
    do k=0,kmx
      theta_z_gem_flx(k)=interpol_spl_1d(theta_z_gem(k),0,spl_thflx)       
    enddo
    call finalize_1d_interpolation(spl_thflx) 
 if(myid==0 .and. ii==150)write(*,*)'point 4.3, GEM_XGC',MPI_WTIME()
    !generate the (x,y,z_xgc) through cubic interpolation
    density_gem_z_xgc=0.0
    do j=0,jmx
      !the gem data in z direction
      do k=0,kmx
        density_gem_z(k)=density_gem(i,j,k)
      enddo
      if(q0<0)then
        do k=0,kmx
           density_gem_z_inv(k)=density_gem_z(kmx-k)
        enddo
        density_gem_z=density_gem_z_inv
      endif
      !cubic interpolation
      call init_1d_interpolation(spl_z,theta_z_gem_flx,density_gem_z,kmx+1,.false.)
      do k=1,surf_len_local
        density_gem_z_xgc(j,k)=interpol_spl_1d(theta_flux(k),0,spl_z)
      enddo
      call finalize_1d_interpolation(spl_z)
    enddo
 if(myid==0 .and. ii==150)write(*,*)'point 4.4, GEM_XGC',MPI_WTIME()

    !generate the xgc data through linear interpolation
    do j=1,surf_len_local
      do k=1,nphi
        !the index for density_xgc and coef is from 1 to cce_all_surface_node_number
        iii=surf_idx_local(j)-node_start+1
        j1=coef%jy1(iii,k)
        j2=coef%jy2(iii,k)
        w1=coef%wy1(iii,k)
        w2=coef%wy2(iii,k)
        density_xgc(iii,k)=w1*density_gem_z_xgc(j1,j)+w2*density_gem_z_xgc(j2,j)
      enddo
    enddo
  if(myid==0 .and. ii==150)write(*,*)'point 4.5, GEM_XGC',MPI_WTIME()
  enddo
  if(myid==0)write(*,*)'point 5, GEM_XGC,',MPI_WTIME()
end subroutine mapping_GEM_XGC_core

#else
!the major part for mapping from GEM to XGC
subroutine mapping_GEM_XGC_core_new(start_num,end_num,grid,coef,nor)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi
  use gem_com, only:imx,jmx,kmx,pi,pi2,lz,myid,pi
  use gem_equil, only:ntheta,delz,thfnz,q0,dth,thflx
  use EZspline_obj
  use EZspline

  integer, intent(in) :: start_num, end_num
  real, intent(in) :: nor
  type(grid_type), intent(in) :: grid
  type(mapping_GEM_XGC_linear_coef_type), intent(in) :: coef
  type(EZspline1_r8) :: spl_z,spl_thflx
  !real, intent(in) :: density_gem(0:imx,0:jmx,0:kmx)
  !real, intent(out) :: density_xgc(1:cce_all_surface_node_number,nphi)
  integer :: i,ii,j,iii,idx_coef,icount,ierr,j1,j2,k,node_start,local_start,i_diag
  integer :: surf_len_local,surf_idx_local(grid%surf_maxlen)
  real :: w1,w2,wz0,wz1,th_tmp,theta_z_gem(0:kmx),density_gem_z(0:kmx), &
    theta_local(grid%surf_maxlen),theta_flux(grid%surf_maxlen),density_gem_z_xgc(0:jmx,grid%surf_maxlen), &
    theta_grid(0:ntheta),z_gem(0:kmx),theta_z_gem_inv(0:kmx),density_gem_z_inv(0:kmx),theta_z_gem_flx(0:kmx)!,density_xgc_tmp(cce_all_surface_node_number,nphi)
  real, external :: interpol_spl_1d

  logical, save :: test_logical=.true. 

  !density_xgc=0.0

  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)then
    write(111,*)'gx_id, before real mapping', gx_id, nphi, cce_all_surface_node_number
    !write(111,*)'density_xgc',density_xgc
    call flush(111)
  endif

  if(gx_id==0)then
    !the start node number
    node_start=grid%surf_idx(1,cce_first_surface)
    local_start=grid%surf_idx(1,cce_surface_start)

    if(myid==0)then
       write(111,*)'point 1, GEM_XGC,',MPI_WTIME()
       call flush(111)
    endif
    !the theta in z direction
    do i=0,kmx
      z_gem(i)=real(i)/real(kmx)*lz
      k=floor(z_gem(i)/delz)
      wz0=real(k+1)-z_gem(i)/delz
      wz1=1-wz0
      theta_z_gem(i)=wz0*thfnz(k)+wz1*thfnz(modulo(k,ntheta))
    enddo
    !write(111,*)'z_gem',z_gem
    !write(111,*)'delz,lz',delz,lz
    !write(111,*)'thfnz',thfnz/pi
    !write(111,*)'theta_z_gem',theta_z_gem/pi
    !call flush(111)
    if(myid==0)then
       write(111,*)'point 2, GEM_XGC,',MPI_WTIME()
       call flush(111)
    endif
    theta_grid(ntheta/2)=0.
    do i=ntheta/2+1,ntheta
      th_tmp=real(i-ntheta/2)*dth
      theta_grid(i)=th_tmp
      theta_grid(ntheta-i)=-th_tmp
    enddo
    if(myid==0)then
      write(111,*)'point 3, GEM_XGC,',MPI_WTIME()
      call flush(111)
    endif
    if(q0<0)then
      do i=0,kmx
        theta_z_gem_inv(i)=theta_z_gem(kmx-i)
      enddo
      theta_z_gem=theta_z_gem_inv
    endif
    if(myid==0)then
      write(111,*)'point 4, GEM_XGC,',MPI_WTIME()
      call flush(111)
    endif
    !generate each surface
    do ii=cce_surface_start,cce_surface_end
      surf_idx_local=grid%surf_idx(:,ii)
      surf_len_local=grid%surf_len(ii)
      theta_local=0.0
      do j=1,surf_len_local
        theta_flux(j)=grid%theta_flx(surf_idx_local(j))
      enddo
      !i is the idex in x direction
     ! do i_diag=0,size(theta_grid)-2
     !      if((theta_grid(i_diag+1)-theta_grid(i_diag))<0)then
     !        write(111,*)'error theta_gridi',myid,ii
     !        write(111,*)size(theta_grid),i_diag,theta_grid(i_diag),theta_grid(i_diag+1)
     !        call flush(111)
     !      endif
     ! enddo

      i=ii-cce_first_surface
      call init_1d_interpolation(spl_thflx,theta_grid,thflx(i,:),ntheta+1,.false.)
      do k=0,kmx
        theta_z_gem_flx(k)=interpol_spl_1d(theta_z_gem(k),0,spl_thflx)       
      enddo
      call finalize_1d_interpolation(spl_thflx) 
      !generate the (x,y,z_xgc) through cubic interpolation
      density_gem_z_xgc=0.0
      do j=0,jmx
        !the gem data in z direction
        do k=0,kmx
          density_gem_z(k)=density_3d(i,j,k)
        enddo
        if(q0<0)then
          do k=0,kmx
            density_gem_z_inv(k)=density_gem_z(kmx-k)
          enddo
          density_gem_z=density_gem_z_inv
        endif
        !cubic interpolation
        !do i_diag=1,size(theta_z_gem_flx)-2
        !   if((theta_z_gem_flx(i_diag+1)-theta_z_gem_flx(i_diag))<0)then
        !     write(111,*)'error theta_z_gem_flx',myid,ii
        !     write(111,*)size(theta_z_gem_flx),i_diag,theta_z_gem_flx(i_diag),theta_z_gem_flx(i_diag+1)
        !     call flush(111)
        !   endif
        !enddo
        !write(111,*)'diag',i,ii
        !write(111,*)'theta_z_gem_flx'
        !write(111,*)theta_z_gem_flx/pi
        !write(111,*)'theta_grid'
        !write(111,*)theta_grid/pi
        !write(111,*)'thflx'
        !write(111,*)thflx(i,:)/pi
        !write(111,*)'theta_z_gem'
        !write(111,*)theta_z_gem/pi
        !write(111,*)'thfnz'
        !write(111,*)thfnz/pi
        !call flush(111)
        call init_1d_interpolation(spl_z,theta_z_gem_flx,density_gem_z,kmx+1,.false.)
        !!$omp parallel do
        do k=1,surf_len_local
          density_gem_z_xgc(j,k)=interpol_spl_1d(theta_flux(k),0,spl_z)
          !if(abs(theta_flux(k))==pi)then
          !  density_gem_z_xgc(j,k)=density_gem_z(0)
          !endif
        enddo
        call finalize_1d_interpolation(spl_z)
      enddo
      !if(ii==50 .and. (test_mapping_logical .eq. .true.))then
      !  test_mapping_logical = .false.
      !  open(124,file='test_GEM_XGC_z.out',status='unknown',position='append')
      !  do k=1,surf_len_local
      !     write(124,*)density_gem_z_xgc(0,k)
      !  enddo
      !  close(124)
      !endif
      !generate the xgc data through linear interpolation
      !!$omp parallel do
      do j=1,surf_len_local
        do k=1,nphi
          !the index for density_xgc and coef is from 1 to cce_all_surface_node_number
          idx_coef=surf_idx_local(j)-local_start+1
          iii=surf_idx_local(j)-node_start+1+start_num
          j1=coef%jy1(idx_coef,k)
          j2=coef%jy2(idx_coef,k)
          w1=coef%wy1(idx_coef,k)
          w2=coef%wy2(idx_coef,k)
          pot_density_3d(iii,k)=w1*density_gem_z_xgc(j1,j)+w2*density_gem_z_xgc(j2,j)*nor
        enddo
      enddo
    enddo
  endif
    
  !if(myid==0)write(*,*)'point 5, GEM_XGC,',MPI_WTIME()
  icount=cce_all_node_number*nphi
  call mpi_allreduce(MPI_IN_PLACE,pot_density_3d,icount,MPI_REAL8,mpi_sum,x_comm,ierr)
  call mpi_bcast(pot_density_3d,icount,MPI_REAL8,0,mpi_comm_world,ierr)
  call mpi_barrier(gx_comm,ierr)
  !if(myid==0)write(*,*)'point 6, GEM_XGC,',MPI_WTIME()
end subroutine mapping_GEM_XGC_core_new

#endif

#ifndef __MAP_PARALLEL
!generate the coefficient for GEM to XGC linear interpolation
subroutine setup_GEM_XGC_coef(coef,grid)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi,nwedge
  use gem_equil, only: sf,r0,q0,q0abs,dth,ntheta,thflx
  use gem_com, only: dy,ly,pi,pi2,myid
  use EZspline_obj
  use EZspline

  type(mapping_GEM_XGC_linear_coef_type), intent(out) :: coef
  type(grid_type), intent(in) :: grid

  type(EZspline1_r8) :: spl_thflx
  real, external :: interpol_spl_1d
  
  integer :: i,j,ii,iii,j1,j2,surf_len_local,node_start,idx_coef,surf_idx_local(grid%surf_maxlen)
  real :: q_local,phi_tmp,y_tmp,th_tmp,w1,w2,theta_grid(0:ntheta),theta_flux(grid%surf_maxlen)

  if(myid==0)write(111,*)'point 0, GEM XGC setup,',MPI_WTIME()
  !nwedge=lymult
  node_start=grid%surf_idx(1,cce_first_surface)
  !claim the memory for coef_GEM_XGC_linear_interpolation
  allocate(coef%jy1(cce_all_surface_node_number,nphi),coef%jy2(cce_all_surface_node_number,nphi),coef%wy1(cce_all_surface_node_number,nphi),coef%wy2(cce_all_surface_node_number,nphi))

  coef%jy1=0
  coef%jy2=0
  coef%wy1=0.0
  coef%wy2=0.0
  
  if(myid==0)write(111,*)'point 1, GEM XGC steup,',MPI_WTIME()
  theta_grid(ntheta/2)=0.
  do i=ntheta/2+1,ntheta
     th_tmp=real(i-ntheta/2)*dth
     theta_grid(i)=th_tmp
     theta_grid(ntheta-i)=-th_tmp
  enddo
  if(myid==0)write(111,*)'point 2, GEM XGC steup,',MPI_WTIME()
  do ii=cce_first_surface,cce_last_surface
    if(myid==0 .and. ii==150)write(*,*)'point 2.1, GEM XGC steup,',MPI_WTIME()
    surf_idx_local=grid%surf_idx(:,ii)
    surf_len_local=grid%surf_len(ii)
    q_local=sf(ii-cce_first_surface) !the local q for the surface
    iii=ii-cce_first_surface ! the corresponding GEM grid

    !write(*,*)'point 2, GEM-XGC','iii=',iii
    call init_1d_interpolation(spl_thflx,theta_grid,thflx(iii,:),ntheta+1,.false.)
    !write(*,*)'point 3, GEM-XGC', 'iii=',iii
    do i=1,surf_len_local
      th_tmp=grid%theta_geo(surf_idx_local(i))
      theta_flux(i)=interpol_spl_1d(th_tmp,0,spl_thflx)  
    enddo

    call finalize_1d_interpolation(spl_thflx)
    if(myid==0 .and. ii==150)write(111,*)'point 2.2, GEM XGC steup,',MPI_WTIME()
    do i=1,surf_len_local
      do j=1,nphi
        phi_tmp=real(j-1)*pi2/real(nphi)/real(nwedge)
        y_tmp=modulo(r0/q0*(q_local*theta_flux(i)-phi_tmp),ly)
        !calcaulte the index and weight for linear interpolation
        call search_y(j1,j2,w1,w2,dy,ly,y_tmp)
        !the idx for GEM to XGC, noting data from 1 to cce_all_surface_node_number
        idx_coef=surf_idx_local(i)-node_start+1
        coef%jy1(idx_coef,j)=j1
        coef%jy2(idx_coef,j)=j2
        coef%wy1(idx_coef,j)=w1
        coef%wy2(idx_coef,j)=w2
      enddo
    enddo
    if(myid==0 .and. ii==150)write(*,*)'point 2.3, GEM XGC steup,',MPI_WTIME()
  enddo
  
  if(myid==0)write(*,*)'point 3, GEM XGC steup,',MPI_WTIME() 
end subroutine setup_GEM_XGC_coef

#else
!generate the coefficient for GEM to XGC linear interpolation
subroutine setup_GEM_XGC_coef_new(coef,grid)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi,nwedge
  use gem_equil, only: sf,r0,q0,q0abs,dth,ntheta,thflx
  use gem_com, only: dy,ly,pi,pi2,myid
  use EZspline_obj
  use EZspline

  type(mapping_GEM_XGC_linear_coef_type), intent(out) :: coef
  type(grid_type), intent(in) :: grid

  type(EZspline1_r8) :: spl_thflx
  real, external :: interpol_spl_1d

  integer :: i,j,ii,iii,ierr,icount,j1,j2,surf_len_local,local_total_len,node_start,local_start,idx_coef,surf_idx_local(grid%surf_maxlen)
  real :: q_local,phi_tmp,y_tmp,th_tmp,w1,w2,theta_grid(0:ntheta),theta_flux(grid%surf_maxlen)

  !if(gx_id==0)then
     !write(*,*)'gx_id==0,myid=',myid
  if(myid==0)then
     write(111,*)'point 0, GEM XGC setup,',MPI_WTIME()
     call flush(111)
  endif
  !nwedge=lymult
  node_start=grid%surf_idx(1,cce_first_surface)
  local_start=grid%surf_idx(1,cce_surface_start)

  local_total_len=0
  do ii=cce_surface_start,cce_surface_end
    surf_len_local=grid%surf_len(ii)
    local_total_len=local_total_len+surf_len_local
  enddo
  !claim the memory for coef_GEM_XGC_linear_interpolation
  allocate(coef%jy1(local_total_len,nphi),coef%jy2(local_total_len,nphi),coef%wy1(local_total_len,nphi),coef%wy2(local_total_len,nphi))

  coef%jy1=0
  coef%jy2=0
  coef%wy1=0.0
  coef%wy2=0.0
  if(gx_id==0)then
    if(myid==0)then
       write(111,*)'point 1, GEM XGC steup,',MPI_WTIME()
       call flush(111)
    endif
    theta_grid(ntheta/2)=0.
    do i=ntheta/2+1,ntheta
      th_tmp=real(i-ntheta/2)*dth
      theta_grid(i)=th_tmp
      theta_grid(ntheta-i)=-th_tmp
    enddo
    if(myid==0)then
      write(111,*)'point 2, GEM XGC steup,',MPI_WTIME()
      call flush(111)
    endif
    do ii=cce_surface_start,cce_surface_end
      !if(ii==150)write(*,*)'point 2.1, GEM XGC steup,',MPI_WTIME()
      surf_idx_local=grid%surf_idx(:,ii)
      surf_len_local=grid%surf_len(ii)
      q_local=sf(ii-cce_first_surface) !the local q for the surface
      iii=ii-cce_first_surface ! the corresponding GEM grid

      !write(*,*)'point 2, GEM-XGC','iii=',iii
      call init_1d_interpolation(spl_thflx,theta_grid,thflx(iii,:),ntheta+1,.false.)
      !write(*,*)'point 3, GEM-XGC', 'iii=',iii
      do i=1,surf_len_local
        th_tmp=grid%theta_geo(surf_idx_local(i))
        theta_flux(i)=interpol_spl_1d(th_tmp,0,spl_thflx)
      enddo

      call finalize_1d_interpolation(spl_thflx)
      !if(ii==150)write(*,*)'point 2.2, GEM XGC steup,',MPI_WTIME()
      if(ii==50)then
        open(123,file='coef_GEM_XGC.out',status='unknown',position='append')
      endif
      !!$omp parallel do
      do i=1,surf_len_local
        do j=1,nphi
          phi_tmp=real(j-1)*pi2/real(nphi)/real(nwedge)
          y_tmp=modulo(r0/q0*(q_local*theta_flux(i)-phi_tmp),ly)
          !calcaulte the index and weight for linear interpolation
          call search_y(j1,j2,w1,w2,dy,ly,y_tmp)
          !the idx for GEM to XGC, noting data from 1 to cce_all_surface_node_number
          idx_coef=surf_idx_local(i)-local_start+1
          coef%jy1(idx_coef,j)=j1
          coef%jy2(idx_coef,j)=j2
          coef%wy1(idx_coef,j)=w1
          coef%wy2(idx_coef,j)=w2
          if(ii==50)then
             write(123,*)i,j,j1,j2,w1,w2,y_tmp
             call flush(123)
          endif
        enddo
      enddo
      !if(ii==150)write(*,*)'point 2.3, GEM XGC steup,',MPI_WTIME()
      if(ii==50)then
        close(123)
      endif
    enddo
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  if(myid==0)then
    write(111,*)'point 3.0, GEM XGC setup,',MPI_WTIME()
    call flush(111)
  endif
  !sum of the coefficient

  !if(gx_id==0)then
    !icount=cce_all_surface_node_number*nphi
    !tmp_jy=0
    !call mpi_allreduce(coef%jy1,tmp_jy,icount,MPI_INTEGER,mpi_sum,gx_comm,ierr)
    !coef%jy1=tmp_jy
    !call mpi_barrier(gx_comm,ierr)
    !tmp_jy=0
    !call mpi_allreduce(coef%jy2,tmp_jy,icount,MPI_INTEGER,mpi_sum,gx_comm,ierr)
    !coef%jy2=tmp_jy
    !call mpi_barrier(gx_comm,ierr)
    !tmp_wy=0
    !call mpi_allreduce(coef%wy1,tmp_wy,icount,MPI_REAL8,mpi_sum,gx_comm,ierr)
    !coef%wy1=tmp_wy
    !call mpi_barrier(gx_comm,ierr)
    !tmp_wy=0
    !call mpi_allreduce(coef%wy2,tmp_wy,icount,MPI_REAL8,mpi_sum,gx_comm,ierr)
    !coef%wy2=tmp_wy
    !call mpi_barrier(gx_comm,ierr)
  !endif
  !if(myid==0)write(*,*)'point 3, GEM XGC steup,',MPI_WTIME()
  !broadcast the index for mapping
  !call mpi_bcast(coef%jy1,cce_all_surface_node_number*nphi,MPI_INTEGER,0,mpi_comm_world,ierr)
  !call mpi_bcast(coef%jy2,cce_all_surface_node_number*nphi,MPI_INTEGER,0,mpi_comm_world,ierr)
  !broadcast the weight for mapping
  !call mpi_bcast(coef%wy1,cce_all_surface_node_number*nphi,MPI_REAL8,0,mpi_comm_world,ierr)
  !call mpi_bcast(coef%wy2,cce_all_surface_node_number*nphi,MPI_REAL8,0,mpi_comm_world,ierr)
  !if(myid==0)write(*,*)'point 4, GEM XGC steup,',MPI_WTIME()
end subroutine setup_GEM_XGC_coef_new
#endif

!generate the coefficient in y direction for each points
subroutine search_y(j1,j2,w1,w2,dlength,length,tmp)
  integer, intent(out) :: j1,j2
  real, intent(out) :: w1,w2
  real, intent(in) :: dlength,length
  real, intent(inout) :: tmp

  if(tmp<0. .or. tmp>=length)then
    tmp=modulo(tmp,length)
  endif

  j1=int(tmp/dlength)
  j2=j1+1
  w2=(tmp-real(j1)*dlength)/dlength
  w1=1.-w2

end subroutine search_y

#ifndef __MAP_PARALLEL
!mapping from XGC to GEM
subroutine mapping_XGC_GEM(field_xgc,field_gem,field_gem_3d,grid,coef)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi,&
  !  cce_first_field,cce_last_field,cce_all_field_node_number,cce_all_node_number
  use gem_com, only:MyId,imx,jmx,kmx,pi,pi2

  integer :: ierr
  real, intent(in) :: field_xgc(cce_all_node_number,nphi)
  real, intent(out) :: field_gem(0:imx,0:jmx,0:1)
  real, intent(out) :: field_gem_3d(0:imx,0:jmx,0:kmx)
  type(grid_type), intent(in) :: grid
  type(mapping_XGC_GEM_linear_coef_type), intent(in) :: coef

  field_gem=0.0
  field_gem_3d=0.0

  !if(MyId==0)write(*,*)'before mapping XGC-GEM',MPI_WTIME()
  if(MyId==0)then
    !mapping
    call mapping_XGC_GEM_core(field_xgc,field_gem_3d,grid,coef)
    !write(*,*)'after mapping core, field_pot_3d=', field_xgc(63388,1),'field_gem_3d=', field_gem_3d(imx/2,jmx/2,0)
  endif

  call mpi_barrier(mpi_comm_world,ierr)

  !if(MyId==0)write(*,*)'after mapping XGC GEM, begin GEM field 3d-2d',MPI_WTIME()
  !generate the field in GEM each z domain
  call generate_field_gem(field_gem,field_gem_3d)
  !if(Myid==0)then
    !write(*,*)'phi=',field_gem(imx/2,jmx/2,0)
  !endif
  call mpi_barrier(mpi_comm_world,ierr)
  if(MyId==0)write(*,*)'after GEM filed 3d-2d',MPI_WTIME()
end subroutine mapping_XGC_GEM
#else
!mapping from XGC to GEM
subroutine mapping_XGC_GEM_new(field_gem,grid,coef,nor)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi,&
  !  cce_first_field,cce_last_field,cce_all_field_node_number,cce_all_node_number
  use gem_com, only:MyId,imx,jmx,kmx,pi,pi2

  integer :: ierr
  !real, intent(inout) :: field_xgc(cce_all_node_number,nphi)
  real, intent(out) :: field_gem(0:imx,0:jmx,0:1)
  real, intent(in) :: nor
  !real, intent(out) :: field_gem_3d(0:imx,0:jmx,0:kmx)
  type(grid_type), intent(in) :: grid
  type(mapping_XGC_GEM_linear_coef_type), intent(in) :: coef

  field_gem=0.0
  density_3d=0.0

  !if(MyId==0)write(*,*)'before mapping XGC-GEM',MPI_WTIME()
  !mapping
  call mapping_XGC_GEM_core_new(grid,coef,nor)
  !write(*,*)'after mapping core, field_pot_3d=', field_xgc(63388,1),'field_gem_3d=', field_gem_3d(imx/2,jmx/2,0)

  call mpi_barrier(mpi_comm_world,ierr)

  !if(MyId==0)write(*,*)'after mapping XGC GEM, begin GEM field 3d-2d',MPI_WTIME()
  !generate the field in GEM each z domain
  call generate_field_gem(field_gem)
  !if(Myid==0)then
    !write(*,*)'phi=',field_gem(imx/2,jmx/2,0)
  !endif
  call mpi_barrier(mpi_comm_world,ierr)
  !if(MyId==0)write(*,*)'after GEM filed 3d-2d',MPI_WTIME()
end subroutine mapping_XGC_GEM_new
#endif

#ifndef __MAP_PARALLEL
!major part for mapping from XGC to GEM
subroutine mapping_XGC_GEM_core(field_xgc,field_gem_3d,grid,coef)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi,&
  !  cce_all_node_number
  use gem_com, only:MyId,imx,jmx,kmx,pi,pi2,lz
  use gem_equil, only:ntheta
#ifdef DEBUG_MAPPING
  use gem_com, only:timestep
#endif
  use gem_equil, only:delz,thfnz
  use EZspline_obj
  use EZspline

  real, intent(in) :: field_xgc(cce_all_node_number,nphi)
  real, intent(out) :: field_gem_3d(0:imx,0:jmx,0:kmx)
  type(grid_type), intent(in) :: grid
  type(mapping_XGC_GEM_linear_coef_type), intent(in) :: coef
  type(EZspline1_r8) :: spl_theta

  real, external :: interpol_spl_1d

  integer :: ierr,ii,i,j,k,kk,z1,z2,jl1,jl2,jr1,jr2
  integer :: surf_len_local,surf_idx_local(grid%surf_maxlen)
  real :: wz1,wz2,wl1,wl2,wr1,wr2,theta_start,z_tmp
  real :: theta_local(grid%surf_maxlen),theta_interpolation(grid%surf_maxlen+1),theta_z_gem(0:kmx)
  real :: field_local(grid%surf_maxlen,nphi),field_gem_z_xgc(0:jmx,grid%surf_maxlen),&
    field_gem_z_xgc_interpolation(0:jmx,grid%surf_maxlen+1),z_gem(0:kmx)
  
#ifdef DEBUG_MAPPING
  integer :: nstep,ndiag,stepid,err
  integer*8 :: buf_id,buf_size,total_size
  character(5) :: stepid_str
  character(512) :: cce_filename

  nstep=timestep-1
  ndiag=2
#endif
 
#ifdef DEBUG_MAPPING 
  if(modulo(nstep,ndiag) .eq. 0)then
    stepid=int(nstep/ndiag)+1
    write(stepid_str,'(I0.5)')stepid
    cce_filename='xgc_gem_1.3d.'//trim(stepid_str)//'.bp'
    call ADIOS_OPEN(buf_id,'xgc_gem_debug',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*5+cce_all_node_number*nphi*8+(imx+1)*(jmx+1)*(kmx+1)*8+100 !100 is buff
   
    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"nphi",nphi,err)
    call ADIOS_WRITE(buf_id,"cce_all_node_number",cce_all_node_number,err)
    call ADIOS_WRITE(buf_id,"imx_1",imx+1,err)
    call ADIOS_WRITE(buf_id,"jmx_1",jmx+1,err)
    call ADIOS_WRITE(buf_id,"kmx_1",kmx+1,err)
    !actrual data
    call ADIOS_WRITE(buf_id,"xgc_3d",field_xgc,err)
    call ADIOS_WRITE(buf_id,"gem_3d",field_gem_3d,err)
    
    call ADIOS_CLOSE(buf_id,err)
  endif
#endif

  if(myid==0)write(*,*)'point 0, XGC GEM,',MPI_WTIME()
  !generate each surface
#ifndef __MAP_PARALLEL  
  do ii=cce_first_surface,cce_last_surface
#else
  do ii=cce_surface_start,cce_surface_end
#endif
    surf_idx_local=grid%surf_idx(:,ii)
    surf_len_local=grid%surf_len(ii)

    !i is the idex in x direction
#ifndef __MAP_PARALLEL
    i=ii-cce_first_surface
#else
    i=ii-cce_surface_start+1
#endif
  if(myid==0 .and. ii==150)write(*,*)'point 0.1, XGC GEM,',MPI_WTIME()
    !local theta and field in XGC
    theta_local=0.0
    field_local=0.0
    do j=1,surf_len_local
      theta_local(j)=grid%theta_geo(surf_idx_local(j))
      do k=1,nphi
        field_local(j,k)=field_xgc(surf_idx_local(j),k)
      enddo
    enddo
  if(myid==0 .and. ii==150)write(*,*)'point 0.2, XGC GEM,',MPI_WTIME()
    !theta in z direction
    do k=0,kmx
      z_gem(k)=lz/real(kmx)*real(k)
      kk=int(z_gem(k)/delz)
      wz1=((kk+1)*delz-z_gem(k))/delz
      wz2=1-wz1
      theta_z_gem(k)=wz1*thfnz(kk)+wz2*thfnz(modulo(kk+1,ntheta))
    enddo
  if(myid==0 .and. ii==150)write(*,*)'point 0.3, XGC GEM,',MPI_WTIME()
    !generate (x,y,z_xgc) through linear interpolation
    field_gem_z_xgc=0.0
    do j=0,jmx
      do k=1,surf_len_local
        z1=coef%z1(i,j,k)
        z2=coef%z2(i,j,k)
        wz1=coef%wz1(i,j,k)
        wz2=coef%wz2(i,j,k)
        
        jl1=coef%jl1(i,j,k)
        jl2=coef%jl2(i,j,k)
        wl1=coef%wl1(i,j,k)
        wl2=coef%wl2(i,j,k)

        jr1=coef%jr1(i,j,k)
        jr2=coef%jr2(i,j,k)
        wr1=coef%wr1(i,j,k)
        wr2=coef%wr2(i,j,k)

        field_gem_z_xgc(j,k)=wz1*(wl1*field_local(jl1,z1)+wl2*field_local(jl2,z1))+wz2*(wr1*field_local(jr1,z2)+wr2*field_local(jr2,z2))        
      enddo
    enddo
  if(myid==0 .and. ii==150)write(*,*)'point 0.4, XGC GEM,',MPI_WTIME()
    !construct the \theta in xgc for cubic interpolation
    theta_start=theta_local(1)
    do j=1,surf_len_local
      theta_interpolation(j)=modulo(theta_local(j)-theta_start,pi2)
    enddo
    theta_interpolation(surf_len_local+1)=pi2
    !extend field_gem_z_xgc for cubic interpolation
    do j=1,surf_len_local
      field_gem_z_xgc_interpolation(:,j)=field_gem_z_xgc(:,j)
    enddo
    field_gem_z_xgc_interpolation(:,surf_len_local+1)=field_gem_z_xgc(:,1) 

  if(myid==0 .and. ii==150)write(*,*)'point 0.5, XGC GEM,',MPI_WTIME()
    !generate (x,y,z) through cubic interpolation
    do j=0,jmx
      call init_1d_interpolation(spl_theta,theta_interpolation(1:(surf_len_local+1)),field_gem_z_xgc_interpolation(j,1:(surf_len_local+1)),surf_len_local+1,.false.)
      do k=0,kmx
        z_tmp=modulo(theta_z_gem(k)-theta_start,pi2)
        field_gem_3d(i,j,k)=interpol_spl_1d(z_tmp,0,spl_theta)        
      enddo
      call finalize_1d_interpolation(spl_theta)
    enddo
  if(myid==0 .and. ii==150)write(*,*)'point 0.6, XGC GEM,',MPI_WTIME()
  enddo
  if(myid==0)write(*,*)'point 1, XGC GEM,',MPI_WTIME()
#ifdef DEBUG_MAPPING
  if(modulo(nstep,ndiag) .eq. 0)then
    stepid=int(nstep/ndiag)+1
    write(stepid_str,'(I0.5)')stepid
    cce_filename='xgc_gem_2.3d.'//trim(stepid_str)//'.bp'
    call ADIOS_OPEN(buf_id,'xgc_gem_debug',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*5+cce_all_node_number*nphi*8+(imx+1)*(jmx+1)*(kmx+1)*8+100 !100 is buff

    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"nphi",nphi,err)
    call ADIOS_WRITE(buf_id,"cce_all_node_number",cce_all_node_number,err)
    call ADIOS_WRITE(buf_id,"imx_1",imx+1,err)
    call ADIOS_WRITE(buf_id,"jmx_1",jmx+1,err)
    call ADIOS_WRITE(buf_id,"kmx_1",kmx+1,err)
    !actrual data
    call ADIOS_WRITE(buf_id,"xgc_3d",field_xgc,err)
    call ADIOS_WRITE(buf_id,"gem_3d",field_gem_3d,err)

    call ADIOS_CLOSE(buf_id,err)
  endif
#endif

end subroutine mapping_XGC_GEM_core

#else
!major part for mapping from XGC to GEM
subroutine mapping_XGC_GEM_core_new(grid,coef,nor)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi,cce_all_node_number
  use gem_com, only:MyId,imx,jmx,kmx,pi,pi2,lz
  use gem_equil, only:delz,thfnz,ntheta
  use EZspline_obj
  use EZspline

  !real, intent(inout) :: field_xgc(cce_all_node_number,nphi)
  !real, intent(out) :: field_gem_3d(0:imx,0:jmx,0:kmx)
  real, intent(in) :: nor
  type(grid_type), intent(in) :: grid
  type(mapping_XGC_GEM_linear_coef_type), intent(in) :: coef
  type(EZspline1_r8) :: spl_theta

  real, external :: interpol_spl_1d

  integer :: ierr,i,ii,iii,icount,j,k,kk,z1,z2,jl1,jl2,jr1,jr2
  integer :: surf_len_local,surf_idx_local(grid%surf_maxlen)
  real :: wz1,wz2,wl1,wl2,wr1,wr2,theta_start,z_tmp
  real :: theta_local(grid%surf_maxlen),theta_interpolation(grid%surf_maxlen+1),theta_z_gem(0:kmx)
  real :: field_local(grid%surf_maxlen,nphi),field_gem_z_xgc(0:jmx,grid%surf_maxlen),&
    field_gem_z_xgc_interpolation(0:jmx,grid%surf_maxlen+1),z_gem(0:kmx),field_gem_3d_tmp(0:imx,0:jmx,0:kmx)
 
  !field_gem_3d=0.0
  density_3d=0.0
  icount=cce_all_node_number*nphi
  call mpi_bcast(pot_density_3d,icount,MPI_REAL8,0,mpi_comm_world,ierr)
  call mpi_barrier(mpi_comm_world,ierr) 
  if(myid==0)then
    write(111,*)'point 0, XGC GEM,',MPI_WTIME()
    call flush(111)
  endif
  if(gx_id==0)then
    !generate each surface
    do ii=cce_surface_start,cce_surface_end
      surf_idx_local=grid%surf_idx(:,ii)
      surf_len_local=grid%surf_len(ii)

      !i is the idex in x direction
      i=ii-cce_surface_start
      iii=ii-cce_first_surface
      if(ii==150)write(*,*)'point 0.1, XGC GEM,',MPI_WTIME()
      !local theta and field in XGC
      theta_local=0.0
      field_local=0.0
      do j=1,surf_len_local
        theta_local(j)=grid%theta_geo(surf_idx_local(j))
        do k=1,nphi
          field_local(j,k)=pot_density_3d(surf_idx_local(j),k)
        enddo
      enddo
      if(ii==150)write(*,*)'point 0.2, XGC GEM,',MPI_WTIME()
      !theta in z direction
      do k=0,kmx
        z_gem(k)=lz/real(kmx)*real(k)
        kk=int(z_gem(k)/delz)
        wz1=((kk+1)*delz-z_gem(k))/delz
        wz2=1-wz1
        theta_z_gem(k)=wz1*thfnz(kk)+wz2*thfnz(modulo(kk+1,ntheta))
      enddo
      if(ii==150)write(*,*)'point 0.3, XGC GEM,',MPI_WTIME()
      !generate (x,y,z_xgc) through linear interpolation
      field_gem_z_xgc=0.0
      do j=0,jmx
        !!$omp parallel do
        do k=1,surf_len_local
          z1=coef%z1(i+1,j,k)
          z2=coef%z2(i+1,j,k)
          wz1=coef%wz1(i+1,j,k)
          wz2=coef%wz2(i+1,j,k)
          
          jl1=coef%jl1(i+1,j,k)
          jl2=coef%jl2(i+1,j,k)
          wl1=coef%wl1(i+1,j,k)
          wl2=coef%wl2(i+1,j,k)

          jr1=coef%jr1(i+1,j,k)
          jr2=coef%jr2(i+1,j,k)
          wr1=coef%wr1(i+1,j,k)
          wr2=coef%wr2(i+1,j,k)
          !if( j==jmx/2 .and. k==surf_len_local/2)write(*,*)'ii,z1,z2,wz1,wz2,jl1,jl2,wl1,wl2,jr1,jr2,wr1,wr2,field_local',ii,z1,z2,wz1,wz2,jl1,jl2,wl1,wl2,jr1,jr2,wr1,wr2,field_local(jl1,z1)
          field_gem_z_xgc(j,k)=wz1*(wl1*field_local(jl1,z1)+wl2*field_local(jl2,z1))+wz2*(wr1*field_local(jr1,z2)+wr2*field_local(jr2,z2))        
        enddo
      enddo
    if(ii==150)write(*,*)'point 0.4, XGC GEM,',MPI_WTIME()
      !construct the \theta in xgc for cubic interpolation
      theta_start=theta_local(1)
      do j=1,surf_len_local
        theta_interpolation(j)=modulo(theta_local(j)-theta_start,pi2)
      enddo
      theta_interpolation(surf_len_local+1)=pi2
      !extend field_gem_z_xgc for cubic interpolation
      !!$omp parallel do
      do j=1,surf_len_local
        field_gem_z_xgc_interpolation(:,j)=field_gem_z_xgc(:,j)
      enddo
      field_gem_z_xgc_interpolation(:,surf_len_local+1)=field_gem_z_xgc(:,1) 

    if(ii==150)write(*,*)'point 0.5, XGC GEM,',MPI_WTIME()
      !generate (x,y,z) through cubic interpolation
      do j=0,jmx
        call init_1d_interpolation(spl_theta,theta_interpolation(1:(surf_len_local+1)),field_gem_z_xgc_interpolation(j,1:(surf_len_local+1)),surf_len_local+1,.false.)
        do k=0,kmx
          z_tmp=modulo(theta_z_gem(k)-theta_start,pi2)
          density_3d(iii,j,k)=interpol_spl_1d(z_tmp,0,spl_theta)*nor        
        enddo
        call finalize_1d_interpolation(spl_theta)
      enddo
    if(ii==150)write(*,*)'point 0.6, XGC GEM,',MPI_WTIME()
    enddo
  endif

  icount=(imx+1)*(jmx+1)*(kmx+1)
  call mpi_allreduce(MPI_IN_PLACE,density_3d,icount,MPI_REAL8,mpi_sum,x_comm,ierr)
  call mpi_barrier(x_comm,ierr)
  
  !if(myid==0)field_gem_3d=field_gem_3d_tmp
  !if(myid==0)write(*,*)'point 1, XGC GEM,',MPI_WTIME()

end subroutine mapping_XGC_GEM_core_new

#endif

!release the GEM 3d data to each domain
subroutine generate_field_gem(field_gem)
  use gem_com, only:imx,jmx,kmx,gclr

  real, intent(out) :: field_gem(0:imx,0:jmx,0:1)
  !real, intent(inout) :: field_gem_3d(0:imx,0:jmx,0:kmx)

  integer :: ierr

  !broadcast the field_gem_3d
  call mpi_bcast(density_3d, (imx+1)*(jmx+1)*(kmx+1), MPI_REAL8, 0, mpi_comm_world, ierr)
  call mpi_barrier(mpi_comm_world,ierr)

  field_gem(:,:,0:1)=density_3d(:,:,gclr:(gclr+1))

end subroutine generate_field_gem

#ifndef __MAP_PARALLEL
!generate the coefficient for XGC to GEM linear interpolation
subroutine setup_XGC_GEM_coef(coef,grid)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi,nwedge
  use gem_equil, only: sf,r0,q0,q0abs,dth,ntheta,thflx
  use gem_com, only: imx,jmx,dy,ly,pi,pi2,MyId
  use EZspline_obj
  use EZspline

  type(mapping_XGC_GEM_linear_coef_type), intent(out):: coef
  type(grid_type), intent(in) :: grid

  type(EZspline1_r8) :: spl_thflx,spl_thflx_inv
  real, external :: interpol_spl_1d

  integer :: ii,i,j,k,z1,z2,jl1,jl2,jr1,jr2,max_len,surf_len_local,surf_idx_local(grid%surf_maxlen)
  real :: q_local,y_tmp,wz1,wz2,wl1,wl2,wr1,wr2,zeta_tmp,zeta_tmp_left,zeta_tmp_right,&
    theta_tmp,theta_left_tmp,theta_right_tmp,th_tmp,theta_grid(0:ntheta),theta_local(grid%surf_maxlen),theta_flux(grid%surf_maxlen)

  if(myid==0)write(*,*)'point 0, XGC GEM,',MPI_WTIME()
  !max length in XGC theta direction 
  max_len=grid%surf_maxlen
  !claim the memory for coef_XGC_GEM_linear_interpolation
  allocate(coef%z1(0:imx,0:jmx,max_len),coef%z2(0:imx,0:jmx,max_len),coef%wz1(0:imx,0:jmx,max_len),&
    coef%wz2(0:imx,0:jmx,max_len),coef%jl1(0:imx,0:jmx,max_len),coef%jl2(0:imx,0:jmx,max_len),&
    coef%wl1(0:imx,0:jmx,max_len),coef%wl2(0:imx,0:jmx,max_len),coef%jr1(0:imx,0:jmx,max_len),&
    coef%jr2(0:imx,0:jmx,max_len),coef%wr1(0:imx,0:jmx,max_len),coef%wr2(0:imx,0:jmx,max_len))

  coef%z1=0
  coef%z2=0
  coef%jl1=0
  coef%jl2=0
  coef%jr1=0
  coef%jr2=0
  coef%wz1=0.0
  coef%wz2=0.0
  coef%wl1=0.0
  coef%wl2=0.0
  coef%wr1=0.0
  coef%wr2=0.0

  if(myid==0)write(*,*)'point 1, XGC GEM,',MPI_WTIME()
  theta_grid(ntheta/2)=0.
  do i=ntheta/2+1,ntheta
     th_tmp=real(i-ntheta/2)*dth
     theta_grid(i)=th_tmp
     theta_grid(ntheta-i)=-th_tmp
  enddo

  if(myid==0)write(*,*)'point 2, XGC GEM,',MPI_WTIME()
  !if(MyId==0)then
  !  write(*,*)theta_grid
  !endif
  !generate the coefficiet
  do ii=cce_first_surface,cce_last_surface
    
    if(myid==0 .and. ii==250)write(*,*)'point 2.0, XGC GEM,',MPI_WTIME() 
    surf_idx_local=grid%surf_idx(:,ii)
    surf_len_local=grid%surf_len(ii)
    
    !x direction start from 0 to imx
    i=ii-cce_first_surface
    !local safety factor
    q_local=sf(ii-cce_first_surface)
 
    !write(*,*)'point1, XGC GEM','ii=',ii
    !if(ii==31)then
    !  write(*,*)'thflx=',thflx(i,:)
    !endif
    call init_1d_interpolation(spl_thflx_inv,thflx(i,:),theta_grid,ntheta+1,.false.)
    call init_1d_interpolation(spl_thflx,theta_grid,thflx(i,:),ntheta+1,.false.)
    !write(*,*)'point2, XGC GEM','ii=',ii
    if(myid==0 .and. ii==250)write(*,*)'point 2.1, XGC GEM,',MPI_WTIME()
    !local theta in XGC 
    theta_local=0.0
    do j=1,surf_len_local
      theta_local(j)=grid%theta_geo(surf_idx_local(j))
      theta_flux(j)=interpol_spl_1d(theta_local(j),0,spl_thflx)
    enddo
    !write(*,*)'point3, XGC GEM', 'ii=',ii
    call finalize_1d_interpolation(spl_thflx)
    if(myid==0 .and. ii==250)write(*,*)'point 2.2, XGC GEM,',MPI_WTIME() 
    do j=1,surf_len_local
      do k=0,jmx
        if(myid==0 .and. ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.0, XGC GEM,',MPI_WTIME()
        !the y in GEM
        y_tmp=real(k)*dy

        !zeta_tmp related to y
        zeta_tmp=modulo(q_local*theta_flux(j)-y_tmp/r0*q0,pi2/real(nwedge))
        !calc the idx and coefficient in zeta direction
        call search_zeta(z1,z2,wz1,wz2,pi2/real(nwedge)/real(nphi),pi2/real(nwedge),nphi,zeta_tmp)
        if(myid==0 .and. ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.1, XGC GEM,',MPI_WTIME()
        !left theta
        zeta_tmp_left=real(z1-1)*pi2/real(nphi)/real(nwedge)
        theta_left_tmp=(q_local*theta_flux(j)-zeta_tmp+zeta_tmp_left)/q_local
        theta_tmp=modulo(theta_left_tmp-thflx(i,0),pi2)+thflx(i,0)
        !theta_tmp=interpol_spl_1d(theta_left_tmp,0,spl_thflx_inv)
        call search_theta(jl1,jl2,wl1,wl2,theta_flux(1:surf_len_local),surf_len_local,theta_tmp)
        if(myid==0 .and. ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.2, XGC GEM,',MPI_WTIME()
        !right theta
        zeta_tmp_right=real(z2-1)*pi2/real(nphi)/real(nwedge)
        if(z2==1)then
          zeta_tmp_right=pi2/real(nwedge)
        endif
        theta_right_tmp=(q_local*theta_flux(j)-zeta_tmp+zeta_tmp_right)/q_local
        theta_tmp=modulo(theta_right_tmp-thflx(i,0),pi2)+thflx(i,0)
        !theta_tmp=interpol_spl_1d(theta_right_tmp,0,spl_thflx_inv)
        call search_theta(jr1,jr2,wr1,wr2,theta_flux(1:surf_len_local),surf_len_local,theta_tmp)
        if(myid==0 .and. ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.3, XGC GEM,',MPI_WTIME()
        coef%z1(i,k,j)=z1
        coef%z2(i,k,j)=z2
        coef%wz1(i,k,j)=wz1
        coef%wz2(i,k,j)=wz2
        coef%jl1(i,k,j)=jl1
        coef%jl2(i,k,j)=jl2
        coef%wl1(i,k,j)=wl1
        coef%wl2(i,k,j)=wl2
        coef%jr1(i,k,j)=jr1
        coef%jr2(i,k,j)=jr2
        coef%wr1(i,k,j)=wr1
        coef%wr2(i,k,j)=wr2
        if(myid==0 .and. ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.4, XGC GEM,',MPI_WTIME()
      enddo
      
      !call finalize_1d_interpolation(spl_thflx_inv) 
    enddo
    call finalize_1d_interpolation(spl_thflx_inv)
    !write(*,*)'point4, XGC GEM', 'ii=',ii
    if(myid==0 .and. ii==250)write(*,*)'point 2.3, XGC GEM,',MPI_WTIME()
  enddo

  if(myid==0)write(*,*)'point 3, XGC GEM,',MPI_WTIME()
end subroutine setup_XGC_GEM_coef

#else
!generate the coefficient for XGC to GEM linear interpolation
subroutine setup_XGC_GEM_coef_new(coef,grid)
  !use coupling_core_edge, only:cce_first_surface,cce_last_surface,cce_all_surface_node_number,nphi,nwedge
  use gem_equil, only: sf,r0,q0,q0abs,dth,ntheta,thflx
  use gem_com, only: imx,jmx,dy,ly,pi,pi2,MyId
  use EZspline_obj
  use EZspline

  type(mapping_XGC_GEM_linear_coef_type), intent(out):: coef
  type(grid_type), intent(in) :: grid

  type(EZspline1_r8) :: spl_thflx,spl_thflx_inv
  real, external :: interpol_spl_1d

  integer :: i,ii,iii,j,k,icount,ierr,z1,z2,jl1,jl2,jr1,jr2,max_len,x_length,surf_len_local,&
    surf_idx_local(grid%surf_maxlen)
  real :: q_local,y_tmp,wz1,wz2,wl1,wl2,wr1,wr2,zeta_tmp,zeta_tmp_left,zeta_tmp_right,&
    theta_tmp,theta_left_tmp,theta_right_tmp,th_tmp,theta_grid(0:ntheta),theta_local(grid%surf_maxlen),theta_flux(grid%surf_maxlen)

  if(myid==0)write(*,*)'point 0, XGC GEM,',MPI_WTIME()
  !max length in XGC theta direction 
  max_len=grid%surf_maxlen
  x_length=cce_surface_end-cce_surface_start+1
  !claim the memory for coef_XGC_GEM_linear_interpolation
  allocate(coef%z1(x_length,0:jmx,max_len),coef%z2(x_length,0:jmx,max_len),coef%wz1(x_length,0:jmx,max_len),&
    coef%wz2(x_length,0:jmx,max_len),coef%jl1(x_length,0:jmx,max_len),coef%jl2(x_length,0:jmx,max_len),&
    coef%wl1(x_length,0:jmx,max_len),coef%wl2(x_length,0:jmx,max_len),coef%jr1(x_length,0:jmx,max_len),&
    coef%jr2(x_length,0:jmx,max_len),coef%wr1(x_length,0:jmx,max_len),coef%wr2(x_length,0:jmx,max_len))

  coef%z1=0
  coef%z2=0
  coef%jl1=0
  coef%jl2=0
  coef%jr1=0
  coef%jr2=0
  coef%wz1=0.0
  coef%wz2=0.0
  coef%wl1=0.0
  coef%wl2=0.0
  coef%wr1=0.0
  coef%wr2=0.0

  if(gx_id==0)then
    if(myid==0)write(*,*)'point 1, XGC GEM,',MPI_WTIME()
    theta_grid(ntheta/2)=0.
    do i=ntheta/2+1,ntheta
      th_tmp=real(i-ntheta/2)*dth
      theta_grid(i)=th_tmp
      theta_grid(ntheta-i)=-th_tmp
    enddo

    if(myid==0)write(*,*)'point 2, XGC GEM,',MPI_WTIME()
    !if(MyId==0)then
    !  write(*,*)theta_grid
    !endif
    !generate the coefficiet
    do ii=cce_surface_start,cce_surface_end

      if(ii==250)write(*,*)'point 2.0, XGC GEM,',MPI_WTIME()
      surf_idx_local=grid%surf_idx(:,ii)
      surf_len_local=grid%surf_len(ii)

      !x direction start from 0 to imx
      i=ii-cce_first_surface
      iii=ii-cce_surface_start+1
      !local safety factor
      q_local=sf(ii-cce_first_surface)

      !write(*,*)'point1, XGC GEM','ii=',ii
      !if(ii==31)then
      !  write(*,*)'thflx=',thflx(i,:)
      !endif
      call init_1d_interpolation(spl_thflx_inv,thflx(i,:),theta_grid,ntheta+1,.false.)
      call init_1d_interpolation(spl_thflx,theta_grid,thflx(i,:),ntheta+1,.false.)
      !write(*,*)'point2, XGC GEM','ii=',ii
      if(ii==250)write(*,*)'point 2.1, XGC GEM,',MPI_WTIME()
      !local theta in XGC 
      theta_local=0.0
      do j=1,surf_len_local
        theta_local(j)=grid%theta_geo(surf_idx_local(j))
        theta_flux(j)=interpol_spl_1d(theta_local(j),0,spl_thflx)
      enddo
      !write(*,*)'point3, XGC GEM', 'ii=',ii
      call finalize_1d_interpolation(spl_thflx)
      if(ii==250)write(*,*)'point 2.2, XGC GEM,',MPI_WTIME()
      do j=1,surf_len_local
        do k=0,jmx
          if(ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.0, XGC GEM,',MPI_WTIME()
          !the y in GEM
          y_tmp=real(k)*dy

          !zeta_tmp related to y
          zeta_tmp=modulo(q_local*theta_flux(j)-y_tmp/r0*q0,pi2/real(nwedge))
          !calc the idx and coefficient in zeta direction
          call search_zeta(z1,z2,wz1,wz2,pi2/real(nwedge)/real(nphi),pi2/real(nwedge),nphi,zeta_tmp)
          if(ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.1, XGC GEM,',MPI_WTIME()
          !left theta
          zeta_tmp_left=real(z1-1)*pi2/real(nphi)/real(nwedge)
          theta_left_tmp=(q_local*theta_flux(j)-zeta_tmp+zeta_tmp_left)/q_local
          theta_tmp=modulo(theta_left_tmp-thflx(i,0),pi2)+thflx(i,0)
          !theta_tmp=interpol_spl_1d(theta_left_tmp,0,spl_thflx_inv)
          call search_theta(jl1,jl2,wl1,wl2,theta_flux(1:surf_len_local),surf_len_local,theta_tmp)
          if(ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.2, XGC GEM,',MPI_WTIME()
          !right theta
          zeta_tmp_right=real(z2-1)*pi2/real(nphi)/real(nwedge)
          if(z2==1)then
            zeta_tmp_right=pi2/real(nwedge)
          endif
          theta_right_tmp=(q_local*theta_flux(j)-zeta_tmp+zeta_tmp_right)/q_local
          theta_tmp=modulo(theta_right_tmp-thflx(i,0),pi2)+thflx(i,0)
          !theta_tmp=interpol_spl_1d(theta_right_tmp,0,spl_thflx_inv)
          call search_theta(jr1,jr2,wr1,wr2,theta_flux(1:surf_len_local),surf_len_local,theta_tmp)
          if(ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.3, XGC GEM,',MPI_WTIME()
          coef%z1(iii,k,j)=z1
          coef%z2(iii,k,j)=z2
          coef%wz1(iii,k,j)=wz1
          coef%wz2(iii,k,j)=wz2
          coef%jl1(iii,k,j)=jl1
          coef%jl2(iii,k,j)=jl2
          coef%wl1(iii,k,j)=wl1
          coef%wl2(iii,k,j)=wl2
          coef%jr1(iii,k,j)=jr1
          coef%jr2(iii,k,j)=jr2
          coef%wr1(iii,k,j)=wr1
          coef%wr2(iii,k,j)=wr2
          if(ii==250 .and. j==1 .and. k==0)write(*,*)'point 2.2.4, XGC GEM,',MPI_WTIME()
        enddo
        !call finalize_1d_interpolation(spl_thflx_inv) 
      enddo
      call finalize_1d_interpolation(spl_thflx_inv)
      if(ii==250)write(*,*)'point 2.3, XGC GEM,',MPI_WTIME()
    enddo
  endif
  call mpi_barrier(mpi_comm_world,ierr)

  if(myid==0)write(*,*)'point 3, XGC GEM,',MPI_WTIME() 
  
end subroutine setup_XGC_GEM_coef_new
#endif

!generate the coefficient in zeta direction for each points
subroutine search_zeta(j1,j2,w1,w2,dlength,length,nlength,tmp)
  integer, intent(out) :: j1,j2
  real, intent(out) :: w1,w2
  integer, intent(in) :: nlength
  real :: dlength,length,tmp

  if(tmp<0. .or. tmp>=length)then
    tmp=modulo(tmp,length)
  endif

  j1=floor(tmp/dlength)+1
  j2=j1+1
  if(j1==nlength)then
    j2=1
  endif

  w2=(tmp-real(j1-1)*dlength)/dlength
  w1=1.0-w2

end subroutine search_zeta

!generate the coefficient in theta direction for each points
subroutine search_theta(j1,j2,w1,w2,theta_local,surf_len,tmp)
  use gem_com, only:pi2

  integer, intent(out) :: j1,j2
  real, intent(out) :: w1,w2
  integer, intent(in) :: surf_len
  real, intent(in) :: theta_local(surf_len),tmp

  integer :: i
  real :: theta_start,tmp1,theta_local_tmp(surf_len+1)

  theta_start=theta_local(1)

  !construct the theta from 0 to pi2
  do i=1,surf_len
    theta_local_tmp(i)=modulo(theta_local(i)-theta_start,pi2)
  enddo

  theta_local_tmp(surf_len+1)=pi2
  
  !change tmp in the range (0,pi2)
  tmp1=modulo(tmp-theta_start,pi2)
  
  do i=1,surf_len
    if( tmp1>=theta_local_tmp(i) .and. tmp1<theta_local_tmp(i+1))then
      j1=i
      j2=i+1
      w2=(tmp1-theta_local_tmp(i))/(theta_local_tmp(i+1)-theta_local_tmp(i))
      w1=1.0-w2
      if(j1==surf_len)then
        j2=1
      endif
    endif
  enddo

end subroutine search_theta

#ifndef __MAP_PARALLEL
!store the near coef GEM-XGC XGC-GEM
subroutine store_coef(coef_GEM_XGC,coef_XGC_GEM,grid)
  !use coupling_core_edge, only:cce_all_surface_node_number,nphi
  use gem_com, only:MyId,imx,jmx
  use gem_equil, only:nr,ntheta,thflx,qhat,yfn,radius,hght,thfnz,zfnth,sf
  integer :: i,j,max_len,err
  integer*8 :: buf_id,buf_size,total_size

  character(512) :: cce_filename

  type(mapping_GEM_XGC_linear_coef_type), intent(in) :: coef_GEM_XGC
  type(mapping_XGC_GEM_linear_coef_type), intent(in) :: coef_XGC_GEM
  type(grid_type), intent(in) :: grid

  max_len=grid%surf_maxlen
 

  if(MyId==0)then
    !output the thflx
#ifdef EQUIL_OUTPUT
    cce_filename='equil.bp'
    call ADIOS_OPEN(buf_id,'equil',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*2+40*(nr+1)*(ntheta+1)+16*(ntheta+1)+100 !last 100 is buffer

    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"nr_1",nr+1,err)
    call ADIOS_WRITE(buf_id,"ntheta_1",ntheta+1,err)
    !actual data
    call ADIOS_WRITE(buf_id,"thflx",thflx,err)
    call ADIOS_WRITE(buf_id,"qhat",qhat,err)
    call ADIOS_WRITE(buf_id,"yfn",yfn,err)
    call ADIOS_WRITE(buf_id,"sf",sf,err)
    call ADIOS_WRITE(buf_id,"radius",radius,err)
    call ADIOS_WRITE(buf_id,"hght",hght,err)
    call ADIOS_WRITE(buf_id,"zfnth",zfnth,err)
    call ADIOS_WRITE(buf_id,"thfnz",thfnz,err)

    call ADIOS_CLOSE(buf_id,err) 
#endif

    cce_filename='coef_GEM_XGC.bp'
    call ADIOS_OPEN(buf_id,'coef_GEM_XGC',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*2+24*nphi*cce_all_surface_node_number+100 !last 100 is buffer

    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"nphi",nphi,err)
    call ADIOS_WRITE(buf_id,"cce_all_surface_node_number",cce_all_surface_node_number,err)
    !actual data
    call ADIOS_WRITE(buf_id,"coef_GEM_XGC%jy1",coef_GEM_XGC%jy1,err)
    call ADIOS_WRITE(buf_id,"coef_GEM_XGC%jy2",coef_GEM_XGC%jy2,err)
    call ADIOS_WRITE(buf_id,"coef_GEM_XGC%wy1",coef_GEM_XGC%wy1,err)
    call ADIOS_WRITE(buf_id,"coef_GEM_XGC%wy2",coef_GEM_XGC%wy2,err)

    call ADIOS_CLOSE(buf_id,err)

    cce_filename='coef_XGC_GEM.bp'
    call ADIOS_OPEN(buf_id,'coef_XGC_GEM',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*3+(imx+1)*(jmx+1)*max_len*72+100 !100 is buff

    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"imx_1",imx+1,err)
    call ADIOS_WRITE(buf_id,"jmx_1",jmx+1,err)
    call ADIOS_WRITE(buf_id,"max_len",max_len,err)
    !actual data
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%z1",coef_XGC_GEM%z1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%z2",coef_XGC_GEM%z2,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%jl1",coef_XGC_GEM%jl1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%jl2",coef_XGC_GEM%jl2,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%jr1",coef_XGC_GEM%jr1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%jr2",coef_XGC_GEM%jr2,err)

    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wz1",coef_XGC_GEM%wz1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wz2",coef_XGC_GEM%wz2,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wl1",coef_XGC_GEM%wl1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wl2",coef_XGC_GEM%wl2,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wr1",coef_XGC_GEM%wr1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wr2",coef_XGC_GEM%wr2,err)

    call ADIOS_CLOSE(buf_id,err)

  endif

end subroutine store_coef

!!!#else

!store the near coef GEM-XGC XGC-GEM
subroutine store_coef_new(coef_GEM_XGC,coef_XGC_GEM,grid)
  !use coupling_core_edge, only:cce_all_surface_node_number,nphi
  use gem_com, only:MyId,imx,jmx
  use gem_equil, only:nr,ntheta,thflx,qhat,yfn,radius,hght,thfnz,zfnth,sf
  integer :: i,j,max_len,err
  integer*8 :: buf_id,buf_size,total_size

  character(512) :: cce_filename

  type(mapping_GEM_XGC_linear_coef_type), intent(in) :: coef_GEM_XGC
  type(mapping_XGC_GEM_linear_coef_type), intent(in) :: coef_XGC_GEM
  type(grid_type), intent(in) :: grid

  max_len=grid%surf_maxlen
 

  if(MyId==0)then
    !output the thflx
#ifdef EQUIL_OUTPUT
    cce_filename='equil.bp'
    call ADIOS_OPEN(buf_id,'equil',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*2+40*(nr+1)*(ntheta+1)+16*(ntheta+1)+100 !last 100 is buffer

    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"nr_1",nr+1,err)
    call ADIOS_WRITE(buf_id,"ntheta_1",ntheta+1,err)
    !actual data
    call ADIOS_WRITE(buf_id,"thflx",thflx,err)
    call ADIOS_WRITE(buf_id,"qhat",qhat,err)
    call ADIOS_WRITE(buf_id,"yfn",yfn,err)
    call ADIOS_WRITE(buf_id,"sf",sf,err)
    call ADIOS_WRITE(buf_id,"radius",radius,err)
    call ADIOS_WRITE(buf_id,"hght",hght,err)
    call ADIOS_WRITE(buf_id,"zfnth",zfnth,err)
    call ADIOS_WRITE(buf_id,"thfnz",thfnz,err)

    call ADIOS_CLOSE(buf_id,err) 
#endif

    !to do real parallize later
    cce_filename='coef_GEM_XGC.bp'
    call ADIOS_OPEN(buf_id,'coef_GEM_XGC',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*2+24*nphi*cce_all_surface_node_number+100 !last 100 is buffer

    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"nphi",nphi,err)
    call ADIOS_WRITE(buf_id,"cce_all_surface_node_number",cce_all_surface_node_number,err)
    !actual data
    call ADIOS_WRITE(buf_id,"coef_GEM_XGC%jy1",coef_GEM_XGC%jy1,err)
    call ADIOS_WRITE(buf_id,"coef_GEM_XGC%jy2",coef_GEM_XGC%jy2,err)
    call ADIOS_WRITE(buf_id,"coef_GEM_XGC%wy1",coef_GEM_XGC%wy1,err)
    call ADIOS_WRITE(buf_id,"coef_GEM_XGC%wy2",coef_GEM_XGC%wy2,err)

    call ADIOS_CLOSE(buf_id,err)

  endif

  if(gx_id==0)then
    cce_filename='coef_XGC_GEM.bp'
    call ADIOS_OPEN(buf_id,'coef_XGC_GEM',cce_filename,'w',MPI_COMM_SELF,err)
    buf_size=4*3+(imx+1)*(jmx+1)*max_len*72+100 !100 is buff

    call ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    call ADIOS_WRITE(buf_id,"imx_1",imx+1,err)
    call ADIOS_WRITE(buf_id,"jmx_1",jmx+1,err)
    call ADIOS_WRITE(buf_id,"max_len",max_len,err)
    !actual data
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%z1",coef_XGC_GEM%z1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%z2",coef_XGC_GEM%z2,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%jl1",coef_XGC_GEM%jl1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%jl2",coef_XGC_GEM%jl2,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%jr1",coef_XGC_GEM%jr1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%jr2",coef_XGC_GEM%jr2,err)

    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wz1",coef_XGC_GEM%wz1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wz2",coef_XGC_GEM%wz2,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wl1",coef_XGC_GEM%wl1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wl2",coef_XGC_GEM%wl2,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wr1",coef_XGC_GEM%wr1,err)
    call ADIOS_WRITE(buf_id,"coef_XGC_GEM%wr2",coef_XGC_GEM%wr2,err)

    call ADIOS_CLOSE(buf_id,err)

  endif

end subroutine store_coef_new

subroutine restore_coef(coef_GEM_XGC,coef_XGC_GEM,grid)
  use ADIOS_READ_MOD
  use mpi
  use gem_com,only:MyId,imx,jmx
  !use coupling_core_edge,only:nphi,cce_all_surface_node_number

  integer :: err,max_len
  integer*8 :: buf_id,sel1=0
  logical :: ex
  character(512) :: cce_filename
  
  type(mapping_GEM_XGC_linear_coef_type), intent(out) :: coef_GEM_XGC
  type(mapping_XGC_GEM_linear_coef_type), intent(out) :: coef_XGC_GEM
  type(grid_type), intent(in) :: grid

  if(MyId==0)then
    !claim the memory for coef_GEM_XGC_linear_interpolation
    allocate(coef_GEM_XGC%jy1(cce_all_surface_node_number,nphi),coef_GEM_XGC%jy2(cce_all_surface_node_number,nphi),coef_GEM_XGC%wy1(cce_all_surface_node_number,nphi),coef_GEM_XGC%wy2(cce_all_surface_node_number,nphi))

    coef_GEM_XGC%jy1=0
    coef_GEM_XGC%jy2=0
    coef_GEM_XGC%wy1=0.0
    coef_GEM_XGC%wy2=0.0

    cce_filename='coef_GEM_XGC.bp'
    call ADIOS_READ_OPEN_FILE(buf_id,cce_filename,0,MPI_COMM_SELF,err)
    if(err/=0) then
      write(*,*)'coupling receive error: could not open file', cce_filename
    endif 

    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_GEM_XGC%jy1",0,1,coef_GEM_XGC%jy1,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_GEM_XGC%jy2",0,1,coef_GEM_XGC%jy2,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_GEM_XGC%wy1",0,1,coef_GEM_XGC%wy1,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_GEM_XGC%wy2",0,1,coef_GEM_XGC%wy2,err) 
    call ADIOS_PERFORM_READS(buf_id,err)
    call ADIOS_READ_CLOSE(buf_id,err)
    write(*,*)'read coef_GEM_XGC completed'

   
    !max length in XGC theta direction 
    max_len=grid%surf_maxlen
    !claim the memory for coef_XGC_GEM_linear_interpolation
    allocate(coef_XGC_GEM%z1(0:imx,0:jmx,max_len),coef_XGC_GEM%z2(0:imx,0:jmx,max_len),coef_XGC_GEM%wz1(0:imx,0:jmx,max_len),&
      coef_XGC_GEM%wz2(0:imx,0:jmx,max_len),coef_XGC_GEM%jl1(0:imx,0:jmx,max_len),coef_XGC_GEM%jl2(0:imx,0:jmx,max_len),&
      coef_XGC_GEM%wl1(0:imx,0:jmx,max_len),coef_XGC_GEM%wl2(0:imx,0:jmx,max_len),coef_XGC_GEM%jr1(0:imx,0:jmx,max_len),&
      coef_XGC_GEM%jr2(0:imx,0:jmx,max_len),coef_XGC_GEM%wr1(0:imx,0:jmx,max_len),coef_XGC_GEM%wr2(0:imx,0:jmx,max_len))

    coef_XGC_GEM%z1=0
    coef_XGC_GEM%z2=0
    coef_XGC_GEM%jl1=0
    coef_XGC_GEM%jl2=0
    coef_XGC_GEM%jr1=0
    coef_XGC_GEM%jr2=0
    coef_XGC_GEM%wz1=0.0
    coef_XGC_GEM%wz2=0.0
    coef_XGC_GEM%wl1=0.0
    coef_XGC_GEM%wl2=0.0
    coef_XGC_GEM%wr1=0.0
    coef_XGC_GEM%wr2=0.0
 
    cce_filename='coef_XGC_GEM.bp'
    call ADIOS_READ_OPEN_FILE(buf_id,cce_filename,0,MPI_COMM_SELF,err)
    if(err/=0) then
      write(*,*)'coupling receive error: could not open file', cce_filename
    endif
    
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%z1",0,1,coef_XGC_GEM%z1,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%z2",0,1,coef_XGC_GEM%z2,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%jl1",0,1,coef_XGC_GEM%jl1,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%jl2",0,1,coef_XGC_GEM%jl2,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%jr1",0,1,coef_XGC_GEM%jr1,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%jr2",0,1,coef_XGC_GEM%jr2,err)

    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%wz1",0,1,coef_XGC_GEM%wz1,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%wz2",0,1,coef_XGC_GEM%wz2,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%wl1",0,1,coef_XGC_GEM%wl1,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%wl2",0,1,coef_XGC_GEM%wl2,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%wr1",0,1,coef_XGC_GEM%wr1,err)
    call ADIOS_SCHEDULE_READ(buf_id,sel1,"coef_XGC_GEM%wr2",0,1,coef_XGC_GEM%wr2,err)
   
    call ADIOS_PERFORM_READS(buf_id,err)
    call ADIOS_READ_CLOSE(buf_id,err)
    write(*,*)'read coef_XGC_GEM completed'
  endif 
 
end subroutine restore_coef
#endif 
end module mapping
!< Initializes a simple 1D cubic spline interpolation for use with
!< interpol_spl_1d. Sample points are
!< denoted (x_i,y_i)
!! @param spline (inout) type(EZspline1_r8), spline data object
!! @param levelsin (in), x-values of sample points
!! @param data_in (in), y-values of sample points
!! @param nlevels (in) integer, number of sample points
!! @param pbc (in) logical, .true. for periodic bd conditions, .false. else
subroutine init_1d_interpolation(spline,levelsin,data_in,nlevels,pbc)
  use gem_com, only: myid
  use EZspline_obj
  use EZspline
  implicit none
  integer :: i,j
  integer, dimension(2) :: BCS1
  integer :: ier
  integer :: nlevels
  real, dimension(nlevels) :: levelsin, data_in
  type(EZspline1_r8) :: spline
  logical :: pbc

  if (pbc) then
    BCS1=(/-1, -1/) ! periodic boundary conditions
  else
    BCS1 =(/0, 0/) ! not a knot
  endif
  !if(myid==0)then
  !  write(111,*)'before EZspline init',nlevels,BCS1
  !  call flush(111)
  !endif 
  call EZspline_init(spline,nlevels,BCS1,ier)
  call EZspline_error(ier)

  spline%x1 = levelsin

  !if(myid==0)then
    !write(111,*)'before EZspline setup',data_in
    !call flush(111)
  !endif
  do i=1,size(levelsin)-1
     if((levelsin(i+1)-levelsin(i))<0.)then
       write(111,*)'error in EZspline setup'
       write(111,*)size(levelsin),i,levelsin(i),levelsin(i+1)
       write(111,*)levelsin
       call flush(111)
     endif 
  enddo
  call EZspline_setup(spline,data_in,ier)
  call EZspline_error(ier)

end subroutine init_1d_interpolation

!< Frees a spline data object.
!! @param spline (inout) type(EZspline1_r8), spline data object
subroutine finalize_1d_interpolation(spline)
  use EZspline_obj
  use EZspline
  implicit none
  integer ier
  type(EZspline1_r8) :: spline

  call EZspline_free(spline,ier)
  call EZspline_error(ier)
end subroutine finalize_1d_interpolation

!< 1D interpolation using a pre-calculated cubic spline
!! @param spline (in) type(EZspline1_r8), spline data object
!! @param ideriv (in) integer, order of derivative
!! @param inval (in) real(8), x-value for interpolation
real function interpol_spl_1d(inval,ideriv,spline)
    use EZspline_obj
    use EZspline
    implicit none
    real, intent(in) :: inval
    integer, intent(in) :: ideriv
    type(EZspline1_r8), intent(in) :: spline
    integer :: ier
    real :: r8value

    if (ideriv==0) then
      call EZspline_interp(spline,inval,r8value,ier)
    else
      call EZspline_derivative(spline,ideriv,inval,r8value,ier)
    endif
    call EZspline_error(ier)
    interpol_spl_1d=r8value
end function interpol_spl_1d

