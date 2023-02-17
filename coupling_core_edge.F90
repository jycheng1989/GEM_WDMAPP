!author Junyi Cheng

#ifdef __GEM_XGC_COUPLING
module coupling_core_edge
  implicit none

  integer :: cce_overlap_first_surface, cce_overlap_last_surface, cce_closed_last_surface
  character(256) :: cce_folder
  real :: cce_in_boundary,cce_out_boundary

  integer :: iocomm, iocolor=0, iostart, iocount=1

contains


subroutine cce_setup_io_comm
  use gem_com, only: gem_comm_world, myid, numprocs
  use mapping, only: procsperwriter, nphi
  use mpi

  implicit none
  integer :: onnodecomm, onnoderank, onnodesize
  integer :: crossnodecomm, crossnoderank
  integer :: newnode=0, nodecount=0, procmin=0, procmax=0
  integer :: ierr, extra=0, each=1

  if (numprocs.lt.nphi) then
	  write(*, *) 'Must run with numproc >= nphi'
	  stop
  end if

  call mpi_comm_split_type(gem_comm_world, MPI_COMM_TYPE_SHARED, myid, MPI_INFO_NULL, onnodecomm, ierr)
  call mpi_comm_rank(onnodecomm, onnoderank, ierr)
  call mpi_comm_size(onnodecomm, onnodesize, ierr)
  call mpi_comm_split(gem_comm_world, onnoderank, myid, crossnodecomm, ierr)

  if (onnoderank.eq.0) then
	  call mpi_reduce(onnodesize, procmin, 1, MPI_INTEGER, MPI_MIN, 0, crossnodecomm, ierr)
	  call mpi_reduce(onnodesize, procmax, 1, MPI_INTEGER, MPI_MAX, 0, crossnodecomm, ierr)
	  if ((procmin.ne.procmax) .and. (onnoderank.eq.0)) then
		  write(*, *) 'Variable number of processes per node not currently supported in XGC coupling'
		  stop
	  end if
  end if

  if (onnoderank.eq.0) then
	  newnode = 1
  end if
  call mpi_allreduce(newnode, nodecount, 1, MPI_INTEGER, MPI_SUM, gem_comm_world, ierr)

  call mpi_comm_rank(crossnodecomm, crossnoderank, ierr)
  if (nodecount.lt.nphi) then
	  each = nphi / nodecount
	  extra = mod(nphi, nodecount)
  end if
  iostart = crossnoderank * each
  
  if ((mod(onnoderank,procsperwriter).eq.0) .and. (onnoderank/procsperwriter.lt.each) .and. (crossnoderank.lt.nphi)) then
	  iocolor = 1
	  iostart = iostart + onnoderank
  end if
  if ((nodecount.lt.nphi) .and. (crossnoderank.lt.extra) .and. (mod(onnoderank,procsperwriter).eq.0) .and. (onnoderank.eq.each)) then
	  iocolor = 1
	  iostart = nodecount * each + crossnoderank
  end if

  call mpi_comm_split(gem_comm_world, iocolor, myid, iocomm, ierr)
end subroutine cce_setup_io_comm


subroutine cce_initialize
  use mpi
  implicit none

  namelist /coupling/ cce_overlap_first_surface,cce_overlap_last_surface,cce_closed_last_surface,cce_folder

  open(unit=20,file='mapping.in', status='old',action='read')
  READ(20, NML=coupling)
  close(unit=20)

  cce_in_boundary=real(cce_overlap_first_surface-1)/real(cce_closed_last_surface-1)
  cce_out_boundary=real(cce_overlap_last_surface-1)/real(cce_closed_last_surface-1)

  call cce_setup_io_comm
  
end subroutine cce_initialize

subroutine cce_send_density_ion(iden)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_GEM_XGC_new, density_3d, pot_density_3d, grid, mapping_GEM_XGC_linear_coef,nphi,cce_all_node_number
  use gem_equil, only : nu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: iden(0:imx,0:jmx,0:1)
  logical, save :: init=.false.
  integer(8), dimension(2) :: gdims, goffset, ldims
  type(adios2_io), save :: io
  type(adios2_engine), save :: engine
  type(adios2_variable) :: varid

  integer :: i
  !real :: iden_xgc(cce_all_node_number,nphi)
  
  call mapping_GEM_XGC_new(iden,grid,mapping_GEM_XGC_linear_coef,nu)

  !do i=2,nphi
  !   iden_xgc(:,i)=pot_density_3d(:,nphi-i+2)*nu
  !enddo
  !pot_density_3d(:,:)=pot_density_3d(:,:)*nu

  if (iocolor==1) then
    if (.not. init) then
      gdims(1)=cce_all_node_number !1M
      gdims(2)=nphi !64
      goffset(1)=0
      goffset(2)=iostart
      ldims(1)=cce_all_node_number
      ldims(2)=iocount
      
      call adios2_declare_io(io,adios2obj,'ion_density',ierr)
      call adios2_define_variable(varid, io, "idensity",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'ion_density.bp', adios2_mode_write, iocomm, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "idensity",  pot_density_3d(:, iostart), ierr)
      call adios2_end_step(engine,ierr) 
  endif

  call mpi_barrier(mpi_comm_world,ierr)
end subroutine cce_send_density_ion

subroutine cce_send_density_electron(eden)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_GEM_XGC_new, density_3d, pot_density_3d, grid, mapping_GEM_XGC_linear_coef,nphi,cce_all_node_number
  use gem_equil, only : nu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: eden(0:imx,0:jmx,0:1)
  logical, save :: init=.false.
  integer(8), dimension(2) :: gdims, goffset, ldims
  type(adios2_io), save :: io
  type(adios2_engine), save :: engine
  type(adios2_variable) :: varid

  integer :: i

  call mapping_GEM_XGC_new(eden,grid,mapping_GEM_XGC_linear_coef,nu)

  !do i=2,nphi
  !   eden_xgc(:,i)=pot_density_3d(:,nphi-i+2)*nu
  !enddo

  if (iocolor==1) then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=iostart
      ldims(1)=cce_all_node_number
      ldims(2)=iocount

      call adios2_declare_io(io,adios2obj,'electron_density',ierr)
      call adios2_define_variable(varid, io, "edensity",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'electron_density.bp', adios2_mode_write, iocomm, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "edensity",  pot_density_3d(:, iostart), ierr)
      call adios2_end_step(engine,ierr)
  endif
end subroutine cce_send_density_electron

subroutine cce_send_current_ion(ijpar)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_GEM_XGC_new, density_3d, pot_density_3d, grid, mapping_GEM_XGC_linear_coef,nphi,cce_all_node_number
  use gem_equil, only : nu,vu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: ijpar(0:imx,0:jmx,0:1)
  logical, save :: init=.false.
  integer(8), dimension(2) :: gdims, goffset, ldims
  type(adios2_io), save :: io
  type(adios2_engine), save :: engine
  type(adios2_variable) :: varid
  real :: e,j_env

  integer :: i
  !real :: ijpar_xgc(cce_all_node_number,nphi)

  e=1.6022e-19
  j_env=nu*vu

  call mapping_GEM_XGC_new(ijpar,grid,mapping_GEM_XGC_linear_coef,j_env)
  !ijpar_xgc(:,:)=pot_density_3d(:,:)*j_env
  !do i=2,nphi
  !   ijpar_xgc(:,i)=pot_density_3d(:,nphi-i+2)*j_env
  !enddo

  if (iocolor==1) then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=iostart
      ldims(1)=cce_all_node_number
      ldims(2)=iocount

      call adios2_declare_io(io,adios2obj,'ion_current',ierr)
      call adios2_define_variable(varid, io, "ijpar",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'ion_current.bp', adios2_mode_write, iocomm, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "ijpar", pot_density_3d(:, iostart), ierr)
      call adios2_end_step(engine,ierr)
  endif

end subroutine cce_send_current_ion

subroutine cce_send_current_electron(ejpar)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_GEM_XGC_new, density_3d, pot_density_3d, grid, mapping_GEM_XGC_linear_coef,nphi,cce_all_node_number
  use gem_equil, only : nu,vu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: ejpar(0:imx,0:jmx,0:1)
  logical, save :: init=.false.
  integer(8), dimension(2) :: gdims, goffset, ldims
  type(adios2_io), save :: io
  type(adios2_engine), save :: engine
  type(adios2_variable) :: varid
  real :: e,j_env
  
  integer :: i
  !real :: ejpar_xgc(cce_all_node_number,nphi)
  e=1.6022e-19
  j_env=nu*vu

  call mapping_GEM_XGC_new(ejpar,grid,mapping_GEM_XGC_linear_coef,j_env)

  !ejpar_xgc(:,:)=pot_density_3d(:,:)*j_env
  !do i=2,nphi
  !   ejpar_xgc(:,i)=pot_density_3d(:,nphi-i+2)*j_env
  !enddo

  if (iocolor==1) then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=iostart
      ldims(1)=cce_all_node_number
      ldims(2)=iocount

      call adios2_declare_io(io,adios2obj,'electron_current',ierr)
      call adios2_define_variable(varid, io, "ejpar",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'electron_current.bp', adios2_mode_write, iocomm, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "ejpar", pot_density_3d(:, iostart), ierr)
      call adios2_end_step(engine,ierr)
  endif

end subroutine cce_send_current_electron

subroutine cce_receive_pot(phi_gem)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_XGC_GEM_new, pot_density_3d, density_3d, grid, mapping_XGC_GEM_linear_coef, nphi, cce_all_node_number
  use gem_equil, only: Tu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: phi_gem(0:imx,0:jmx,0:1)
  integer(8), dimension(2) :: starts,counts
  type(adios2_io), save :: io_field, io_field_electron
  type(adios2_engine), save :: engine, engine_electron
  type(adios2_variable) :: varid
  logical, save :: init=.false.
  real :: e,phi_env

  integer :: i
  !real :: dpot_xgc(cce_all_node_number,nphi)
  real :: readarr(cce_all_node_number, iocount)

  e=1.6022e-19
  phi_env=e/Tu

  if (iocolor==1) then
    starts(1)=0
    counts(1)=cce_all_node_number
    starts(2)=iostart
    counts(2)=iocount
    if(.not. init)then
      call adios2_declare_io(io_field,adios2obj,'pot',ierr)
      call adios2_open(engine, io_field, trim(cce_folder) // '/' // 'pot.bp', adios2_mode_read, iocomm, ierr)
      init=.true.
    endif
    call adios2_begin_step(engine, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(varid, io_field, "pot", ierr)
    call adios2_set_selection(varid, 2, starts, counts, ierr)
    call adios2_get(engine, varid, readarr, adios2_mode_deferred, ierr)
    call adios2_end_step(engine, ierr)
    !pot_density_3d=pot_density_3d*phi_env
    !idpot_xgc(:,:)=pot_density_3d(:,:)
    !do i=2,nphi
    !   dpot_xgc(:,i)=pot_density_3d(:,nphi-i+2)
    !enddo

	call mpi_gather( &
		readarr,        iocount*cce_all_node_number, MPI_REAL8, &
		pot_density_3d, iocount*cce_all_node_number, MPI_REAL8, &
		0, iocomm, ierr &
		)
  endif

  call mapping_XGC_GEM_new(phi_gem,grid,mapping_XGC_GEM_linear_coef,phi_env)

end subroutine cce_receive_pot

subroutine cce_receive_as(as_gem)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_XGC_GEM_new, pot_density_3d, density_3d, grid, mapping_XGC_GEM_linear_coef, nphi, cce_all_node_number
  use gem_equil, only: Tu,vu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: as_gem(0:imx,0:jmx,0:1)
  integer(8), dimension(2) :: starts,counts
  type(adios2_io), save :: io_field, io_field_electron
  type(adios2_engine), save :: engine, engine_electron
  type(adios2_variable) :: varid
  logical, save :: init=.false.
  real :: e,A_env
  
  integer :: i
  !real :: As_xgc(cce_all_node_number,nphi)
  real :: readarr(cce_all_node_number, iocount)

  e=1.6022e-19
  A_env=e*vu/Tu

  if (iocolor==1) then
    starts(1)=0
    counts(1)=cce_all_node_number
    starts(2)=iostart
    counts(2)=iocount
    if(.not. init)then
      call adios2_declare_io(io_field,adios2obj,'As',ierr)
      call adios2_open(engine, io_field, trim(cce_folder) // '/' // 'as.bp', adios2_mode_read, iocomm, ierr)
      init=.true.
    endif
    call adios2_begin_step(engine, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(varid, io_field, "As", ierr)
    call adios2_set_selection(varid, 2, starts, counts, ierr)
    call adios2_get(engine, varid, readarr, adios2_mode_deferred, ierr)
    call adios2_end_step(engine, ierr)
    !pot_density_3d=pot_density_3d*A_env
    !As_xgc(:,:)=-pot_density_3d(:,:)
    !do i=2,nphi
    !   As_xgc(:,i)=pot_density_3d(:,nphi-i+2)
    !enddo

	call mpi_gather( &
		readarr,        iocount*cce_all_node_number, MPI_REAL8, &
		pot_density_3d, iocount*cce_all_node_number, MPI_REAL8, &
		0, iocomm, ierr &
		)
  endif

  call mapping_XGC_GEM_new(as_gem,grid,mapping_XGC_GEM_linear_coef,A_env)
end subroutine cce_receive_as
subroutine cce_receive_ah(ah_gem)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_XGC_GEM_new, pot_density_3d, density_3d, grid, mapping_XGC_GEM_linear_coef, nphi, cce_all_node_number
  use gem_equil, only: Tu,vu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: ah_gem(0:imx,0:jmx,0:1)
  integer(8), dimension(2) :: starts,counts
  type(adios2_io), save :: io_field, io_field_electron
  type(adios2_engine), save :: engine, engine_electron
  type(adios2_variable) :: varid
  logical, save :: init=.false.
  real :: e,A_env

  integer :: i
  !real :: Ah_xgc(cce_all_node_number,nphi)
  real :: readarr(cce_all_node_number, iocount)
  e=1.6022e-19
  A_env=e*vu/Tu

  if (iocolor==1) then
    starts(1)=0
    counts(1)=cce_all_node_number
    starts(2)=iostart
    counts(2)=iocount
    if(.not. init)then
      call adios2_declare_io(io_field,adios2obj,'Ah',ierr)
      call adios2_open(engine, io_field, trim(cce_folder) // '/' // 'ah.bp', adios2_mode_read, iocomm, ierr)
      init=.true.
    endif
    call adios2_begin_step(engine, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(varid, io_field, "Ah", ierr)
    call adios2_set_selection(varid, 2, starts, counts, ierr)
    call adios2_get(engine, varid, readarr, adios2_mode_deferred, ierr)
    call adios2_end_step(engine, ierr)
    !pot_density_3d=pot_density_3d*A_env
    !Ah_xgc(:,:)=-pot_density_3d(:,:)
    !do i=2,nphi
    !   Ah_xgc(:,i)=pot_density_3d(:,nphi-i+2)
    !enddo

	call mpi_gather( &
		readarr,        iocount*cce_all_node_number, MPI_REAL8, &
		pot_density_3d, iocount*cce_all_node_number, MPI_REAL8, &
		0, iocomm, ierr &
		)
  endif

  call mapping_XGC_GEM_new(ah_gem,grid,mapping_XGC_GEM_linear_coef,A_env)
end subroutine cce_receive_ah
subroutine cce_send_pot(ejpar)
  use gem_com, only : imx,jmx,ntube
  use mapping, only : mapping_GEM_XGC_new, density_3d, pot_density_3d, grid, mapping_GEM_XGC_linear_coef,nphi,cce_all_node_number
  use gem_equil, only : nu,vu,Tu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: ejpar(0:imx,0:jmx,0:1)
  logical, save :: init=.false.
  integer(8), dimension(2) :: gdims, goffset, ldims
  type(adios2_io), save :: io
  type(adios2_engine), save :: engine
  type(adios2_variable) :: varid
  real :: e,j_env

  integer :: i
  real :: ejpar_xgc(cce_all_node_number,nphi)
  e=1.6022e-19
  j_env=e/Tu*real(ntube)

  call mapping_GEM_XGC_new(ejpar,grid,mapping_GEM_XGC_linear_coef,j_env)

  ejpar_xgc(:,:)=pot_density_3d(:,:)
  !do i=2,nphi
  !   ejpar_xgc(:,i)=pot_density_3d(:,nphi-i+2)*j_env
  !enddo

  if(myid==0)then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=0
      ldims(1)=cce_all_node_number
      ldims(2)=nphi

      call adios2_declare_io(io,adios2obj,'dpotrev',ierr)
      call adios2_define_variable(varid, io, "dpotrev",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'dpotrev.bp', adios2_mode_write, mpi_comm_self, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "dpotrev",  ejpar_xgc, ierr)
      call adios2_end_step(engine,ierr)
  endif

end subroutine cce_send_pot
subroutine cce_send_as(ejpar)
  use gem_com, only : imx,jmx,ntube
  use mapping, only : mapping_GEM_XGC_new, density_3d, pot_density_3d, grid, mapping_GEM_XGC_linear_coef,nphi,cce_all_node_number
  use gem_equil, only : nu,vu,Tu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: ejpar(0:imx,0:jmx,0:1)
  logical, save :: init=.false.
  integer(8), dimension(2) :: gdims, goffset, ldims
  type(adios2_io), save :: io
  type(adios2_engine), save :: engine
  type(adios2_variable) :: varid
  real :: e,j_env

  integer :: i
  real :: ejpar_xgc(cce_all_node_number,nphi)
  e=1.6022e-19
  j_env=e*vu/Tu*real(ntube)

  call mapping_GEM_XGC_new(ejpar,grid,mapping_GEM_XGC_linear_coef,j_env)

  ejpar_xgc(:,:)=pot_density_3d(:,:)
  !do i=2,nphi
  !   ejpar_xgc(:,i)=pot_density_3d(:,nphi-i+2)*j_env
  !enddo

  if(myid==0)then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=0
      ldims(1)=cce_all_node_number
      ldims(2)=nphi

      call adios2_declare_io(io,adios2obj,'asrev',ierr)
      call adios2_define_variable(varid, io, "asrev",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'asrev.bp', adios2_mode_write, mpi_comm_self, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "asrev",  -ejpar_xgc, ierr)
      call adios2_end_step(engine,ierr)
  endif

end subroutine cce_send_as
subroutine cce_send_ah(ejpar)
  use gem_com, only : imx,jmx,ntube
  use mapping, only : mapping_GEM_XGC_new, density_3d, pot_density_3d, grid, mapping_GEM_XGC_linear_coef,nphi,cce_all_node_number
  use gem_equil, only : nu,vu,Tu
  use adios2_comm_module
  implicit none

  integer :: ierr
  real, intent(inout) :: ejpar(0:imx,0:jmx,0:1)
  logical, save :: init=.false.
  integer(8), dimension(2) :: gdims, goffset, ldims
  type(adios2_io), save :: io
  type(adios2_engine), save :: engine
  type(adios2_variable) :: varid
  real :: e,j_env

  integer :: i
  real :: ejpar_xgc(cce_all_node_number,nphi)
  e=1.6022e-19
  j_env=e*vu/Tu*real(ntube)

  call mapping_GEM_XGC_new(ejpar,grid,mapping_GEM_XGC_linear_coef,j_env)

  ejpar_xgc(:,:)=pot_density_3d(:,:)
  !do i=2,nphi
  !   ejpar_xgc(:,i)=pot_density_3d(:,nphi-i+2)*j_env
  !enddo

  if(myid==0)then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=0
      ldims(1)=cce_all_node_number
      ldims(2)=nphi

      call adios2_declare_io(io,adios2obj,'ahrev',ierr)
      call adios2_define_variable(varid, io, "ahrev",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'ahrev.bp', adios2_mode_write, mpi_comm_self, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "ahrev",  -ejpar_xgc, ierr)
      call adios2_end_step(engine,ierr)
  endif

end subroutine cce_send_ah
end module coupling_core_edge
#endif
