!author Junyi Cheng

#ifdef __GEM_XGC_COUPLING
module coupling_core_edge
  implicit none

  !cce_overlap_first_surface, cce_overlap_last_surface: first/last overlap surface
  !cce_closed_last_surface: the total surface number in the closed-flux region
  integer :: cce_overlap_first_surface, cce_overlap_last_surface, cce_closed_last_surface
  !cce_folder : folder for data exchange
  character(256) :: cce_folder
  !cce_in_boundary, cce_out_boundary: the overlap weight
  real :: cce_in_boundary,cce_out_boundary

contains

subroutine cce_initialize
  use mpi
  implicit none

  namelist /coupling/ cce_overlap_first_surface,cce_overlap_last_surface,cce_closed_last_surface,cce_folder

  open(unit=20,file='mapping.in', status='old',action='read')
  READ(20, NML=coupling)
  close(unit=20)

  cce_in_boundary=real(cce_overlap_first_surface-1)/real(cce_closed_last_surface-1)
  cce_out_boundary=real(cce_overlap_last_surface-1)/real(cce_closed_last_surface-1)
  
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
  real :: iden_xgc(cce_all_node_number,nphi)
  
  call mapping_GEM_XGC_new(density_3d,pot_density_3d,iden,grid,mapping_GEM_XGC_linear_coef)
  iden_xgc(:,:)=pot_density_3d(:,:)*nu

  !do i=2,nphi
  !   iden_xgc(:,i)=pot_density_3d(:,nphi-i+2)*nu
  !enddo
  if(myid==0)then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=0
      ldims(1)=cce_all_node_number
      ldims(2)=nphi
      
      call adios2_declare_io(io,adios2obj,'ion_density',ierr)
      call adios2_define_variable(varid, io, "idensity",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'ion_density.bp', adios2_mode_write, mpi_comm_self, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "idensity",  iden_xgc, ierr)
      call adios2_end_step(engine,ierr) 
  endif
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
  real :: eden_xgc(cce_all_node_number,nphi)

  call mapping_GEM_XGC_new(density_3d,pot_density_3d,eden,grid,mapping_GEM_XGC_linear_coef)

  eden_xgc(:,:)=pot_density_3d(:,:)*nu
  !do i=2,nphi
  !   eden_xgc(:,i)=pot_density_3d(:,nphi-i+2)*nu
  !enddo

  if(myid==0)then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=0
      ldims(1)=cce_all_node_number
      ldims(2)=nphi

      call adios2_declare_io(io,adios2obj,'electron_density',ierr)
      call adios2_define_variable(varid, io, "edensity",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'electron_density.bp', adios2_mode_write, mpi_comm_self, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "edensity",  eden_xgc, ierr)
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
  real :: ijpar_xgc(cce_all_node_number,nphi)

  e=1.6022e-19
  j_env=nu*vu

  call mapping_GEM_XGC_new(density_3d,pot_density_3d,ijpar,grid,mapping_GEM_XGC_linear_coef)
  ijpar_xgc(:,:)=pot_density_3d(:,:)*j_env
  !do i=2,nphi
  !   ijpar_xgc(:,i)=pot_density_3d(:,nphi-i+2)*j_env
  !enddo

  if(myid==0)then
    if(.not. init)then
      gdims(1)=cce_all_node_number
      gdims(2)=nphi
      goffset(1)=0
      goffset(2)=0
      ldims(1)=cce_all_node_number
      ldims(2)=nphi

      call adios2_declare_io(io,adios2obj,'ion_current',ierr)
      call adios2_define_variable(varid, io, "ijpar",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'ion_current.bp', adios2_mode_write, mpi_comm_self, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "ijpar", ijpar_xgc, ierr)
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
  real :: ejpar_xgc(cce_all_node_number,nphi)
  e=1.6022e-19
  j_env=nu*vu

  call mapping_GEM_XGC_new(density_3d,pot_density_3d,ejpar,grid,mapping_GEM_XGC_linear_coef)

  ejpar_xgc(:,:)=pot_density_3d(:,:)*j_env
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

      call adios2_declare_io(io,adios2obj,'electron_current',ierr)
      call adios2_define_variable(varid, io, "ejpar",  adios2_type_dp, 2, gdims, goffset, ldims, adios2_constant_dims, ierr)
      call adios2_open(engine, io, trim(cce_folder) // '/' //'electron_current.bp', adios2_mode_write, mpi_comm_self, ierr)
      call adios2_comm_engine_push(engine)
      init=.true.
    endif

      call adios2_begin_step(engine, adios2_step_mode_append, ierr)
      call adios2_put(engine, "ejpar",  ejpar_xgc, ierr)
      call adios2_end_step(engine,ierr)
  endif

end subroutine cce_send_current_electron

subroutine cce_receive_pot(phi_gem)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_XGC_GEM_new, pot_field_3d, field_3d, grid, mapping_XGC_GEM_linear_coef, nphi, cce_all_node_number
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
  real :: dpot_xgc(cce_all_node_number,nphi)

  e=1.6022e-19
  phi_env=e/Tu

  if(MyId==0)then
    starts(1)=0
    counts(1)=cce_all_node_number
    starts(2)=0
    counts(2)=nphi
    if(.not. init)then
      call adios2_declare_io(io_field,adios2obj,'pot',ierr)
      call adios2_open(engine, io_field, trim(cce_folder) // '/' // 'pot.bp', adios2_mode_read, mpi_comm_self, ierr)
      init=.true.
    endif
    call adios2_begin_step(engine, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(varid, io_field, "pot", ierr)
    call adios2_set_selection(varid, 2, starts, counts, ierr)
    call adios2_get(engine, varid, pot_field_3d, adios2_mode_deferred, ierr)
    call adios2_end_step(engine, ierr)
    pot_field_3d=pot_field_3d*phi_env
    dpot_xgc(:,:)=pot_field_3d(:,:)
    !do i=2,nphi
    !   dpot_xgc(:,i)=pot_field_3d(:,nphi-i+2)
    !enddo
  endif

  call mapping_XGC_GEM_new(dpot_xgc,phi_gem,field_3d,grid,mapping_XGC_GEM_linear_coef)

end subroutine cce_receive_pot
subroutine cce_receive_as(as_gem)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_XGC_GEM_new, pot_field_3d, field_3d, grid, mapping_XGC_GEM_linear_coef, nphi, cce_all_node_number
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
  real :: As_xgc(cce_all_node_number,nphi)

  e=1.6022e-19
  A_env=e*vu/Tu

  if(MyId==0)then
    starts(1)=0
    counts(1)=cce_all_node_number
    starts(2)=0
    counts(2)=nphi
    if(.not. init)then
      call adios2_declare_io(io_field,adios2obj,'As',ierr)
      call adios2_open(engine, io_field, trim(cce_folder) // '/' // 'as.bp', adios2_mode_read, mpi_comm_self, ierr)
      init=.true.
    endif
    call adios2_begin_step(engine, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(varid, io_field, "As", ierr)
    call adios2_set_selection(varid, 2, starts, counts, ierr)
    call adios2_get(engine, varid, pot_field_3d, adios2_mode_deferred, ierr)
    call adios2_end_step(engine, ierr)
    pot_field_3d=pot_field_3d*A_env
    As_xgc(:,:)=-pot_field_3d(:,:)
    !do i=2,nphi
    !   As_xgc(:,i)=pot_field_3d(:,nphi-i+2)
    !enddo
  endif

  call mapping_XGC_GEM_new(As_xgc,as_gem,field_3d,grid,mapping_XGC_GEM_linear_coef)
end subroutine cce_receive_as
subroutine cce_receive_ah(ah_gem)
  use gem_com, only : imx,jmx
  use mapping, only : mapping_XGC_GEM_new, pot_field_3d, field_3d, grid, mapping_XGC_GEM_linear_coef, nphi, cce_all_node_number
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
  real :: Ah_xgc(cce_all_node_number,nphi)
  e=1.6022e-19
  A_env=e*vu/Tu

  if(MyId==0)then
    starts(1)=0
    counts(1)=cce_all_node_number
    starts(2)=0
    counts(2)=nphi
    if(.not. init)then
      call adios2_declare_io(io_field,adios2obj,'Ah',ierr)
      call adios2_open(engine, io_field, trim(cce_folder) // '/' // 'ah.bp', adios2_mode_read, mpi_comm_self, ierr)
      init=.true.
    endif
    call adios2_begin_step(engine, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(varid, io_field, "Ah", ierr)
    call adios2_set_selection(varid, 2, starts, counts, ierr)
    call adios2_get(engine, varid, pot_field_3d, adios2_mode_deferred, ierr)
    call adios2_end_step(engine, ierr)
    pot_field_3d=pot_field_3d*A_env
    Ah_xgc(:,:)=-pot_field_3d(:,:)
    !do i=2,nphi
    !   Ah_xgc(:,i)=pot_field_3d(:,nphi-i+2)
    !enddo
  endif

  call mapping_XGC_GEM_new(Ah_xgc,ah_gem,field_3d,grid,mapping_XGC_GEM_linear_coef)
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

  call mapping_GEM_XGC_new(density_3d,pot_density_3d,ejpar,grid,mapping_GEM_XGC_linear_coef)

  ejpar_xgc(:,:)=pot_density_3d(:,:)/j_env
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

  call mapping_GEM_XGC_new(density_3d,pot_density_3d,ejpar,grid,mapping_GEM_XGC_linear_coef)

  ejpar_xgc(:,:)=pot_density_3d(:,:)/j_env
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

  call mapping_GEM_XGC_new(density_3d,pot_density_3d,ejpar,grid,mapping_GEM_XGC_linear_coef)

  ejpar_xgc(:,:)=pot_density_3d(:,:)/j_env
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
