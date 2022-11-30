
module adios2_comm_module
    use gem_com,  only: myid
    use mpi
#ifdef __ADIOS2
    use adios2
#endif
    implicit none
#ifdef __ADIOS2
    type(adios2_adios) :: adios2obj
    type(adios2_engine), allocatable :: list_engines(:)
    integer :: n_engines

    !! very simple single timer
    character(len=128) :: timer_name
    integer :: timer_comm
    integer :: timer_index
    real(kind=8) :: t_start

contains
    subroutine adios2_comm_init(initfile)
        implicit none
        character(len=*), intent(in) :: initfile
        integer :: ierr
		
		call adios2_init(adios2obj, initfile, mpi_comm_world, .true., ierr)
        allocate(list_engines(16))
        n_engines = 0
    end subroutine adios2_comm_init

    subroutine adios2_comm_finalize()
        implicit none
        integer :: ierr
        integer :: i

        do i = 1, n_engines
            if (myid.eq.0) print *, 'ADIOS2: close output ', trim(list_engines(i)%name)
            call adios2_close(list_engines(i), ierr)
        enddo
        call adios2_finalize(adios2obj, ierr)
    end subroutine adios2_comm_finalize

    subroutine adios2_comm_engine_push(engine)
        implicit none
        type(adios2_engine), intent(in) :: engine
        type(adios2_engine), allocatable :: tmp(:)

        if (n_engines.ge.size(list_engines)) then
            if (myid.eq.0) print *, 'Increasing the size of list_engines to ', size(list_engines)*2
            allocate(tmp(size(list_engines)*2))
            tmp(1:n_engines) = list_engines(1:n_engines)
            deallocate(list_engines)
            call move_alloc(tmp,list_engines)
        endif
        n_engines = n_engines+1
        list_engines(n_engines) = engine
        if (myid.eq.0) print *, 'ADIOS2: push to close on finalizing ', trim(list_engines(n_engines)%name)
    end subroutine 
#endif
end module adios2_comm_module
