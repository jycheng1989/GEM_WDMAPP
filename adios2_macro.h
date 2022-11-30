/* Adios2 macros */
#define ADIOS2_DECLARE(io,grp,err) \
  call adios2_declare_io(io,adios2obj,trim(grp)//char(0),err)

#if defined(__PGI) && (__PGIC__<15) || defined(__GFORTRAN__)
#define ADIOS2_DEFINE(variable,io,var,err) \
  call adios2_comm_define_variable(variable,io,'var'//char(0),var,err)
#else
#define ADIOS2_DEFINE(variable,io,var,err) \
  call adios2_comm_define_variable(variable,io,#var//char(0),var,err)
#endif
#define ADIOS2_DEFINE_LBL(variable,io,label,var,err) \
  call adios2_comm_define_variable(variable,io,trim(label)//char(0),var,err)

#define ADIOS2_OPEN(engine,io,file,mode,grp_comm,err) \
  call adios2_open(engine,io,trim(file)//char(0),mode,grp_comm,err)

#if defined(__PGI) && (__PGIC__<15) || defined(__GFORTRAN__)
#define ADIOS2_WRITE(engine,var,err) \
  call adios2_put(engine,'var'//char(0),var,err)
#else
#define ADIOS2_WRITE(engine,var,err) \
  call adios2_put(engine,#var//char(0),var,err)
#endif
#define ADIOS2_WRITE_LBL(engine,label,var,err) \
  call adios2_put(engine,trim(label)//char(0),var,err)
#define ADIOS2_WRITE_LBL_SYNC(engine,label,var,err) \
  call adios2_put(engine,trim(label)//char(0),var,adios2_mode_sync,err)

#define ADIOS2_CLOSE(engine,err) \
  call adios2_close(engine,err)
