# =========================================================================
# If NO_WMBASIC_SRCDIRS is defined on the command line, all f90 files are
# assumed to be in the current directory; otherwise, they are assumed to be
# in subdirectories "general", "hydro", "occup", etc. under the parent dir.
# To enable compile time and run time debug information, use "make debug=1"
# To compile the 3D version of occup, use "make rt3d=1"; this can also
# be used in combination with debug.


# =========================================================================
# optional Math Kernel Library settings

ifndef MKL_SHLIB
  MKL_SHLIB = -mkl=parallel
  MKL_FLAGS = -mkl=parallel
endif

# =========================================================================
# operating system conventions

ifdef COMSPEC # Windows
  ds           = /
  exeext       = .exe
  objext       = .obj
  modext       = .mod
  ifdef rt3d # should be set in "submit"
    m3d        = /module:3d
  else
    m3d        =
  endif
  Fortran      = ifort
  compile      = /compile_only
  outexe       = /exe:
  outobj       = /object:
  generflags   = /Od /real_size:64 /fpe:0 /align:dcommons /Qsave
  ifdef debug # Compile with all checks enabled
    generflags = += /traceback /check:bounds /Qfpstkchk /warn:all /check:all /check:noarg_temp_created /nogen-interfaces /Qtrapuv
    optimflags = $(generflags)
  else
    optimflags = $(generflags) /O3
  endif
  occupflags   = $(optimflags)
else # Linux
  ds           = /
  exeext       =
  objext       = .o
  modext       = .mod
  ifdef rt3d # should be set in "submit"
    m3d        = -module 3d
  else
    m3d        =
  endif
  Fortran      = ifort
  compile      = -c
  outexe       = -o
  outobj       = -o
  generflags   = -V -O0  -r8 -fpe0 -align dcommons -save
  ifdef debug # Compile with all checks enabled
    generflags += -traceback -check bounds -fpstkchk -warn all -fltconsistency -check all -check noarg_temp_created -nogen-interfaces -ftrapuv
    optimflags = $(generflags)
  else
    optimflags = $(generflags) -O3
  endif
  ifdef perfana # Compile with performance analysis
    occupflags = $(optimflags) -openmp $(MKL_FLAGS) $(m3d) -pg
  else
    occupflags = $(optimflags) -openmp $(MKL_FLAGS) $(m3d)
  endif
endif

# =========================================================================
# source and destination paths

ifdef NO_WMBASIC_SRCDIRS # should be set in "submit"
  general_dir  =
  gtempd_dir   =
  hydro_dir    =
  occup_dir    =
  formal_dir   =
  occup_su_dir =
else # standard
  general_dir  = ..$(ds)general$(ds)
  gtempd_dir   = ..$(ds)gtemp$(ds)
  hydro_dir    = ..$(ds)hydro$(ds)
  occup_dir    = ..$(ds)occup$(ds)
  formal_dir   = ..$(ds)formal$(ds)
  occup_su_dir = ..$(ds)occup$(ds)supply$(ds)
  atom_dir     = ..$(ds)atom$(ds)
  mcgamma_dir  = ..$(ds)mcgamma$(ds)
endif

  Datadefs_f90 = $(general_dir)Datadefs.f90
  Modgen1d_f90 = $(general_dir)Modgen.f90
  Modnum_f90   = $(general_dir)Modnum.f90
  fsplit_f90   = $(general_dir)fsplit.f90
  abund_f90    = $(general_dir)abund.f90
  gtempd_f90   = $(gtempd_dir)gtempd.f90
  stwindg_f90  = $(hydro_dir)stwindg.f90
  SNmodel_f90  = $(hydro_dir)SNmodel.f90
  modhydro_f90 = $(hydro_dir)modhydro.f90
  Modocc_f90   = $(occup_dir)Modocc.f90
  f0_cray_f90  = $(occup_dir)f0_cray.f90
  fcino_f90    = $(occup_dir)fcino.f90
  fcinoa_f90   = $(occup_dir)fcinoa.f90
  fdatm_f90    = $(occup_dir)fdatm.f90
  file1_f90    = $(occup_dir)file1.f90
  file2_f90    = $(occup_dir)file2.f90
  file2b_f90   = $(occup_dir)file2b.f90
  file2c_f90   = $(occup_dir)file2c.f90
  zzf2_f90     = $(occup_dir)zzf2.f90
  zzf3_f90     = $(occup_dir)zzf3.f90
  formal_f90   = $(formal_dir)formal.f90
  formals_f90  = $(formal_dir)formal_sup.f90
  formalm_f90  = $(formal_dir)formal_mls.f90
  fpara_f90    = $(occup_su_dir)fpara.f90
  mlinreg_f90  = $(occup_su_dir)mlinreg.f90
  atom_f90     = $(atom_dir)mkatom.f90
  mcgamma_f90  = $(mcgamma_dir)mcgamma.f90

ifdef rt3d # should be set in "submit"
  d3d          = 3d$(ds)
  Modgen_f90   = $(general_dir)Modgen3d.f90
else
  Modgen_f90   = $(general_dir)Modgen.f90
  d3d          = .$(ds)
endif

  Datadefs_obj = Datadefs$(objext)
  Modgen_obj   = $(d3d)Modgen$(objext)
  Modgen1d_obj = Modgen$(objext)
  Modnum_obj   = $(d3d)Modnum$(objext)
  Modnum1d_obj = Modnum$(objext)
  fsplit_obj   = fsplit$(objext)
  abund_obj    = abund$(objext)
  gtempd_obj   = gtempd$(objext)
  stwindg_obj  = stwindg$(objext)
  SNmodel_obj  = SNmodel$(objext)
  modhydro_obj = modhydro$(objext)
  Modocc_obj   = $(d3d)Modocc$(objext)
  f0_cray_obj  = $(d3d)f0_cray$(objext)
  fcino_obj    = $(d3d)fcino$(objext)
  fcinoa_obj   = $(d3d)fcinoa$(objext)
  fdatm_obj    = $(d3d)fdatm$(objext)
  file1_obj    = $(d3d)file1$(objext)
  file2_obj    = $(d3d)file2$(objext)
  file2b_obj   = $(d3d)file2b$(objext)
  file2c_obj   = $(d3d)file2c$(objext)
  zzf2_obj     = $(d3d)zzf2$(objext)
  zzf3_obj     = $(d3d)zzf3$(objext)
  formal_obj   = $(d3d)formal$(objext)
  formals_obj  = $(d3d)formals$(objext)
  formalm_obj  = $(d3d)formalm$(objext)
  fpara_obj    = fpara$(objext)
  mlinreg_obj  = mlinreg$(objext)
  atom_obj     = mkatom$(objext)
  mcgamma_obj  = mcgamma$(objext)


  split_exe    = fsplit$(exeext)
  abund_exe    = abund$(exeext)
  gtemp_exe    = gtemp$(exeext)
  hydro_exe    = hydro$(exeext)
  snhyd_exe    = snhyd$(exeext)
  occup_exe    = $(d3d)occup$(exeext)
  fpara_exe    = fpara$(exeext)
  atom_exe     = mkatom$(exeext)
  mcgamma_exe  = mcgamma$(exeext)

# =========================================================================
# individual compiler settings

ifdef COMSPEC # Windows
  splitflags   = $(generflags)
  abundflags   = $(generflags)
  hydroflags   = $(generflags)
  snhydflags   = $(generflags)
  gtempflags   = $(generflags)
  fparaflags   = $(generflags)
  atomflags    = $(generflags)
  generlibs    = $(MKL_SHLIB)
  abundlibs    = $(generlibs)
  splitlibs    = $(generlibs)
  gtemplibs    = $(generlibs)
  snhydlibs    = $(generlibs)
  hydrolibs    = $(generlibs)
  fparalibs    = $(generlibs)
  occuplibs    = $(generlibs) /link /stack:0x3000000
  atomlibs     =
  mcgamlibs    =
else # Linux
  splitflags   = $(generflags)
  abundflags   = $(generflags)
  hydroflags   = $(generflags)
  snhydflags   = $(generflags)
  gtempflags   = $(generflags)
  fparaflags   = $(generflags)
  atomflags    = $(generflags)
  generlibs    = -static
  usemklibs    = $(MKL_SHLIB)
  abundlibs    = $(generlibs)
  splitlibs    = $(generlibs)
  gtemplibs    = $(generlibs) $(usemklibs)
  snhydlibs    = $(generlibs) $(usemklibs)
  hydrolibs    = $(generlibs) $(usemklibs)
  fparalibs    = $(generlibs) $(usemklibs)
  ifdef perfana # link with performance analysis
    occuplibs  = $(generlibs) -pg
  else
    occuplibs  = $(generlibs) $(usemklibs)
  endif
  atomlibs     =
  mcgamlibs    =
endif

# =========================================================================
# dependencies

ifdef rt3d # should be set in "submit"
  programs = $(abund_exe)   \
             $(occup_exe)
else
  programs = $(split_exe)   \
             $(abund_exe)   \
             $(gtemp_exe)   \
             $(hydro_exe)   \
             $(snhyd_exe)   \
             $(occup_exe)   \
             $(fpara_exe)
endif

  all: $(programs)

  atom: $(atom_exe)

  mcgamma: $(mcgamma_exe)

  clean:
	-rm -r 3d *$(objext) *$(modext) $(programs) $(atom_exe) $(mcgamma_exe)

  genermods   = $(Datadefs_obj) $(Modgen_obj)   $(Modnum_obj)
  splitmods   = $(Datadefs_obj)
  abundmods   = $(Datadefs_obj) $(Modgen1d_obj)
  gtempmods   = $(genermods)
  hydromods   = $(genermods) $(modhydro_obj)
  snhydmods   = $(genermods) $(modhydro_obj)
  occupmods   = $(genermods) $(Modocc_obj) \
                $(file2b_obj) $(zzf2_obj) $(formals_obj) $(fdatm_obj)
  fparamods   = $(genermods)
  atommods    = $(Datadefs_obj) $(Modgen1d_obj)
  mcgammamods = $(Datadefs_obj) $(Modgen1d_obj)

  occupclean:
	-rm \
      $(occupmods) \
      $(f0_cray_obj) $(fcinoa_obj) $(fcino_obj) $(fdatm_obj) \
      $(file1_obj) $(file2_obj) $(file2b_obj) $(file2c_obj) \
      $(zzf3_obj) $(formals_obj) $(formalm_obj) $(formal_obj)


  $(Datadefs_obj): $(Datadefs_f90)
	$(Fortran) $(compile) \
      $(Datadefs_f90) $(outobj)$(Datadefs_obj)

  $(Modgen_obj): $(Modgen_f90) $(Datadefs_obj)
	-mkdir -p $(d3d)
	$(Fortran) $(compile) $(generflags) \
      $(Modgen_f90) $(outobj)$(Modgen_obj)

  $(Modgen1d_obj): $(Modgen1d_f90) $(Datadefs_obj)
	$(Fortran) $(compile) $(generflags) \
      $(Modgen1d_f90) $(outobj)$(Modgen1d_obj)

  $(Modnum_obj): $(Modnum_f90) $(Datadefs_obj) $(Modgen_obj)
	$(Fortran) $(compile) $(generflags) \
      $(Modnum_f90) $(outobj)$(Modnum_obj)


  $(Modnum1d_obj): $(Modnum_f90) $(Datadefs_obj) $(Modgen1d_obj)
	$(Fortran) $(compile) $(generflags) \
      $(Modnum_f90) $(outobj)$(Modnum1d_obj)


  $(split_exe): $(fsplit_obj) $(splitmods)
	$(Fortran) \
      $(fsplit_obj) $(splitlibs) $(outexe)$(split_exe)

  $(fsplit_obj): $(fsplit_f90) $(splitmods)
	$(Fortran) $(compile) $(splitflags) $(splitlibs) \
      $(fsplit_f90) $(outobj)$(fsplit_obj)

  $(abund_exe): $(abund_obj)
	$(Fortran) $(abund_obj) $(abundmods) $(abundlibs) \
               $(outexe)$(abund_exe)

  $(abund_obj): $(abund_f90) $(abundmods)
	$(Fortran) $(compile) $(abundflags) \
      $(abund_f90) $(outobj)$(abund_obj)

  $(gtemp_exe): $(gtempd_obj)
	$(Fortran) \
      $(gtempd_obj) $(gtempmods) $(gtemplibs) $(outexe)$(gtemp_exe)

  $(gtempd_obj): $(gtempd_f90) $(gtempmods)
	$(Fortran) $(compile) $(gtempflags) \
      $(gtempd_f90) $(outobj)$(gtempd_obj)

  $(hydro_exe): $(stwindg_obj)
	$(Fortran) \
      $(stwindg_obj) $(hydromods) $(hydrolibs) $(outexe)$(hydro_exe)

  $(stwindg_obj): $(stwindg_f90) $(hydromods)
	$(Fortran) $(compile) $(hydroflags) \
      $(stwindg_f90) $(outobj)$(stwindg_obj)

  $(snhyd_exe): $(SNmodel_obj)
	$(Fortran) \
    $(SNmodel_obj) $(snhydmods) $(snhydlibs) $(outexe)$(snhyd_exe)

  $(SNmodel_obj): $(SNmodel_f90) $(snhydmods)
	$(Fortran) $(compile) $(snhydflags) \
      $(SNmodel_f90) $(outobj)$(SNmodel_obj)

  $(modhydro_obj): $(modhydro_f90) $(genermods)
	$(Fortran) $(compile) $(hydroflags) \
      $(modhydro_f90) $(outobj)$(modhydro_obj)

  $(occup_exe): $(occupmods) \
      $(f0_cray_obj) $(file2c_obj) $(fcinoa_obj) $(fcino_obj) \
      $(zzf3_obj) $(file1_obj) $(file2_obj) \
      $(formalm_obj) $(formal_obj)
	$(Fortran) \
      $(f0_cray_obj) $(fcinoa_obj) $(fcino_obj) \
      $(file1_obj) $(file2_obj) $(file2c_obj) \
      $(zzf3_obj) $(formalm_obj) $(formal_obj) \
      $(occupmods) \
      $(outexe)$(occup_exe) $(occuplibs)


  $(zzf2_obj): $(zzf2_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(zzf2_f90) $(outobj)$(zzf2_obj)

  $(Modocc_obj): $(Modocc_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(Modocc_f90) $(outobj)$(Modocc_obj)

  $(f0_cray_obj): $(f0_cray_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(f0_cray_f90) $(outobj)$(f0_cray_obj)

  $(fcinoa_obj): $(fcinoa_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(fcinoa_f90) $(outobj)$(fcinoa_obj)

  $(fcino_obj): $(fcino_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(fcino_f90) $(outobj)$(fcino_obj)

  $(fdatm_obj): $(fdatm_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(fdatm_f90) $(outobj)$(fdatm_obj)

  $(file1_obj): $(file1_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(file1_f90) $(outobj)$(file1_obj)

  $(file2b_obj): $(file2b_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(file2b_f90) $(outobj)$(file2b_obj)

  $(file2_obj): $(file2_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(file2_f90) $(outobj)$(file2_obj)

  $(file2c_obj): $(file2c_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(file2c_f90) $(outobj)$(file2c_obj)

  $(zzf3_obj): $(zzf3_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(zzf3_f90) $(outobj)$(zzf3_obj)

  $(formals_obj): $(formals_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(formals_f90) $(outobj)$(formals_obj)

  $(formalm_obj): $(formalm_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(formalm_f90) $(outobj)$(formalm_obj)

  $(formal_obj): $(formal_f90) $(genermods)
	$(Fortran) $(compile) $(occupflags) \
      $(formal_f90) $(outobj)$(formal_obj)

  $(fpara_exe): $(fpara_obj) $(mlinreg_obj)
	$(Fortran) \
      $(fpara_obj) $(mlinreg_obj) $(fparamods) $(fparalibs) \
      $(outexe)$(fpara_exe)

  $(fpara_obj): $(fpara_f90) $(fparamods)
	$(Fortran) $(compile) $(fparaflags) \
      $(fpara_f90) $(outobj)$(fpara_obj)

  $(mlinreg_obj): $(mlinreg_f90) $(fparamods)
	$(Fortran) $(compile) $(fparaflags) \
      $(mlinreg_f90) $(outobj)$(mlinreg_obj)

  $(atom_exe): $(atom_obj)
	$(Fortran) \
      $(atom_obj) $(atommods) $(atomlibs) \
      $(outexe)$(atom_exe)

  $(atom_obj): $(atom_f90) $(atommods)
	$(Fortran) $(compile) $(atomflags) \
      $(atom_f90) $(outobj)$(atom_obj)

  $(mcgamma_exe): $(mcgamma_obj)
	$(Fortran) \
      $(mcgamma_obj) $(mcgammamods) $(mcgamlibs) \
      $(outexe)$(mcgamma_exe)

  $(mcgamma_obj): $(mcgamma_f90) $(mcgammamods)
	$(Fortran) $(compile) $(optimflags) \
      $(mcgamma_f90) $(outobj)$(mcgamma_obj)

# ========================================================================
# vi: set tabstop=4:
# ========================================================================
