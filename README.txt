to build, you need:
  -> cblas, boost
  -> haders from https://github.com/IwM341/GridObject
    * in CMAKE set variable GO_PATH to GridObject/include

executables:
  -> annihilation : creates annihilation matrix from Lgrid and Hgrid
  -> annihilation_rd : gives radial distribution of annihilation of distrib between captured and captured

  -> annihilation_vector : creates annihilation vector from grid
  -> annihilation_rd_vector : gives radial distribution of annihilation of distrib betwee captured and halo
  -> r_distrib : gives radial distribution of captured wimps
	
  -> matrix_mult : evolution of inelastic DM
  -> matrix_mult_elastic_mult : evolution of elastic DM

  -> conversions : convert histogramm (grid + values) to density function
  all parametrs in executables are listed in form: -[parametr_name] [parametr_value] or
  -[parametr_name]=[parametr_value] or use config file in json format, writing -config [path/to/config/file],
  in that case all filename pathes would be relative to config file.

form factors:
  -> in folder factors the file factors.hpp contain typedefs for using form factors:
    * dF_Nuc - nuclear form factor struct
    * Phi_Fac_S - phi factor of scatter 
    * Phi_Fac_Electron_S - phi factor, using to scatter with hydrogen, interacing with electron
    * Phi_Fac_Ann - annihilation phi factor
  -> file __factors.hpp contain implementation of this factors



