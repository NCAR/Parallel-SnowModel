---   
# MODULES AND ENVIRONMENTS  
---     
  
## **MINIMAL REQUIREMENTS**  
Before executing the Parallel SnowModel the appropriate compilers and modules need to be downloaded. Below are some examples of a few compilers and modules that can be used to run Parallel SnowModel and test scripts to confirm that the modules are working appropriately.   
***Examples below can be used as a template but may vary slightly depending on the computing environment and versions of the downloaded modules. Additionally, there are other fortran compilers and coarray implementations that can be used but are not exampled here.*** 

### **Required Compilers / Modules** 
The following modules and compilers can be used execute Parallel SnowModel in serial or parallel computing environments.

  - Fortran 
    - recommended: gfortran
      - https://gcc.gnu.org/ 
  - OpenMPI
    - https://www.open-mpi.org/software/ompi/v4.1/
  - OpenCoarrays
    - http://www.opencoarrays.org/  
    
#### **Testing Fortran + OpenCoarrays**  
1). Load gfortran, OpenMPI, OpenCoarray modules using linux terminal on supercomputing environment.    
    - ***Example:***    
               `module load gnu/9.1.0`    
               `module load openmpi/4.0.5`  
               `module load opencoarrays/2.9.0`    
2). Compile hello_world.f90  
    `caf -o hello_world hello_world.f90`    
3). Run executable with four processes   
    `cafrun -np 4 hello_world`    
4). Confirm output... *should be similar to below*  
    `Hello world! I am process number: 0 out of 4`  
    `Hello world! I am process number: 1 out of 4`  
    `Hello world! I am process number: 2 out of 4`  
    `Hello world! I am process number: 3 out of 4`  


### **Additional Compilers / Modules** 
In addition to the required compilers / modules, a python environment is recommended for pre- and post-processing analysis and a netcdf module allows for forcing Parallel SnowModel with netcdf files. 
      
    
#### **NetCDF** 
  1). Download NetCDF Fortran Libraries    
      - https://downloads.unidata.ucar.edu/netcdf/    
  2). Link NetCDF and Fortran Libraries   
  
##### **Testing Fortran + OpenCoarrays + NetCDF**  
  1). Load gfortran, OpenMPI, OpenCoarray, and NetCDF modules using linux terminal on supercomputing environment.        
      - ***Example:***        
                 `module load gnu/9.1.0`      
                 `module load openmpi/4.0.5`    
                 `module load opencoarrays/2.9.0`  
                 `module load netcdf/4.7.3`  
  2). Identify absolute paths for netcdf include and lib directories  
      - ***Example:*** (creating the following environment variables with the appropriate paths):  
      `NCPATH=/glade/u/apps/ch/opt/netcdf/4.7.3/gnu/9.1.0`        
      `NC_INCLUDE=${NCPATH}/include`        
      `NC_LIB=${NCPATH}/lib`        
  2). Compile hello_world.f90    
      `caf -o hello_world_nc -I${NC_INCLUDE} -L${NC_LIB} -lnetcdf hello_world_nc.f90`      
  3). Run executable with four processes     
      `cafrun -np 4 hello_world_nc`      
  4). Confirm output *should be similar to below*    
      `NetCDF worked! I am process number: 0 out of 4`    
      `NetCDF worked! I am process number: 1 out of 4`    
      `NetCDF worked! I am process number: 2 out of 4`    
      `NetCDF worked! I am process number: 3 out of 4`    

---    
  

