---  
# PARALLEL SNOWMODEL CODE DIRECTORY 
---    
This directory contains the fortran files and compile script required to run Parallel SnowModel using either parallel or serial compilers.   
## **Description**    
The `./compile_snowmodel.script` contains several examples of compile statements using different serial and parallel compilers, as well as additional modules (e.g. NETCDF) on lines `91 - 121`. Alterations to the script may be needed to accommodate other fortran compilers. Below are the steps to compile the code: 
1. Load the appropriate module(s) using the following terminal command `module load {...}`.  
2. Update the `./compile_snowmodel.script` to reflect the desired module combinations.  
    - Uncomment the desired compile line and comment the undesired compile statements with the `#` symbol.  
    - If necessary, change the file paths of the linked and included modules to reflect the appropriate locations within the computing environment.  
3. Compile the code by running `./compile_snowmodel.script` in the terminal.
4. Confirm the executable was created (or updated) in `../` directory.