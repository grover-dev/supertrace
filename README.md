# Raytracer

## Building
Lines 24-35 of "stl_raytracer.cpp" contains macros that can be edited to change the input file and camera parameters. To adjust the rotation of the object, refer to the rotate_stl functions (lines 223,224,225). These are relative, meaning that they will adjust the rotation relative to the current orientation. Various STL files have been included in this repo. 
Run the startup script to configure your node: 
```
./build_script
```
If the build script does not work, try manually loading openmpi (module add openmpi).
To change the number of active nodes, edit the number of nodes in the run_on_mpi_cpu script. This value needs to be greater than or equal to 2 (1 node reserved for receiving+saving the data, the others for  rendering)
## Running
Run the following command
```
sbatch run_on_mpi_cpu
```

## Evaluating Results
Performance will be a printed in the output file with the naming convention "cpu_job.xx.out". To evaluate the resuling images, scp the output directory to your host machine from the Discovery Cluster and view them with a .ppm viewer.  On linux, imagemagick can be used to generate a gif out of these output files. A few are linked in the output folder to demonstrate what may be visible.