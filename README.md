# Raytracer

## Building
Lines 10-25 of "stl_raytracer.cu" contain configuration macros, including file selection, image dimensioning, etc. that can be edited to change the input file and camera parameters. To adjust the rotation of the object, refer to the rotate_stl functions (lines 184-185, 196-198). These are relative, meaning that they will adjust the rotation relative to the current orientation. Various STL files have been included in this repo.
Run the following sequence of commands to configure your node: 
```
module load cuda/10.2
module load openmpi
./build_script
```
Running this on hardware other than Pascal GPUs requries an adjustment of the Makefile (line 8) specifying the GPU architecture. Refer to the link in the Makefile for the value for your respective GPU.

## Running
Run the following command
```
sbatch run_on_pascal
```

## Evaluating Results
Performance will be a printed in the output file with the naming convention "pascal_job.xx.out". To evaluate the resuling images, scp the output directory to your host machine from the Discovery Cluster and view them with a .ppm viewer. On linux, imagemagick can be used to generate a gif out of these output files. A few are linked in the output folder to demonstrate what may be visible. Some stuttering is visible in the gifs, however, this is an issue coming from imagemagick, as the jumps are not present in the generated files.