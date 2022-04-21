# Raytracer

## Building
Lines 78-90 of "stl_raytracer.cu" contains macros that can be edited to change the input file and camera parameters. Various STL files have been included in this repo.
Run the startup script to configure your node: 
```
./build_script
```

## Running
Run the following command
```
sbatch run_on_pascal
```

## Evaluating Results
Performance will be a printed in the output file with the naming convention "pascal_job.xx.out". To evaluate the resuling images, scp the output directory to your host machine from the Discovery Cluster and view them with a .ppm viewer