# Tutorial

## Getting Started
First thing that you need to do to implement this workflow is to download an up-to-date version of conda. From there you can read in the environment.yml file to build the conda environment that you will be using. \
Once the appropriate environment has been built and activated, the workflow can be run in the following two ways: 

1. Run each of the python scripts individually and pass the output from one to the next: 
```
python3 geom.py
python3 flow.py
python3 render.py
```
2. Run the workflow altogether using the snakefile:
```
snakemake
```
The second method is not only less work for you, but it will automatically clean up the files so that only the final output animation is saved. \
This method also works if you would like to rerun the simulation with different flow variable, but on the same mesh or render the same velocity gradient differently without rerunning all parts. Snakemake will avoid rerunning the aspects that do not need to be ran again.
 
## Trial Run
Go ahead and try running the workflow with the provided geometry file:
```
aorta.geo
```
and hopefully you get the following output: 
