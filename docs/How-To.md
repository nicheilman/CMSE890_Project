## How-To 

# FAQ
Q: My solution values are blowing up. \
A: Your discretization might be too fine or coarse for your application. Try to keep dx/dt approximately equal to the maximum velocity that you expect to see. \
 \
Q: How do I change the animation recording window? \
A: This is done within the flow solver. At the end of the linear system solving loop, you can edit when and how often the snapshots are recorded. \
 \
Q: I do not have enough RAM to handle the big memory requirements for my simulation. \
A: Try using the provided batch file on an HPC that utilizes SLURM scheduling. 
```
sbatch batch.sb
```

