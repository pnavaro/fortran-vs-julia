# maxwell_yee_2d_with_mpi
Maxwell Solver in 2D using Yee numerical scheme and MPI topology


## Build

```
make
````

## Run

```
mpirun -np 4 ./maxyee_par 
```

## Plot the magnetic field

```
gnuplot Bz.gnu
```

##  Note

A sequential version is also built named `maxyee`
