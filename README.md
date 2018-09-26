# Maxwell parallel solver in 2D

We use the Yee numerical scheme FDTD: [Finite-Difference Time-Domain method](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) and MPI topology

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
