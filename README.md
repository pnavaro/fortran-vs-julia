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

## Julia version

Test your [MPI.jl](https://juliaparallel.github.io/MPI.jl/stable/installation/) installation with 

```
$ mpirun --oversubscribe -np 4 julia --project hello_mpi.jl
Hello world, I am 0 of 4
Hello world, I am 3 of 4
Hello world, I am 1 of 4
Hello world, I am 2 of 4
```

Run the maxwell simulation

```
mpirun -np 4 julia --project main.jl
```

![](bz_field.gif)
