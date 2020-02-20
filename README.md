# Maxwell parallel solver in 2D

We use the Yee numerical scheme FDTD: [Finite-Difference Time-Domain method](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) and MPI topology

## Build and run
```
make
mpirun -np 4 ./maxyee_par 
```
You can modify `Makefile` and `input_data`.

NB: A sequential version is also built named `maxyee`

## Plot the magnetic field

```
gnuplot Bz.gnu
```

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

# Performances

1000 x 1000 on 4 processors 

- `gfortran -O3` : 22 seconds 
- `gfortran -O3` serial : 40 seconds 
- `julia -O3 --check-bounds=no` serial : 16 seconds
- `julia -O3 --check-bounds=no MPI` : 5 seconds
