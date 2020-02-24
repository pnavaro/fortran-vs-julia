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

Uncomment the plot_fields call in Julia programs or change idiag value in input_data for fortran.

```
gnuplot Bz.gnu
```
![](bz_field.gif)

## Julia version

Test your [MPI.jl](https://juliaparallel.github.io/MPI.jl/stable/installation/) installation with 

```
$ mpirun -np 4 julia --project hello_mpi.jl
Hello world, I am 0 of 4
Hello world, I am 3 of 4
Hello world, I am 1 of 4
Hello world, I am 2 of 4
```


# Performances (without disk IO)

On small program like this Julia is really fast.

## Serial computation

### 1200 x 1200 and 1000 iterations.

- `julia -O3 --check-bounds=no maxwell_serial.jl` : 14 seconds
- `make && time ./maxwell_serial_fortran` : 31 seconds 

### 1200 x 1200 on 9 processors and 1000 iterations

- `make && time mpirun -np 9 ./maxwell_mpi_fortran` : 7 seconds 
- `mpirun -np 9 julia --project -O3 --check-bounds=no ` : 5 seconds
