# Julia Syntax: Comparison with Fortran

This is a simple cheatsheet and some performance comparison for scientific programmers who are interested in discover Julia.
It is not an exhaustive list. This page is inspired from [A Cheatsheet for Fortran 2008 Syntax: Comparison with Python 3](https://github.com/wusunlab/fortran-vs-python/).

<table>
    <tr>
        <td></td>
        <td>Fortran</td>
        <td>Julia</td>
    </tr>
    <tr>
        <td>Top-level constructs</td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>the main program</td>
        <td>
<pre lang="fortran">
program my_program
    ...
end program
</pre>
        </td>
        <td><pre lang="julia">
function my_program()
    ...
end
my_program()
</pre>
Not required but it is recommended to use a function for your main program  </td>
    </tr>
    <tr>
        <td>modules</td>
        <td>
<pre lang="fortran">
module my_module
    ...
end module my_module
</pre>
        </td>
        <td>
      <pre lang="julia">
module MyModule
    ...
end
</pre>
      </td>
    </tr>
    <tr>
        <td>subroutines</td>
        <td>
<pre lang="fortran">
subroutine my_subroutine
    ...
end subroutine my_subroutine
</pre>
        </td>
        <td><pre lang="fortran">
function my_subroutine!
    ...
end
</pre>
      The bang in the function name is a convention if a function mutates one or more of its arguments. The convention is that the modified arguments should (if possible) come first.   
      </td>
    </tr>
    <tr>
        <td>functions</td>
        <td>
<pre lang="fortran">
function f(x) result(res)
    res = ...
end function f
</pre>
        </td>
        <td>
<pre lang="julia">
function my_function(x)
    ...
    return res
end
</pre>
        </td>
    </tr>
    <tr>
        <td>submodules</td>
        <td>
<pre lang="fortran">
module main_module
    ...
contains
&lt;<i>submodule statements</i>&gt;
end module main_module
</pre>
        </td>
        <td><pre lang="julia">
module MainModule
    ...
    module SubModule 
        ...
    end
end
</pre></td>
    </tr>
    <tr>
        <td>import statement</td>
        <td>
<pre lang="fortran">
use my_module
use my_module, only : fun1, var1
</pre>
        </td>
        <td>
<pre lang="julia">
using MyModule
import MyModule: fun1, var1
</pre>
        </td>
    </tr>
    <tr>
        <td>call subroutines and functions</td>
        <td>
<pre lang="fortran">
call my_subroutine(<i>args</i>)
my_function(<i>args</i>)
</pre>
        </td>
        <td>
<pre lang="julia">
my_function(<i>args</i>)
</pre>
        </td>
    </tr>
    <tr>
        <td>abort a program</td>
        <td>
<pre lang="fortran">
stop
</pre>
        </td>
        <td>
<pre lang="julia">
exit()
</pre>
        </td>
    </tr>
    <tr>
        <td>inline comments</td>
        <td>
<pre lang="fortran">
! This is a comment
</pre>
        </td>
        <td>
<pre lang="julia">
# This is a comment
</pre>
        </td>
    </tr>
    <tr>
        <td>include external source files</td>
        <td>
<pre lang="fortran">
include <i>'source_file_name'</i>
</pre>
        </td>
        <td><pre lang="julia">
include(<i>"source_file_name"</i>)
</pre></td>
    </tr>
    <tr>
        <td>Control flow patterns</td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td><code>if</code> construct</td>
        <td>
<pre lang="fortran">
if <i>&lt;logical expr&gt;</i> then
    ...
else if <i>&lt;logical expr&gt;</i> then
    ...
else
    ...
end if
</pre>
        </td>
        <td>
<pre lang="julia">
if <i>&lt;logical expr&gt;</i>
    ...
elseif <i>&lt;logical expr&gt;</i>
    ...
else
    ...
end
</pre>
        </td>
    </tr>
    <tr>
        <td><code>case</code> construct</td>
        <td>
<pre lang="fortran">
select case <i>&lt;expr&gt;</i>
    case <i>&lt;value&gt;</i>
        ...
    case <i>&lt;value&gt;</i>
        ...
    case default
        ...
end select
</pre>
        </td>
        <td>Not supported.
        </td>
    </tr>
    <tr>
        <td><code>do</code> construct</td>
        <td>
<pre lang="fortran">
do <i>i = start_value</i>, <i>end_value</i>, <i>step</i>
    ...
end do
</pre>
        </td>
        <td>
<pre lang="julia">
for i in start:step:end
    ...
end
</pre>
        </td>
    </tr>
    <tr>
        <td><code>do while</code> construct</td>
        <td>
<pre lang="fortran">
do while <i>&lt;logical expr&gt;</i>
    ...
end do
</pre>
        </td>
        <td>
<pre lang="julia">
while <i>&lt;logical expr&gt;</i>
    ...
end
</pre>
        </td>
    </tr>
    <tr>
        <td>break from a loop</td>
        <td>
<pre lang="fortran">
exit
</pre>
        </td>
        <td>
<pre lang="julia">
break
</pre>
        </td>
    </tr>
    <tr>
        <td>leave this iteration and continue to the next iteration</td>
        <td>
<pre lang="fortran">
cycle
</pre>
        </td>
        <td>
<pre lang="julia">
continue
</pre>
        </td>
    </tr>
    <tr>
        <td>Data types</td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>declaration</td>
        <td>
<pre lang="fortran">
integer(kind=kind_number) :: n = 0
real(kind=kind_number) :: x = 0.
</pre>
        </td>
        <td>
<pre lang="julia">
n = 0
n :: Int64 = 0
x = 0.
x :: Float64 = 0.
</pre>
        </td>
    </tr>
    <tr>
        <td>named constants</td>
        <td>
<pre lang="fortran">
integer, parameter :: answer = 42
real(8), parameter :: pi = 4d0 * atan(1d0)
</pre>
        </td>
        <td>
      <pre lang="julia">
const answer = 42
</pre>
            <i>pi</i> is a named constant in Julia standard.
      </td>
    </tr>
    <tr>
        <td>complex number</td>
        <td>
<pre lang="fortran">
complex :: z = (1., -1.)
</pre>
        </td>
        <td>
<pre lang="julia">
z = 1 - 1im
</pre>
        </td>
    </tr>
    <tr>
        <td>string</td>
        <td>
<pre lang="fortran">
character(len=10) :: str_fixed_length
character(len=:), allocatable :: str_var_length
</pre>
        </td>
        <td>
<pre lang="julia">
string = "this is a string"
</pre>
        </td>
    </tr>
    <tr>
        <td>pointer</td>
        <td>
<pre lang="fortran">
real, pointer :: p
real, target :: r
p => r
</pre>
        </td>
        <td>
          <pre lang="julia">
p = Ref(r)
</pre>
</pre>
      </td>
    </tr>
    <tr>
        <td>boolean</td>
        <td>
<pre lang="fortran">
.true.
.false.
</pre>
        </td>
        <td>
<pre lang="julia">
true
false
</pre>
        </td>
    </tr>
    <tr>
        <td>logical operators</td>
        <td>
<pre lang="fortran">
.not.
.and.
.or.
.eqv.
.neqv.
</pre>
        </td>
        <td>
<pre lang="julia">
!
&&
||
</pre>
Other logical operators do not have built-in support.
        </td>
    </tr>
    <tr>
        <td>equal to</td>
        <td>
<pre lang="fortran">
==, .eq.
</pre>
        </td>
        <td>
<pre lang="julia">
==
</pre>
        </td>
    </tr>
    <tr>
        <td>not equal to</td>
        <td>
<pre lang="fortran">
/=, .ne.
</pre>
        </td>
        <td>
<pre lang="julia">
!==
</pre>
        </td>
    </tr>
    <tr>
        <td>greater than</td>
        <td>
<pre lang="fortran">
>, .gt.
</pre>
        </td>
        <td>
<pre lang="julia">
>
</pre>
        </td>
    </tr>
    <tr>
        <td>less than</td>
        <td>
<pre lang="fortran">
<, .lt.
</pre>
        </td>
        <td>
<pre lang="julia">
<
</pre>
        </td>
    </tr>
    <tr>
        <td>greater than or equal to</td>
        <td>
<pre lang="fortran">
>=, .ge.
</pre>
        </td>
        <td>
<pre lang="julia">
>=
</pre>
        </td>
    </tr>
    <tr>
        <td>less than or equal to</td>
        <td>
<pre lang="fortran">
<=, .ge.
</pre>
        </td>
        <td>
<pre lang="julia">
<=
</pre>
        </td>
    </tr>
    <tr>
        <td>array declaration</td>
        <td>
<pre lang="fortran">
real(8), dimension(3) :: a = [1., 2., 3.]
</pre>
        </td>
        <td>
<pre lang="julia">
a = [1., 2., 3.]
</pre>
        </td>
    </tr>
    <tr>
        <td>string array declaration</td>
        <td>
<pre lang="fortran">
character(len=20), dimension(3, 4) :: char_arr
</pre>
        </td>
        <td>
<pre lang="julia">
char_arr = String[]
push!(char_arr, new_string)
</pre>
There is no easy way to preallocate space for strings.
        </td>
    </tr>
    <tr>
        <td>elementwise array operations</td>
        <td>
<pre lang="fortran">
a <i>op</i> b
</pre>
<i><code>op</code></i> can be <code>+, -, *, /, **, =, ==</code>, etc.<br>
This is supported since the Fortran 90 standard.
        </td>
        <td>Supported by using the broadcast operator `.`and the `f.(x)` syntax</td>
    </tr>
    <tr>
        <td>first element</td>
        <td>
<pre lang="fortran">
a(1)
</pre>
        </td>
        <td>
<pre lang="julia">
a[1]
</pre>
        </td>
    </tr>
    <tr>
        <td>slicing</td>
        <td>
<pre lang="fortran">
a(1:5)
</pre>
This slice includes <code>a(5)</code>.
        </td>
        <td>
<pre lang="julia">
a[1:5]
</pre>
This slice includes <code>a[5]</code>.
        </td>
    </tr>
    <tr>
        <td>slicing with steps</td>
        <td>
<pre lang="fortran">
a(1:100:2)
</pre>
        </td>
        <td>
<pre lang="julia">
a[0:2:100]
</pre>
        </td>
    </tr>
    <tr>
        <td>size</td>
        <td>
<pre lang="fortran">
size(a)
</pre>
        </td>
        <td>
<pre lang="julia">
length(a)
</pre>
        </td>
    </tr>
    <tr>
        <td>shape</td>
        <td>
<pre lang="fortran">
shape(a)
</pre>
        </td>
        <td>
<pre lang="julia">
size(a)
</pre>
        </td>
    </tr>
    <tr>
        <td>shape along a dimension</td>
        <td>
<pre lang="fortran">
size(a, dim)
</pre>
        </td>
        <td>
<pre lang="julia">
size(a)[dim]
</pre>
        </td>
    </tr>
    <tr>
        <td>Type conversion</td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>to integer by truncation</td>
        <td>
<pre lang="fortran">
int(x)
</pre>
        </td>
        <td>
<pre lang="julia">
trunc(Int, x )
</pre>
        </td>
    </tr>
    <tr>
        <td>to integer by rounding</td>
        <td>
<pre lang="fortran">
nint()
</pre>
        </td>
        <td>
<pre lang="julia">
round()
</pre>
        </td>
    </tr>
    <tr>
        <td>integer to float</td>
        <td>
<pre lang="fortran">
real(a[, kind])
</pre>
        </td>
        <td>
<pre lang="julia">
float()
</pre>
        </td>
    </tr>
    <tr>
        <td>complex to real</td>
        <td>
<pre lang="fortran">
real(z[, kind])
</pre>
        </td>
        <td>
<pre lang="julia">
real()
</pre>
        </td>
    </tr>
    <tr>
        <td>to complex</td>
        <td>
<pre lang="fortran">
cmplx(x [, y [, kind]])
</pre>
        </td>
        <td>
<pre lang="julia">
complex()
</pre>
        </td>
    </tr>
    <tr>
        <td>to boolean</td>
        <td>
<pre lang="fortran">
logical()
</pre>
        </td>
        <td>
<pre lang="julia">
Bool()
</pre>
        </td>
    </tr>
    <tr>
        <td>Derived data types</td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>definition</td>
        <td>
<pre lang="fortran">
type Point
    real(8) :: x, y
end type Point
</pre>
        </td>
        <td>
<pre lang="julia">
struct Point
    x :: Float64
    y :: Float64
end
</pre>
        </td>
    </tr>
    <tr>
        <td>instantiation</td>
        <td>
<pre lang="fortran">
type(Point) :: point1 = Point(-1., 1.)
</pre>
        </td>
        <td>
<pre lang="julia">
point1 = Point(-1., 1.)
</pre>
        </td>
    </tr>
    <tr>
        <td>get attributes</td>
        <td>
<pre lang="fortran">
point1%x
point1%y
</pre>
        </td>
        <td>
<pre lang="julia">
point1.x
point1.y
</pre>
        </td>
    </tr>
    <tr>
        <td>array of derived type</td>
        <td>
<pre lang="fortran">
type(Point), dimension(:), allocatable :: point_arr
</pre>
        </td>
        <td><pre lang="julia">
point_arr = Vector{Point}
</pre></td>
    </tr>
    <tr>
        <td>type bound procedures (aka class method)</td>
        <td>Assume that <code>Circle</code> has a type bound procedure (subroutine) <code>print_area</code>.
<pre lang="fortran">
type(Circle) :: c
call c%print_area
</pre>
        </td>
        <td>Assume that <code>Circle</code> has a method <code>print_area(c :: Circle)</code>.
<pre lang="python">
c = Circle()
print_area(c)
</pre>
        </td>
    </tr>
    <tr>
    <td>Built-in mathematical functions</td>
    <td></td>
    <td></td>
    </tr>
    <tr>
        <td>functions with the same names</td>
        <td>
<pre lang="fortran">
abs(), cos(), cosh(), exp(), floor(), log(),
log10(), max(), min(), sin(), sinh(), sqrt(),
sum(), tan(), tanh(), acos(), asin(), atan()
</pre>
        </td>
        <td>Have the same name in Julia.
        </td>
    </tr>
    <tr>
        <td>functions with different names</td>
        <td>
<pre lang="fortran">
aimag()
atan2(x, y)
ceiling()
conjg(z)
modulo()
call random_number()
</pre>
        </td>
        <td>
<pre lang="julia">
imag()
atan(x, y)
ceil()
conj()
mod(), %
Random.rand()
</pre>
        </td>
    </tr>
    <tr>
        <td>Built-in string functions</td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>string length</td>
        <td>
<pre lang="fortran">
len()
</pre>
        </td>
        <td>
<pre lang="julia">
length()
</pre>
        </td>
    </tr>
    <tr>
        <td>string to ASCII code</td>
        <td>
<pre lang="fortran">
iachar()
</pre>
        </td>
        <td>
<pre lang="julia">
Int()
</pre>
        </td>
    </tr>
    <tr>
        <td>ASCII code to string</td>
        <td>
<pre lang="fortran">
achar()
</pre>
        </td>
        <td>
<pre lang="julia">
String()
</pre>
        </td>
    </tr>
    <tr>
        <td>string slicing</td>
        <td>Same as 1D array slicing.</td>
        <td>Same as 1D array slicing.</td>
    </tr>
    <tr>
        <td>find the position of a substring</td>
        <td>
<pre lang="fortran">
index(string, substring)
</pre>
        </td>
        <td>
<pre lang="julia">
findfirst(substring, string)
</pre>
        </td>
    </tr>
    <tr>
        <td>string concatenation</td>
        <td>
<pre lang="fortran">
"hello" // "world"
</pre>
        </td>
        <td>
<pre lang="julia">
"hello" * "world"
string("hello", "world")
</pre>
        </td>
    </tr>
    <tr>
        <td>Array constructs</td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td><code>where</code> construct</td>
        <td>
<pre lang="fortran">
where a > 0
    b = 0
elsewhere
    b = 1
end where
</pre>
        </td>
        <td>
<pre lang="julia">
b[a .> 0] .= 0
b[a .<= 0] .= 1
</pre>
function <code>filter()</code> or comprehension can be used also.
        </td>
    </tr>
    <tr>
        <td><code>forall</code> construct</td>
        <td>
<pre lang="fortran">
real, dimension(10, 10) :: a = 0
int :: i, j
...
forall (i = 1:10, j = 1:10, i <= j)
    a(i, j) = i + j
end forall
</pre>
        </td>
        <td>
<pre lang="julia">
a = zeros(Float32, 10, 10)
for i in 1:10, j in 1:10
    if i <= j
        a[i, j] = i + j
    end
end
</pre>
        </td>
    </tr>
    <tr>
        <td>CPU time</td>
        <td>
<pre lang="fortran">
call cpu_time()
</pre>
        </td>
        <td>
<pre lang="julia">
time = @elapsed begin
...
end
</pre>
        </td>
    </tr>
    <tr>
        <td>command line arguments</td>
        <td>
<pre lang="fortran">
call command_argument_count()
call get_command()
call get_command_argument()
</pre>
        </td>
        <td>For basic parsing, use <code>ARGS</code></td>
    </tr>
    <tr>
        <td>Input/output</td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>print</td>
        <td>
<pre lang="fortran">
print fmt, <i>&lt;output list&gt;</i>
</pre>
        </td>
        <td>
<pre lang="julia">
println()
</pre>
        </td>
    </tr>
    <tr>
        <td>read from the command prompt</td>
        <td>
<pre lang="fortran">
read fmt, <i>&lt;input list&gt;</i>
</pre>
        </td>
        <td>
<pre lang="julia">
readline()
</pre>
        </td>
    </tr>
    <tr>
        <td>open a file</td>
        <td>
<pre lang="fortran">
open(unit, file, ...)
</pre>
        </td>
        <td>
<pre lang="julia">
f = open(file, mode='r', ...)
</pre>
        </td>
    </tr>
    <tr>
        <td>read from a file</td>
        <td>
<pre lang="fortran">
read(unit, fmt, ...) <i>&lt;input list&gt;</i>
</pre>
        </td>
        <td>
<pre lang="julia">
read(f)
readlines(f)
</pre>
        </td>
    </tr>
    <tr>
        <td>write to a file</td>
        <td>
<pre lang="fortran">
write(unit, fmt, ...) <i>&lt;output list&gt;</i>
</pre>
        </td>
        <td>
<pre lang="julia">
write(f, ...)
</pre>
        </td>
    </tr>
    <tr>
        <td>close a file</td>
        <td>
<pre lang="fortran">
close(unit, ...)
</pre>
        </td>
        <td>
<pre lang="python">
close(f)
</pre>
        </td>
    </tr>
    <tr>
        <td>file inquiry</td>
        <td>
<pre lang="fortran">
inquire(unit, ...)
</pre>
        </td>
        <td><pre lang="julia">
isfile("my_file.txt")
</pre></td>
    </tr>
    <tr>
        <td>backspace in a file</td>
        <td>
<pre lang="fortran">
backspace(unit, ...)
</pre>
        </td>
        <td>
<pre lang="julia">
skip(f, -1)
</pre>
        </td>
    </tr>
    <tr>
        <td>end of file (EOF)</td>
        <td>
<pre lang="fortran">
endfile(unit, ...)
</pre>
        </td>
        <td>
<pre lang="julia">
seekend(f)
</pre>
        </td>
    </tr>
    <tr>
        <td>return to the start of a file</td>
        <td>
<pre lang="fortran">
rewind(unit, ...)
</pre>
        </td>
        <td>
<pre lang="julia">
seekstart(f)
</pre>
        </td>
    </tr>
</table>


## Maxwell parallel solver in 2D

Here an example of a Fortran to Julia translation. We use the Yee numerical scheme FDTD: [Finite-Difference Time-Domain method](https://en.wikipedia.org/wiki/Finite-difference_time-domain_method) and MPI topology.
You can find a serial version and a parallel version using MPI library.

Test your [MPI.jl](https://juliaparallel.github.io/MPI.jl/stable/installation/) installation with 

```
$ mpirun -np 4 julia --project hello_mpi.jl
Hello world, I am 0 of 4
Hello world, I am 3 of 4
Hello world, I am 1 of 4
Hello world, I am 2 of 4
```
### Performances (without disk IO)

On small program like this in Julia is really fast.

### Serial computation

#### 1200 x 1200 and 1000 iterations.

- `julia -O3 --check-bounds=no maxwell_serial.jl` : 14 seconds
- `make && time ./maxwell_serial_fortran` : 31 seconds 

#### 1200 x 1200 on 9 processors and 1000 iterations

- `make && time mpirun -np 9 ./maxwell_mpi_fortran` : 7 seconds 
- `mpirun -np 9 julia --project -O3 --check-bounds=no ` : 5 seconds

### Plot the magnetic field

Uncomment the plot_fields call in Julia programs or change idiag value in input_data for fortran.

```
gnuplot bz.gnu
```
![](bz_field.gif)
