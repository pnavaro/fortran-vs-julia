using Hwloc
Hwloc.num_physical_cores()

using Distributed

nprocs()

nworkers()

addprocs(4)

workers()

@fetch begin
    println(myid())
    rand(2, 2)
end

@sync begin
    pids = workers()
    @spawnat pids[1] (sleep(2); println("Today is reverse day!"))
    @spawnat pids[2] (sleep(1); println(" class!"))
    @spawnat pids[3] println("Hello")
end;
println("Done!")

@everywhere begin # execute this block on all workers
    using Random

    function complicated_calculation()
        sleep(1)
        randexp(5) # lives in Random
    end
end

@fetch complicated_calculation()


ch = Channel{Int}(5)

isready(ch)


put!(ch, 3)

isready(ch)

take!(ch)

put!(ch, 4)

fetch(ch)

take!(ch)

isready(ch)


const mychannel = RemoteChannel(() -> Channel{Int}(10), workers()[2])

function whohas(s::String)
    @everywhere begin
        var = Symbol($s)
        if isdefined(Main, var)
            println("$var exists.")
        else
            println("Doesn't exist.")
        end
    end
    return nothing
end

whohas("mychannel")


@everywhere const mychannel = $mychannel

# +
whohas("mychannel")


# +
function do_something()
    rc = RemoteChannel(() -> Channel{Int}(10)) # lives on the master
    @sync for p in workers()
        @spawnat p put!(rc, myid())
    end
    return rc
end

r = do_something()
# -

using Distributed, BenchmarkTools;
rmprocs(workers());
addprocs(4);
nworkers();

# +
# serial version - count heads in a series of coin tosses
function add_serial(n)
    c = 0
    for i in 1:n
        c += rand(Bool)
    end
    return c
end

@btime add_serial(200_000_000);

# +
# distributed version
function add_distributed(n)
    c = @distributed (+) for i in 1:n
        Int(rand(Bool))
    end
    return c
end

@btime add_distributed(200_000_000);

# +
# verbose distributed version
function add_distributed(n)
    c = @distributed (+) for i in 1:n
        x = Int(rand(Bool))
        println(x)
        x
    end
    return c
end

add_distributed(8);
# -

@everywhere using SharedArrays # must be loaded everywhere

A = rand(2, 3)

S = SharedArray(A)

# +
function fill_shared_problematic(N)
    S = SharedMatrix{Int64}(N, N)
    @sync @distributed for i in 1:length(S) # added @sync here
        S[i] = i
    end
    return S
end

S = fill_shared_problematic(100)
minimum(S)
# -


minimum(S)
