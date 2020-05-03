using Hwloc
Hwloc.num_physical_cores()

using Distributed

nprocs()

nworkers()

addprocs(4)

workers()

@fetch begin
    println(myid());
    rand(2,2)
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



