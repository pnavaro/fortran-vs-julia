using MPI

function main()

    # Initialize the MPI environment
    MPI.Init()

    tag = 1111
    comm = MPI.COMM_WORLD

    # Get the number of processes
    world_size = MPI.Comm_size(comm)

    # Get the rank of the process
    world_rank = MPI.Comm_rank(comm)

    # Print off a hello world message
    println("Processor $world_rank out of $world_size processors")

    next_proc = mod(world_rank + 1, world_size)
    last_proc = mod(world_size + world_rank - 1, world_size)

    if world_rank == 0
        token = 10 # Set the token's value if you are process 0
        MPI.Send(token, next_proc, tag, comm)
        (token, _) = MPI.Recv(Int, last_proc, tag, comm)
    else
        (token, _) = MPI.Recv(Int, last_proc, tag, comm)
        MPI.Send(token + 1, next_proc, tag, comm)
    end

    println("Process $world_rank received token $token from process $(last_proc)")

    # Finalize the MPI environment.
    MPI.Finalize()

end

main()
