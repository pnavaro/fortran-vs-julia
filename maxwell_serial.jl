const c = 1.0 # speed of light 
const csq = c * c

include("mesh.jl")

struct MeshFields

   mesh :: Mesh
   ex :: Array{Float64, 2}
   ey :: Array{Float64, 2}
   bz :: Array{Float64, 2}

    function MeshFields( mesh )

        mx, my = mesh.mx, mesh.my
 
        ex = zeros(Float64, (mx,my+1))
        ey = zeros(Float64, (mx+1,my))
        bz = zeros(Float64, (mx,my))
        new( mesh, ex, ey, bz)
    end

end 

include("fdtd.jl")

function main( nstep )


    cfl    = 0.2 	# Courant-Friedrich-Levy
    tfinal = 1.	    # final time
    nstepmax = 500  # max steps
    md = 2		    # md : wave number x (initial condition)
    nd = 2		    # nd : wave number y (initial condition)
    mx = 1000	    # x number of points
    my = 1000	    # y number of points
    dimx = 1.0	    # width
    dimy = 1.0	    # height

    dx = dimx / mx
    dy = dimy / my

    mesh = Mesh( dimx, mx, dx, dimy, my, dy)
    
    dt = cfl / sqrt(1/dx^2+1/dy^2) / c
    
    nstep  = min( nstepmax, nstep)
    
    fields = MeshFields(mesh)
    
    omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)
    for j=1:my, i=1:mx
        fields.bz[i,j] = (- cos(md*pi*((i-0.5)*dx/dimx)) 
                      * cos(nd*pi*((j-0.5)*dy/dimy))
                      * cos(omega*(-0.5*dt)) )
    end  
    
    tag = 1111
    
    for istep = 1:nstep #*** Loop over time
    
       # E(n) [1:mx]*[1:my] --> B(n+1/2) [1:mx-1]*[1:my-1]
       faraday!(fields, 1, mx, 1, my, dt)   
    
       periodic_bc!(fields, 1, mx, 1, my, dt)

       # Bz(n+1/2) [1:mx]*[1:my] --> Ex(n+1) [1:mx]*[2:my]
       # Bz(n+1/2) [1:mx]*[1:my] --> Ey(n+1) [2:mx]*[1:my]
       ampere_maxwell!(fields, 1, mx, 1, my, dt) 
    
       err_l2 = 0.0
       time = (istep-0.5)*dt
       for j = 1:my, i = 1:mx
           th_bz = (- cos(md*pi*((i-0.5)*dx/dimx))
                    * cos(nd*pi*((j-0.5)*dy/dimy))
                    * cos(omega*time))
           err_l2 += (fields.bz[i,j] - th_bz)^2
       end
    
       println(err_l2)
    
    end # next time step
    
end

main( 1 ) # trigger building

@time main( 500 )
