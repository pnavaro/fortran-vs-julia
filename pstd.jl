using LinearAlgebra
using FFTW


struct PSTD

    kx 
    ky 
    dx :: Array{ComplexF64, 2}
    dy :: Array{ComplexF64, 2}

   function PSTD( mesh )

        nx, ny = mesh.nx, mesh.ny

        kx = 2π ./ mesh.dimx .* vcat(0:nx÷2-1,-nx÷2:-1)
        ky = 2π ./ mesh.dimy .* vcat(0:ny÷2-1,-ny÷2:-1)

        dx = zeros(ComplexF64, (nx, ny))
        dy = zeros(ComplexF64, (ny, nx))

        new( kx, ky, dx, dy)

    end

end

function faraday!( bz, pstd, ex, ey, dt )

    transpose!(pstd.dy, ex)
    fft!(pstd.dy, 1)
    pstd.dy .*= 1im .* pstd.ky
    ifft!(pstd.dy, 1)

    pstd.dx .= ey
    fft!(pstd.dx, 1)
    pstd.dx .*= 1im .* pstd.kx
    ifft!(pstd.dx, 1)

    nx, ny = size(bz)

    for j in 1:ny, i in 1:nx
        @inbounds bz[i,j] += dt * (pstd.dy[j,i].re - pstd.dx[i,j].re)
    end

end

function ampere_maxwell!( ex, ey, pstd, bz, dt )

    transpose!(pstd.dy, bz)
    fft!(pstd.dy, 1)
    pstd.dy .*= 1im .* pstd.ky
    ifft!(pstd.dy, 1)
    
    nx, ny = size(bz)

    for j in 1:ny, i in 1:nx
        @inbounds ex[i,j] += dt * pstd.dy[j,i].re
    end

    pstd.dx .= bz
    fft!(pstd.dx, 1)
    pstd.dx .*= 1im .* pstd.kx
    ifft!(pstd.dx, 1)

    for i in eachindex(ey)
        @inbounds ey[i] -= dt * pstd.dx[i].re
    end

end 

