using FFTW

function faraday!( bz, mesh :: Mesh, ex, ey, dt )

    nx, ny = mesh.nx, mesh.ny
    lx, ly = mesh.dimx, mesh.dimy
    kx = 2π ./ lx .* vcat(0:nx÷2-1,-nx÷2:-1)
    ky = 2π ./ ly .* vcat(0:ny÷2-1,-ny÷2:-1)

    dex_dy = ifft( 1im .* ky .* fft( ex, 2), 2)
    dey_dx = ifft( 1im .* kx .* fft( ey, 1), 1)
    bz .= bz .+ dt .* (dex_dy .- dey_dx)

end

function ampere_maxwell!( ex, ey, mesh :: Mesh, bz, dt )

    nx, ny = mesh.nx, mesh.ny
    lx, ly = mesh.dimx, mesh.dimy
    kx = 2π ./ lx .* vcat(0:nx÷2-1,-nx÷2:-1)
    ky = 2π ./ ly .* vcat(0:ny÷2-1,-ny÷2:-1)

    dbz_dy = ifft( 1im .* ky .* fft( bz, 2), 2)
    ex .= ex .+ dt .* csq .* dbz_dy 

    dbz_dx = ifft( 1im .* kx .* fft( bz, 1), 1)
    ey .= ey .- dt .* csq .* dbz_dx 

end 
