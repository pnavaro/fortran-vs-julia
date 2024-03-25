function faraday!(fields, ix, jx, iy, jy, dt)

    dx, dy = fields.mesh.dx, fields.mesh.dy

    for j = iy:jy, i = ix:jx
        dex_dy = (fields.ex[i, j+1] - fields.ex[i, j]) / dy
        dey_dx = (fields.ey[i+1, j] - fields.ey[i, j]) / dx
        fields.bz[i, j] = fields.bz[i, j] + dt * (dex_dy - dey_dx)
    end

end

function ampere_maxwell!(fields, ix, jx, iy, jy, dt)

    dx, dy = fields.mesh.dx, fields.mesh.dy

    for j = iy+1:jy, i = ix:jx
        dbz_dy = (fields.bz[i, j] - fields.bz[i, j-1]) / dy
        fields.ex[i, j] = fields.ex[i, j] + dt * csq * dbz_dy
    end

    for j = iy:jy, i = ix+1:jx
        dbz_dx = (fields.bz[i, j] - fields.bz[i-1, j]) / dx
        fields.ey[i, j] = fields.ey[i, j] - dt * csq * dbz_dx
    end

end

function periodic_bc!(fields, ix, jx, iy, jy, dt)

    for i = ix:jx
        dbz_dy = (fields.bz[i, iy] - fields.bz[i, jy]) / fields.mesh.dy
        fields.ex[i, iy] = fields.ex[i, iy] + dt * csq * dbz_dy
        fields.ex[i, jy+1] = fields.ex[i, iy]
    end

    for j = iy:jy
        dbz_dx = (fields.bz[ix, j] - fields.bz[jx, j]) / fields.mesh.dx
        fields.ey[ix, j] = fields.ey[ix, j] - dt * csq * dbz_dx
        fields.ey[jx+1, j] = fields.ey[ix, j]
    end

    fields.bz[:, jy+1] .= fields.bz[:, iy]
    fields.bz[jx+1, :] .= fields.bz[ix, :]

end
