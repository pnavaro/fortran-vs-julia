function plot_fields(mesh, rank, proc, field, xp, yp, iplot )

    dx, dy = mesh.dx, mesh.dy
    ix, jx = 1, mesh.mx
    iy, jy = 1, mesh.my

    if iplot == 1
        mkpath("data/$rank")
    end

    io = open("data/$(rank)/$(iplot)", "w")
    for j=iy:jy
        for i=ix:jx
            @printf( io, "%f %f %f \n", xp+(i-0.5)*dx, yp+(j-1)*dy, field[i,j])
        end
        @printf( io, "\n")
    end
    close(io)
   
    # write master file

    if rank == 0

      if iplot == 1 
         io = open( "field.gnu", "w" )
         write(io, "set xr[-0.1:1.1]\n")
         write(io, "set yr[-0.1:1.1]\n")
         write(io, "set zr[-1.1:1.1]\n")
         write(io, "set surf\n")
      else
         io = open( "field.gnu", "a" )
      end
      write(io, "set title '$(iplot)' \n")
      write(io, "splot 'data/$(rank)/$(iplot)' u 1:2:3 w lines")
   
      for p = 1:proc-1
         write(io, ", 'data/$(p)/$(iplot)' u 1:2:3 w lines")
      end
      write( io, "\n")
      write( io, "set term gif \n")
      write( io, "set output 'image$(lpad(iplot,3,"0")).gif'\n")
      write( io, "replot\n")

      close(io)

    end

end 

