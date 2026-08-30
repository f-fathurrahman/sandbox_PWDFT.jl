function integrate_inward!(grid, f, y, c, el, idx_r)
     # prepare inward integration
     # charlotte froese can j phys 41,1895(1963)
     #
     # start at min( rmax, 10*rmatch )
     #
     Nrmesh = grid.Nrmesh
     ir_start = Nrmesh     
     rstart = 10.0*grid.r[idx_r]
     if rstart < grid.r[Nrmesh]
          for ir in idx_r:Nrmesh
               ir_start = ir
               if grid.r[ir] >= rstart
                    break
               end
          end
          # make ir_start odd number
          if ir_start%2 == 0
               ir_start += 1
          end
     end
     #
     #  set up a, l, and c vectors
     #
     ir = idx_r+1
     el[ir] = 10.0*f[ir] - 12.0
     c[ir] = -f[idx_r] * y[idx_r]
     for ir in (idx_r+2):ir_start
          di = 10.0*f[ir] - 12.0
          el[ir] = di - f[ir]*f[ir-1]/el[ir-1]
          c[ir] = -c[ir-1]*f[ir-1]/el[ir-1]
     end
     #
     # start inward integration by the froese's tail procedure
     ir = ir_start - 1
     expn = exp( -sqrt( 12*abs(1.0 - f[ir]) ) )
     y[ir] = c[ir]/( el[ir] + f[ir_start]*expn )
     y[ir_start] = expn*y[ir]
     #
     # and integrate inward
     for ir in range(ir_start-2, stop = idx_r + 1, step = -1)
          y[ir] = ( c[ir] - f[ir+1]*y[ir+1] )/el[ir]
     end
     return ir_start
end

