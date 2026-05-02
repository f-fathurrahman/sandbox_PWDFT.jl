from string import Template
import subprocess

# Define the template
inp_template = Template("""
tasks
  0

! use LAPW nxoapwlo=1
!nxoapwlo
! 1

rgkmax
 $rgkmax

! Libxc Slater exchange + PW92 correlation
xctype
 100 1 12

tforce
 .true.

! use Broyden mixing
mixtype
  3

avec
  5.13  5.13  0.00
  5.13  0.00  5.13
  0.00  5.13  5.13

! this is the relative path to the species files
sppath
  './'

atoms
  1                            : nspecies
  'Si.in'                     : spfname
  2                            : natoms; atposl below
  0.0  0.0  0.0
  0.25 0.25 0.25

ngridk
  2 2 2

vkloff
  0.0  0.0  0.0
""")


import os
etot_list = []
rgkmax_list = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
for rgkmax in rgkmax_list:
    #
    inp = inp_template.substitute({"rgkmax": rgkmax})
    with open("elk.in", "w") as file:
        file.write(inp)
    subprocess.run(["elk-6.3.2-gfortran"], capture_output=False, text=False)
    subprocess.run(["mv", "INFO.OUT", f"INFO_rgkmax_{rgkmax}.OUT"], capture_output=False, text=False)
    #
    #
    #outstr = subprocess.getoutput(f"grep 'total energy' INFO_rgkmax_{rgkmax}.OUT | tail -1")
    #etot_list.append(float(outstr.split()[-1]))

#import matplotlib.pyplot as plt
#plt.plot(rgkmax_list, etot_list)
#plt.grid(True)
#plt.savefig("IMG_rgkmax.png", dpi=150)
#plt.show()