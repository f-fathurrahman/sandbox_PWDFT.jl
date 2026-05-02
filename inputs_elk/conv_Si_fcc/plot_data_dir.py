import subprocess
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.style.use("dark_background")

# Global var
rgkmax_list = np.array([5.0, 5.5, 6.0, 6.5, 7.0, 7.5])

def process_dir(dir_name):
    etot_list = []
    for rgkmax in rgkmax_list:
        outstr = subprocess.getoutput(f"grep 'total energy' {dir_name}/INFO_rgkmax_{rgkmax}.OUT | tail -1")
        etot_list.append(float(outstr.split()[-1]))
    return np.array(etot_list)


for d in ["01", "02", "03", "04"]:
    etot_list = process_dir("TEMP_data_" + d)
    plt.plot(rgkmax_list, etot_list, label=d, marker="o")
    #diff_e = np.abs(etot_list - etot_list[-1])
    #plt.semilogy(rgkmax_list, diff_e, label=d, marker="o")

plt.grid(True)
plt.legend()
plt.savefig("IMG_rgkmax.png", dpi=150)
plt.show()
