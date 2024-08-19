import numpy as np
import matplotlib.pyplot as plt

f1i = np.reshape( np.loadtxt("fort.102"), (4,273) )
f2i = np.reshape( np.loadtxt("fort.103"), (4,273) )

f1o = np.reshape( np.loadtxt("fort.104"), (49,124) )
f2o = np.reshape( np.loadtxt("fort.105"), (49,124) )

import matplotlib.pyplot as pl
f1is = 