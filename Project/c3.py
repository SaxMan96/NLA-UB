import numpy as np
from c2 import c2
import time
import matplotlib.pyplot as plt
#%matplotlib inline
def c3():
    n_list = []
    t_list = []
    for n in range(2,100):
        n_list.append(n)
        t0 = time.clock()
        c2(n)
        t_list.append(time.clock()-t0)
    return n_list, t_list


n_list, t_list  = c3()

plt.plot(n_list,t_list)
plt.show()