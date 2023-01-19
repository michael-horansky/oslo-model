from oslo_lattice_class import *

import matplotlib.pyplot as plt

def L_scaling_test(L_space, method, *args):
    res_array = []
    print(f"--------- Commencing testing on {len(L_space)} instances --------")
    for L_val in L_space:
        cur_lattice = oslo_lattice(L_val)
        print(f"Current L = {L_val}")
        res = method(cur_lattice, *args)
        print(f"  Result = {res}")
        res_array.append(res)
    return(res_array)


"""
ALGORITHM TIME BENCHMARK TEST

import time

naive_start = time.time()
kim = oslo_lattice(400)
kim.simulate(1e3, 'never', 'naive')
print("Naive: --- %s seconds ---" % (time.time() - naive_start))

ceilidh_start = time.time()
kim = oslo_lattice(400)
kim.simulate(1e3, 'never', 'ceilidh')
print("Ceilidh: --- %s seconds ---" % (time.time() - ceilidh_start))"""

"""
kim = oslo_lattice(10)


my_h, my_z_th = kim.h.copy(), kim.z_th.copy()

print(kim)

random.seed(11)
kim.simulate(100, 'never')
print(kim)

for i in range(1):
    random.seed(10)
    cur_lat = oslo_lattice(10, 0.5, my_h, my_z_th)
    cur_lat.simulate(10, 'never', 'naive')
    print(cur_lat)"""

L_space = np.linspace(10, 100, 10, dtype=int)
res_space = L_scaling_test(L_space, oslo_lattice.average_t_c, 1e9, 10)
plt.plot(L_space, res_space)
plt.show()

"""
kim = oslo_lattice(16)
chris = oslo_lattice(32)
#kim.simulate(1000)
#chris.simulate(5000)
print(kim.average_t_c(5000, 10))
print(chris.average_t_c(5000, 10))

plt.subplot(2, 1, 1)
plt.plot(kim.evolution_array)
plt.subplot(2, 1, 2)
plt.plot(chris.evolution_array)
plt.show()"""



