from oslo_lattice_class import *



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

kim = oslo_lattice(16)
chris = oslo_lattice(32)
print(kim.aggregate_measurement(500, 10))
print(chris.aggregate_measurement(5000, 10))

