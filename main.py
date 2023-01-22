from oslo_lattice_class import *

import matplotlib.pyplot as plt

def L_scaling_test(L_space, method, N_m=1, *args):
    # L_space = array of L values to test on
    # method = method of the oslo lattice class that returns a measurement value
    # N_m = number of repetitions per L value datapoint
    # *args = arguments to pass to oslo_lattice.method
    res_array = []
    print(f"--------- Commencing testing on {len(L_space)} instances --------")
    for L_val in L_space:
        print(f"Current L = {L_val}")
        cur_lattice = oslo_lattice(L_val)
        cur_res_sum = method(cur_lattice, *args)
        for i in range(N_m - 1):
            cur_lattice = oslo_lattice(L_val)
            cur_res_sum += method(cur_lattice, *args)
        if N_m > 1:
            cur_res_sum /= N_m
        print(f"  Result = {cur_res_sum}")
        res_array.append(cur_res_sum)
    return(res_array)


# --------------------------------------------------------
# ------------------- TASK functions ---------------------
# --------------------------------------------------------

# ----------- TASK 1 --------------

def task1_matching_seed():
    kim = oslo_lattice(10)
    my_h, my_z_th = kim.h.copy(), kim.z_th.copy()

    print(kim)

    random.seed(11)
    kim.simulate(100, 'never')
    print(kim)

    for i in range(10):
        random.seed(10)
        cur_lat = oslo_lattice(10, 0.5, my_h, my_z_th)
        cur_lat.simulate(100, 'never', 'naive')
        print(cur_lat)

def task1_benchmark_test():
    import time

    naive_start = time.time()
    kim = oslo_lattice(400)
    kim.simulate(1e3, 'never', 'naive')
    print("Naive: --- %s seconds ---" % (time.time() - naive_start))

    ceilidh_start = time.time()
    kim = oslo_lattice(400)
    kim.simulate(1e3, 'never', 'ceilidh')
    print("Ceilidh: --- %s seconds ---" % (time.time() - ceilidh_start))


# ----------- TASK 2 --------------

def task2a():
    L_space = [4 , 8  , 16 , 32  , 64 , 128, 256 , 512 ]
    t_space = [50, 200, 600, 2500, 1e4, 4e4, 16e4, 64e4]
    x_plotsize, y_plotsize = subplot_dimensions(len(L_space))
    for n in range(len(L_space)):
        plt.subplot(x_plotsize, y_plotsize, n+1)
        print(f"Simulating lattice with L = {L_space[n]} until t = {int(t_space[n])}")
        cur_lattice = oslo_lattice(L_space[n])
        cur_evolution = cur_lattice.simulate(t_space[n])
        plt.title(f"L = {L_space[n]}")
        plt.ticklabel_format(style='sci', scilimits=(-3, 3))
        plt.xlabel('t')
        plt.ylabel('h')
        plt.axhline(y=L_space[n], linestyle='dotted', color='red', label='bound.')
        plt.axhline(y=2*L_space[n], linestyle='dotted', color='red')
        plt.plot(cur_evolution, label='h(t;L)')
        plt.legend()
    plt.tight_layout()
    plt.show()

def task2b():

    L_space = [4 , 8  , 16 , 32  , 64 , 128]
    res_space = L_scaling_test(L_space, oslo_lattice.get_t_c, 10)
    plt.title("$\\langle t_c\\rangle(L)$ (10 measurements per datapoint)")
    plt.xlabel('L')
    plt.ylabel('$\\langle t_c\\rangle$')
    plt.plot(L_space, res_space)
    plt.show()


#task2b()


# Proof check
lulu = oslo_lattice(32)
trotel, _ = lulu.simulate(1000)
print(trotel)

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



