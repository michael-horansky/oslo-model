# -------------------------------------------------------------
#  This file provides the class oslo_lattice, which implements
#  the main features of the lattice algorithm.
#
#  Created by Michal Horansky on 18/01/2023
#  Imperial College London
# -------------------------------------------------------------


import numpy as np
import random


class oslo_lattice():
    
    # -------------- Helping functions ---------------
    
    def get_random_slope_threshold(self):
        # z_th = 1 if p, = 2 if 1-p
        if random.random() < self.p:
            return(1)
        else:
            return(2)
    
    # --------- Constructor and descriptors ----------
    
    def reset(self):
        self.h = [0] * self.L
        self.evolution_array = [0]
        self.z_th = []
        for i in range(self.L):
            self.z_th.append(self.get_random_slope_threshold())
        self.t_c = None
    
    def __init__(self, L, p=0.5, heights = False, threshold_slopes = False):
        
        self.L = L
        self.p = p
        if heights == False:
            self.h = [0] * self.L
        else:
            self.h = heights.copy()
        self.evolution_array = [self.h[0]]
        
        if threshold_slopes == False:
            self.z_th = []
            for i in range(self.L):
                self.z_th.append(self.get_random_slope_threshold())
        else:
            self.z_th = threshold_slopes.copy()
        self.t_c = None
    
    def print_lattice_state(self):
        return(f"heights = [{', '.join(str(item) for item in self.h)}]\nslopes  = [{', '.join(str(item) for item in self.z_array())}]\ns thres = [{', '.join(str(item) for item in self.z_th)}]")
    
    def __repr__(self):
        return(f"L={self.L},p={self.p},h=[{','.join(str(item) for item in self.h)}],z_th=[{','.join(str(item) for item in self.z_th)}]")
    def __str__(self):
        return(f"Oslo-model lattice; L = {self.L}, p = {self.p}\n{self.print_lattice_state()}")
    
    # ------------- Iterative functions --------------
    
    def z(self, index):
        # returns the slope at index i
        # we could save z as an array but since it contains the same information as h[]
        # this reduces the need for constant syncing of the two arrays
        if index == self.L - 1:
            return(self.h[index])
        return(self.h[index] - self.h[index+1])
    
    def z_array(self):
        # returns an array of z vals. To be used only for descriptors and such!
        slopes = []
        for i in range(self.L):
            slopes.append(self.z(i))
        return(slopes)
    
    def drive(self):
        self.h[0] += 1
        
    def relaxation(self):
        # returns a tuple (size_of_avalanche, did_crossover_occur)
        s = 0
        did_crossover_occur = False
        # The "ceilidh algorithm" is implemented here
        cur_i = 0
        max_i_affected = cur_i
        
        while(cur_i < self.L):
            # check if site i is supercritical
            if self.z(cur_i) > self.z_th[cur_i]:
                # increment avalanche size
                s += 1
                # relax current site
                self.h[cur_i] -= 1
                if cur_i != self.L - 1:
                    self.h[cur_i + 1] += 1
                else:
                    did_crossover_occur = True
                # Change the threshold slope of the relaxed site
                self.z_th[cur_i] = self.get_random_slope_threshold()
                # sites up to cur_i + 1 can be supercritical now. Update max_i_affected if needed
                if max_i_affected < cur_i + 1:
                    max_i_affected = cur_i + 1
                # step back
                cur_i = max(cur_i - 1, 0)
            else:
                # step forward
                cur_i += 1
                # if cur_i > max_i_affected, we are outside of the range of possible supercriticality. exit.
                if cur_i > max_i_affected:
                    return(s, did_crossover_occur)
        return(s, did_crossover_occur)
        
    def relaxation_naive(self):
        s = 0
        can_exit = False
        while(not can_exit):
            can_exit = True
            for i in range(self.L):
                # check if avalanche continues
                if self.z(i) > self.z_th[i]:
                    # avalanche continues. Increment the avalanche size
                    s += 1
                    # Relax site i
                    self.h[i] -= 1
                    if i != self.L - 1:
                        self.h[i + 1] += 1
                    # Change the threshold slope of the relaxed site
                    self.z_th[i] = self.get_random_slope_threshold()
                    # you have to do the loop again after this one finishes
                    can_exit = False
        return(s)
    
    def step(self, step_alg = 'ceilidh'):
        # Wrapper function for all the steps in one iteration of the algorithm
        # returns the avalanche size s associated with this step and whether a crossover occured on that step
        self.drive()
        if step_alg == 'ceilidh':
            s, did_crossover_occur = self.relaxation()
        if step_alg == 'naive':
            s, did_crossover_occur = self.relaxation_naive()
        return(s, did_crossover_occur)
    
    # ------------ Operative functions -------------
    
    def simulate(self, N_steps = 100, print_mode = 'never', step_alg = 'ceilidh', terminate_on_crossover = False):
        # The main loop for a basic simulation
        # updates the hiegh of the pile as an array
        t = 0
        printed_last_it = False
        while(t < N_steps):
            # Check whether to print in anticipation of avalanche
            if print_mode == 'avalanche':
                if self.z(0) == self.z_th[0]:
                    if printed_last_it == False:
                        print(f"------- t = {t} -------------")
                        print(self.print_lattice_state())
                    printed_last_it = True
                else:
                    printed_last_it = False
            cur_s, did_crossover_occur = self.step(step_alg)
            t += 1
            self.evolution_array.append(self.h[0])
            if print_mode == 'avalanche' and printed_last_it:
                print(f"------- t = {t}, s = {cur_s} --------")
                print(self.print_lattice_state())
            if print_mode == 'always':
                print(f"------- t = {t} -------------")
                print(self.print_lattice_state())
            
            if did_crossover_occur:
                # save the crossover time
                self.t_c = t
                if terminate_on_crossover:
                    if print_mode in ['always', 'avalanche', 'events']:
                        print(f"Simulation terminated early on t = {t} due to cross-over occuring.")
                    break
    
    def average_pile_height(self, N_steps = 100, N_repetitions = 10):
        sum_pile_height = 0
        for rep in range(N_repetitions):
            self.reset()
            self.simulate(N_steps, 'never')
            sum_pile_height += self.h[0]
        return(sum_pile_height / N_repetitions)
    
    def average_t_c(self, N_steps = 100, N_repetitions = 10):
        sum_t_c = 0
        for rep in range(N_repetitions):
            self.reset()
            self.simulate(N_steps, terminate_on_crossover = True)
            sum_t_c += self.t_c
        return(sum_t_c / N_repetitions)


