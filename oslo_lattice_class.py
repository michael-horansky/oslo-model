# -------------------------------------------------------------
#  This file provides the class oslo_lattice, which implements
#  the main features of the lattice algorithm.
#
#  Created by Michal Horansky on 18/01/2023
#  Imperial College London
# -------------------------------------------------------------


import numpy
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
    
    def __init__(self, L, p=0.5, heights = False, threshold_slopes = False):
        
        self.L = L
        self.p = p
        if heights == False:
            self.h = [0] * self.L
        else:
            self.h = heights.copy()
        
        if threshold_slopes == False:
            self.z_th = []
            for i in range(self.L):
                self.z_th.append(self.get_random_slope_threshold())
        else:
            self.z_th = threshold_slopes.copy()
    
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
        s = 0
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
            else:
                # If site i wasn't relaxed, sites with j>i won't be relaxed under
                # any circumstances. Hence we can optimise by exiting the loop early
                break
        return(s)
    
    def step(self):
        # Wrapper function for all the steps in one iteration of the algorithm
        # returns the avalanche size s associated with this step
        self.drive()
        s = self.relaxation()
        return(s)
    
    # ------------ Operative functions -------------
    
    def simulate(self, N_steps = 100, print_mode = 'avalanche'):
        # The main loop for a basic simulation
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
            cur_s = self.step()
            t += 1
            if print_mode == 'avalanche' and printed_last_it:
                print(f"------- t = {t}, s = {cur_s} --------")
                print(self.print_lattice_state())
            if print_mode == 'always':
                print(f"------- t = {t} -------------")
                print(self.print_lattice_state())


