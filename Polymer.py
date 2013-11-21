__author__ = 'Ronen Abravanel'

import numpy as np
import random


class Polymer:


    def __init__(self, D, L, IS_SELF_AVOIDING=False, reflect=False):
        """
        Initialize Polymer constructor..
        @param D: The dimensionality of the polymer.
        @param L: The length of each dimension of the Polymer's lattice.
        @param IS_SELF_AVOIDING: Bool, Toggle creation as simple Random Walk or as Self-avoiding walk (by loop-removal)
        @param reflect: Bool, Use reflecting boundary conditions.
        """
        self.D = D
        self.L = L
        self.IS_SELF_AVOIDING = IS_SELF_AVOIDING
        self.reflect = reflect

        self.dims = np.array([L] * D, dtype=np.int)


        self.start_point = self.dims/2

        self.zero_step = np.array([0] * D)

        #Initilize list of all possible directions.
        self.possible_steps = [np.array(a,dtype=np.int) for a in zip(*[(0,0)*i + (1,-1) + (0,0) * (D-1-i) for i in range(D)])]

    def create_polymer(self, steps):
        """
        Creates one realization of a polymer with 'steps' monomers.
        @param steps: The number of monomers
        @return: (final_point, step_number, grid)
         final point - Cartesian coordinates of R_N, the last monomer.
         step_number - The actual number of monomers.
         grid - A grid contains the current realization.
        """
        if steps>2**16:
            print "Too many steps"
            return

        grid = np.zeros(self.dims, dtype=np.int16)

        prev_step = self.zero_step
        cur_point = self.start_point

        #A while loop, as the number of steps can varies in SAW.
        i = 0
        while i<steps:
            i += 1

            #Choose Random direction.
            step = random.choice(self.possible_steps)

            #Make sure we do not go backwards.  (Should be only in SAW?)
            while (prev_step + step == self.zero_step).all():
                step = random.choice(self.possible_steps)

            cur_point = cur_point + step

            #Boundary check.
            if (cur_point >= self.L).any() or (cur_point<0).any():
                if self.reflect:
                    #In order to reflect, revert the last step, and make another step backwards.
                    cur_point = cur_point - step*2
                else:
                    print "curdir = ", cur_point , " brakeing. "
                    print "Did only %d steps!!!!" % i
                    break

            prev_step = step
            point_tup = tuple(cur_point)

            #IF SAW, check for loop and remove them.
            if self.IS_SELF_AVOIDING:
                if grid[point_tup] != 0:
                    start_collision = grid[point_tup]
                    end_collision = i

                    #Clean leftovers.
                    #print "clean from  %d to %d " % (start_collision, end_collision)
                    grid[np.logical_and(grid>start_collision , grid<end_collision) ] = 0
                    i = start_collision

            grid[point_tup] = i

        return cur_point, i, grid