from Item import Item
from MIPmodelSolver import MIPSolver
from packingVisualisation import draw_packing_mass

import random

if __name__ == '__main__':

    # set number if items
    N = 10
    # set maximum possible deviation of gravity center from the geometrical center of the knapsack
    x_tol = 15
    y_tol = 15
    # generate random items
    items = [Item(random.randint(15, 20), random.randint(15, 20), random.randint(20, 50)) for _ in range(N)]
    # set knapsack parameters
    W = 100
    H = 100
    # set time limit for solver
    time_limit = 240000

    positions = MIPSolver(items, W, H, x_tol, y_tol, time_limit)

    draw_packing_mass(items, positions, W, H, x_tol, y_tol)

