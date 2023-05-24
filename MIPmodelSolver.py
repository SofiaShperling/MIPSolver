from ortools.linear_solver import pywraplp
from Item import Position


def MIPSolver(items, pallet_width, pallet_height, x_tolerance, y_tolerance, time_limit):
    n = len(items)
    # [START solver]
    # Create the mip solver
    solver = pywraplp.Solver.CreateSolver('SCIP')
    # [END solver]
    solver.SetTimeLimit(time_limit)

    # [START variables]
    infinity = solver.infinity()
    h = list(map(lambda i: i.height, items))
    w = list(map(lambda i: i.width, items))
    big_num = sum(i.width + i.height for i in items)
    x = [solver.NumVar(0, infinity, 'x') for i in range(n)]
    y = [solver.NumVar(0, infinity, 'y') for i in range(n)]
    z = [solver.BoolVar('z') for i in range(n)]
    r = [solver.BoolVar('r') for i in range(n)]
    l = [[solver.BoolVar('l') for j in range(n)] for i in range(n)]
    b = [[solver.BoolVar('b') for j in range(n)] for i in range(n)]
    cx = [solver.NumVar(0, solver.infinity(), f'cx{i}') for i in range(n)]
    cy = [solver.NumVar(0, solver.infinity(), f'cy{i}') for i in range(n)]

    print('Number of variables =', solver.NumVariables())
    # [END variables]

    # [START constraints]
    for i in range(n):
        solver.Add(cx[i] <= big_num * z[i])
        solver.Add(cy[i] <= big_num * z[i])
        solver.Add(cx[i] <= x[i] + w[i] / 2 * (1 - r[i]) + h[i] / 2 * r[i])
        solver.Add(cy[i] <= y[i] + h[i] / 2 * (1 - r[i]) + w[i] / 2 * r[i])
        solver.Add(cx[i] >= x[i] + w[i] / 2 * (1 - r[i]) + h[i] / 2 * r[i] - big_num * (1 - z[i]))
        solver.Add(cy[i] >= y[i] + h[i] / 2 * (1 - r[i]) + w[i] / 2 * r[i] - big_num * (1 - z[i]))
        solver.Add(x[i] + w[i] * (1 - r[i]) + h[i] * r[i] <= pallet_width + (1 - z[i]) * big_num)
        solver.Add(y[i] + h[i] * (1 - r[i]) + w[i] * r[i] <= pallet_height + (1 - z[i]) * big_num)
        for j in range(n):
            if i != j:
                solver.Add(x[i] + w[i] * (1 - r[i]) + h[i] * r[i] <= x[j] + (1 - l[i][j]) * big_num)
                solver.Add(y[i] + h[i] * (1 - r[i]) + w[i] * r[i] <= y[j] + (1 - b[i][j]) * big_num)
                solver.Add(l[i][j] + l[j][i] + b[i][j] + b[j][i] >= 1)
                solver.Add(sum((items[i].mass * cx[i] if z[i] else 0) for i in range(n)) <= (pallet_width / 2 + x_tolerance) * sum(
                        items[i].mass * z[i] for i in range(n)))
                solver.Add(sum((items[i].mass * cx[i] if z[i] else 0) for i in range(n)) >= (pallet_width / 2 - x_tolerance) * sum(
                        items[i].mass * z[i] for i in range(n)))
                solver.Add(sum((items[i].mass * cy[i] if z[i] else 0) for i in range(n)) <= (pallet_height / 2 + y_tolerance) * sum(
                        items[i].mass * z[i] for i in range(n)))
                solver.Add(sum((items[i].mass * cy[i] if z[i] else 0) for i in range(n)) >= (pallet_height / 2 - y_tolerance) * sum(
                        items[i].mass * z[i] for i in range(n)))

    print('Number of constraints =', solver.NumConstraints())
    # [END constraints]

    # [START objective]
    solver.Minimize(pallet_height * pallet_width - sum(h[i] * w[i] * z[i] for i in range(n)))
    # [END objective]

    # [START solve]
    status = solver.Solve()
    # [END solve]

    # [START print_solution]
    if status == pywraplp.Solver.OPTIMAL:
        print('Solution:')
        print('Objective value =', solver.Objective().Value())
    else:
        print('The problem does not have an optimal solution.')
    # [END print_solution]
    print('Solution:', solver.Objective().Value())
    # [START advanced]
    print('\nAdvanced information:')
    print('Problem solved in %f milliseconds' % solver.wall_time())
    print('Problem solved in %d iterations' % solver.iterations())
    print('Problem have a bound', solver.Objective().BestBound())
    # [END advanced]
    zv = [bool(t.solution_value()) for t in z]
    xv = [t.solution_value() for t in x]
    yv = [t.solution_value() for t in y]
    rv = [bool(t.solution_value()) for t in r]
    rez = [Position(xv[i], yv[i], rv[i]) if zv[i] else None for i in range(n)]
    return rez