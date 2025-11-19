from juliadraftsolver import NonlinearPoissonSolver

print("precompiling & solving...")
solver=NonlinearPoissonSolver()
x,u=solver.solve(m=2)
print(x)
print(u)
