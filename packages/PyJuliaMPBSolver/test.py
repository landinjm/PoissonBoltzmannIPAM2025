from PyJuliaMPBSolver import JuliaMPBSolver

print("precompiling & solving...")
solver = JuliaMPBSolver()
x, c0, cp, cm = solver.mpbpsolve()
print(x)
print(c0)
print(cm)
print(cp)

solver = JuliaMPBSolver()
x, c0, cp, cm = solver.icmpbpsolve()
print(x)
print(c0)
print(cm)
print(cp)
