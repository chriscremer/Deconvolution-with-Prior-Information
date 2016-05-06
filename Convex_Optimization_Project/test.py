

# # from matlab import engine
# import matlab.engine
# eng = matlab.engine.start_matlab()
# tf = eng.isprime(37)
# print(tf)


from cvxopt import matrix, solvers
c = matrix([-4., -5.])
G = matrix([[2., 1., -1., 0.], [1., 2., 0., -1.]])
h = matrix([3., 3., 0., 0.])

solvers.options['show_progress'] = False
sol = solvers.lp(c, G, h)

print(sol['x'])

