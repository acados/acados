import casadi as cs
import numpy as np

d = cs.SX.sym('d')
s = cs.SX.sym('s')
x = cs.vertcat(d, s)

# y = -1
y = -0.001

## Nominal QP
# f = 0.5*1e-3*d**2 -d
# g = cs.vertcat((y-1) + d, (4-y**2)-2*y*d)
# nlp2 = {'x':d, 'f':f, 'g':g}

# solver2 = cs.nlpsol('solver2', 'ipopt',nlp2)
# lbg = cs.vertcat(-cs.inf,-cs.inf)
# ubg = cs.vertcat(0.0, 0.0)
# res = solver2(lbg=lbg, ubg=ubg)
# print(res['x'])

# print(res2['x'])

# f = s
f = 0.5*1e-4*d**2 + s
g = cs.vertcat((y-1) + d, (4-y**2)-2*y*d - s, s)
lbg = cs.vertcat(-cs.inf, -cs.inf, 0.0)
ubg = cs.vertcat(0.0, 0.0, cs.inf)
nlp = {'x':x, 'f':f, 'g':g}
solver = cs.nlpsol('solver', 'ipopt',nlp)
res = solver(lbg=lbg, ubg=ubg)
print(res['x'])

lbg2 = cs.vertcat(-cs.inf,-cs.inf)
ubg2 = cs.vertcat(0, res['x'][1])
f2 = 0.5*1e-3*d**2 - d
g2 = cs.vertcat((y-1) + d, (4-y**2)-2*y*d)
nlp2 = {'x':d, 'f':f2, 'g':g2}
solver2 = cs.nlpsol('solver2', 'ipopt',nlp2)
res2 = solver2(lbg=lbg2, ubg=ubg2)
print(res2['x'])
