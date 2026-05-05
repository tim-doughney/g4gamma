import sys; sys.path.insert(0, '.')
import numpy as np, g4gamma as g
edges = np.linspace(0, 3000, 3001) * g.units.keV
res = g.build_spectrum(g.IsotopeKey(19, 40, 0), -1, edges)
print('peak at 661 keV:', np.array(res.counts)[:])