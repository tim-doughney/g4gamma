import sys
import os

# Adjust this if your build dir is elsewhere
if not any('g4gamma' in os.path.basename(p) for p in sys.path):
    for cand in ('build', '../build', '.'):
        if os.path.isdir(cand) and any(f.startswith('g4gamma') and f.endswith('.so')
                                       for f in os.listdir(cand)):
            sys.path.insert(0, cand)
            break
import g4gamma as g
import numpy as np
edges = np.linspace(0, 3000, 3001) * g.units.keV
res = g.build_spectrum(g.IsotopeKey(55, 137, 0), -1, edges)
print('peak at 661 keV:', np.array(res.counts)[661])