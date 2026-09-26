"""Map a sliced G-code with the fast mapper using an existing deformed pickle; compare to a baseline."""
import argparse
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np  # noqa: E402

from s4 import fast_map, meshio_s4  # noqa: E402
from s4.timing import TIMER  # noqa: E402
from tools.compare_gcode import compare  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--model", required=True)
ap.add_argument("--pkl", required=True)
ap.add_argument("--sliced", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--baseline")
a = ap.parse_args()

with TIMER("tetgen + offset"):
    inp = meshio_s4.load_and_tetrahedralize(a.model)
    in_cells = meshio_s4.cells_of(inp)
    in_centers = np.asarray(inp.cell_centers().points)
with TIMER("load pickle"):
    d = pickle.load(open(a.pkl, "rb"))
    d = meshio_s4.deformed_grid_from_points(inp, np.asarray(d.points))
    d_centers = np.asarray(d.cell_centers().points)
with TIMER("MAP total (fast)"):
    pts, stats = fast_map.map_gcode(in_cells, np.asarray(inp.points), in_centers, np.asarray(d.points), d_centers,
                                    meshio_s4.make_find_cells(d), a.sliced)
    with TIMER("write gcode"):
        fast_map.write_gcode(pts, a.out)
print(stats)
print(TIMER.report("fast map timings"))
if a.baseline:
    ok, _ = compare(a.baseline, a.out)
    os.system(f'cmp "{a.baseline}" "{a.out}" && echo BYTE-IDENTICAL')
