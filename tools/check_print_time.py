"""
Parity check: the pipeline's print-time estimate (s4/print_time.py) must equal R-Theta Sim's (sim/r-theta-simulator.html).

    venv\\Scripts\\python tools\\check_print_time.py

The expected numbers are the sim's own parseGcode + plan totals for the fixtures in tools/fixtures/, computed in the
browser with the settings in SIM_SETTINGS (the sim's defaults when they were taken). After changing either planner,
recompute them in the sim (browser console on the sim page):

    const t = await fetch('FIXTURE_URL').then(r => r.text()); const S2 = {...DEFAULTS};
    const Q = await parseGcode(S2, t.split(/\\r?\\n/), null); computeActual(S2, Q); plan(S2, Q); [Q.total, Q.n]

Exit code 0 = PASS.
"""
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from s4 import print_time  # noqa: E402

# sim setting -> value used for the expected totals, as pipeline map names
SIM_SETTINGS = dict(
    MAX_SPEED_X=150, MAX_SPEED_Z=50, MAX_SPEED_B=180, MAX_SPEED_C=360, MAX_SPEED_E=60,
    MAX_ACCEL_X=2000, MAX_ACCEL_Z=200, MAX_ACCEL_B=5000, MAX_ACCEL_C=3000, MAX_ACCEL_E=3000,
    CORNER_SPEED=8, HOME_X=140, HOME_Z=200, HOME_B=0, HOME_SPEED=40,
)  # feedMode 'linear', g0Inverse 'inverse' (the only modes the port implements)

EXPECTED = {  # fixture: (sim total seconds, sim move count)
    "print_time_cases.gcode": (31.73139114280078, 25),
    "pi_head.gcode": (163.2995216703945, 5535),
}
TOL = 1e-6  # relative; both sides are deterministic float64 (feed/dwell stored as float32 on both)


def main():
    ok = True
    for name, (want, want_n) in EXPECTED.items():
        path = os.path.join(HERE, "tools", "fixtures", name)
        with open(path, errors="replace") as fh:
            lines = fh.read().split("\n")
        m = {**print_time.MACHINE_DEFAULTS, **SIM_SETTINGS}
        n = len(print_time.parse(lines, m)[0])
        got = print_time.estimate(lines, SIM_SETTINGS)
        rel = abs(got - want) / want
        good = rel <= TOL and n == want_n
        ok &= good
        print(f"{'PASS' if good else 'FAIL'}  {name:24s} port {got:.9f} s, sim {want:.9f} s (rel {rel:.1e}); "
              f"moves {n} vs {want_n}")
    print("PRINT-TIME PARITY PASS" if ok else "PRINT-TIME PARITY FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
