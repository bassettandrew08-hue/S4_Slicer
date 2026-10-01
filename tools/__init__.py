"""
Developer tools for the S4 pipeline (run as scripts: venv\\Scripts\\python tools\\<tool>.py ...).

    check_equivalence.py  gate: the fast pipeline must reproduce the reference (notebook port) byte for byte
    check_print_time.py   gate: s4/print_time.py must give R-Theta Sim's print time on the fixtures in tools/fixtures/
    compare_gcode.py      compare two 4-axis G-code files line by line within tolerances
    model_suite.py        slice all test models in parallel and flag quality regressions against a baseline run
    run_reference.py      run the reference port's deform or map stage on its own (e.g. to map a notebook pickle)
"""
