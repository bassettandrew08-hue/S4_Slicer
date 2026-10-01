"""
S4 non-planar slicer pipeline: STL in, 4-axis (C/X/Z/B/E) G-code out. Entry point: s4_slice.py -> pipeline.run.

Data flow (one run):

    profile / params      settings: built-in defaults < params/<model>.json < --set < --cura-* / --notebook-exact
      -> fast_deform      tetgen mesh + the notebook's rotation (tilt) field; the deformation itself is
                          island_free (default) or the notebook's least-squares solve
      -> cura             headless CuraEngine slices the deformed STL flat (planar G-code)
      -> fast_map         maps the planar G-code back through the deformation to the 4-axis machine
                          (reference.py is the line-for-line notebook port it must match)
      -> feed_limits      slows moves that would drive an axis past its speed limit
      -> support_check    checks: extrusion printed in mid-air (support_check), poles and surface quality
         quality          (quality), the print-time estimate R-Theta Sim will show (print_time)
         print_time
      -> sim_header       the settings comment block at the top of the final G-code

Modules:

    pipeline            run(): the stages above, with per-stage timings
    profile             build profiles (deform / map / cura sections): defaults, loading, --set, JSON output
    params              deformation defaults (the profile's "deform" section) and the "iterations" expansion
    fast_deform         deformation stage: mesh attributes, rotation field optimisation, deformation dispatch
    island_free         fold-free, island-free deformation (fit + lifting of floating height minima)
    island_free_solver  its numba kernels and L-BFGS solvers (FitProblem, PenaltyProblem)
    fast_lsq            faster Jacobian products inside scipy's TRF least-squares (same results)
    geometry            shared numpy geometry: tangential axes, rotation matrices, nozzle-tip position
    meshio_s4           mesh loading / tetrahedralisation and VTK point-location helpers
    cura                Cura project (.3mf) setting resolution and headless CuraEngine slicing
    fast_map            planar G-code -> 4-axis G-code (notebook cells 16-18, vectorised)
    reference           faithful port of main.ipynb (deform and map), the golden reference
    feed_limits         per-axis speed limits applied to the final G-code (G93 inverse-time feeds)
    support_check       print-order support check: floating / ungrounded extrusion
    quality             final 4-axis G-code checks: poles, extrusion along the nozzle axis, path roughness
    print_time          print-time estimate, a port of R-Theta Sim's planner
    sim_header          settings comments (Cura header, ;SETTING_3, "; s4:" lines) for R-Theta Sim
    timing              stage timer
"""
