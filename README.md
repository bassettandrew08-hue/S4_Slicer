# S4 Slicer
A generic non-planar slicer, that can print almost any part without support.

> **This fork** adds a one-command headless pipeline for the Core R-Theta printer: STL in, 4-axis G-code out,
> with CuraEngine called directly (no Cura window). It has per-model settings files and a deformation that avoids
> floating islands.
>
> **Start here: [S4_PIPELINE.md](S4_PIPELINE.md)** (setup, slicing, the settings tutorial).
>
> Quick start, after the one-time setup in S4_PIPELINE.md:
>
> ```
> venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl"
> ```
>
> This writes `output_gcode/pi 3mm.gcode`. Open it in [R-Theta Sim](sim/README.md) to look at it.
>
> **Where things are**
>
> | path | what it is |
> |---|---|
> | `s4_slice.py` | the command-line entry point |
> | `s4/` | the pipeline's code (deform, slice, map, checks) |
> | `params/` | per-model settings profiles, picked up automatically by model name |
> | `tools/` | the checks for code changes (equivalence gate, six-model suite, print-time check) |
> | `sim/` | R-Theta Sim, a single-file G-code viewer |
> | `input_models/` | the STL test models |
> | `output_gcode/` | the pipeline's results |
> | `build/` | the pipeline's work files (git-ignored) |
> | `main.ipynb` | Joshua Bird's original notebook; it alone uses `input_gcode/`, `output_models/`, `pickle_files/` and `gifs/` |
> | `docs/` | [SETTINGS.md](docs/SETTINGS.md) (every setting), [DEVELOPING.md](docs/DEVELOPING.md) (checks, rules, modules) |
>
> The pipeline writes only `output_gcode/` and `build/`. The other data folders belong to the notebook.
>
> What changed from the original notebook, and why: [CHANGELOG.md](CHANGELOG.md).
>
> Upstream: [jyjblrd/S4_Slicer](https://github.com/jyjblrd/S4_Slicer). Everything below is the original README.

Please use the [dicussions tab](https://github.com/jyjblrd/S4_Slicer/discussions) to ask questions and help others.

[Try it now](https://colab.research.google.com/github/jyjblrd/S4_Slicer) on Google Colab! (note: colab free tier is only powerful enough to slice very simple models)

[![Watch the video](https://github.com/jyjblrd/S4_Slicer/blob/main/thumnail.jpeg?raw=true)](https://www.youtube.com/watch?v=M51bMMVWbC8)

Check out my [YouTube video](https://youtu.be/M51bMMVWbC8?si=pfud7bHgjYDnO2_z) for more details!

Thank you to JLCCNC for helping create the extruder mount and build plate for my [4 Axis Core R-Theta Printer](https://github.com/jyjblrd/Core-R-Theta-4-Axis-Printer).



Bibtex Citation:
```
@software{Bird_S4_Slicer,
author = {Bird, Joshua},
license = {GPL-3.0},
title = {{S4 Slicer}},
url = {https://github.com/jyjblrd/S4_Slicer}
}
```
