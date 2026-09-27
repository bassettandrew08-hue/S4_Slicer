# S4 Slicer
A generic non-planar slicer, that can print almost any part without support.

> **This fork** adds a one-command headless pipeline for the Core R-Theta printer: STL in, 4-axis G-code out,
> with CuraEngine called directly (no Cura window). It has per-model settings files and a deformation that avoids
> floating islands.
>
> **Start here: [S4_PIPELINE.md](S4_PIPELINE.md)** (setup, tutorial, settings reference).
>
> What changed from the original notebook, and why: [CHANGELOG.md](CHANGELOG.md).
>
> ```
> venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl"
> ```
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
