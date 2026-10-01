# R-Theta Sim

A single-file viewer for the 4-axis G-code this pipeline writes: `r-theta-simulator.html`. It shows the print as the
printer would draw it, with the nozzle tilt, the plastic flow and the path quality coloured on the beads, and it
compares the sliced file with your machine's settings. The pipeline's print-time estimate is a port of this sim's
planner, so the two agree ([how it is kept in step](../docs/DEVELOPING.md#other-tools)).

## Opening a print

Open `sim/r-theta-simulator.html` in a browser (double-click it). It needs internet once per load, for three.js from
cdnjs; without it the page says so. Then drop a G-code file anywhere on the view, or press **Open G-code**. **Load
demo** shows a sample cup.

## The File panel

- **Print time** is the sim's own time. Next to it is the **Slicer estimate** from the file's
  [settings block](../S4_PIPELINE.md#the-settings-block-at-the-top-of-the-g-code). With the same machine settings
  the two agree to within rounding.
- **Sliced for different machine settings** appears when the file's machine values (`NOZZLE_OFFSET`, the B range,
  `MAX_SPEED_*`, `MAX_ACCEL_*`, `CORNER_SPEED`, `HOME_*`, the layer height) differ from the sim's settings.
  **Use the file's values** copies them in. The sim never does that on its own, because its settings describe your
  real printer. With default settings the table always lists the nozzle offset and usually the layer height, until
  you settle them.

## Colour modes

The colour menu changes what the beads show: nozzle tilt (B), plastic per mm, the support check, height, path
roughness, or plain filament.

**Path roughness** colours each bead by how sharply the path turns at it. A turn only counts when a neighbouring
segment also turns sharply, so a part's real corners stay neutral and zig-zag patches light up (for example a rough
bridge underside). The same measure, as a count, is in `tools/model_suite.py`.

## Settings are remembered

The sim stores its settings in the browser, and only the values you changed from the defaults, so later default fixes
still reach you. Chrome shares that storage between all local HTML files; Firefox may forget the settings if the file
moves. **Reset to defaults**, then Apply, clears them.

## On the real printer

Two sim settings describe your printer and cannot be known from a file. Confirm them once:

- **C direction.** The sim assumes the S4 convention. Jog C+ on the printer and check which way the bed turns
  against the view; if it is the other way round, flip the setting.
- **Nozzle offset.** Measure the distance from the B pivot to the nozzle tip and set it in the sim and in the
  pipeline (`NOZZLE_OFFSET`, [see its note](../docs/SETTINGS.md#map-the-machine-and-the-final-4-axis-moves)). The two
  defaults differ until you do, which is why the File panel lists it.
