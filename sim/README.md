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

## What the view shows

- **The real tool path.** The printer moves C, X, Z and B together, each in a straight line in joint space, so between
  two G-code points the nozzle tip follows a curve (a spiral when C turns), not the straight chord. The beads, the
  upcoming path, the travel moves and the collision checks all follow that curve. The tip sits at
  X + sin(B)·L, Z + (1 − cos(B))·L from the G-code values (L = nozzle offset), the same as `s4/geometry.py`.
- **Beads shaped by the plastic.** Each bead is as thick as a layer, measured along the nozzle, and as wide as the
  E value fills (Cura's bead model: a rectangle with round sides). It hangs below the tip, so the tip rides on its top
  face. The thickness is the file's layer height (or **Bead thickness** in Settings); a bead on the bed is as thick as
  the gap under the tip. Primes keep the nominal width. **Bead shape: Round** gives the old fixed-width tubes. Click a
  bead to see its width and thickness.
- **The bead being laid** grows up to the nozzle tip while the move plays.
- **Travel and retractions** (View panel, off by default): every move that lays no plastic, up to the current move,
  with a red dot where filament is pulled back and a green dot where it is pushed in without the tip moving.
- The view only redraws when something changes, so a large print costs nothing while it sits still. The part casts a
  shadow on the bed up to 400,000 bead segments.

## Colour modes

The colour menu changes what the beads show: nozzle tilt (B), plastic per mm, the support check, height, path
roughness, or plain filament.

**Path roughness** colours each bead by how sharply the path turns at it. A turn only counts when a neighbouring
segment also turns sharply, so a part's real corners stay neutral and zig-zag patches light up (for example a rough
bridge underside). The same measure, as a count, is in `tools/model_suite.py`.

## Preflight check

**Run preflight check** looks at the whole file: axis limits, speed, hotend flow, over- and under-extrusion, primes in
mid-air, floating plastic, and the nozzle against the bed and the printed part. The nozzle, its shaft and the heater
block are checked at several poses along each move (every 5° of C or B, or 4 mm of tip travel), not only where the
move ends, because a big turn can swing the head through the part between two G-code points. **Travel through
plastic** flags a move that lays no plastic while its tip goes about half a layer deep into a bead printed earlier
(skimming the top of the layer, as combing does, is not flagged). It uses the support-check layer height in Settings,
so set that to the file's layer height (the File panel offers it).

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
