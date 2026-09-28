; Parity fixture: s4/print_time.py vs R-Theta Sim's parseGcode + plan (tools/check_print_time.py).
; Each line exercises one parser or planner path.
G90
M82
G28 ; home all axes
G92 E0
G1 X10 Z5 F3000
G1 X20 E1.5 ; sticky F
X25 Z6 E2.0 ; axis words without a G: the last motion command
G1 B-20 F1200 ; B only (kind 2, length in degrees)
G1 C45 ; C only
G1 C90 B-10
G1 X30 C100 E2.5 ; X/Z again after B/C: kind change, junction speed 0
G1 X30 ; no change: skipped
G4 S0.5 ; dwell in seconds
G4 P250 ; dwell in milliseconds
G91
G1 X-5 Z1 E0.2 ; relative moves
G1 C-30
G90
M83
G1 E-1 F2400 ; E only (kind 3)
G1 E1
G93
G1 X40 Z10 B-30 C200 E0.5 F120 ; inverse time: 0.5 s
G0 X45 Z12 F600 ; G0 in G93 counts as inverse time
G1 E-0.8 F300 ; E only in G93
G1 E0.8 F300
G1 C220 B-35 F200 ; B/C only in G93
(a comment in parentheses) G1 X46 F150
G94
G28 X
G28 Z B
g1 x50 z20 f1e4 ; lower case, exponent
G0 C0 X0 Z20 B0
