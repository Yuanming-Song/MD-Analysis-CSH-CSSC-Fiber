# Clear previous drawings
draw delete all

# Select atoms C1 and C20
set seltext "resid 10 and resname CSSC"
set C1 [atomselect top "$seltext and name C1"]
set C20 [atomselect top "$seltext and name C20"]

# Get coordinates of C1 and C20
set C1_pos [lindex [$C1 get {x y z}] 0]
set C20_pos [lindex [$C20 get {x y z}] 0]

# Calculate average z-coordinate (aveZ)
# Calculate average x, y, and z coordinates
set aveX [expr {([lindex $C1_pos 0] + [lindex $C20_pos 0]) / 2.0}]
set aveY [expr {([lindex $C1_pos 1] + [lindex $C20_pos 1]) / 2.0}]
set aveZ [expr {([lindex $C1_pos 2] + [lindex $C20_pos 2]) / 2.0}]

set ave_pos [list $aveX $aveY $aveZ]
# Define start and end points of the cylinder with full 3D vectors
set start_point [list 0.0 0.0 [expr $aveZ - 1.0]]
set end_point [list 0.0 0.0 [expr $aveZ + 1.0]]

# Draw a cylinder centered at (0, 0, aveZ) on the xy-plane with radius 27 and height 2
draw color yellow
draw material Transparent
graphics top cylinder $start_point $end_point radius 15 resolution 30
# Set C1 and C20 positions with aveZ as the z-coordinate
set C1_flat [list [lindex $C1_pos 0] [lindex $C1_pos 1] $aveZ]
set C20_flat [list [lindex $C20_pos 0] [lindex $C20_pos 1] $aveZ]
# Calculate the direction vector from C1_flat to C20_flat
set direction_vector [vecsub $C20_flat $C1_flat]

# Normalize the direction vector to a unit vector and scale it to length 5
set half_vector [vecscale [expr {5.0 / [veclength $direction_vector]}] $direction_vector]

# Define start and end points for the arrow centered around ave_pos
set C1_start [vecsub $ave_pos $half_vector]
set C1_end [vecadd $ave_pos $half_vector]

# Draw the arrow with fixed length 10, centered around ave_pos
vmd_draw_arrow top $C1_start $C1_end purple 0.5

# Define average position and origin
set origin [list 0 0 $aveZ]

# Calculate the mirrored position of ave_pos across origin
set mirrored_pos [list [expr {-$aveX}] [expr {-$aveY}] $aveZ]

# Draw arrow from ave_pos to mirrored_pos
vmd_draw_arrow top $ave_pos $mirrored_pos black 0.5

# Draw a solid sphere at ave_pos with radius 0.7
draw color green
draw material Opaque
draw sphere $ave_pos radius 0.7

# Draw a solid sphere at ave_pos with radius 0.7
draw color black
draw material Opaque
draw sphere $origin radius 1
