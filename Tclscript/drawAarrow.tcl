# Clear all previous drawings
draw delete all

# Define the arrow-drawing procedure with custom radius, color, and transparency
proc vmd_draw_arrow {mol start end color radius} {
    # Set transparency material for the arrow
    graphics $mol material Transparent
    # Calculate middle point for the cylinder part of the arrow
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    # Draw cylinder for the arrow shaft
    graphics $mol color $color
    graphics $mol cylinder $start $middle radius $radius
    # Draw cone for the arrowhead
    graphics $mol cone $middle $end radius [expr $radius * 1.5]
}

# Select atoms with resid 20 and resname CSSC
set seltext "resid 10 and resname CSSC"
set C1 [atomselect top "$seltext and name C1"]
set C20 [atomselect top "$seltext and name C20"]

# Get the coordinates of C1 and C20
set C1_pos [lindex [$C1 get {x y z}] 0]
set C20_pos [lindex [$C20 get {x y z}] 0]

# Draw an arrow from C1 to C20 in purple with radius 0.5
vmd_draw_arrow top $C1_pos $C20_pos purple 0.5

# Calculate the z-direction unit vector (C120dir) for C1 to C20
set dirz [expr [lindex $C20_pos 2] - [lindex $C1_pos 2]]
set C120dir [expr $dirz / abs($dirz)]

# Calculate the start and end points for the new arrow in the z-direction
set C1_start [list [lindex $C1_pos 0] [lindex $C1_pos 1] [expr [lindex $C1_pos 2] - 10 * $C120dir]]
set C1_end [list [lindex $C1_pos 0] [lindex $C1_pos 1] [expr [lindex $C1_pos 2] + 10 * $C120dir]]

# Draw the arrow in the z-direction with length 20, color black, and radius 0.5
vmd_draw_arrow top $C1_start $C1_end black 0.5