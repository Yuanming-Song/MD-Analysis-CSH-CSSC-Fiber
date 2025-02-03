
set outfile [open cylxy.dat w]

set molwork [mol new ../CSH_CSSC_3to7_cyl_water_ions.psf] 
 package require pbctools

 proc iota n {
        set res {}
             for {set i 0} {$i<$n} {incr i} {
                     lappend res $i
                         }
                             set res
                              }

proc transposeMatrix m {
    set cols [iota [llength [lindex $m 0]]]
    foreach row $m {
        foreach element $row   *col $cols {
            lappend ${*col} $element
        }
    }
    eval list $[join $cols " $"]
 }

animate read dcd ../npt09n.dcd beg 0 end -1 skip 500 waitfor all $molwork
pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all
set reslist "CSH CSSC"
set coormin 0
set coormax 0
set totframes [molinfo top get numframes]
for {set j [expr $totframes -1]} {$j < $totframes} {incr j} {

animate goto $j
		set sel1 [atomselect top "resname CSH or resname CSSC and (index 1 to 10)" ]
			set coor [$sel1 get {x y}]
			set $coor [transposeMatrix $coor]
			puts $outfile "$coor"
	}



exit

