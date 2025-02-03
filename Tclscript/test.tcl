
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

set coor {4.161861896514893 -12.648401260375977} {5.2132039070129395 -11.84918212890625} {5.2174153327941895 -10.67877197265625}

transposeMatrix $coor

