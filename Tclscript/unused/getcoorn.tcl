
source /dfs2/tw/yuanmis1/statistics.tcl 
package require math::statistics 
set outfile [open data/cylsizen1.dat w]

set molwork [mol new data/run2_redo.psf]
 package require pbctools

animate read dcd data/run2n1_sum09to12.dcd beg 0 end -1 skip 1 waitfor all $molwork
#pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all
set reslist "CSH CSSC"
set coormin 0
set coormax 0
set totframes [molinfo top get numframes]
for {set j 0} {$j < $totframes} {incr j} {

animate goto $j
#foreach resn $reslist {
#	for {set i 1} {$i < 146} {incr i} {
#		set sel1 [atomselect top "resname $resn and resid $i" ]
#		set atomlist [$sel1 get name]
#		foreach atomn $atomlist {
##			set sel1 [atomselect top "resname $resn and resid $i and name $atomn" ] 
#			set coor [$sel1 get {z}]
#			if {[lindex $coor 0] < $coormin} {
#				set coormin [lindex $coor 0]
#			}
#			if {[lindex $coor 0] > $coormax} {
#                                set coormax [lindex $coor 0]
#                        }
#		}
#	}
#}
set sel1 [atomselect top "resname CSH or resname CSSC"]
    set coor [$sel1 get {z}]
if {0} {
set coormin 0
set coormax 0
	set res [lindex $coor 0]
    foreach element $coor {
        if {$element >= $coormax} {
		set coormax $element
	} 
	if {$element <= $coormin } {
                set coormin $element
	}
	
    }
}
if {1} {
proc QuantilesRawData { data confidence } {
    variable TOOFEWDATA
    variable OUTOFRANGE



    set sorted_data [lsort -real -increasing $data]

    set result      {}
    set number_data [llength $sorted_data]
    foreach cond $confidence {
        set elem [expr {round($number_data*$cond)-1}]
        if { $elem < 0 } {
            set elem 0
        }
        lappend result [lindex $sorted_data $elem]
    }

    return $result
}
proc quantiles { arg1 arg2 {arg3 {}} } {
    variable TOOFEWDATA

    if { [catch {
        if { $arg3 == {} } {
            set result [QuantilesRawData $arg1 $arg2]
        } else {
            set result [QuantilesHistogram $arg1 $arg2 $arg3]
        }
    } msg] } {
    }
    return $result
}
}
set CI 0.98
set coormax [quantiles $coor  $CI ]
set coormin [quantiles $coor [expr 1- $CI] ]
set l [expr $coormax - $coormin]
set rho [expr 414/(15*15*3.14159*$l)]
puts $outfile "$j $coormax $coormin $l $rho "
}



exit

