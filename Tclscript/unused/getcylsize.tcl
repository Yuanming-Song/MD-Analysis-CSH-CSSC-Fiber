
 package require pbctools
source /dfs2/tw/yuanmis1/statistics.tcl 
set CI 0.98

set ini 10
set fini 26
set molwork [mol new ../CSH_CSSC_3to7_cyl_water_ions.psf]
set dcdlist {../npt09n.dcd}
for {set i $ini} {$i <= $fini } {incr i} {
 lappend dcdlist ../npt${i}n.dcd
}

set outfile [open data/cylsizerun9to${fini}_$CI.dat w]

foreach dcdfile $dcdlist {
        set otf [molinfo top get numframes]
	if {[regexp -all -inline -- {[0-9]+} $dcdfile] >= 16} {
                animate read dcd $dcdfile beg 0 end -1 skip 1 waitfor all $molwork
        } else {
                animate read dcd $dcdfile beg 0 end -1 skip 10 waitfor all $molwork
        }
        puts "was [expr $otf/10.] steps, with $dcdfile, now [expr [molinfo top get numframes]/10.] "
}
animate write dcd temp.dcd  beg 0 end -1 skip 10 waitfor all top
animate delete all
animate read dcd temp.dcd beg 0 end -1 skip 1 waitfor all $molwork

#pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all
set reslist "CSH CSSC"
set coormin 0
set coormax 0
set totframes [molinfo top get numframes]
for {set j 0} {$j < $totframes} {incr j} {

animate goto $j
set sel1 [atomselect top "resname CSH or resname CSSC and not name H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22"]
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
set coormax [quantiles $coor  $CI ]
set coormin [quantiles $coor [expr 1- $CI] ]
set l [expr $coormax - $coormin]
set rho [expr 414/(15*15*3.14159*$l)]
puts $outfile "$j $coormax $coormin $l $rho "
}

#animate write dcd data/run2n_sum09to14_every50.dcd sel $sel1 beg 0 end -1 waitfor all top

molinfo top get numframes

exit

