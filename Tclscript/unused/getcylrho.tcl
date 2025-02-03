
set begframe 132
set endframe 150
 package require pbctools
source /dfs2/tw/yuanmis1/statistics.tcl 
set ini 17
set fini 25
set molwork [mol new ../CSH_CSSC_3to7_cyl_water_ions.psf]
if {0} {
set dcdlist {}
for {set i $ini} {$i <= $fini } {incr i} {
 lappend dcdlist ../npt${i}n.dcd
}

foreach dcdfile $dcdlist {
        if {$i > 16} {
                animate read dcd $dcdfile beg 0 end -1 skip 1 waitfor all $molwork
        } else {
                animate read dcd $dcdfile beg 0 end -1 skip 10 waitfor all $molwork
        }
        molinfo top get numframes
}
animate write dcd data/temp.dcd  beg 0 end -1 skip 10 waitfor all top
}

#set round [expr [molinfo top get numframes] / (10*($endframe - $begframe + 1))]

set round [expr 2965 / (10*($endframe - $begframe + 1))]

#set lastframe [expr [molinfo top get numframes] -1]
set lastframe [expr 2965/10 - 1]
for {set m 0} {$m<$round} {incr m} {
        animate delete all
        if {[expr ($m + 1)*($endframe - $begframe ) ] <  $lastframe} {
                animate read dcd data/temp.dcd beg $begframe end $endframe skip 1 waitfor all $molwork
        } else {
                animate read dcd data/temp.dcd beg $begframe end -1 skip 1 waitfor all $molwork
        }
        pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all


	for {set i 1} {$i<124} {incr i} {
	        set sel1 [atomselect top "resname CSH and resid $i"]
	        $sel1 set resid [expr 145 + $i]
	}
#         pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all

#         for {set i 1} {$i<145} {incr i} {
#                 set sel1 [atomselect top "resname CSH and resid $i"]
#                         $sel1 set resid [expr 124 + $i]
#                         }

#                        set sel1 [atomselect top "not resname TIP3"]


#pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all
set coormin 0
set coormax 0
set reslist {CSH CSSC}
	set totframes [molinfo top get numframes]
	puts "check"
	for {set j 0} {$j < $totframes} {incr j} {
	
	foreach resn $reslist {
		set rlist {}
		set zlist {}

        set outfile [open data/cylrho${ini}to${fini}${resn}frame[expr $begframe + $j ].dat w]
	animate goto $j
	set sel1 [atomselect top "resname $resn and not name H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22"]

	set atomlist [$sel1 get index]

	foreach atomn $atomlist {
		set sel2 [atomselect top "index $atomn"]
		set x [$sel2 get x]
		set y [$sel2 get y]
		set z [$sel2 get z]
		set r [expr ($x**2+$y**2)**0.5]	
		lappend rlist $r
		lappend zlist $z	
	}
        		for {set i 0} {$i < [llength ${rlist}]} {incr i} {
        			puts $outfile "[lindex ${rlist} $i] [lindex ${zlist} $i]"


			}
			close $outfile

		}
		puts "round $m frame [expr $begframe + $j ]"
	}
	set endframe [expr $endframe + $endframe - $begframe + 1]
	set begframe [expr $begframe + $j]

}

#animate write dcd data/run2n_sum09to14_every50.dcd sel $sel1 beg 0 end -1 waitfor all top

#molinfo top get numframes

exit

if {0} {
    set coor [$sel1 get {z}]

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

set coormax [quantiles $coor  $CI ]
set coormin [quantiles $coor [expr 1- $CI] ]
set l [expr $coormax - $coormin]
set rho [expr 414/(15*15*3.14159*$l)]
puts $outfile "$j $coormax $coormin $l $rho "
}
