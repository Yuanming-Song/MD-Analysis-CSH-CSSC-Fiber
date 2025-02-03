
set ini 17
set fini 19
set begframe 0
set endframe 20

set ref [list 0.0 0.0 1.0] 
 package require pbctools
source /dfs2/tw/yuanmis1/statistics.tcl 


set molwork [mol new ../CSH_CSSC_3to7_cyl_water_ions.psf]
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
        puts [molinfo top get numframes]
}
set selt {name C5 C6 C7 C8 C9 C10}
set selt2 {name C12 C13  C16 C17 C14 C18 C19}

animate write dcd data/tempori.dcd  beg 0 end -1 skip 10 waitfor all top
set round [expr [molinfo top get numframes] / (10*($endframe - $begframe + 1))]
set lastframe [expr [molinfo top get numframes] -1]
for {set m 0} {$m<$round} {incr m} {
	animate delete all
	if {[expr ($m + 1)*($endframe - $begframe ) ] <  $lastframe} {
		animate read dcd data/tempori.dcd beg $begframe end $endframe skip 1 waitfor all $molwork
	} else {
		animate read dcd data/tempori.dcd beg $begframe end -1 skip 1 waitfor all $molwork
	}
	pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all


         for {set i 1} {$i<145} {incr i} {
                 set sel1 [atomselect top "resname CSH and resid $i"]
                         $sel1 set resid [expr 145 + $i]
                         }

                         set sel1 [atomselect top "not resname TIP3"]


	set reslist "CSH CSSC"
	set coormin 0
	set coormax 0
	set totframes [molinfo top get numframes]
for {set j 0} {$j < $totframes} {incr j} {
	animate goto $j
	set outfile [open data/${ini}to${fini}rdf[expr $begframe + ${j}].dat w]

	for {set i 1} {$i <= 269 } {incr i} {
		set sel1 [atomselect top "resid $i and $selt"]
		set res1 [lindex [$sel1 get resname] 0 ]
		set com1 [measure center $sel1 ] 
		for {set k [expr $i + 1 ]} {$k <= 269 } {incr k} {	
			set sel2 [atomselect top "resid $k and $selt"]
			set com2 [measure center $sel2 ]
			set dis [veclength [vecsub $com1 $com2 ] ] 
			set res2 [lindex [$sel2 get resname] 0]
			puts $outfile "$res1 $res2 $dis"

			if {$k < 146 } {
       				set sel2 [atomselect top "resid $k and $selt2"]
				set com2 [measure center $sel2 ]
                        	set dis [veclength [vecsub $com1 $com2 ] ] 
                        	set res2 [lindex [$sel2 get resname] 0]
                        	puts $outfile "$res1 $res2 $dis"
			} 
		}
		if {$i < 146 } {
			set sel1 [atomselect top "resid $i and $selt2"]
                	set res1 [lindex [$sel1 get resname] 0 ]
                	set com1 [measure center $sel1 ] 
	                for {set k [expr $i + 1 ]} {$k <= 269 } {incr k} {
        	                set sel2 [atomselect top "resid $k and $selt"]
                	        set com2 [measure center $sel2 ]
                        	set dis [veclength [vecsub $com1 $com2 ] ] 
	                        set res2 [lindex [$sel2 get resname] 0]
        	                puts $outfile "$res1 $res2 $dis"

        	                if {$k < 146 } { 
                	                set sel2 [atomselect top "resid $k and $selt2"]
                	                set com2 [measure center $sel2 ] 
                	                set dis [veclength [vecsub $com1 $com2 ] ]
                	                set res2 [lindex [$sel2 get resname] 0]
                	                puts $outfile "$res1 $res2 $dis"
                	        }
                	}
		}
	}
	puts "frame $j"	
	close $outfile
}
set endframe [expr $endframe + $endframe - $begframe + 1]
set begframe [expr $begframe + $j+1] 

}
molinfo top get numframes

exit

