
 package require pbctools
source /dfs2/tw/yuanmis1/statistics.tcl 

set begframe 0

set ini 17
set fini 19
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
        molinfo top get numframes
}

animate write dcd data/tempori.dcd  beg 0 end -1 skip 300 waitfor all top
animate delete all
animate read dcd data/tempori.dcd beg $begframe end -1 skip 1 waitfor all $molwork
 package require pbctools
pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all


         for {set i 1} {$i<145} {incr i} {
                 set sel1 [atomselect top "resname CSH and resid $i"]
                         $sel1 set resid [expr 145 + $i]
                         }

                         set sel1 [atomselect top "not resname TIP3"]


#pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all
set reslist "CSH CSSC"
set coormin 0
set coormax 0
set totframes [molinfo top get numframes]
for {set j 0} {$j < $totframes} {incr j} {
	animate goto $j
set outfile [open data/${ini}to${fini}pheorient[expr [$begframe + ${j}].dat w]

	for {set i 1} {$i <= 269 } {incr i} {
	set sel1 [atomselect top "resid $i and name C8"]
	set sel2 [atomselect top "resid $i and name C5"]
	set coor1 [lindex [$sel1 get {x y z} ] 0]
	set coor2 [lindex [$sel2 get {x y z} ] 0]
	set coorv [vecsub $coor1 $coor2]	
if { ([llength $coor1] != 3) || ([llength $coor2 ]!= 3) ||([llength $coorv]!=3) } {
	puts "$j $i $coor1 $coor2 $coorv"
}
        set costheta [expr [vecdot $coorv $ref] / ([veclength $coorv] * [veclength $ref])]
	#set theta [expr acos([vecdot $coorv $ref] / ([veclength $coorv] * [veclength $ref]))]
	set theta [expr acos($costheta)]
	set dir [expr [lindex $coor2 2] - [lindex $coor1 2]]
	if {$dir < 0 } {
		set theta [expr 180 - $theta]
	} 
	set resn [$sel1 get resname]


	puts $outfile "$resn $theta"
        #puts $outfile "$resn $dir $costheta"

	}
	puts "frame $j"	
close $outfile
}


molinfo top get numframes

#exit

