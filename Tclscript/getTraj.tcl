set inilist [list 6]
set fini 15
set endframe -1
set steps 25
set molwork [mol new ../setup/CSH_CSSC_mix_27A_cyl_water_ions.psf]
set dcdlist {}
foreach inin $inilist {
for {set i $inin} {$i <= $fini } {incr i} {
	set j [expr $i + 1 ]
	if {$j<10} {
		set dcdname ../nvt0${j}.dcd
	} else {
		set dcdname ../nvt${j}.dcd
		#set dcdname ../nvt${j}.dcd
	}
	if {[file exists $dcdname] == 1 } {
		if {$i < 10} {
			lappend dcdlist ../nvt0${i}.dcd
		} else {
			lappend dcdlist ../nvt${i}.dcd
			#lappend dcdlist ../nvt${i}.dcd
		}
	} else {
		set fini [expr $i - 1 ]
		break
	}
	
} 



animate delete all
set logfile [open traj.log a]
foreach dcdfile $dcdlist {
        set otf [molinfo top get numframes]
        animate read dcd $dcdfile beg 0 end -1 skip $steps waitfor all $molwork
        puts $logfile "$dcdfile $steps $otf"
	#puts "was [expr $otf/$steps.] steps, with $dcdfile, now [expr [molinfo top get numframes]/$steps.] "
}
close $logfile

#animate write dcd data/mix_27A_sum${ini}to${fini}.dcd  beg 0 end -1 skip $steps waitfor all $molwork
#animate delete all
#animate read dcd data/mix_27A_sum${ini}to${fini}.dcd beg 0 end -1 skip 1 waitfor all $molwork
 package require pbctools
pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all

for {set i 1} {$i<353} {incr i} {
	set sel1 [atomselect top "resname CSH and resid $i"]
	$sel1 set resid [expr 411 + $i]
}

animate write dcd data/mix_27A_sum${inin}to${fini}_all_every_${steps}.dcd beg 0 end $endframe skip 1 waitfor all $molwork
#animate write dcd data/mix_27A_sum${inin}to${fini}_all_every_${steps}_lowf.dcd beg 0 end $endframe skip 1 waitfor all $molwork

set sel1 [atomselect top "resname CSH or resname CSSC"]

animate write dcd data/mix_27A_sum${inin}to${fini}_every_${steps}.dcd sel $sel1 beg 0 end $endframe skip 1 waitfor all $molwork
#animate write dcd data/mix_27A_sum${inin}to${fini}_every_${steps}_lowf.dcd sel $sel1 beg 0 end $endframe skip 1 waitfor all $molwork
}
#$sel1 writepdb data/mix_27A.pdb

exit

