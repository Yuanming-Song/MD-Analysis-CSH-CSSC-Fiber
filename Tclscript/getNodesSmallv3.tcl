
set myPSF data/CSH_CSSC_mix_27A_cyl_lattice_redo.psf
set resname [list CSH CSSC]


#-------------------------------------------------------------

set molwork [mol new $myPSF waitfor all]

set myTypes(CSSC) {}
set myAnames(CSSC) {}
set myMoids(CSSC) {}
set myNames(CSSC) {}

lappend myTypes(CSSC) NOP
lappend myAnames(CSSC) "C5 C6 C7 C8 C9 C10"
lappend myMoids(CSSC) 1
lappend myNames(CSSC) PHE

lappend myTypes(CSSC)  DIP
lappend myAnames(CSSC) "O2 C4 N2"
lappend myMoids(CSSC) 1
lappend myNames(CSSC) AM1

lappend myTypes(CSSC)  DIP
lappend myAnames(CSSC) "O1 C1 N1"
lappend myMoids(CSSC) 1
lappend myNames(CSSC)  AM2

lappend myTypes(CSSC) THIO
lappend myAnames(CSSC) "C3 S1 C2"
lappend myMoids(CSSC) 1
lappend myNames(CSSC) CYS

lappend myTypes(CSSC)  NOP
lappend myAnames(CSSC) "C12 C13 C14 C16 C18 C19"
lappend myMoids(CSSC) 2
lappend myNames(CSSC) PHE

lappend myTypes(CSSC) DIP
lappend myAnames(CSSC) "O3 C17 N3"
lappend myMoids(CSSC) 2
lappend myNames(CSSC) AM1

lappend myTypes(CSSC)   DIP
lappend myAnames(CSSC) "N4 C20 O4"
lappend myMoids(CSSC) 2
lappend myNames(CSSC) AM2

lappend myTypes(CSSC) THIO
lappend myAnames(CSSC) "C11 S2 C15"
lappend myMoids(CSSC) 2
lappend myNames(CSSC) CYS

set myTypes(CSH) {}
set myAnames(CSH) {}
set myMoids(CSH) {}
set myNames(CSH) {}

lappend myTypes(CSH) NOP
lappend myAnames(CSH) "C5 C6 C7 C8 C9 C10"
lappend myMoids(CSH) 1
lappend myNames(CSH) PHE

lappend myTypes(CSH)  DIP
lappend myAnames(CSH) "O2 C4 N2"
lappend myMoids(CSH) 1
lappend myNames(CSH) AM1

lappend myTypes(CSH)  DIP
lappend myAnames(CSH) "O1 C1 N1"
lappend myMoids(CSH) 1
lappend myNames(CSH)  AM2

lappend myTypes(CSH) THIO
lappend myAnames(CSH) "C3 S C2"
lappend myMoids(CSH) 1
lappend myNames(CSH) CYS




#--------------------------------------------------------------
set idf [open ${myPSF}.nodes w]


set count 1
foreach res $resname {
	set sel [atomselect $molwork "resname $res"]
	set resids [lsort -unique -integer [$sel get resid]]
	$sel delete
	foreach resid $resids {
		foreach thing $myAnames($res) ntype $myTypes($res) nname $myNames($res) moid $myMoids($res) {
			set sel [atomselect $molwork "resname $res and resid $resid and name $thing"]
			puts $idf "$count $res $resid $moid $nname $ntype [$sel get index]"
			$sel delete
			incr count
		}
	}	
}


close $idf

exit
