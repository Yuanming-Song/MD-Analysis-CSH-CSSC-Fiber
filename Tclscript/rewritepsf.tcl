
#mol new /dfs2/tw/yuanmis1/mrsec/cyl1/mix_27A/setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb
#set DIR /dfs2/tw/yuanmis1/mrsec/
package require psfgen
set myPSF       [list /dfs2/tw/yuanmis1/mrsec/cssc/CSSC.psf /dfs2/tw/yuanmis1/mrsec/csh/CSH.psf ]

set topoList {/dfs2/tw/yuanmis1/mrsec/ff/top_all36_cgenff.rtf /dfs2/tw/yuanmis1/mrsec/ff/CSSC.str /dfs2/tw/yuanmis1/mrsec/ff/CSH.str}

psfgen_logfile "psf.log"

foreach topo $topoList {
        topology $topo
}
set myseg [list CSSC CSH]
#set myseg [list A]
foreach seg $myseg {
segment $seg {
	pdb data/mix_27A.pdb

}
}                
coordpdb data/mix_27A.pdb

#for {set i 1} {$i<411} {incr i} {
#        set sel1 [atomselect top "resname CSH and resid $i"]
#        $sel1 set resid [expr 411 + $i]
#}

guesscoord
writepsf CSH_CSSC_mix_27A_cyl_lattice_redo.psf


