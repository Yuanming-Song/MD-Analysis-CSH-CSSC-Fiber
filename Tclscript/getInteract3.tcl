#
#
#
# getInteract.tcl
#
# simple script to summarize the information in the .edgeatt files
#
# this is a tclsh script
#

#pairs to be summarized
set mytypes {PHE-PHE AM-AM AM1-PHE AM2-PHE DISU STER}

#where are .edgeatt files
set dataDir "edge/" 
set step 25
set moiety "_CSH"
set analysis nottip${moiety}
#prefix of the .edgeatt files 
#also the output file prefix
set infile mix_27A_3to7_cyl_${analysis}_every_${step}

#total number of .edgeatt files
set totframes 200


############# END OF UI ########################

set idf2 [open data/${infile}_intertypes_evol.dat w]


for {set frame 0} {$frame < $totframes} {incr frame} {
            	incr nframes
                if {$nframes < 10} {
                        set framenum 00${nframes}
                } elseif {$nframes < 100} {
                        set framenum 0${nframes}
                } else {
                        set framenum $nframes
                }
		set idf [open ${dataDir}${infile}_${framenum}.edgeatt r]

	foreach type $mytypes {
		set intertype($type) 0
	}
	set count 0
	while {[gets $idf line]>=0} {
		incr count
		incr intertype([lindex $line 2])
	}
	close $idf
	set ans {}
	foreach type $mytypes {
		if {$count > 0} {
			lappend ans [expr 1.0 * $intertype($type)/$count]
		} else {
			lappend ans 0
		}
		
	}
	puts $idf2 "$frame $ans"
}

	close $idf2
