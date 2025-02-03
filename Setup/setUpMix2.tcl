#Set up solution with the solutes in a lattice configuration

#Input: single-molecule (solute) PSF and DCD 
#	lattice and cell
#	output name

#Final Output: three files
# $outname_water_ions.psf
# $outname_water_ions.pdb
# $outname_water_ions.cons.pdb
#
# where $myPSF is the input name

set myPSF 	[list ../../../cssc/CSSC.psf ../../../csh/CSH.psf ]
set myDCD	[list ../../../cssc/nve_cssc.dcd ../../../csh/nve_csh.dcd ]
set myseg [list CSSC CSH]

set outname CSH_CSSC_mix_27A_cyl

#lattice points along each direction
set nx 5
set ny 4
set nz 39
#set nz 31

#number of molecs to add
set nmolmax [list 411 352]
#set nmolmax [list 651 850]

#cell
set a 94.
set b 94.
set c 400.
#set c 377.
#ionic strength in M
set ionics .0

set topoList {../../../ff/top_all36_cgenff.rtf ../../../ff/CSSC.str ../../../ff/CSH.str}


#----------- END OF USER INTERFACE -------------------
lappend auto_path /data/homezvol0/yuanmis1/tempotools/libs
package require tempoUserVMD
package require tempoUtils
package require solvate
package require autoionize
package require psfgen

set myminx [expr -0.5 * $a]
set mymaxx [expr 0.5 * $a]
set myminy [expr -0.5 * $b]
set mymaxy [expr 0.5 * $b]
set myminz [expr -0.5 * $c]
set mymaxz [expr 0.5 * $c]
#grid spacing (in A)
set sizex [expr ($a - 24  )/$nx]
set sizey [expr ($b - 24)/$ny]
set sizez [expr ($c - 24)/$nz]
set origx [expr -0.5*$a + 18.]
set origy [expr -0.5*$b + 18.]
set origz [expr -0.5*$c + 14.]

set count 0
set lista {}

for {set x 0} {$x < $nx} {incr x } {
        for {set y 0} {$y < $ny} {incr y } {
                for {set z 0} {$z < $nz} {incr z } {
                        set lattice($count) [list $x $y $z]
                        lappend lista $count
                        incr count
                }
        }
}


set acum 0

foreach seg $myseg psf $myPSF dcd $myDCD mynmolmax $nmolmax  {

	#load the single-molecule configs
	set mol($seg) [mol new $psf waitfor all]
	animate read dcd $dcd waitfor all $mol($seg)
	#center 
	centering -mol $mol($seg) 
	#place random configs on lattice points output single-molecule pdb
	set nmax [expr [molinfo $mol($seg) get numframes] -1]
	set myfrac [expr 1.0 * $mynmolmax/(($nx * $ny * $nz) - $acum)]
	set mylista [tempoUtils::randomSel $lista $myfrac]
	foreach thing $mylista {
		set ind [lsearch -integer $lista $thing]
		set lista [lreplace $lista $ind $ind]
	}
	set count 1

	puts "$seg $lista **** $mylista #### $mynmolmax %%%% $myfrac"
	#do the deed
	foreach thing $mylista {
	set point $lattice($thing)
        foreach {x y z} $point {}
                set xc [expr $origx + $x * $sizex ]
                set yc [expr $origy + $y * $sizey ]
                set zc [expr $origz + $z * $sizez ]
				set ind [expr round(rand()* $nmax)]
				set sel [atomselect $mol($seg) all frame $ind]
				puts "xxxx $ind oooo $thing **** [list $xc $yc  $zc] ***** "
				$sel moveby [list $xc $yc $zc]
				if {$count < 10} {
					set outfile ${outname}_${seg}_0${count}.pdb
				} else {
					set outfile ${outname}_${seg}_${count}.pdb
				}
				$sel set segname $seg
				$sel set resid $count
				$sel writepdb $outfile
				$sel moveby [vecscale -1.0 [list $xc $yc $zc]]
				$sel delete
				incr count
	}
	set acum [expr $count + $acum]
}















#generate lattice psf/pdb


resetpsf

foreach topo $topoList {
        topology $topo
}

#set mynmolmax 16
#set seg CSH

#foreach seg {CSH CSSC} mynmolmax {16 8} 

foreach seg $myseg mynmolmax $nmolmax {

segment $seg {
        for {set count 1} {$count <= $mynmolmax} {incr count} {
                if {$count < 10} {
                                set outfile ${outname}_${seg}_0$count.pdb
                        } else {
                                set outfile ${outname}_${seg}_${count}.pdb
                        }
                pdb $outfile
        }
}
        for {set count 1} {$count <= $mynmolmax} {incr count} {
                if {$count < 10} {
                                set outfile ${outname}_${seg}_0$count.pdb
                        } else {
                                set outfile ${outname}_${seg}_${count}.pdb
                        }
                coordpdb $outfile
        }

}

guesscoord
writepsf ${outname}_lattice.psf
writepdb ${outname}_lattice.pdb

#solvate 

solvate ${outname}_lattice.psf ${outname}_lattice.pdb -o ${outname}_solvated -s WAT -minmax [list [list $myminx $myminy $myminz] [list $mymaxx $mymaxy $mymaxz]]
#add ions

autoionize -psf ${outname}_solvated.psf -pdb ${outname}_solvated.pdb -sc $ionics -o ${outname}_water_ions -seg IONS

#load the final configuration to generate the constraints file

mol new ${outname}_water_ions.psf
mol addfile ${outname}_water_ions.pdb waitfor all

[atomselect top all] set beta 0.0
[atomselect top "segname CSH CSSC"] set beta 20.0
[atomselect top all] writepdb ${outname}_water_ions.cons.pdb

puts "check check ########################################## check check"
puts $myminx 
puts $myminy 
puts $myminz
puts $mymaxx 
puts $mymaxy 
puts $mymaxz
exit
