#
# getEdges.tcl
#
# JAF and ECW
#
# jfreites@uci.edu
#
# version 1 08/14
#
# use a psf to generate network nodes according to
# Benson and Daggett scheme
#
# Input: PSF/DCDs and psfile.nodes psfile.bonds (node list and bonded list)
# Output: two files
#       {outputfilename}_{framenumber}.edges one edge per line with total number of contacts in 3d column
#       {outputfilename}_{framenumber}.edgatt one edge per line with nature of interaction in 3d column
#
# Uses C-C and C-X cutoff distances as in Benson and Daggett

# System/Trajectory parameters
#------------------------------------------------------------------------
# dataDir: this is concatenated to the myPSF and file names in theFiles and myReference
# workDir: this is concatenated to the output files
# myPSF: your topology file
# trajFileType: file type for your trajectory file
# step: step length used when reading the trajectory files

set myPSF data/mix_27A.pdb
set workDir "edge/"
set myNodeList data/CSH_CSSC_mix_27A_cyl_lattice_redo.psf.nodes
set theFiles {}


set step 25
set trajFileType dcd
set firstrun 6
set totalrun 10
set begframe 0
set dcdfile data/mix_27A_sum6to10_every_${step}.dcd
#-----------------------------------------------------------------------------
# Output file name

set outfile mix_27A_3to7_cyl_tip_CSH_every_$step

# theFiles:
# Provide a TCL list of trajectory file names or use your TCL skills to build it

#we'll ignore the dcds with zero frames after round(numframe/1000)*1000=0
#see deamidated worksheet for details
#the current workDir is set to the "second part" first part adapted from it


#for {set i 1} {$i <= 9} {incr i} {
#        lappend theFiles nvt0${i}.dcd
#}

#for {set i $firstrun} {$i <= $totalrun} {incr i} {
#	if {$i<10} {
#		lappend theFiles nvt0${i}.dcd
#	} else {
#        	lappend theFiles nvt${i}.dcd
#	}
#}

# theFileRange:
# Provide a TCL list with the first and last frame number to be analyzed in each
# trajectory file.
# Leave theFileRange empty (set to "") if you want all the frames of all the files
# to be analyzed.

#set theFileRange [list first1 last1 first2 last2 ...]
#

set theFileRange ""

#------------------------------------------------------------------------
# main selection

set  mySelection "resname CSH"
#CSH "

#
#------------------------------------------------------------------------


#------------------------------------------------------------------------
# Cutoff Parameters 
#
#
set carbonCutoff 5.4
set Cutoff 4.6
set sulfurCutoff 6.3

#
#------------------------------------------------------------------------
#
#
# Initial frame count

#---------------------END of USER INTERFACE ---------------------------------------------------
#
#
#------------------------------------------------------------------------


#PBCed distance function to measure contacts

proc myMeasureContacts {sel ref cutoff box } {
  set selpos [$sel get {x y z}]
  set refpos [$ref get {x y z}]
  set selind [$sel get index]
  set refind [$ref get index]
  set selres [$sel get resid]
  set refres [$ref get resid]
  set selList {}
  set refList {}
  foreach vec1 $selpos ind1 $selind res1 $selres {
    foreach vec2 $refpos ind2 $refind res2 $refres {
      set dist {}
      if { $res1 != $res2} {
        set dif [vecsub $vec1 $vec2]
        foreach x $dif  {
          #lappend dist [expr {$x - $box*round($x/$box)}]
          lappend dist $x
        }
        if {[veclength $dist] <= $cutoff} {
          lappend selList $ind1
          lappend refList $ind2
        }
      }
    }
  }
  return [list $selList $refList]
  
}


#Create a dictionary based on Eric's node list

set idf [open $myNodeList r]
while {[gets $idf line] >= 0} {
  set nindex [lindex $line 0]
  dict set network $nindex resname [lindex $line 1]
  dict set network $nindex resid [lindex $line 2]
  dict set network $nindex moid [lindex $line 3]
  dict set network $nindex nname [lindex $line 4]
  dict set network $nindex ntype [lindex $line 5]
  foreach thing [lrange $line 6 end] {
    dict set indices $thing $nindex
  }
}
close $idf

set fid [open "tips_index_6_to_10_every_25_.dat" r]
set tiplines [split [read $fid] "\n"]
close $fid


#------------------------------------------------------------------------
#
#

set molwork [mol new $myPSF]

animate delete all $molwork


set nframes $begframe




# Process beg/end trajectory files lists

#if {$theFileRange == ""} {
#        foreach dcdfile $theFiles {
#                append theFileRange "0 -1 "
#        }
#} else {
#        if {[llength $theFileRange] != [expr 2 * [llength $theFiles]]} {
#                puts "the File Range list inconsistent with Files list"
#                exit
#        }
#}

# Loop over files / get one frame at a time /Loop over selections
#------------------------------------------------------------------------
#for {set i $firstrun} {$i <= $totalrun} {incr i} {
set logfile [open edge.log a]
#        if {$i<10} {
#                set dcdfile nvt0${i}.dcd
#        } else {
#                set dcdfile  nvt${i}.dcd
#        }
#	if {$i<16} {
#		set setpi [expr $step *10]
#	} else {
#		set stepi $step
#	}
animate delete all $molwork
animate read $trajFileType $dcdfile beg 0 end -1 skip 1 waitfor all $molwork
set totframes [molinfo $molwork get numframes]
for {set frame 0} {$frame < $totframes} {incr frame} {
  incr nframes
  if {$nframes < 10} {
    set framenum 00${nframes}
  } elseif {$nframes < 100} {
    set framenum 0${nframes}
  } else {
    set framenum $nframes
  }
  set idf [open ${workDir}${outfile}_${framenum}.edges w]
  set idf2 [open ${workDir}${outfile}_${framenum}.edgeatt w]
  animate goto $frame
  set selresind [lindex $tiplines $frame]
  set nohs [atomselect $molwork "($mySelection) and noh and index $selresind"]
  set nocarbons [atomselect $molwork "($mySelection) and not (hydrogen or carbon or sulfur) and index $selresind"]
  set carbonAtoms [atomselect $molwork "($mySelection) and carbon and index $selresind"]
  set sulfurAtoms [atomselect $molwork "($mySelection) and sulfur and index $selresind"]
  #               prepareFrame $molwork
  $carbonAtoms update
  $nocarbons update
  $nohs update
  set box [molinfo $molwork get a]
  set carbonContacts [myMeasureContacts $carbonAtoms $carbonAtoms $carbonCutoff $box]
  set sulfurContacts [myMeasureContacts $sulfurAtoms $sulfurAtoms $sulfurCutoff $box]
  set Contacts [myMeasureContacts $nocarbons $nohs $Cutoff $box]
  set Ai [concat [lindex $carbonContacts 0] [lindex $Contacts 0] [lindex $sulfurContacts 0]]
  set Aj [concat [lindex $carbonContacts 1] [lindex $Contacts 1] [lindex $sulfurContacts 1]]
  #		write something
  set weight {}
  set interactType {}
  foreach thing1 $Ai thing2 $Aj {
    #			puts -nonewline "**** atoms are $thing1 $thing2 "
    set Ni [dict get $indices $thing1]
    set Nj [dict get $indices $thing2]
    if {$Ni != $Nj} {
      set pair [lsort -integer [list $Ni $Nj]]
      #				puts "nodes are $pair"
      dict incr weight $pair
      if {![dict exists $interactType $pair]} {
        set interaction [list [dict get $network [lindex $pair 0] nname] [dict get $network [lindex $pair 1] nname]]
        switch -- $interaction {
          {PHE PHE} {
            dict set interactType $pair PHE-PHE
          }
          {AM1 AM2} -
            {AM2 AM1} -
            {AM1 AM1} -
            {AM2 AM2} {
              dict set interactType $pair AM-AM
            }
          {CYS CYS} {
            dict set interactType $pair DISU
          }
          {AM1 PHE} -
            {PHE AM1} {
              dict set interactType $pair AM1-PHE
            }
          {AM2 PHE} -
            {PHE AM2} {
              dict set interactType $pair AM2-PHE
            }
          default {
            dict set interactType $pair STER
          }
        }
      }
    }
  }
  dict for {key val} $weight {
    puts $idf "$key $val"
  }
  dict for {key val} $interactType {
    puts $idf2 "$key $val"
  }
  close $idf
  close $idf2
  puts "frame $nframes"
}
puts $logfile "$dcdfile $nframes"
close $logfile
#}
#------------------------------------------------------------------------

exit
