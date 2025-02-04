# we will use the following variable to calculate and print the average 
# number of ions found outside the cylinder at each step

#set avgNumIons 0

wrapmode cell

##################################################

proc calcforces {step unique} {

  global cylinderCenter cylinderRadius maxForce forceCoef avgNumIons begID lastID

  if { $step > 0 && $step % 500 == 0 } { 
	  
    #set avgNumIons [expr $avgNumIons / 100.]
    #print "Step $step, average number of ions outside the cylinder: $avgNumIons"
    #set avgNumIons 0
	  
    cleardrops 
	  
  }

  while {[nextatom]} { 

    if { [getid] < $begID || [getid] > $lastID } {
      dropatom ;# not a CSSC/CSH atom, forget about it
      continue
    }
	
    # vector between the atom and the cylinder's center

    set coor [getcoord]
   # #set xycoor [list [lindex $coor 0] [lindex $coor 1] 0.0]
   # #set relativePosition [vecsub $xycoor $cylinderCenter] 
   # #set rho              [veclength $relativePosition]
#	#set zcoor [list 0.0 0.0 [lindex $coor 2]]
#	#set relativePositionz [vecsub $zcoor $cylinderCenter]
#	#set rhoz [veclength $relativePositionz]
	set xyzcoor [list [lindex $coor 0] [lindex $coor 1] [lindex $coor 2]]
	set relativePositionxyz [vecsub $xyzcoor $cylinderCenter]
	set rhoxyz              [veclength $relativePositionxyz]
	if { $rhoxyz > 40 } {
		addforce [vecscale $relativePositionxyz [expr -$maxForce/$rhoxyz]]
	}
  }
}
