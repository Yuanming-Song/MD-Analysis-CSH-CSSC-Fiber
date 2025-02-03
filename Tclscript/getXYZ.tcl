set outputdir "/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/analysis/data/"
set NewAna 0 ;# FALSE in R is equivalent to 0 in Tcl
set midfix ""

if {0} {
  set dcddir "/dfs9/tw/yuanmis1/mrsec/cyl1/mix_27A/"
  set outputpre "mix_27A_fibre_snapshot"
  set mypdb [file join $dcddir "setup/CSH_CSSC_mix_27A_cyl_water_ions.psf"] 
  
  if {1} { 
    set midfix "cylcssc."
    set outputpre "mix_27A_cylcssc_fibre_snapshot"
  }
} else {
  set dcddir "/dfs9/tw/yuanmis1/mrsec/cyl1/cssc_27A/"
  set outputpre "pure_27A_fibre_snapshot"
  set mypdb [file join $dcddir "setup/CSSC_pure_27A_cyl_water_ions.psf"] 
}
file mkdir $outputpre

mol new $mypdb
animate delete all
# Further processing can be added based on the values of variables like $outputdir, $midfix, $outputpre, $maxCSH, etc.

set FrameNeeded 1000
set dcdprefix "nvt"
set firstrun 12
set lastrun 100
set step 50

proc getdcdlist {firstrun midfix dcddir} {
  set dcdname_list [list]
  for {set runindex $firstrun} {$runindex <= 100} {incr runindex} {
    set dcdname [format "%s%02d.%sdcd" $::dcdprefix $runindex $midfix]
    set dcdfile [file join $dcddir $dcdname]
    if {[file exists $dcdfile]} {
      lappend dcdname_list $dcdname
    } else {
      break
    }
  }
  return $dcdname_list
}

set dcdlist [getdcdlist $firstrun $midfix $dcddir]
set num_dcd_files [llength $dcdlist]
set min_frames_per_file [expr {floor($FrameNeeded / $num_dcd_files)}]
set max_frames_per_file [expr {ceil((1+$FrameNeeded - $min_frames_per_file * $num_dcd_files) / $num_dcd_files)}]

proc generate_random_frame {min max} {
  #return [expr {int(rand() * ($max - $min + 1)) + $min}]
  return [expr int(floor(int(rand()*($max +1 ))))]
}

set randomdcdlistframe [list]
for {set i 0} {$i < $num_dcd_files} {incr i} {
  lappend randomdcdlistframe [expr {$min_frames_per_file}]
}

# Generate the initial list of available files
set available_files [list]
for {set i 0} {$i < $num_dcd_files} {incr i} {
  lappend available_files $i
}
set remainder [expr {$FrameNeeded % $num_dcd_files}]

# Distribute the remainder frames randomly to files
while {$remainder > 0} {
  set selected_file [lindex $available_files [expr {int(rand() * [llength $available_files])}]]
  set current_value [lindex $randomdcdlistframe $selected_file]
  set current_value_numeric [expr {$current_value}]
  lset randomdcdlistframe $selected_file [expr {$current_value_numeric + 1}]
  
  set index [lsearch -exact $available_files $selected_file]
  set available_files [lreplace $available_files $index $index]
  
  incr remainder -1
  
}





set totframe 1000
for {set i 0} {$i < $num_dcd_files} {incr i} {
  set frames_to_extract [lindex $randomdcdlistframe $i]
  if {[expr $frames_to_extract < 1]} {
    continue
  }
  set dcdname [lindex $dcdlist $i]
  puts $dcdname
  set random_frames {}
  for {set j 0} {$j < $frames_to_extract} {incr j} {
    lappend random_frames [expr {int(rand() * $totframe) }]
  }
  foreach frame $random_frames {
    set extractedframe [molinfo top get numframes]
    animate read dcd ${dcddir}$dcdname beg  $frame  end  $frame 
    # if $totframe equals [molinfo top get numframes], add 1 to i+1 term of $randomdcdlistframe
    if {$extractedframe == [molinfo top get numframes]} {
      # Increase the next term by 1
      set randomdcdlistframe [lset randomdcdlistframe $i+1 [expr {[lindex $randomdcdlistframe $i+1] + 1}]]
    }
  }
  
}
package require pbctools
package require topotools 1.5

pbc wrap -centersel " resname CSH CSSC" -center com -compound resid -all
for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame} {
  animate goto $frame
  set movesel [atomselect 0 "all" frame $frame]
  set fibersel [atomselect 0 "resname CSH CSSC" frame $frame]
  set fibcom [measure center $fibersel]
  set movelist {}
  foreach element $fibcom {
    lappend movelist [expr {$element * -1}]
  }
  $movesel moveby   $movelist
  set fibcom [measure center $fibersel]
  set movelist {}
  foreach element $fibcom {
    lappend movelist [expr {$element * -1}]
  }
  
  $movesel moveby   $movelist

}
set sel [atomselect 0 "all"]
animate write dcd ${outputpre}/${outputpre}.dcd beg 0 end -1 


for {set frame 0} {$frame < [molinfo top get numframes]} {incr frame} {
  animate goto $frame 
  # load a molecule
  # do the magic
  set newmol [::TopoTools::replicatemol 0 1 6 3 ]
  [atomselect top all] writexyz ${outputpre}/${outputpre}_${frame}.xyz
  mol delete top
}
exit


  