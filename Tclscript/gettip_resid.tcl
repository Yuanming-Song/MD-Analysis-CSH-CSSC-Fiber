mol new ../setup/CSH_CSSC_mix_27A_cyl_water_ions.pdb
#read index
set fid [open "tips_index_6_to_10_every_25_.dat" r]
#get each line
set lines [split [read $fid] "\n"]
#how many lines do we have
set last_line [expr {[llength $lines] - 1}]
#close index file
close $fid
# open the output file for writing
set fid [open "tips_id_6_to_10_every_25.dat" w]
# loop through frames from the first to the last line of tips_index.txt
for {set frame 0} {$frame < $last_line} {incr frame} {
  # read index for this frame
  set selind [lindex $lines $frame]
  # select these ind
  if {[llength $selind] == 0} {
    puts $fid $selind
  } else {
  set sel [atomselect top "index $selind"]
  # get unique number from [$sel get resid]
  set resids [lsort -unique [$sel get resid]]
  # write the unique resids to the output file
  puts $fid $resids
  }
}

# close the output file
close $fid
