# set step number
set step_num 50
# set the bin width for the histogram
set bin_width 5.0
# set the number of parallel processes to use
set num_procs 19
# Set the name of the PSF file
set psf_file "/dfs2/tw/yuanmis1/mrsec/cyl1/mix_27A/setup/CSH_CSSC_mix_27A_cyl_water_ions.psf"
# set namelist for measure dihed angle
set namelist {C3 S1 S2 C11}
# set dcd file prefix
set dcd_prefix "/dfs2/tw/yuanmis1/mrsec/cyl1/mix_27A/nvt"
# create an empty list for dcd file names
set dcd_files {}
# loop through possible dcd file names and add existing ones to list
for {set i 1} {$i <= 99} {incr i} {
  if {$i < 10} {
    set dcd_name "${dcd_prefix}0$i.dcd"
  } else {
    set dcdname "${dcd_prefix}$i.dcd"
  }
  if {[file exists $dcd_name]} {
    lappend dcd_files $dcd_name
  } else {
    break
  }
} 
set dcd_files {/dfs2/tw/yuanmis1/mrsec/cyl1/mix_27A/nvt01.dcd}
# function to update histogram table with dihedral angle measurements
# create an empty histogram table
set hist_table_inside {}
set hist_table_surface {}
# load each dcd file with step set earlier
mol new $psf_file
package require pbctools
proc get_histogram_inside {start end bin_width} {
  # create an empty histogram table
  set hist_table_inside {}
  # loop through each frame
  for {set currframe $start} {$currframe <= $end} {incr currframe} {
    animate goto $currframe
    #set sel_surface [atomselect top "resname CSSC and within 5 of water"]
    set sel_inside [atomselect top "resname CSSC and not within 5 of water"]
    
    # Loop through each residue in the surface selection
    # foreach resid [lsort -unique [$sel_surface get resid]] {
    #   set sel_res [atomselect top "resname CSSC and resid $resid and name $namelist"]
    #   set sel_atoms [$sel_res get index]
    #   # Make sure we have 4 atoms to calculate dihedral angle
    #   set dihed [measure dihed $sel_atoms]
    #   set bin [expr int($dihed/$bin_width)]
    #   if {[dict exists $hist_table_surface $bin]} {
    #     dict incr hist_table_surface $bin
    #   } else {
    #     dict set hist_table_surface $bin 1
    #   }
    # }
    
    # Loop through each residue in the inside selection
    foreach resid [lsort -unique [$sel_inside get resid]] {
      set sel_res [atomselect top  "resname CSSC and resid $resid and name $namelist"]
      set sel_atoms [$sel_res get index]
      set dihed [measure dihed $sel_atoms]
      set bin [expr int($dihed/$bin_width)]
      if {[dict exists $hist_table_inside $bin]} {
        dict incr hist_table_inside $bin
      } else {
        dict set hist_table_inside $bin 1
      }
    }
    return $hist_table_inside
  }
}
set hist_table_inside_all {}
foreach dcd_file $dcd_files {
  animate delete all
  animate read dcd $dcd_file skip $step_num waitfor all top
  pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all
  set totframes [molinfo top get numframes]
  # set chunk size
  set chunk_size [expr int($totframes/$num_procs)+1]
  # split the frames into chunks
  set chunks {}
  for {set i 0} {$i < $num_procs} {incr i} {
    set start [expr {$i * $chunk_size}]
    set end [expr {($i + 1) * $chunk_size - 1}]
    lappend chunks [list $start $end]
  }
  foreach chunk $chunks {
  set start_frame [lindex $chunk 0]
  set end_frame [lindex $chunk 1]
  set chunk_hist_table [parallel get_histogram_inside $start_frame $end_frame $bin_width]
}
puts $chunk_hist_table
# open two output files
set outfile_surface [open "data/dihed_surface.dat" w]
set outfile_inside [open "data/dihed_inside.dat" w]
foreach {key value} [dict get $hist_table_surface] {
  puts $outfile_surface "$key $value"
}
foreach {key value} [dict get $hist_table_inside] {
  puts $outfile_inside "$key $value"
}
# close output files
close $outfile_surface
close $outfile_inside


