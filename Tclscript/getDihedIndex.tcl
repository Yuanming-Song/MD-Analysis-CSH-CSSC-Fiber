# set step number
set step_num 20
# set the bin width for the histogram
set bin_width 5.0
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
    set dcd_name "${dcd_prefix}$i.dcd"
  }
  if {[file exists $dcd_name]} {
    lappend dcd_files $dcd_name
  } else {
    break
  }
} 
# function to update histogram table with dihedral angle measurements
# load each dcd file with step set earlier
mol new $psf_file
package require pbctools
set outfile_surface [open "data/dihed_surface.ind" w]
set outfile_inside [open "data/dihed_inside.ind" w]
foreach dcd_file $dcd_files {
  animate delete all
  animate read dcd $dcd_file skip $step_num waitfor all top
  pbc wrap -centersel "resname CSSC or resname CSH" -center com -compound residue -all
  set totframes [molinfo top get numframes]
  # loop through each frame
  for {set currframe 0} {$currframe < $totframes} {incr currframe} {
    animate goto $currframe
    set sel_surface [atomselect top "resname CSSC and name S1 S2 and within 5 of water"]
    set index_surface [lsort -unique [$sel_surface get resid]]
    puts $outfile_surface $index_surface
    set sel_inside [atomselect top "resname CSSC and name S1 S2 and not within 5 of water"]
    set index_inside [lsort -unique [$sel_inside get resid]]
    puts $outfile_inside $index_inside
  }
}
close $outfile_surface
close $outfile_inside
exit
