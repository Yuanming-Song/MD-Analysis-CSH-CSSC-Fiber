# VMD Tcl Script

# Get the data file from arguments
set data_file [lindex $argv 0]

# Extract the number at the end before .dat
regexp {.*_(\d+)\.dat$} $data_file match number

# Define the directory
set dcddir "/dfs9/tw/yuanmis1/mrsec/CSH-CSSC/cyl1/mix_27A/"

# Build the dcd file path
set dcd_file "${dcddir}/nvt${number}.dcd"

# Find the psf file containing 'ions' in the name
set psf_files [glob "${dcddir}/setup/*ions*.psf"]

# Load the psf file
mol new $psf_files type psf waitfor all

# Load the dcd file
mol addfile $dcd_file type dcd waitfor all

# Load the pbcwrap plugin
package require pbctools

# Perform pbc wrap
pbc wrap -sel "resname CSH CSSC" -compound fragment -all

# Open the data file
set fp [open $data_file r]

# Read all lines
set lines [split [read $fp] "\n"]

# Close the file
close $fp

# Iterate over lines, skip first line
foreach line [lrange $lines 1 end] {
    if {[string length $line] == 0} continue
    set fields [split $line " "]
    
    # Assign fields
    set donor_resname [lindex $fields 0]
    set donor_AM [lindex $fields 1]
    set acceptor_resname [lindex $fields 2]
    set acceptor_AM [lindex $fields 3]
    set donor_index [lindex $fields 4]
    set acceptor_index [lindex $fields 5] 
    set frame [expr {[lindex $fields 9] - 1}]
    
    # Go to frame
    animate goto $frame
    
    # Select donor and acceptor residues
    set res_donor [atomselect top "same residue as index [expr {$donor_index} - 1]"]
    set res_acceptor [atomselect top "same residue as index [expr {$acceptor_index} - 1]"]
    
    # Combine selections
    set sel_combined [atomselect top "same residue as index [expr {$donor_index} + 1] or same residue as index [expr {$acceptor_index} + 1]"]
    
    # Define output pdb file name
    set pdb_name "data/hbond_visual/${donor_resname}_${donor_AM}_${acceptor_resname}_${acceptor_AM}.pdb"
    
    # Write pdb file
    $sel_combined writepdb $pdb_name
}