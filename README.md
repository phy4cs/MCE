# MCE
MCE and CCS simulation - Chris Symonds
This program allows a wide array of simulations using the MCE and CCS equations
To run, use the run.sh file in the run folder, with the following arguments:

$1 - Total number of repeat calculations desired
$2 - Number of parallel threads per folder
$3 - Number of folders to split the job into.

For example, to run 128 jobs using 4 parallel cores, with the load split into 4 folders,
you would use the command 

                ./run.sh 128 4 4
               
which would have a total of 16 parallel threads running simultaneously. This allows openmp
execution to be carried out over many more cores than would be present on a single node, which
on arc1 has a maximum of 8, on arc 2 a maximum of 16.

If restarts are needed due to server time limits, use the restart.sh script which uses the same 
arguments as the run.sh script. This file should be called from the running folder of the most recent
partial run. 

To ensure that you are using the latest version of the code, make sure that you download it from 
github, using the url

            https://github.com/phy4cs/MCE
            
If you are having trouble accessing it, contact me (email : phy4cs@leeds.ac.uk)

The run.sh script creates a second script, called result.sh which when run calls the collate.sh script
which combines the results from all the completed runs. This script deletes the original output files
so use with care!

If multiple partial runs are used, the combine.sh script should be used which requires the file 
<<folderlist.dat>>, a file containing a list of the required folders which should be in order. This file is 
automatically created by the restart.sh script, however care should be taken to ensure that no confusion 
occurs when there are multiple simulations happening at the same time. This script calls the collate.sh 
file so again be aware that the original output files will be deleted.

Any problems, contact me.

