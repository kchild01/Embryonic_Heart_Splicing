#Make sample list file
ls b*.txt | uniq > sample_list.txt
#Make script executeable (if not already)
chmod +x make_rMATS_slurms.sh
#Run MakeFileConcatenate_rbg.sh (see attached file below)
./make_rMATS_slurms.sh
#Make script to launch files
ls *rMATS_prep*.slurm | awk '{print "sbatch "$1}' > launch_rMATS_prep.sh
#Make script executeable
chmod +x launch_rMATS_prep.sh
#LAUNCH SCRIPT!!
./launch_rMATS_prep.sh
