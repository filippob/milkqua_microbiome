In order to use SLURM on INDACO, you need to generate a script containing the following informations:

#!/bin/bash
#SBATCH -o /gpfs/home/users/nome.cognome/test/myjob.%j.%N.out (%j codify the job ID in output file, while %N code for job master node and needs to be added as ID could repeat themselves in differente clusters).
#SBATCH -D /gpfs/home/users/nome.cognome (job wd)
#SBATCH -J Prova (job name, max 10 typos)
#SBATCH --get-user-env
#SBATCH -p light 
#SBATCH --nodes=1
#SBATCH -c 1 # Number of cores per task
#SBATCH --mem-per-cpu=1000 # Mem for each core
#SBATCH --mail-type=all (you can choose, begin or end or fail or all)
#SBATCH --mail-user=nome.cognome@unimi.it
#SBATCH --time=0:1:00 (max is 8hours)
#SBATCH --account=nomeprogetto

module load sss_example (to upload requested modules in the environment)

cd /gpfs/home/users/nome.cognome/test
./my_prog_seriale (program to be executed)


Once your script is ready, you can submit it with the command sbatch.

sbatch myjob.slurm.sh

You can verify your job status with the followings:

squeue –clusters=all: verify job status
squeue -u <username>: list jobs for a user
squeue -u <username> - t RUNNING: list ongoing jobs for a user
squeue -u <username> -t PENDING: list pending jobs for a user
scancel <jobid>: delete a job
scancel -u <username> : delete jobs by user
scancel -t PENDING -u <username>:  delete pending jobs for a user
scancel –name myJobName: delete a job using its name
scontrol hold <jobid>: pause a job
scontrol resume <jobid>: restart a paused job
scontrol requeue <jobid><username>: cancel and rerun a job



