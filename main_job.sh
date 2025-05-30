#!/bin/bash

declare -a threads=(1 5 10 20 40 50 80 100)

declare -a stages=(1) #(1 2 4 5)

declare -a sweeps=(5 10 20 25 40 50)

mkdir -p slurm

output_file="SA_ConvergenceEfficiency_and_TimePerAnneal.txt"
>"$output_file"

for THREADS in "${threads[@]}"
do
    for STAGES in "${stages[@]}"
    do
        for SWEEPS in "${sweeps[@]}"
        do
            export OMP_PROC_BIND=TRUE       # Enables binding threads to CPUs
            export OMP_PLACES=cores         # Maps each thread to a distinct core
            export OMP_NUM_THREADS=${THREADS}  # Set number of threads to your cores count
            SAMPLES_PER_THREAD=$((10000 / THREADS))
        
            echo "Running with THREADS=${THREADS}, STAGES=${STAGES}, SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}, SWEEPS=${SWEEPS}"
            
            SCRIPT_DIR="job_scripts"
            mkdir -p "$SCRIPT_DIR"
            
            JOB_FILE="${SCRIPT_DIR}/job_${THREADS}_${STAGES}_${SWEEPS}.sh"
        
            cat > "$JOB_FILE" << EOF
#!/bin/bash
#SBATCH --job-name=run_${THREADS}
#SBATCH --nodes=1
#SBATCH --ntasks=${THREADS}
#SBATCH --time=04:00:00

export SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}

make clean
make && make run THREADS=${THREADS} STAGES=${STAGES} SWEEPS=${SWEEPS} SAMPLES_PER_THREAD=\$SAMPLES_PER_THREAD
EOF
    
            chmod +x "$JOB_FILE"
            
            # Submit job and capture Job ID
            jobid=$(sbatch "$JOB_FILE" | awk '{print $4}')
            echo "Submitted job slurm-$jobid for THREADS=$THREADS , STAGES=$STAGES , SWEEPS=$SWEEPS"
            
            echo ""
            
            scontrol show job "$jobid"
            
            echo ""
            
            # Wait for the job to complete before continuing
            while true; do
                state=$(sacct -j "$jobid" --format=State --noheader | awk '{print $1}' | head -n 1)
                
                squeue -j "$jobid"
                #echo "Job slurm-$jobid status: $state"
                
                if [[ "$state" == "COMPLETED" || "$state" == "FAILED" || "$state" == "CANCELLED" || "$state" == "TIMEOUT" ]]; then
                    mv slurm-${jobid}.out slurm/
                    break
                fi
            
                sleep 10
            done
            
            echo ""
            
            echo "Job slurm-$jobid ($state) completed. Proceeding to next."
            
            echo ""
        done
    done
done
