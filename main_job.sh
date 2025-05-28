#!/bin/bash

declare -a threads=(1) #(1 2 4 5)

declare -a stages=(1) #(1 2 4 5)

for THREADS in "${threads[@]}"
do
    for STAGES in "${stages[@]}"
    do
        export OMP_NUM_THREADS=${THREADS}
        SAMPLES_PER_THREAD=$((10000 / THREADS))
    
        echo "Running with THREADS=${THREADS}, STAGES=${STAGES}, SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}"
        
        SCRIPT_DIR="job_scripts"
        mkdir -p "$SCRIPT_DIR"
        
        JOB_FILE="${SCRIPT_DIR}/job_${THREADS}_${STAGES}.sh"
    
        cat > "$JOB_FILE" << EOF
#!/bin/bash
#SBATCH --job-name=run_${THREADS}
#SBATCH --nodes=1
#SBATCH --ntasks=${THREADS}
#SBATCH --time=04:00:00

export SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}

make clean
make && make run THREADS=${THREADS} STAGES=${STAGES} SAMPLES_PER_THREAD=\$SAMPLES_PER_THREAD
EOF
    
        chmod +x "$JOB_FILE"
        
        # Submit job and capture Job ID
        jobid=$(sbatch "$JOB_FILE" | awk '{print $4}')
        echo "Submitted job slurm-$jobid for THREADS=$THREADS"
    
        scontrol show job "$jobid"
        
        # Wait for the job to complete before continuing
        while squeue -j "$jobid"; do
            echo "Job slurm-$jobid is still running..."
            # Print detailed status of the job
            squeue -j "$jobid"   # This prints the job's status
            sleep 10
        done
    
        echo "Job slurm-$jobid completed. Proceeding to next."
    done
done
