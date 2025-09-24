#!/bin/bash

tuples=("(5 15 14 12)" "(4 20 13 14)" "(4 16 17 16)" "(3 15 15 12)" "(2 28 17 28)" "(2 16 21 21)" "(2 10 40 21)")

# Iterate through each tuple by name
for t in "${tuples[@]}"; do
    
    t="${t//[\(\)]/}"  # Remove parentheses
    read -r VERTICES TRACKS SWEEPS DA_SWEEPS <<< "$t"

    declare -a threads=(1) #20 50)
    
    declare -a stages=(1) #5 10)
    
    declare -a sweeps=(10 50 100 150) #(${SWEEPS})
    
    mkdir -p slurm
    
    for THREADS in "${threads[@]}"
    do
        for STAGES in "${stages[@]}"
        do
            for SWEEPS in "${sweeps[@]}"
            do
                if (( SWEEPS % STAGES == 0 && SWEEPS != STAGES)); then
                
                    DA_SWEEPS=${SWEEPS}
                    #DA_SWEEPS=${DA_SWEEPS}
                    
                    output_file="ConvergenceEfficiency_and_TimePerAnneal.txt"
                    >"$output_file"
    
                    output_file2="ConvergenceEfficiency_and_DunnIndex_Binning.txt"
                    >"$output_file2"
                
                    export OMP_PROC_BIND=TRUE       # Enables binding threads to CPUs
                    export OMP_PLACES=cores         # Maps each thread to a distinct core
                    export OMP_NUM_THREADS=${THREADS}  # Set number of threads to your cores count
                    SAMPLES_PER_THREAD=$((10000 / THREADS))
                
                    echo "Running with ${VERTICES} VERTICES and ${TRACKS} TRACKS, THREADS=${THREADS}, STAGES=${STAGES}, SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}, SWEEPS=${SWEEPS}, DA_SWEEPS=${DA_SWEEPS}"
                    
                    SCRIPT_DIR="job_scripts"
                    mkdir -p "$SCRIPT_DIR"
                    
                    JOB_FILE="${SCRIPT_DIR}/job_${VERTICES}_${TRACKS}_${THREADS}_${STAGES}_${SWEEPS}_${DA_SWEEPS}.sh"
                
                    cat > "$JOB_FILE" << EOF
#!/bin/bash
#SBATCH --account=physics
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --job-name=run_${THREADS}
#SBATCH --nodes=1
#SBATCH --ntasks=${THREADS}
#SBATCH --time=02:00:00

export SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}

make clean
make && make run VERTICES=${VERTICES} TRACKS=${TRACKS} THREADS=${THREADS} STAGES=${STAGES} SWEEPS=${SWEEPS} DA_SWEEPS=${DA_SWEEPS} SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}
EOF
        
                    chmod +x "$JOB_FILE"
                    
                    # Submit job and capture Job ID
                    jobid=$(sbatch "$JOB_FILE" | awk '{print $4}')
                    echo "Submitted job slurm-$jobid for VERTICES=$VERTICES, TRACKS=$TRACKS, THREADS=$THREADS , STAGES=$STAGES , SWEEPS=$SWEEPS, DA_SWEEPS=$DA_SWEEPS"
                    
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
                fi
            done
        done
    done
    
    unset threads
    
    unset stages
    
    unset sweeps
done
