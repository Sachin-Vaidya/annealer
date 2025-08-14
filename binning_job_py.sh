#!/bin/bash

module load conda
module load use.own
module load conda-env/mypackages-py3.12.5
module load gnuplot

tuples=("(3 15)") #("(2 10)" "(2 16)" "(2 28)" "(3 15)" "(4 16)" "(4 20)" "(5 15)")

# Iterate through each tuple by name
for t in "${tuples[@]}"; do
    
    t="${t//[\(\)]/}"  # Remove parentheses
    read -r VERTICES TRACKS <<< "$t"

    declare -a threads=(1) #25 50 100)
    declare -a stages=(1) #2 5 10)
    declare -a sweeps=(27) #2 5 10 20)
    
    # Define the directory where AltMakefile and source files are located
    PROJECT_DIR=$(pwd) # Adjust this if AltMakefile is in a different directory
    MAKEFILE="${PROJECT_DIR}/AltMakefile"
    SRC_FILE="${PROJECT_DIR}/binning_the_efficiency_results.cpp"
    
    # Verify that AltMakefile and source file exist
    if [[ ! -f "$MAKEFILE" ]]; then
        echo "Error: AltMakefile not found at $MAKEFILE"
        exit 1
    fi
    if [[ ! -f "$SRC_FILE" ]]; then
        echo "Error: Source file not found at $SRC_FILE"
        exit 1
    fi
    
    mkdir -p slurm
    
    for THREADS in "${threads[@]}"
    do
        for STAGES in "${stages[@]}"
        do
            for SWEEPS in "${sweeps[@]}"
            do
                if (( SWEEPS % STAGES == 0 && SWEEPS != STAGES)); then
                    output_file="ConvergenceEfficiency_and_TimePerAnneal.txt"
                    >"$output_file"
    
                    output_file2="ConvergenceEfficiency_and_DunnIndex_Binning.txt"
                    >"$output_file2"
                
                    export OMP_PROC_BIND=TRUE
                    export OMP_PLACES=cores
                    export OMP_NUM_THREADS=${THREADS}
                    SAMPLES_PER_THREAD=$((10000 / THREADS))
                
                    echo "Running with ${VERTICES} VERTICES and ${TRACKS} TRACKS, THREADS=${THREADS}, STAGES=${STAGES}, SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}, SWEEPS=${SWEEPS}, DA_SWEEPS=${SWEEPS}"
                    
                    SCRIPT_DIR="job_scripts"
                    mkdir -p "$SCRIPT_DIR"
                    
                    JOB_FILE="${SCRIPT_DIR}/job_${VERTICES}_${TRACKS}_${THREADS}_${STAGES}_${SWEEPS}_${SWEEPS}.sh"
                
                    cat > "$JOB_FILE" << EOF
#!/bin/bash
#SBATCH --account=physics
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --job-name=run_${THREADS}
#SBATCH --nodes=1
#SBATCH --ntasks=${THREADS}
#SBATCH --time=00:10:00
#SBATCH --chdir=${PROJECT_DIR}  # Set the working directory for the Slurm job

export SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores
export OMP_NUM_THREADS=${THREADS}

makefile="AltMakefile_py"

# Verify Makefile exists
if [[ ! -f "\${makefile}" ]]; then
    echo "Error: \${makefile} not found in ${PROJECT_DIR}"
    exit 1
fi

# Run make commands with explicit directory
make -C ${PROJECT_DIR} -f \${makefile} clean
make -C ${PROJECT_DIR} -f \${makefile} && make -C ${PROJECT_DIR} -f \${makefile} run VERTICES=${VERTICES} TRACKS=${TRACKS} THREADS=${THREADS} STAGES=${STAGES} SWEEPS=${SWEEPS} DA_SWEEPS=${SWEEPS} SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}
EOF
        
                    chmod +x "$JOB_FILE"
                    
                    # Submit job and capture Job ID
                    jobid=$(sbatch "$JOB_FILE" | awk '{print $4}')
                    echo "Submitted job slurm-$jobid for VERTICES=$VERTICES, TRACKS=$TRACKS, THREADS=$THREADS, STAGES=$STAGES, SWEEPS=$SWEEPS, DA_SWEEPS=$SWEEPS"
                    
                    echo ""
                    
                    scontrol show job "$jobid"
                    
                    echo ""
                    
                    # Wait for the job to complete
                    while true; do
                        state=$(sacct -j "$jobid" --format=State --noheader | awk '{print $1}' | head -n 1)
                        
                        squeue -j "$jobid"
                        
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