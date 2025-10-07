#!/bin/bash

module load conda
module load use.own
module load conda-env/mypackages-py3.12.5
conda activate mypackages
module load gnuplot

tuples=("(5 15 14 12)" "(4 20 13 15)" "(4 16 17 16)" "(3 15 15 12)" "(3 15 27 25)" "(2 28 17 28)" "(2 16 21 21)" "(2 10 40 21)")

myfile="fit_parameters_and_errors.txt"
mkdir -p "$(dirname "$myfile")"
> "$myfile"

# Iterate through each tuple by name
for t in "${tuples[@]}"; do
    
    t="${t//[\(\)]/}"  # Remove parentheses
    read -r VERTICES TRACKS SWEEPS DA_SWEEPS <<< "$t"

    declare -a threads=(1) #25 50 100)
    declare -a stages=(1) #2 5 10)
    declare -a sweeps=(${SWEEPS}) #(100 150) #(10 50 100 150) 
    
    # Define the directory where AltMakefile and source files are located
    PROJECT_DIR=$(pwd) # Adjust this if AltMakefile is in a different directory
    
    mkdir -p slurm
    
    for THREADS in "${threads[@]}"
    do
        for STAGES in "${stages[@]}"
        do
            for SWEEPS in "${sweeps[@]}"
            do
                if (( SWEEPS % STAGES == 0 && SWEEPS != STAGES)); then
                    
                    #DA_SWEEPS=${SWEEPS}
                    DA_SWEEPS=${DA_SWEEPS}
                
                    export OMP_PROC_BIND=TRUE
                    export OMP_PLACES=cores
                    export OMP_NUM_THREADS=${THREADS}
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
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
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
make -C ${PROJECT_DIR} -f \${makefile} && make -C ${PROJECT_DIR} -f \${makefile} run VERTICES=${VERTICES} TRACKS=${TRACKS} THREADS=${THREADS} STAGES=${STAGES} SWEEPS=${SWEEPS} DA_SWEEPS=${DA_SWEEPS} SAMPLES_PER_THREAD=${SAMPLES_PER_THREAD}
EOF
        
                    chmod +x "$JOB_FILE"
                    
                    # Submit job and capture Job ID
                    jobid=$(sbatch "$JOB_FILE" | awk '{print $4}')
                    echo "Submitted job slurm-$jobid for VERTICES=$VERTICES, TRACKS=$TRACKS, THREADS=$THREADS, STAGES=$STAGES, SWEEPS=$SWEEPS, DA_SWEEPS=$SDA_SWEEPS"
                    
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

rm -f "slurm-${SLURM_JOB_ID}.out" >/dev/null 2>&1

conda deactivate
