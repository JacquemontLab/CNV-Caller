#!/usr/bin/env bash 
    
#======================================================================
#Batch parallel CNV caller with cpu managment. Relies on cnv_caling.sh 
#======================================================================
usage() {
  echo "Usage: $0 --batch_list FILE --plink_data FILE --pfb FILE --gcmodel FILE --gcdir DIR"
  echo ""
  echo "Required options:"
  echo "  --batch_list   FILE    Text file with one file path per line"
  echo "  --plink_data   FILE    .tsv file with individual information"
  echo "  --pfb          FILE    .tsv file with variant-level scores"
  echo "  --gcmodel      FILE    Genome file with GC frequency"
  echo "  --gcdir        DIR     Directory with per-chromosome GC frequency"
  echo "  --help                 Show this help message and exit"
  exit 1
}

#Parse Options
while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch_list)
      batch_list="$2"
      shift 2
      ;;
    --plink_data)
      plink_data="$2"
      shift 2
      ;;
    --pfb)
      pfb="$2"
      shift 2
      ;;
    --gcmodel)
      gcmodel="$2"
      shift 2
      ;;
    --gcdir)
      gcdir="$2"
      shift 2
      ;;
    --help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

# Check required args
if [[ -z "$batch_list" || -z "$plink_data" || -z "$pfb" || -z "$gcmodel" || -z "$gcdir" ]]; then
  echo "Error: All required options must be specified."
  usage
fi



# Get number of available cores
ncores=$(nproc)

# Track core usage and PIDs
declare -A core_pid

# Function to check for free cores
wait_for_core() {
    while true; do
        for ((c=0; c<ncores; c++)); do
            pid="${core_pid[$c]}" # Get the PID currently assigned to core $c
            # Check if the core is free:
            # -z "$pid"                → true if there's no PID assigned yet (i.e., this core hasn't been used)
            # ! -e /proc/$pid          → true if the PID is not running anymore (no such process in /proc)
            if [[ -z "$pid" || ! -e /proc/$pid ]]; then
                echo "$c" # This core is available — return it
                return
            fi
        done
        sleep 0.5
    done
}



# Default parameters can be set here
chr="1:23"
hmm_file="/usr/local/bin/hhall.hmm"
levels="/usr/local/quantisnp/bin/config/levels.dat"
config="/usr/local/quantisnp/bin/config/params.dat"


    while IFS= read -r line; do
        core=$(wait_for_core)
        sample_id=$(basename "$line" | cut -d "." -f1)
        echo "Assigning IID "$sample_id" to CPU core $core"

        # Bind the process to a specific core
        taskset -c $core cnv_calling.sh \
        --sample_id "$sample_id"  \
        --BAF_LRR_Probes "$line" \
        --sexfile $plink_data \
        --pfb $pfb \
        --hmm "$hmm_file" \
        --gcmodel $gcmodel \
        --chr "$chr" \
        --levels "$levels" \
        --config "$config" \
        --gcdir $gcdir &

        # Record the PID for this core
        core_pid["$core"]=$!
    done < $batch_list
    
wait 