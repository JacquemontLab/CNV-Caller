#!/bin/bash

#======================================================================
# Batch parallel CNV caller with cpu managment. Relies on cnv_caling.sh 
#======================================================================
usage() {
  echo "Usage: $0 --batch_list FILE --sexfile FILE --pfb FILE --gcmodel FILE --gcdir DIR --hmm_file FILE --levels FILE --config FILE --chr CHR [--mode MODE] [--cpus INT] --autosome_only"
  echo ""
  echo "Required options:"
  echo "  --batch_list    FILE    Text file with one file path per line"
  echo "  --sexfile       FILE    .tsv file with individual information"
  echo "  --pfb           FILE    .tsv file with variant-level scores"
  echo "  --gcmodel       FILE    Genome file with GC frequency"
  echo "  --gcdir         DIR     Directory with per-chromosome GC frequency"
  echo "  --hmm_file      FILE    Hidden Markov Model parameters file"
  echo "  --levels        FILE    Levels file for analysis"
  echo "  --config        FILE    Configuration file"
  echo "  --chr           CHR     Chromosome to analyze (e.g. 1, 2, X)"
  echo ""
  echo "Optional:"
  echo "  --mode          MODE    'taskset' (default) or 'parallel'"
  echo "  --cpus          INT     Number of CPUs to use (overrides SLURM_CPUS_ON_NODE or nproc)"
  echo "  --help                 Show this help message and exit"
  echo "  --autosome_only FLAG  Only analyze autosomes (true/false)"
  exit 1
}

# Parse Options
cpus=""
mode="taskset"  # default

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch_list) batch_list="$2"; shift 2 ;;
    --autosome_only) autosome_only=1; shift ;;
    --sexfile) sexfile="$2"; shift 2 ;;
    --pfb) pfb="$2"; shift 2 ;;
    --gcmodel) gcmodel="$2"; shift 2 ;;
    --gcdir) gcdir="$2"; shift 2 ;;
    --hmm_file) hmm_file="$2"; shift 2 ;;
    --levels) levels="$2"; shift 2 ;;
    --config) config="$2"; shift 2 ;;
    --chr) chr="$2"; shift 2 ;;
    --mode) mode="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Check required parameters
if [[ -z "$batch_list" || -z "$sexfile" || -z "$pfb" || -z "$gcmodel" || -z "$gcdir" || -z "$hmm_file" || -z "$levels" || -z "$config" || -z "$chr" ]]; then
  echo "❌ Error: All required options must be specified."
  usage
fi


mode="${mode:-taskset}"

# Detect CPU count
if [[ -z "$cpus" ]]; then
  cpus="${SLURM_CPUS_ON_NODE:-$(nproc)}"
fi

echo "💻 Running with $cpus cores"

# Track core usage and PIDs
declare -A core_pid

# Function to check for free cores
wait_for_core() {
    while true; do
        for ((c=0; c<cpus; c++)); do
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



if [[ "$mode" == "taskset" ]]; then
  echo "⚙️ Running with taskset (manual CPU binding)"

  while IFS= read -r file; do
      core=$(wait_for_core)
      sample_id=$(basename "$file" | cut -d "." -f1)
      echo "💻 Assigning SampleID $sample_id to CPU core $core"

      taskset -c $core cnv_calling.sh \
        $([ "$autosome_only" -eq 1 ] && echo "--autosome_only") \
        --sample_id "$sample_id" \
        --BAF_LRR_Probes "$file" \
        --sexfile "$sexfile" \
        --pfb "$pfb" \
        --hmm "$hmm_file" \
        --gcmodel "$gcmodel" \
        --chr "$chr" \
        --levels "$levels" \
        --config "$config" \
        --gcdir "$gcdir" &

      core_pid["$core"]=$!
  done < "$batch_list"

  wait

elif [[ "$mode" == "parallel" ]]; then
  echo "⚙️ Running with GNU parallel using $cpus threads"

  export sexfile pfb gcmodel chr hmm_file levels config gcdir PATH

  cat "$batch_list" | parallel -j "$cpus" --eta --line-buffer '
    file={}
    sample_id=$(basename "$file" | cut -d "." -f1)
    echo "🔄 Processing $sample_id"

    cnv_calling.sh '"$([ "$autosome_only" -eq 1 ] && echo "--autosome_only")"'  \
      --sample_id "$sample_id" \
      --BAF_LRR_Probes "$file" \
      --sexfile "$sexfile" \
      --pfb "$pfb" \
      --hmm "$hmm_file" \
      --gcmodel "$gcmodel" \
      --chr "$chr" \
      --levels "$levels" \
      --config "$config" \
      --gcdir "$gcdir"
  ' 2> parallel.progress.log

else
  echo "❌ Unknown mode: $mode"
  usage
fi
