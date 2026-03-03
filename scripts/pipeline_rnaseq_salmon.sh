#!/bin/bash

# ==============================================================================
# PARALLEL SALMON-BASED RNA-SEQ PIPELINE
# ==============================================================================
# Hardware: 48 cores, 252GB RAM
# Runs 3 samples concurrently with individual logs
# ==============================================================================

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Input Files
GENOME_FA="/home/matheus/R570_reference/sugarcane/assembly/SofficinarumxspontaneumR570_771_v2.0.fa.gz"
GENOME_GTF="/home/matheus/R570_reference/sugarcane/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf"
IDS_FILE="/home/matheus/SB100/sugarcane/rna/ids.txt"

# Destination Paths
PRIMARY_DEST="/mnt/hd02/sb100/sugarcane/rna/results"
BACKUP_DEST="/mnt/hd01/sb100/sugarcane/rna/results"
INDEX_STORE="/mnt/hd02/sb100/sugarcane/rna/indices"

# Working Directory
WORK_DIR="/mnt/hd02/sb100/sugarcane/rna/work"

# Log Files
LOG_DIR="${WORK_DIR}/logs"
LOG_SUCCESS="${LOG_DIR}/EXECUTED_SUCCESS.txt"
LOG_FAIL="${LOG_DIR}/EXECUTED_FAILED.txt"
LOG_RUNNING="${LOG_DIR}/CURRENTLY_RUNNING.txt"

# Parallel Configuration
MAX_PARALLEL=3
CPUS_PER_JOB=44         # Burn!!!
MEMORY_PER_JOB="240GB"  # Burn!!

# Create directories
echo "Creating directories..."
mkdir -p "$PRIMARY_DEST" || { echo "ERROR: Cannot create PRIMARY_DEST"; exit 1; }
mkdir -p "$BACKUP_DEST"  || { echo "ERROR: Cannot create BACKUP_DEST"; exit 1; }
mkdir -p "$INDEX_STORE"  || { echo "ERROR: Cannot create INDEX_STORE"; exit 1; }
mkdir -p "$WORK_DIR"     || { echo "ERROR: Cannot create WORK_DIR"; exit 1; }
mkdir -p "$LOG_DIR"      || { echo "ERROR: Cannot create LOG_DIR"; exit 1; }

# Initialize log files
touch "$LOG_SUCCESS" "$LOG_FAIL" "$LOG_RUNNING"

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

check_disk_space() {
    local path=$1
    local required_gb=$2
    local available=$(df -BG "$path" 2>/dev/null | awk 'NR==2 {print $4}' | sed 's/G//')
    
    if [ -z "$available" ]; then
        echo "ERROR: Cannot check disk space for $path"
        return 1
    fi
    
    if [ "$available" -lt "$required_gb" ]; then
        echo "ERROR: Insufficient disk space on $path. Required: ${required_gb}GB, Available: ${available}GB"
        return 1
    fi
    echo "Disk space OK on $path: ${available}GB available"
    return 0
}

cleanup_failed_run() {
    local sra_id=$1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleaning up failed run for $sra_id..." >> "${LOG_DIR}/${sra_id}.log"
    
    # Log failure
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $sra_id" >> "$LOG_FAIL"
    
    # Remove from running log
    sed -i "/^${sra_id}$/d" "$LOG_RUNNING"

    # Cleanup temporary files
    rm -f "${WORK_DIR}/temp_input_${sra_id}.csv"
    rm -rf "${WORK_DIR}/temp_work_${sra_id}"
    rm -rf "${WORK_DIR}/results_${sra_id}"
    rm -rf "${WORK_DIR}/fetchngs_work_${sra_id}"
    rm -rf "${WORK_DIR}/rnaseq_work_${sra_id}"
    rm -rf "${WORK_DIR}/nf_home_${sra_id}"
}

# ==============================================================================
# PROCESS SINGLE SAMPLE FUNCTION
# ==============================================================================

process_sample() {
    local SRA_ID=$1
    local SAMPLE_NUM=$2
    local TOTAL=$3
    local LOG_FILE="${LOG_DIR}/${SRA_ID}.log"
    
    # Create isolated Nextflow home directory for this sample
    local NF_HOME="${WORK_DIR}/nf_home_${SRA_ID}"
    mkdir -p "$NF_HOME"
    
    # Export environment variables for this process
    export NXF_HOME="$NF_HOME"
    export NXF_WORK="${WORK_DIR}/nxf_work_${SRA_ID}"
    export NXF_TEMP="${WORK_DIR}/nxf_temp_${SRA_ID}"
    mkdir -p "$NXF_WORK" "$NXF_TEMP"
    
    # Redirect all output to log file
    exec > >(tee -a "$LOG_FILE") 2>&1
    
    echo "================================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] PROCESSING SAMPLE $SAMPLE_NUM/$TOTAL: $SRA_ID"
    echo "================================================================"
    echo "Nextflow home: $NF_HOME"
    echo "Work directory: $NXF_WORK"
    
    # Add to running log
    echo "$SRA_ID" >> "$LOG_RUNNING"
    
    # Skip if already processed
    if [ -d "$PRIMARY_DEST/$SRA_ID" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample $SRA_ID already processed. Skipping."
        sed -i "/^${SRA_ID}$/d" "$LOG_RUNNING"
        return 0
    fi
    
    # -------------------------------------------------------------------------
    # STEP 1: INDEX MANAGEMENT
    # -------------------------------------------------------------------------
    if [ -d "${INDEX_STORE}/salmon" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using existing Salmon index"
        INDEX_FLAGS="--salmon_index ${INDEX_STORE}/salmon"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Will create Salmon index"
        INDEX_FLAGS="--save_reference"
    fi
    
    # -------------------------------------------------------------------------
    # STEP 2: FETCH SRA DATA
    # -------------------------------------------------------------------------
    TEMP_DIR="${WORK_DIR}/temp_work_${SRA_ID}"
    TEMP_ID_FILE="${WORK_DIR}/temp_input_${SRA_ID}.csv"
    echo "$SRA_ID" > "$TEMP_ID_FILE"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading SRA data..."
    
    set +e
    nextflow run nf-core/fetchngs \
        -profile docker \
        -w "${WORK_DIR}/fetchngs_work_${SRA_ID}" \
        --input "$TEMP_ID_FILE" \
        --outdir "$TEMP_DIR" \
        --nf_core_pipeline rnaseq \
        --max_cpus "$CPUS_PER_JOB" \
        --max_memory "$MEMORY_PER_JOB" \
        -name "fetchngs_${SRA_ID}_$(date +%s)"
    
    FETCHNGS_EXIT=$?
    set -e
    
    rm -f "$TEMP_ID_FILE"

    if [ $FETCHNGS_EXIT -ne 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: fetchngs failed with exit code $FETCHNGS_EXIT"
        cleanup_failed_run "$SRA_ID"
        return 1
    fi
    
    SAMPLESHEET="${TEMP_DIR}/samplesheet/samplesheet.csv"
    if [ ! -f "$SAMPLESHEET" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Samplesheet not created"
        cleanup_failed_run "$SRA_ID"
        return 1
    fi
    
    # -------------------------------------------------------------------------
    # STEP 3: RUN SALMON PIPELINE
    # -------------------------------------------------------------------------
    RESULTS_DIR="${WORK_DIR}/results_${SRA_ID}"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Salmon quantification..."
    
    set +e
    nextflow run nf-core/rnaseq \
        --input "$SAMPLESHEET" \
        --outdir "$RESULTS_DIR" \
        --fasta "$GENOME_FA" \
        --gtf "$GENOME_GTF" \
        --skip_alignment \
        --pseudo_aligner salmon \
        --extra_salmon_quant_args '--seqBias --gcBias --numGibbsSamples 30 --validateMappings' \
        --skip_qc false \
        --skip_fastqc false \
        --skip_rseqc false \
        --skip_multiqc false \
        --max_memory "$MEMORY_PER_JOB" \
        --max_cpus "$CPUS_PER_JOB" \
        $INDEX_FLAGS \
        -w "${WORK_DIR}/rnaseq_work_${SRA_ID}" \
        -profile docker \
        -name "rnaseq_${SRA_ID}_$(date +%s)"
        
    RNASEQ_EXIT=$?
    set -e
    
    if [ $RNASEQ_EXIT -ne 0 ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: RNA-Seq pipeline failed with exit code $RNASEQ_EXIT"
        cleanup_failed_run "$SRA_ID"
        return 1
    fi
    
    # -------------------------------------------------------------------------
    # STEP 4: SAVE INDEX & ARCHIVE
    # -------------------------------------------------------------------------
    if [[ "$INDEX_FLAGS" == "--save_reference" ]] && [ -d "${RESULTS_DIR}/genome/index/salmon" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Saving Salmon index..."
        mv "${RESULTS_DIR}/genome/index/salmon" "${INDEX_STORE}/"
    fi
    
    SALMON_OUTPUT="${RESULTS_DIR}/salmon/salmon.merged.gene_counts.tsv"
    
    if [ -f "$SALMON_OUTPUT" ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline successful! Archiving..."
        
        # Log success
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $SRA_ID" >> "$LOG_SUCCESS"
        
        # Remove from running log
        sed -i "/^${SRA_ID}$/d" "$LOG_RUNNING"

        # Move to primary destination
        mv "$RESULTS_DIR" "$PRIMARY_DEST/$SRA_ID"
        
        # Asynchronous archiving
        (
            sleep 120
            tar -czf "${PRIMARY_DEST}/${SRA_ID}.tar.gz" -C "$PRIMARY_DEST" "${SRA_ID}"
            cp "${PRIMARY_DEST}/${SRA_ID}.tar.gz" "$BACKUP_DEST/"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Background archive completed" >> "$LOG_FILE"
        ) &
        
        # Cleanup
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleaning up temporary files..."
        rm -rf "$TEMP_DIR"
        rm -rf "${WORK_DIR}/fetchngs_work_${SRA_ID}"
        rm -rf "${WORK_DIR}/rnaseq_work_${SRA_ID}"
        rm -rf "$NF_HOME"
        rm -rf "$NXF_WORK"
        rm -rf "$NXF_TEMP"
        rm -f "$TEMP_ID_FILE"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: $SRA_ID completed"
        return 0
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Output file not found"
        cleanup_failed_run "$SRA_ID"
        return 1
    fi
}

# ==============================================================================
# PRE-FLIGHT CHECKS
# ==============================================================================

echo "================================================================"
echo "PARALLEL RNA-Seq Pipeline with Salmon"
echo "================================================================"
echo "System: $(uname -a)"
echo "Available CPUs: $(nproc)"
echo "Parallel jobs: $MAX_PARALLEL (${CPUS_PER_JOB} CPUs each)"
echo "Memory per job: $MEMORY_PER_JOB"
echo "================================================================"

# Check files
if [ ! -f "$GENOME_FA" ]; then echo "ERROR: GENOME_FA not found"; exit 1; fi
if [ ! -f "$GENOME_GTF" ]; then echo "ERROR: GENOME_GTF not found"; exit 1; fi
if [ ! -f "$IDS_FILE" ]; then echo "ERROR: IDS_FILE not found"; exit 1; fi
if [ ! -s "$IDS_FILE" ]; then echo "ERROR: IDS_FILE is empty"; exit 1; fi

# Check disk space & tools
check_disk_space "$WORK_DIR" 300 || exit 1

if ! docker info > /dev/null 2>&1; then echo "ERROR: Docker is not running"; exit 1; fi
if ! command -v nextflow &> /dev/null; then echo "ERROR: Nextflow not found"; exit 1; fi

echo "All pre-flight checks passed. Starting pipeline..."

# ==============================================================================
# MAIN PIPELINE LOOP - PARALLEL EXECUTION
# ==============================================================================

TOTAL_SAMPLES=$(wc -l < "$IDS_FILE")
CURRENT_SAMPLE=0
ACTIVE_JOBS=0

# Read SRA IDs into array
mapfile -t SRA_IDS < "$IDS_FILE"

for SRA_ID in "${SRA_IDS[@]}"; do
    # Skip empty lines
    [ -z "$SRA_ID" ] && continue
    
    ((CURRENT_SAMPLE++))
    
    # Wait if we have max parallel jobs running
    while [ $ACTIVE_JOBS -ge $MAX_PARALLEL ]; do
        sleep 10
        # Count running jobs
        ACTIVE_JOBS=$(jobs -r | wc -l)
    done
    
    echo "================================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] LAUNCHING: $SRA_ID ($CURRENT_SAMPLE/$TOTAL_SAMPLES)"
    echo "Active jobs: $ACTIVE_JOBS/$MAX_PARALLEL"
    echo "Log file: ${LOG_DIR}/${SRA_ID}.log"
    echo "================================================================"
    
    # Launch in background with isolated environment
    (process_sample "$SRA_ID" "$CURRENT_SAMPLE" "$TOTAL_SAMPLES") &
    
    ((ACTIVE_JOBS++))
    
    # Small delay to avoid race conditions
    sleep 5
done

# Wait for all jobs to complete
echo "================================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] All samples launched. Waiting for completion..."
echo "================================================================"
wait

echo "================================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] PIPELINE COMPLETED"
echo "================================================================"
echo "Success log: $LOG_SUCCESS"
echo "Failed log: $LOG_FAIL"
echo "Individual logs: ${LOG_DIR}/*.log"
echo "================================================================"

# Generate summary
echo ""
echo "SUMMARY:"
echo "--------"
echo "Total samples: $TOTAL_SAMPLES"
echo "Successful: $(wc -l < "$LOG_SUCCESS" 2>/dev/null || echo 0)"
echo "Failed: $(wc -l < "$LOG_FAIL" 2>/dev/null || echo 0)"
echo ""
echo "Background compression jobs may still be running."
