#!/bin/bash

# --- Configuration ---
RAW_LSCPU_DIR="raw_lscpu"
SUMMARY_FILE="node_summary.csv"
REPORT_FILE="hpc_cluster_report.txt"
LOG_FILE="hpc_script.log" # For detailed messages

# List of nodes to skip (by number, e.g., "1" for compute-0-1)
SKIP_NODES=("7" "14" "15" "20" "25")

# --- Setup ---
mkdir -p "$RAW_LSCPU_DIR" || { echo "ERROR: Could not create directory $RAW_LSCPU_DIR. Exiting." | tee -a "$LOG_FILE"; exit 1; }
# Initialize summary file with header
# CORRECTED: Swapped CPU_Model and Clock_Speed_MHz in the header
echo "Node,Cores,Clock_Speed_MHz,CPU_Model,Uptime_Days,Load_Average,Longest_Process_Hours" > "$SUMMARY_FILE"
echo "Starting HPC Cluster analysis script..." | tee "$LOG_FILE"
echo "Detailed logs will be written to $LOG_FILE"

# --- Functions ---

# Function to check if a node should be skipped
should_skip() {
    local node_num=$1
    for skip in "${SKIP_NODES[@]}"; do
        if [ "$skip" = "$node_num" ]; then
            return 0 # Skip
        fi
    done
    return 1 # Don't skip
}

# Function to convert etime (e.g., "2-12:34:56" or "12:34:56") to total hours (float)
# Using awk for robustness
convert_etime_to_hours() {
    local etime_str="$1"
    echo "$etime_str" | awk -F'[-:]' '
        {
            days = 0; hrs = 0; mins = 0; secs = 0;
            if (NF == 4) { # DD-HH:MM:SS
                days = $1; hrs = $2; mins = $3; secs = $4;
            } else if (NF == 3) { # HH:MM:SS
                hrs = $1; mins = $2; secs = $3;
            } else if (NF == 2) { # MM:SS (less common for long processes)
                mins = $1; secs = $2;
            }
            total_seconds = (days * 24 * 3600) + (hrs * 3600) + (mins * 60) + secs;
            printf "%.2f\n", total_seconds / 3600;
        }'
}

# --- Main Loop ---
for i in $(seq 1 28); do
    node="compute-0-$i"
    node_log_prefix="[$(date '+%Y-%m-%d %H:%M:%S')] $node:"

    # Skip specified nodes
    if should_skip "$i"; then
        echo "$node_log_prefix Skipping (excluded)." | tee -a "$LOG_FILE"
        continue
    fi

    echo "$node_log_prefix Processing..." | tee -a "$LOG_FILE"

    # Attempt SSH with timeout (5 seconds) and save raw output
    # Using a here-document for cleaner multi-line remote command
    ssh -o ConnectTimeout=5 -o BatchMode=yes "$node" bash << 'EOF_REMOTE' > "${RAW_LSCPU_DIR}/${node}_raw.txt" 2>/dev/null
lscpu
uptime
ps -eo etime,comm --sort=-etime | head -n 2 # Get header and the longest process
EOF_REMOTE

    # Check if SSH was successful and raw file was created
    if [ $? -ne 0 ] || [ ! -s "${RAW_LSCPU_DIR}/${node}_raw.txt" ]; then
        echo "$node_log_prefix Failed to connect or retrieve data. Skipping." | tee -a "$LOG_FILE"
        # Remove empty or incomplete raw file
        rm -f "${RAW_LSCPU_DIR}/${node}_raw.txt"
        continue
    fi

    # Parse lscpu
    raw_output=$(cat "${RAW_LSCPU_DIR}/${node}_raw.txt")
    cores=$(echo "$raw_output" | grep "CPU(s):" | head -1 | awk '{print $2}' || echo "N/A")
    # Parse clock speed FIRST as it appears before Model name in the new order
    clock_speed=$(echo "$raw_output" | grep "CPU MHz:" | head -1 | awk '{print $3}' || echo "N/A")
    cpu_model=$(echo "$raw_output" | grep "Model name:" | awk -F: '{print $2}' | xargs || echo "N/A")


    # Parse uptime and load average robustly
    uptime_output=$(echo "$raw_output" | grep " up ") # Ensure it's the uptime line
    
    # Extract uptime days (handles cases with or without 'days' explicitly mentioned)
    uptime_days=$(echo "$uptime_output" | grep -oP '\s*up\s+\K[0-9]+\s+days' | awk '{print $1}' || echo "0")
    
    # Extract load average (1-minute average)
    load_avg=$(echo "$uptime_output" | awk -F'load average: ' '{print $2}' | awk -F', ' '{print $1}' || echo "0.00")

    # Parse longest-running process
    ps_output_line=$(echo "$raw_output" | grep -A1 "ELAPSED" | tail -n 1) # Get the process line after header
    etime=$(echo "$ps_output_line" | awk '{print $1}')
    hours=$(convert_etime_to_hours "$etime")
    hours=${hours:-0.00} # Default to 0.00 if conversion fails

    # Append to summary file
    # CORRECTED: Swapped clock_speed and cpu_model variables to match new header order
    echo "$node,$cores,$clock_speed,$cpu_model,$uptime_days,$load_avg,$hours" >> "$SUMMARY_FILE"
    echo "$node_log_prefix Data extracted and added to summary." | tee -a "$LOG_FILE"
done

echo "Node processing complete. Generating report..." | tee -a "$LOG_FILE"

# --- Generate Report ---
{
    echo "HPC Cluster Analysis Report"
    echo "=========================="
    echo
    echo "Report generated on: $(date)"
    echo "-----------------------------------"
    echo

    echo "1. Raw lscpu Outputs"
    echo "-------------------"
    echo "Raw lscpu, uptime, and process data for each processed node are stored in the '$RAW_LSCPU_DIR/' directory."
    echo "You can inspect them individually, e.g., 'cat $RAW_LSCPU_DIR/compute-0-1_raw.txt'."
    echo

    echo "2. Node Summary (node_summary.csv)"
    echo "--------------------------------"
    if [ -s "$SUMMARY_FILE" ]; then
        column -s, -t "$SUMMARY_FILE"
    else
        echo "No summary data available. '$SUMMARY_FILE' is empty or missing."
    fi
    echo

    echo "2.1 Cluster Overview Statistics"
    echo "-------------------------------"
    if [ -s "$SUMMARY_FILE" ] && [ "$(wc -l < "$SUMMARY_FILE")" -gt 1 ]; then
        total_nodes=$(awk 'NR>1 {print $1}' "$SUMMARY_FILE" | wc -l)
        total_cores=$(awk -F, 'NR>1 {sum+=$2} END {print sum}' "$SUMMARY_FILE")
        avg_uptime=$(awk -F, 'NR>1 {sum+=$5} END {if (NR>1) printf "%.1f", sum/(NR-1); else print 0}' "$SUMMARY_FILE")
        
        # Calculate average load, handling "N/A" or empty strings for robustness
        avg_load=$(awk -F, '
            NR>1 {
                if ($6 ~ /^[0-9.]+$/) { # Check if load_avg is a number (now column 6)
                    sum_load += $6; count_load++;
                }
            } 
            END {
                if (count_load > 0) printf "%.2f", sum_load/count_load; 
                else print 0
            }' "$SUMMARY_FILE")
        
        echo "Total Nodes Processed: $total_nodes"
        echo "Total Cores Across Processed Nodes: $total_cores"
        echo "Average Uptime (Days): $avg_uptime"
        echo "Average Load Average (1-min): $avg_load"
    else
        echo "Insufficient data for cluster statistics."
    fi
    echo

    echo "3. Availability Ranking (Highest Uptime)"
    echo "--------------------------------------"
    echo "Nodes ranked by their uptime in days (highest first)."
    if [ -s "$SUMMARY_FILE" ] && [ "$(wc -l < "$SUMMARY_FILE")" -gt 1 ]; then
        awk -F, 'NR>1 {
            uptime_val = ($5 ~ /^[0-9.]+$/) ? $5 : 0; # Uptime_Days is still column 5
            printf "%s,%.0f\n", $1, uptime_val; # Print node and uptime (rounded)
        }' "$SUMMARY_FILE" | sort -t, -k2 -nr | \
        awk -F, '{printf "%d. %s (Uptime: %.0f days)\n", NR, $1, $2}'
    else
        echo "Insufficient data for this ranking."
    fi
    echo

    echo "4. Performance Ranking (Highest Clock Speed)"
    echo "------------------------------------------"
    echo "Nodes ranked by their CPU clock speed in MHz (highest first)."
    if [ -s "$SUMMARY_FILE" ] && [ "$(wc -l < "$SUMMARY_FILE")" -gt 1 ]; then
        awk -F, 'NR>1 {
            clock_speed = ($3 ~ /^[0-9.]+$/) ? $3 : 0; # Clock_Speed_MHz is now column 3
            printf "%s,%.0f\n", $1, clock_speed; # Print node and clock speed (rounded)
        }' "$SUMMARY_FILE" | sort -t, -k2 -nr | \
        awk -F, '{printf "%d. %s (Clock Speed: %.0f MHz)\n", NR, $1, $2}'
    else
        echo "Insufficient data for this ranking."
    fi
    echo

    echo "5. Long-Running Simulations Ranking"
    echo "---------------------------------"
    echo "Rank based on longest process runtime (hours) - higher is more utilized/long-running."
    if [ -s "$SUMMARY_FILE" ] && [ "$(wc -l < "$SUMMARY_FILE")" -gt 1 ]; then
        awk -F, 'NR>1 {
            hours = ($7 ~ /^[0-9.]+$/) ? $7 : 0; # Longest_Process_Hours is still column 7
            printf "%s,%.2f\n", $1, hours;
        }' "$SUMMARY_FILE" | sort -t, -k2 -nr | \
        awk -F, '{printf "%d. %s (Longest Process: %.2f hours)\n", NR, $1, $2}'
    else
        echo "Insufficient data for this ranking."
    fi
    echo

} > "$REPORT_FILE"

echo "Done! Results saved in $SUMMARY_FILE and $REPORT_FILE" | tee -a "$LOG_FILE"
echo "Raw lscpu outputs are in $RAW_LSCPU_DIR/" | tee -a "$LOG_FILE"