#!/bin/bash

############################################
# Usage: ./run_optimization.sh [1|2|3|4] [recipient_email(Mac Only)]
#   1: Run mu=1..24 sequentially (single group)
#   2: Run two groups (odd, even) in parallel
#   3: Run three groups (mu % 3 == 0,1,2) in parallel
#   4: Run four groups (mu % 4 == 0,1,2,3) in parallel
#
# Example:
#   ./run_optimization.sh 2 example@gmail.com
############################################

# 1) Check arguments
if [ $# -lt 1 ]; then
  echo "Usage: $0 [1|2|3|4] [recipient_email]"
  exit 1
fi

MODE="$1"
RECIPIENT="$2"  # Optional second argument

# Validate MODE
if [[ ! "$MODE" =~ ^[1-4]$ ]]; then
  echo "Error: first argument must be one of 1,2,3,4"
  exit 1
fi

# 2) Detect OS (for mail sending)
OS_TYPE=$(uname -s)

# 3) Set environment (script path, permissions, etc.)
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$SCRIPT_DIR"

# Give execution permission to the Julia script
chmod +x main_moo.jl

############################################
# Function: send_mail()
#   $1: Mail subject
#   $2: Mail body
#   $3: Attachment (optional)
############################################
send_mail() {
  local subject="$1"
  local body="$2"
  local attachment="$3"

  # If not macOS, skip sending mail
  if [[ "$OS_TYPE" != "Darwin" ]]; then
    echo "[INFO] Not macOS. Skip sending mail."
    return 0
  fi

  # If recipient email not provided, skip
  if [ -z "$RECIPIENT" ]; then
    echo "[INFO] No recipient email specified. Skip sending mail."
    return 0
  fi

  # On macOS, use osascript (Apple Mail)
  /usr/bin/osascript <<EOD
tell application "Mail"
    set newMessage to make new outgoing message with properties {subject:"$subject", content:"$body", visible:true}
    tell newMessage
        make new to recipient at end of to recipients with properties {address:"$RECIPIENT"}
        if ("$attachment" is not "") then
            make new attachment with properties {file name:POSIX file "$attachment"} at after the last paragraph
        end if
        send
    end tell
end tell
EOD
}

############################################
# Function: run_group()
#   $1: Group name (string)
#   $2: List of mu values (e.g. "1 5 9 13 17 21")
############################################
run_group() {
  local group_name="$1"
  local mu_list="$2"

  # Each group has its own log file
  local log_file="output_${group_name}.txt"
  rm -f "$log_file"

  echo "===== Start group [$group_name] =====" | tee -a "$log_file"
  for mu in $mu_list
  do
    echo "Running optimization for mu = $mu (group=$group_name)..." | tee -a "$log_file"
    # Changed to >> log_file to avoid pipe/tee issues
    julia -t10 main_moo.jl --mu "$mu" --no-circ >> "$log_file" 2>&1
    echo "Completed optimization for mu = $mu (group=$group_name)" | tee -a "$log_file"
  done
  echo "===== All [${group_name}] optimizations completed. =====" | tee -a "$log_file"

  # When this group finishes, send an email alert (if conditions allow)
  send_mail \
    "Group [${group_name}] Optimization Completed" \
    "The optimization for group '${group_name}' has been completed." \
    "$SCRIPT_DIR/$log_file"
}

############################################
# Decide how to group mu=1..24 based on MODE
############################################

if [ "$MODE" -eq 1 ]; then
  # MODE=1: Run mu=1..24 sequentially (single group)
  run_group "all_1_to_24" "$(seq 1 24)"

elif [ "$MODE" -eq 2 ]; then
  # MODE=2: Two groups in parallel (odd, even)
  odd_mus=""
  even_mus=""
  for i in $(seq 1 24); do
    if [ $((i % 2)) -eq 1 ]; then
      odd_mus="$odd_mus $i"
    else
      even_mus="$even_mus $i"
    fi
  done

  run_group "odd"  "$odd_mus"  &
  run_group "even" "$even_mus" &
  wait

elif [ "$MODE" -eq 3 ]; then
  # MODE=3: Three groups in parallel (mu%3=0, mu%3=1, mu%3=2)
  group0=""  # mu % 3 = 0 => 3,6,9,12,15,18,21,24
  group1=""  # mu % 3 = 1 => 1,4,7,10,13,16,19,22
  group2=""  # mu % 3 = 2 => 2,5,8,11,14,17,20,23

  for i in $(seq 1 24); do
    r=$(( i % 3 ))
    if [ "$r" -eq 0 ]; then
      group0="$group0 $i"
    elif [ "$r" -eq 1 ]; then
      group1="$group1 $i"
    else
      group2="$group2 $i"
    fi
  done

  run_group "mod3_0" "$group0" &
  run_group "mod3_1" "$group1" &
  run_group "mod3_2" "$group2" &
  wait

else
  # MODE=4: Four groups in parallel (mu%4=0,1,2,3)
  group0=""  # mu % 4 = 0 => 4,8,12,16,20,24
  group1=""  # mu % 4 = 1 => 1,5,9,13,17,21
  group2=""  # mu % 4 = 2 => 2,6,10,14,18,22
  group3=""  # mu % 4 = 3 => 3,7,11,15,19,23

  for i in $(seq 1 24); do
    r=$(( i % 4 ))
    case "$r" in
      0) group0="$group0 $i";;
      1) group1="$group1 $i";;
      2) group2="$group2 $i";;
      3) group3="$group3 $i";;
    esac
  done

  run_group "mod4_0" "$group0" &
  run_group "mod4_1" "$group1" &
  run_group "mod4_2" "$group2" &
  run_group "mod4_3" "$group3" &
  wait
fi

echo "=== All tasks in MODE=$MODE have completed. ==="
exit 0
