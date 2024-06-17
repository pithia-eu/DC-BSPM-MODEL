#!/bin/bash

# Set the environment variables
#export IRI0_IRILIB64PATH=irilib64.so
#export LD_PRELOAD=iri0.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/Libs/iri2016
# Define the app name
app_name="BSPM_May2024.py"
script_dir="$(pwd)"
# Define the parent directory of the script directory
parent_dir=$(dirname "$script_dir")

# Default values for year, month, and day
year=""
month=""
day=""
executionid=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --year)
      year="$2"
      shift 2
      ;;
    --month)
      month="$2"
      shift 2
      ;;
    --day)
      day="$2"
      shift 2
      ;;
    --executionid)
      executionid="$2"
      shift 2
      ;;
    --rerun)
      rerun="$2"
      shift 2
      ;;
    *)
      echo "{\"code\": -1, \"msg\": \"Unknown option: $1\"}"
      exit 1
      ;;
  esac
done

# Check if year, month, and day are provided
if [ -z "$year" ] || [ -z "$month" ] || [ -z "$day" ]; then
  echo "{\"code\": -1, \"msg\": \"Year, month, and day are required.\"}"
  exit 1
fi

# If executionid is not provided, set it to year-month-day
if [ -z "$executionid" ]; then
  executionid="$(printf "%04d" "$year")-$(printf "%02d" "$month")-$(printf "%02d" "$day")"
fi

# Define the output folder path
output_folder="$parent_dir/out/bspm/$USER/$executionid"
python_command="python3.9 $script_dir/../bspm/$app_name --year $year --month $month --day $day --executionid $executionid"
#echo "$output_folder"

# Check if rerun is provided, default is false, if true, remove the output folder
if [ "$rerun" = "true" ]; then
  if [ -d "$output_folder" ]; then
    rm -rf "$output_folder"
  fi
fi

# Check if the output folder exists; create it if it doesn't
if [ ! -d "$output_folder" ]; then
  mkdir -p "$output_folder"
  # echo "Created output folder: $output_folder"
  #IRI0_IRILIB64PATH=irilib64.so LD_PRELOAD=iri0.so python3.9 "$app_name" --year "$year" --month "$month" --day "$day" --executionid "$executionid" &
  nohup $python_command > $output_folder/app.log 2>&1 &
  echo "{\"code\": 0, \"msg\": \"Started a new execution of $app_name by date $executionid\",\"date\":\"$executionid\",\"status\":\"start\"}"
else
  if pgrep -f "$app_name.*--executionid $executionid" > /dev/null; then
    echo "{\"code\": 0, \"msg\": \"App $app_name on date $executionid is already running.\",\"date\":\"$executionid\",\"status\":\"progressing\"}"
  else
    echo "{\"code\": 1, \"msg\": \"App $app_name on date $executionid is completed.\",\"date\":\"$executionid\",\"status\":\"completed\"}"
  fi
fi
