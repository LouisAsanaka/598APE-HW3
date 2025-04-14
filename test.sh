#!/bin/bash
set -e

BODIES=${1}
TIMESTEPS=${2}

# CORRECT_ICON="\U2705"
# WRONG_ICON="\U274C"

# Run sim and save output
./main.exe "$BODIES" "$TIMESTEPS" > "output/$BODIES-$TIMESTEPS.txt"

# Compare with correct output
# if python diff.py "output/$BODIES-$TIMESTEPS.txt" "output/${BODIES}-${TIMESTEPS}_correct.txt"; then
#     echo -e "$CORRECT_ICON N=$BODIES T=$TIMESTEPS passed"
# else
#     echo -e "$WRONG_ICON N=$BODIES T=$TIMESTEPS failed"
# fi