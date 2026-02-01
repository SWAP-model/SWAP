#!/bin/bash

# Run SWAP on Linux
if [ "$1" = "original" ]; then
    echo "Running with the original SWAP executable..."
    ../../../swap_org/linux/swap420
else
    echo "Running with the development version of the SWAP executable..."
    ../../../builddir/swap
fi
# For automation: simulate pressing Enter so the script doesn't block
echo "Press Enter to continue... (automated)"
printf '\n'

rm swap.swp swap_swap.log swap.ok result.* reruns.log *.tmp

echo "SWAP run completed and temporary files removed."