#!/bin/bash

# Run SWAP on Linux
../../../builddir/swap

# For automation: simulate pressing Enter so the script doesn't block
echo "Press Enter to continue... (automated)"
printf '\n'

rm swap.swp swap_swap.log swap.ok result.* reruns.log

echo "SWAP run completed and temporary files removed."