#!/bin/bash
echo "=== CLEANING START! ==="
cd /scratch3/dflusova/afterburner/slurm
rm *
pwd
ls
cd /scratch3/dflusova/afterburner/out
rm *
pwd
ls
echo "=== CLEANING DONE! ==="
