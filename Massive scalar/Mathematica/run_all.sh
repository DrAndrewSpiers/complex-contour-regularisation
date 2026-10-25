#!/bin/zsh
# Runs every stage of the massive-scalar test and writes logs to ../results/logs/.
# Usage: zsh run_all.sh [smoke]      (smoke: reduced l, n and l'max for a quick check, about 10 minutes)
W=/Applications/Wolfram.app/Contents/MacOS/wolframscript
cd "$(dirname "$0")"
mkdir -p ../results/logs ../figures
flag=$1
for s in massive_stage1_first_order massive_stage2_reference massive_stage3_contour massive_stage4_convergence; do
  echo "=== $s $flag ($(date)) ==="
  "$W" -file $s.wls $flag 2>&1 | tee ../results/logs/$s$flag.log | grep "SUMMARY\|wall time\|FAIL"
done
echo "=== note ==="
(cd ../note && latexmk -pdf -interaction=nonstopmode massive-scalar-test.tex > /dev/null 2>&1 && echo "note built") || echo "note build failed"
