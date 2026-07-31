#!/bin/bash
# Master launcher — runs ALL co-evolution scripts with all 1,299 sequences
# Uses numba GPU acceleration + 24 CPU cores

cd /store/shuvam/E-motioner-X-SBS/datasets/co-evolution
export PYTHONUNBUFFERED=1

LOGDIR=/tmp/coevo_logs
mkdir -p $LOGDIR

echo "=== Starting ALL co-evolution scripts with 1,299 sequences ==="
echo "Logs: $LOGDIR"

# 1. run_kmap_analysis.py (already edited for 1299)
python3 -u run_kmap_analysis.py > $LOGDIR/run_kmap.log 2>&1 &
PID1=$!
echo "run_kmap_analysis.py: PID $PID1"

# 2. boolean_co-evolution.py (edited for 1299 + multiprocessing)
python3 -u boolean_co-evolution.py > $LOGDIR/boolean.log 2>&1 &
PID2=$!
echo "boolean_co-evolution.py: PID $PID2"

# 3. nary_kmap_co-evolution.py (edited for 1299)
python3 -u nary_kmap_co-evolution.py > $LOGDIR/nary.log 2>&1 &
PID3=$!
echo "nary_kmap_co-evolution.py: PID $PID3"

# 4. position_kmap_coevolution.py (edited for 1299)
python3 -u position_kmap_coevolution.py > $LOGDIR/position.log 2>&1 &
PID4=$!
echo "position_kmap_coevolution.py: PID $PID4"

# 5. full_length_analysis.py (edited for 1299)
python3 -u full_length_analysis.py > $LOGDIR/full_length.log 2>&1 &
PID5=$!
echo "full_length_analysis.py: PID $PID5"

# 6. master_boolean.py (already uses all 1299)
python3 -u master_boolean.py > $LOGDIR/master_bool.log 2>&1 &
PID6=$!
echo "master_boolean.py: PID $PID6"

# 7. allseq_constraint_function.py
python3 -u allseq_constraint_function.py > $LOGDIR/allseq_cons.log 2>&1 &
PID7=$!
echo "allseq_constraint_function.py: PID $PID7"

# 8. predictive_constraint_function.py
python3 -u predictive_constraint_function.py > $LOGDIR/pred_cons.log 2>&1 &
PID8=$!
echo "predictive_constraint_function.py: PID $PID8"

# 9. perplexity_coevolution.py
python3 -u perplexity_coevolution.py > $LOGDIR/perplexity.log 2>&1 &
PID9=$!
echo "perplexity_coevolution.py: PID $PID9"

# 10. advanced_co-evolution_analysis.py
python3 -u advanced_co-evolution_analysis.py > $LOGDIR/advanced.log 2>&1 &
PID10=$!
echo "advanced_co-evolution_analysis.py: PID $PID10"

# 11. variable_position_coevolution.py
python3 -u variable_position_coevolution.py > $LOGDIR/var_pos.log 2>&1 &
PID11=$!
echo "variable_position_coevolution.py: PID $PID11"

# 12. dca_boolean_coevolution.py
python3 -u dca_boolean_coevolution.py > $LOGDIR/dca_bool.log 2>&1 &
PID12=$!
echo "dca_boolean_coevolution.py: PID $PID12"

# 13. flipped_boolean_coevolution.py
python3 -u flipped_boolean_coevolution.py > $LOGDIR/flipped_bool.log 2>&1 &
PID13=$!
echo "flipped_boolean_coevolution.py: PID $PID13"

# 14. create_mi_heatmap.py
python3 -u create_mi_heatmap.py > $LOGDIR/mi_heat.log 2>&1 &
PID14=$!
echo "create_mi_heatmap.py: PID $PID14"

# 15. kmap_boolean_coevolution.py
python3 -u kmap_boolean_coevolution.py > $LOGDIR/kmap_bool_coevo.log 2>&1 &
PID15=$!
echo "kmap_boolean_coevolution.py: PID $PID15"

# 16. run_allseq_analysis.py
python3 -u run_allseq_analysis.py > $LOGDIR/run_allseq.log 2>&1 &
PID16=$!
echo "run_allseq_analysis.py: PID $PID16"

# 17. GPU-accelerated full analysis
python3 -u gpu_full_analysis.py > $LOGDIR/gpu_full.log 2>&1 &
PID17=$!
echo "gpu_full_analysis.py: PID $PID17"

echo ""
echo "All 17 scripts launched! Waiting..."
wait

echo ""
echo "=== ALL SCRIPTS COMPLETED ==="
for f in $LOGDIR/*.log; do
    echo "$(basename $f): $(wc -l < $f) lines, last: $(tail -1 $f 2>/dev/null | cut -c1-80)"
done
