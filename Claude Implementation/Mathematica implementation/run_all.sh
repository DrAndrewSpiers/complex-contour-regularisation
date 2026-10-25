#!/bin/zsh
# Runs every stage and writes logs to results/logs/.  Usage:  zsh run_all.sh
W=/Applications/Wolfram.app/Contents/MacOS/wolframscript
cd "${0:A:h}"
mkdir -p results/logs
for s in stage1_singular_mode stage2_regular_field stage3a_source_modes; do
  echo "== $s"; $W -file $s.wls > results/logs/$s.log 2>&1; grep -E "SUMMARY" results/logs/$s.log
done
echo "== stage3b_second_order_solve (lmax = 6)"; $W -file stage3b_second_order_solve.wls 6 > results/logs/stage3b_lmax6.log 2>&1; grep -E "SUMMARY" results/logs/stage3b_lmax6.log
echo "== stage3b_second_order_solve (lmax = 12)"; $W -file stage3b_second_order_solve.wls 12 > results/logs/stage3b_lmax12.log 2>&1; grep -E "SUMMARY" results/logs/stage3b_lmax12.log
echo "== stage4_deferred_shadow_test"; $W -file stage4_deferred_shadow_test.wls > results/logs/stage4_deferred_shadow_test.log 2>&1; grep -E "SUMMARY" results/logs/stage4_deferred_shadow_test.log
echo "== stage5a_gradient_sources"; $W -file stage5a_gradient_sources.wls > results/logs/stage5a_gradient_sources.log 2>&1; grep -E "SUMMARY" results/logs/stage5a_gradient_sources.log
echo "== stage5b_gradient_solve (lmax = 6)"; $W -file stage5b_gradient_solve.wls 6 > results/logs/stage5b_lmax6.log 2>&1; grep -E "SUMMARY" results/logs/stage5b_lmax6.log
echo "== stage5b_gradient_solve (lmax = 12)"; $W -file stage5b_gradient_solve.wls 12 > results/logs/stage5b_lmax12.log 2>&1; grep -E "SUMMARY" results/logs/stage5b_lmax12.log
echo "== stage6_3d_puncture_test"; $W -file stage6_3d_puncture_test.wls > results/logs/stage6_3d_puncture_test.log 2>&1; grep -E "SUMMARY" results/logs/stage6_3d_puncture_test.log
echo "== stage7_lprime_optimum"; $W -file stage7_lprime_optimum.wls > results/logs/stage7_lprime_optimum.log 2>&1; grep -E "SUMMARY" results/logs/stage7_lprime_optimum.log
echo "== stage8_hybrid_split"; $W -file stage8_hybrid_split.wls > results/logs/stage8_hybrid_split.log 2>&1; grep -E "SUMMARY" results/logs/stage8_hybrid_split.log
echo "== stage9_dzz"; $W -file stage9_dzz.wls > results/logs/stage9_dzz.log 2>&1; grep -E "SUMMARY" results/logs/stage9_dzz.log
echo "== stage10_teukolsky_benchmark"; $W -file stage10_teukolsky_benchmark.wls > results/logs/stage10_teukolsky_benchmark.log 2>&1; grep -E "SUMMARY" results/logs/stage10_teukolsky_benchmark.log
echo "== stage11_TQ"; $W -file stage11_TQ.wls > results/logs/stage11_TQ.log 2>&1; grep -E "SUMMARY" results/logs/stage11_TQ.log
echo "== stage12_ss_reuse_convergence"; $W -file stage12_ss_reuse_convergence.wls > results/logs/stage12_ss_reuse_convergence.log 2>&1; grep -E "SUMMARY" results/logs/stage12_ss_reuse_convergence.log
echo "== stage12b_physical_bc"; $W -file stage12b_physical_bc.wls > results/logs/stage12b_physical_bc.log 2>&1; grep -E "SUMMARY" results/logs/stage12b_physical_bc.log
echo "== stage7b_physical_bc"; $W -file stage7b_physical_bc.wls > results/logs/stage7b_physical_bc.log 2>&1; grep -E "SUMMARY" results/logs/stage7b_physical_bc.log
echo "== stage12 (phys)"; $W -file stage12_ss_reuse_convergence.wls phys > results/logs/stage12phys.log 2>&1; grep -E "SUMMARY" results/logs/stage12phys.log
echo "== stage11 (phys)"; $W -file stage11_TQ.wls phys > results/logs/stage11phys.log 2>&1; grep -E "SUMMARY" results/logs/stage11phys.log
echo "== make_prl_figures"; $W -file make_prl_figures.wls > results/logs/make_prl_figures.log 2>&1; grep -E "SUMMARY" results/logs/make_prl_figures.log
echo "== stage14_cost"; $W -file stage14_cost.wls > results/logs/stage14_cost.log 2>&1; grep -E "SUMMARY" results/logs/stage14_cost.log
echo "== stage14b_radial"; $W -file stage14b_radial.wls > results/logs/stage14b_radial.log 2>&1; grep -E "SUMMARY" results/logs/stage14b_radial.log
echo "== stage14b (phys)"; $W -file stage14b_radial.wls phys > results/logs/stage14bphys.log 2>&1; grep -E "SUMMARY" results/logs/stage14bphys.log
echo "== make_notebook"; $W -file make_notebook.wls
