#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 4: deferred shadow test.  The first-order field inside the tube comes from a TRUNCATED puncture
   (Taylor series of the corrected phiS_l about r0 through order n, canonical logs: the analogue of the notebook's
   phiSLocal).  The residual field phi~res = phi~R + kappa (phiS - phiP_n) carries a shadow at the truncation order and is
   used directly in the SR source inside the tube.  How large is the exterior error, and how does it scale with n?
   Run:  wolframscript -file stage4_deferred_shadow_test.wls                                                    *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
Needs["NumericalDifferentialEquationAnalysis`"];   (* GaussianQuadratureWeights *)
resDir = FileNameJoin[{dir, "results"}]; Quiet[CreateDirectory[resDir]];
base = <|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>;
SetCSParameters[base];
Q = "Product"; lmax = 6; nList = {2, 3, 4, 6, 8, 10, 12}; methods = {"ThetaToH", "CorrectedSingular"};
mLabel = <|"ThetaToH" -> "Theta->H", "CorrectedSingular" -> "corrected-phiS"|>;
r0 = csPar["r0"]; q = csPar["q"]; rMin = csPar["rMin"]; rMax = csPar["rMax"];
csResetChecks[];
t0 = AbsoluteTime[];
tick[] := Round[AbsoluteTime[] - t0, 0.1];

Print["=== Stage 4: deferred shadow test (truncated puncture, residual field used in SR) ==="];
Print["    n = ", nList, ", l = 0..", lmax, ", K continuations: ", Values[mLabel]];

(* ---------------- helpers ---------------- *)
(* standard exterior grid (radius 1/2) and the SR scale used for all scaled errors *)
extStd = Join[Subdivide[rMin, 9/2, 35], Subdivide[11/2, rMax, 45]];
extR = Subdivide[11/2, rMax, 45];          (* region 3 of the standard tube *)
solSR = Association@Table[l -> pieceAnalyticSolution[Q, "SR", l], {l, 0, 12}];
scaleSR = Association@Table[l -> Max[Abs[analyticEval[solSR[l], #] & /@ extStd]], {l, 0, lmax}];

dphiExpr[lp_, n_] := dphiExpr[lp, n] = kappaL[lp] (singularModeCorrected[lp][zVar] - singularModeTaylor[lp, n][zVar]);
ktExpr[l_, lp_, method_] := ktExpr[l, lp, method] = KCouplingContinued[l, lp, zVar, method];
evalAt[e_, s_] := N[e /. zVar -> SetPrecision[s, 90], 80];      (* real s: canonical logs give the value from above *)

(* exterior error of a homogeneous solution with data (u, u') at aR, maximised over a grid, scaled *)
homErr[l_, {u_, du_}, aRv_, grid_] := Module[{AB},
   AB = LinearSolve[{{aRv^l, aRv^(-l - 1)}, {l aRv^(l - 1), -(l + 1) aRv^(-l - 2)}}, {u, du}];
   Max[Abs[AB[[1]] #^l + AB[[2]] #^(-l - 1)] & /@ grid]/scaleSR[l]];

(* Lemma-1 decomposition by Gauss-Legendre quadrature on [aL, r0] and [r0, aR] (split at the log singularity):
   response at aR, with zero data at aL, to (a) 2q Sum K Re(dphi) and (b) -2q Sum Im K~ Im(dphi), dphi = kappa(phiS - phiP_n). *)
nodesFor[np_] := Join[GaussianQuadratureWeights[np, tubeL[], r0, 60], GaussianQuadratureWeights[np, r0, tubeR[], 60]];
contrib[l_, lp_, n_, method_, nodes_] := Module[{s = nodes[[All, 1]], w = nodes[[All, 2]], aRv = tubeR[], g0, g1, dp, kt, fa, fb},
   dp = evalAt[dphiExpr[lp, n], #] & /@ s;
   kt = evalAt[ktExpr[l, lp, method], #] & /@ s;
   g0 = (s^2/(2 l + 1)) (aRv^l s^(-l - 1) - s^l aRv^(-l - 1));
   g1 = (s^2/(2 l + 1)) (l aRv^(l - 1) s^(-l - 1) + (l + 1) s^l aRv^(-l - 2));
   fa = 2 q Re[kt] Re[dp]; fb = -2 q Im[kt] Im[dp];
   <|"a" -> {Total[w g0 fa], Total[w g1 fa]}, "b" -> {Total[w g0 fb], Total[w g1 fb]}|>];

(* first-order residual fields (notebook section-1 matching, correct branch) and SR three-region solve *)
residuals[n_, lps_] := Association@Table[lp -> firstOrderRegular[lp, Automatic, n], {lp, lps}];
srSolve[l_, res1_Association, method_] := Module[{ext, rr},
   ext = Association@KeyValueMap[Function[{lp, rs}, lp -> Module[{f}, f[r_?NumericQ] := exteriorValue[rs, r]; f]], KeySelect[res1, # <= 1 &]];
   rr = threeRegionSolve[l, srSourceReal[Q, l, ext], srSourceContour[Q, l, #["sol2"] & /@ res1, method],
     {analyticEval[solSR[l], rMin], analyticEvalD[solSR[l], rMin]}];
   <|"res" -> rr, "delta" -> {Re[rr["valR"]] - analyticEval[solSR[l], tubeR[]], Re[rr["derR"]] - analyticEvalD[solSR[l], tubeR[]]},
     "ext" -> With[{g = Select[extStd, # <= tubeL[] || # >= tubeR[] &]},
       Max[Abs[(exteriorValue[rr, #] - analyticEval[solSR[l], #]) & /@ g]]/scaleSR[l]]|>];
step1Check[res1_, n_, label_] := KeyValueMap[Function[{lp, rs},
    Module[{aG = Subdivide[tubeL[], tubeR[], 40], zs, ref},
     zs = N[yC[SetPrecision[#, 60]], 60] & /@ aG;
     ref = (N[phiRAnalytic[lp][#] + (dphiExpr[lp, n] /. zVar -> #), 40] &) /@ zs;
     csCheck[label <> " step 1: residual l'=" <> ToString[lp] <> ", n=" <> ToString[n] <> " = phiR + kappa(phiS - phiP_n) on the contour (scaled)",
      scaledError[rs["sol2"] /@ aG, ref], 10^-12]]], res1];
consTol[l_, d_] := 10^-6 Norm[d] + 10^-16 scaleSR[l];

(* ---------------- baseline: exact phiS (n = infinity) ---------------- *)
resExact = Association@Table[lp -> firstOrderRegular[lp], {lp, 0, 1}];
noise = Association@Flatten@Table[{l, m} -> srSolve[l, resExact, m]["ext"], {l, 0, lmax}, {m, methods}];
Print["    baseline (exact phiS) scaled exterior SR error, max over l and K: ", sci[Max[Values[noise]]], "   [", tick[], " s]"];

(* ================= steps 1-3: main scan ================= *)
Print["\n--- Steps 1-3: scan over n (l' in {0,1}) ---"];
rows = {};          (* {n, l, method, total(step 2), a, b, total(a+b), consistency} *)
quadConv = {};
Do[
  Module[{res1 = residuals[n, {0, 1}], nodes = nodesFor[60], nodes2 = nodesFor[120]},
   step1Check[res1, n, "amp=1"];
   Do[
    Module[{cs, cs2, ua, ub, s2},
     cs = Table[contrib[l, lp, n, m, nodes], {lp, 0, 1}];
     ua = Total[#["a"] & /@ cs]; ub = Total[#["b"] & /@ cs];
     If[MemberQ[{2, 8, 12}, n] && MemberQ[{0, 3, 6}, l],
      cs2 = Table[contrib[l, lp, n, m, nodes2], {lp, 0, 1}];
      AppendTo[quadConv, Max[Abs[Join[ua - Total[#["a"] & /@ cs2], ub - Total[#["b"] & /@ cs2]]]]/Max[Abs[Join[ua, ub]]]]];
     s2 = srSolve[l, res1, m];
     csCheck["n=" <> ToString[n] <> " l=" <> ToString[l] <> " " <> mLabel[m] <> ": step-2 delta(aR) = (a) + (b) from the Lemma-1 decomposition",
      Max[Abs[s2["delta"] - (ua + ub)]], consTol[l, s2["delta"]]];
     AppendTo[rows, {n, l, m, s2["ext"], homErr[l, ua, tubeR[], extR], homErr[l, ub, tubeR[], extR],
       homErr[l, ua + ub, tubeR[], extR], Max[Abs[s2["delta"] - (ua + ub)]]}]],
    {l, 0, lmax}, {m, methods}];
   Print["    n = ", n, " done  [", tick[], " s]"]],
  {n, nList}];
csCheck["quadrature converged (60 vs 120 Gauss nodes per side, relative)", Max[quadConv], 10^-12];

(* ---------------- tables for the main scan ---------------- *)
get[n_, l_, m_, k_] := First@Cases[rows, {n, l, m, x___} :> {n, l, m, x}[[k]]];
maxOverL[n_, m_, k_] := Max[Table[get[n, l, m, k], {l, 0, lmax}]];
Export[FileNameJoin[{resDir, "stage4_scan.csv"}],
  Prepend[{#[[1]], #[[2]], mLabel[#[[3]]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]], sci[If[#[[5]] == 0, 0, #[[6]]/#[[5]]]], sci[#[[8]]]} & /@ rows,
   {"n", "l", "K continuation", "total_ext_err_step2(scaled)", "a_real_remainder(scaled)", "b_shadow(scaled)", "a+b(scaled)", "b/a", "|delta-(a+b)|"}], "CSV"];

fitRate[pts_] := Module[{p = Select[pts, #[[2]] > 10^3 Max[Values[noise]] &], m},
   If[Length[p] < 3, Missing[], m = Fit[{#[[1]], Log10[#[[2]]]} & /@ p, {1, x}, x]; 10^Coefficient[m, x]]];
summaryRows = Table[
   {n, maxOverL[n, "ThetaToH", 5], maxOverL[n, "ThetaToH", 6], maxOverL[n, "ThetaToH", 4], maxOverL[n, "ThetaToH", 6]/maxOverL[n, "ThetaToH", 5],
    maxOverL[n, "CorrectedSingular", 6], maxOverL[n, "CorrectedSingular", 4], maxOverL[n, "CorrectedSingular", 6]/maxOverL[n, "CorrectedSingular", 5],
    N[(csPar["radius"]/r0)^n]}, {n, nList}];
rates = Table[{k, m, fitRate[Table[{n, maxOverL[n, m, k]}, {n, nList}]]}, {m, methods}, {k, {4, 5, 6}}];
ratesPerL = Flatten[Table[{l, mLabel[m], fitRate[Table[{n, get[n, l, m, 4]}, {n, nList}]], fitRate[Table[{n, get[n, l, m, 5]}, {n, nList}]],
     fitRate[Table[{n, get[n, l, m, 6]}, {n, nList}]]}, {m, methods}, {l, 0, lmax}], 1];
Export[FileNameJoin[{resDir, "stage4_fits.csv"}], Prepend[Map[If[NumericQ[#] && ! IntegerQ[#], ToString[NumberForm[#, 3]], #] &, ratesPerL, {2}] /. _Missing -> "-",
   {"l", "K", "decay factor per order: total", "(a)", "(b)"}], "CSV"];

md1 = StringRiffle[Join[{
     "| n | (a) real remainder | (b) shadow, Theta->H | total, Theta->H | (b)/(a), Theta->H | (b) shadow, corrected-phiS | total, corrected-phiS | (b)/(a), corrected-phiS | (radius/r0)^n |",
     "|---|---|---|---|---|---|---|---|---|"},
    ("| " <> StringRiffle[Join[{ToString[#[[1]]]}, sci /@ Rest[#]], " | "] <> " |") & /@ summaryRows], "\n"];
md2 = StringRiffle[Join[{"| quantity | K continuation | fitted decay factor per order (max over l) |", "|---|---|---|"},
    ("| " <> {4 -> "total", 5 -> "(a)", 6 -> "(b)"}[[#[[1]] - 3, 2]] <> " | " <> mLabel[#[[2]]] <> " | " <> If[MissingQ[#[[3]]], "-", ToString[NumberForm[#[[3]], 3]]] <> " |") & /@ Flatten[rates, 1]], "\n"];
Print["\n    Max over l = 0..6 of the scaled exterior SR error (scale = max|SR analytic| over the exterior)\n"];
Print[md1]; Print[""]; Print[md2];

(* ================= step 4: sensitivity ================= *)
Print["\n--- Step 4: amplitude of the contour (Lemma 1: exterior unchanged) ---"];
ampRows = {};
deltaAmp1 = <||>;
Do[
  Module[{res1 = residuals[n, {0, 1}]},
   Print["    amp = 1, n = ", n, ": max |Im(phi~res - phi~R)| on the contour = ",
    sci[Max[Table[Max[Abs[Im[(res1[lp]["sol2"][#] - N[phiRAnalytic[lp][yC[SetPrecision[#, 40]]], 40]) & /@ Subdivide[tubeL[], tubeR[], 40]]]], {lp, 0, 1}]]]];
   Do[deltaAmp1[{n, l, m}] = srSolve[l, res1, m]["delta"], {l, 0, lmax}, {m, methods}]], {n, {2, 4, 8}}];
Do[
  SetCSParameters[Join[base, <|"amp" -> amp|>]];
  Do[
   Module[{res1 = residuals[n, {0, 1}], imMax},
    step1Check[res1, n, "amp=" <> ToString[amp, InputForm]];
    (* shadow of the residual on the contour, i.e. Im(phi~res - phi~R) = Im kappa(phiS - phiP_n) along y(a)
       (Im phi~res itself is dominated by the harmless Im of the continued phiR_1 = c1 z Sqrt[4Pi/3]) *)
    imMax = Max[Table[Max[Abs[Im[(res1[lp]["sol2"][#] - N[phiRAnalytic[lp][yC[SetPrecision[#, 40]]], 40]) & /@ Subdivide[tubeL[], tubeR[], 40]]]], {lp, 0, 1}]];
    Do[
     Module[{d = srSolve[l, res1, m]["delta"], d1 = deltaAmp1[{n, l, m}]},
      csCheck["amp=" <> ToString[amp, InputForm] <> " n=" <> ToString[n] <> " l=" <> ToString[l] <> " " <> mLabel[m] <> ": delta(aR) equals the amp=1 value",
       Max[Abs[d - d1]], consTol[l, d1]];
      AppendTo[ampRows, {amp, n, l, mLabel[m], imMax, Norm[d1], Max[Abs[d - d1]]}]], {l, 0, lmax}, {m, methods}]],
   {n, {2, 4, 8}}];
  Print["    amp = ", amp, " done; max |Im(phi~res - phi~R)| on the contour for n = 2, 4, 8: ",
   sci /@ Table[Max[Cases[ampRows, row : {amp, n, __} :> row[[5]]]], {n, {2, 4, 8}}], "  [", tick[], " s]"], {amp, {1/2, 2}}];
SetCSParameters[base];
Export[FileNameJoin[{resDir, "stage4_amplitude.csv"}],
  Prepend[{ToString[#[[1]], InputForm], #[[2]], #[[3]], #[[4]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]]} & /@ ampRows,
   {"amp", "n", "l", "K", "max|Im(phiRes - phiR)| on contour", "|delta(aR)| at amp=1", "|delta(amp) - delta(amp=1)|"}], "CSV"];

Print["\n--- Step 4b: tube radius 1/4, 1/2, 1 (decomposition; three-region confirmation at n = 4) ---"];
radRows = {};
gridCommon = Subdivide[6, rMax, 40];       (* common exterior region for all radii *)
Do[
  SetCSParameters[Join[base, <|"radius" -> rad|>]];
  Module[{nodes = nodesFor[60]},
   Do[
    Module[{cs = Table[contrib[l, lp, n, m, nodes], {lp, 0, 1}], ua, ub},
     ua = Total[#["a"] & /@ cs]; ub = Total[#["b"] & /@ cs];
     AppendTo[radRows, {rad, n, l, m, homErr[l, ua, tubeR[], gridCommon], homErr[l, ub, tubeR[], gridCommon]}];
     If[n == 4 && rad != 1/2,
      Module[{s2 = srSolve[l, residuals[4, {0, 1}], m]},
       csCheck["radius=" <> ToString[rad, InputForm] <> " n=4 l=" <> ToString[l] <> " " <> mLabel[m] <> ": step-2 delta(aR) = (a) + (b)",
        Max[Abs[s2["delta"] - (ua + ub)]], consTol[l, s2["delta"]]]]]],
    {n, nList}, {l, 0, lmax}, {m, methods}]];
  Print["    radius = ", rad, " done  [", tick[], " s]"], {rad, {1/4, 1/2, 1}}];
SetCSParameters[base];
radMax[rad_, n_, m_, k_] := Max[Cases[radRows, {rad, n, _, m, a_, b_} :> {a, b}[[k]]]];
radTable = Table[{n, Sequence @@ Flatten@Table[{radMax[rad, n, "ThetaToH", 1], radMax[rad, n, "ThetaToH", 2], radMax[rad, n, "CorrectedSingular", 2]}, {rad, {1/4, 1/2, 1}}]}, {n, nList}];
Export[FileNameJoin[{resDir, "stage4_radius.csv"}],
  Prepend[Map[If[IntegerQ[#], #, sci[#]] &, radTable, {2}],
   {"n", "(a) R=1/4", "(b) ThetaH R=1/4", "(b) corrK R=1/4", "(a) R=1/2", "(b) ThetaH R=1/2", "(b) corrK R=1/2", "(a) R=1", "(b) ThetaH R=1", "(b) corrK R=1"}], "CSV"];

(* ================= step 5 (secondary): l' >= 2 residual modes ================= *)
Print["\n--- Step 5 (secondary): include the l' >= 2 residual modes (pure truncation remainder) ---"];
lpMaxList = {1, 2, 4, 8, 12};
lpRows = {};     (* {n, l, method, lp, |a contribution|, |b contribution|} per l' (scaled) *)
cumRows = {};    (* {n, l, method, lpMax, quadrature total, three-region total} *)
Do[
  Module[{res1 = residuals[n, Range[0, 12]], nodes = nodesFor[60]},
   Do[If[lp >= 2, step1Check[KeyTake[res1, {lp}], n, "l'>=2"]], {lp, 0, 12}];
   Do[
    Module[{cs = Association@Table[lp -> contrib[l, lp, n, m, nodes], {lp, 0, 12}]},
     Do[AppendTo[lpRows, {n, l, m, lp, homErr[l, cs[lp]["a"], tubeR[], extR], homErr[l, cs[lp]["b"], tubeR[], extR]}], {lp, 0, 12}];
     Do[
      Module[{uq = Total[Table[cs[lp]["a"] + cs[lp]["b"], {lp, 0, lpm}]], s2},
       s2 = If[lpm >= 2, srSolve[l, KeyTake[res1, Range[0, lpm]], m], Missing[]];
       If[! MissingQ[s2],
        csCheck["l'max=" <> ToString[lpm] <> " n=" <> ToString[n] <> " l=" <> ToString[l] <> " " <> mLabel[m] <> ": step-2 delta(aR) = sum of per-l' (a)+(b)",
         Max[Abs[s2["delta"] - uq]], consTol[l, s2["delta"]]]];
       AppendTo[cumRows, {n, l, m, lpm, homErr[l, uq, tubeR[], extR], If[MissingQ[s2], Missing[], s2["ext"]]}]], {lpm, lpMaxList}]],
    {l, 0, lmax}, {m, methods}];
   Print["    n = ", n, " done  [", tick[], " s]"]], {n, {4, 8}}];
Export[FileNameJoin[{resDir, "stage4_lprime_contributions.csv"}],
  Prepend[{#[[1]], #[[2]], mLabel[#[[3]]], #[[4]], sci[#[[5]]], sci[#[[6]]]} & /@ lpRows, {"n", "l", "K", "l'", "(a) contribution (scaled)", "(b) contribution (scaled)"}], "CSV"];
Export[FileNameJoin[{resDir, "stage4_lprime_cumulative.csv"}],
  Prepend[{#[[1]], #[[2]], mLabel[#[[3]]], #[[4]], sci[#[[5]]], If[MissingQ[#[[6]]], "-", sci[#[[6]]]]} & /@ cumRows, {"n", "l", "K", "l'max", "total via decomposition (scaled)", "total via three-region solve (scaled)"}], "CSV"];
lpTable = Table[{n, lp, Max[Cases[lpRows, {n, _, "ThetaToH", lp, a_, _} :> a]], Max[Cases[lpRows, {n, _, "ThetaToH", lp, _, b_} :> b]],
    Max[Cases[lpRows, {n, _, "CorrectedSingular", lp, _, b_} :> b]]}, {n, {4, 8}}, {lp, 0, 12}];
cumTable = Table[{n, lpm, Max[Cases[cumRows, {n, _, "ThetaToH", lpm, t_, _} :> t]], Max[Cases[cumRows, {n, _, "CorrectedSingular", lpm, t_, _} :> t]]}, {n, {4, 8}}, {lpm, lpMaxList}];
md3 = StringRiffle[Join[{"| n | l' | (a) contribution | (b) contribution, Theta->H | (b) contribution, corrected-phiS |", "|---|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;;]], " | "] <> " |") & /@ Flatten[lpTable, 1]], "\n"];
md4 = StringRiffle[Join[{"| n | l'max | total, Theta->H | total, corrected-phiS |", "|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;;]], " | "] <> " |") & /@ Flatten[cumTable, 1]], "\n"];
Print["\n    Per-l' contribution (max over l = 0..6, scaled)\n"]; Print[md3];
Print["\n    Cumulative total error vs l'max (max over l, scaled)\n"]; Print[md4];

(* ---------------- Markdown summary and plots ---------------- *)
Export[FileNameJoin[{resDir, "stage4_summary.md"}],
  "# Stage 4: deferred shadow test\n\nScaled exterior SR error = max over r in [1,4.5] U [5.5,10] of |error| / max|SR analytic|, then max over l = 0..6. Cauchy data at r = 1. Baseline with the exact phiS: " <>
   sci[Max[Values[noise]]] <> ".\n\n## Error vs puncture order n (l' in {0,1})\n\n" <> md1 <> "\n\n" <> md2 <>
   "\n\n## Secondary: l' >= 2 residual modes\n\n" <> md3 <> "\n\n" <> md4 <> "\n", "Text"];

plotData[m_, k_] := Table[{n, maxOverL[n, m, k]}, {n, nList}];
p1 = ListLogPlot[{plotData["ThetaToH", 5], plotData["ThetaToH", 6], plotData["CorrectedSingular", 6], plotData["ThetaToH", 4],
    Table[{n, plotData["ThetaToH", 5][[1, 2]] (1/10)^(n - 2)}, {n, nList}]},
   Joined -> True, PlotMarkers -> Automatic, PlotRange -> All, ImageSize -> 520,
   PlotLegends -> {"(a) real truncation remainder", "(b) shadow, \[CapitalTheta]\[RightArrow]H", "(b) shadow, corrected \!\(\*SuperscriptBox[\(\[Phi]\), \(S\)]\)", "total, \[CapitalTheta]\[RightArrow]H", "\[Proportional] (radius/\!\(\*SubscriptBox[\(r\), \(0\)]\)\!\(\*SuperscriptBox[\()\), \(n\)]\)"},
   PlotStyle -> {Black, Blue, Red, {Gray, Dashed}, {Orange, Dotted}}, AxesLabel -> {"n", "scaled error"},
   PlotLabel -> "SR exterior error from the truncated-puncture residual field (max over \[ScriptL] = 0..6)",
   Epilog -> {Gray, Dashed, Line[{{2, Log[Max[Values[noise]]]}, {12, Log[Max[Values[noise]]]}}]}];
p2 = ListLogPlot[Table[{lp, Max[Cases[lpRows, {n, _, "ThetaToH", lp, a_, b_} :> a + b]]}, {n, {4, 8}}, {lp, 0, 12}],
   Joined -> True, PlotMarkers -> Automatic, PlotLegends -> {"n = 4", "n = 8"}, AxesLabel -> {"\[ScriptL]'", "scaled error"},
   PlotLabel -> "secondary: contribution of residual mode \[ScriptL]' to the SR exterior error (\[CapitalTheta]\[RightArrow]H, max over \[ScriptL])", ImageSize -> 520];
Export[FileNameJoin[{resDir, "stage4_error_vs_n.pdf"}], Column[{p1, p2}]];
Print["\n   wrote results/stage4_*.csv, stage4_summary.md, stage4_error_vs_n.pdf"];

Print["\nStage 4 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage4_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
