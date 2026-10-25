#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 2: first-order regular field with the corrected exact singular mode (S_eff = 0).
   Run:  wolframscript -file stage2_regular_field.wls                                        *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]];
resDir = FileNameJoin[{dir, "results"}]; Quiet[CreateDirectory[resDir]];
SetCSParameters[<|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>];
r0 = csPar["r0"]; rMin = csPar["rMin"]; rMax = csPar["rMax"]; aL = tubeL[]; aR = tubeR[];
csResetChecks[];
t0 = AbsoluteTime[];
aGrid = Subdivide[aL, aR, 100];
extL = Subdivide[rMin, aL, 35]; extR = Subdivide[aR, rMax, 45];
z = zVar;

Print["=== Stage 2: first-order regular field, PhiR = c0 + c1 (r cos theta), c0 = ", csPar["c0"], ", c1 = ", csPar["c1"], " ==="];
Print["    q = ", csPar["q"], ", r0 = ", r0, ", tube [", aL, ", ", aR, "], domain [", rMin, ", ", rMax, "], WorkingPrecision ", csPar["wp"]];

rows = {};
testCase[label_, l_, data_, exact_] := Module[{res, full, seffMax, dCont, dAR, imAR, dExt, dFull, dFullAR, imFullAR, line, dLine, phiS = singularModeCorrected[l]},
   res = firstOrderRegular[l, data];
   seffMax = Max[Abs[seffFunction[l] /@ N[Most[Rest[aGrid]], 40]]];
   csCheck[label <> " max |S_eff| on the contour (notebook formula, exact phiS)", seffMax, 10^-20];
   dCont = scaledError[res["sol2"] /@ aGrid, N[exact[yC[#]], 40] & /@ aGrid];
   csCheck[label <> " contour phiR(y(a)) vs analytic continuation (101 pts, scaled)", dCont, 10^-12];
   dAR = scaledError[{res["sol2"][aR], res["sol2"]'[aR]}, N[{exact[aR], exact'[aR]}, 40]];
   csCheck[label <> " contour solution at aR (value + derivative) vs exact (scaled)", dAR, 10^-15];
   imAR = (Abs[res["imAR"][[1]]] + Abs[res["imAR"][[2]]])/Max[1, Abs[N[exact[aR]]]];
   csCheck[label <> " (|Im phiR(aR)| + |Im phiR'(aR)|)/max(1,|phiR(aR)|) (discarded by Re[valR])", imAR, 10^-15];
   dExt = scaledError[Join[res["sol1"] /@ extL, res["sol3"] /@ extR], N[exact /@ Join[extL, extR], 40]];
   csCheck[label <> " exterior phiR vs exact on [1,aL] and [aR,10] (scaled)", dExt, 10^-12];
   (* full field on the contour, source kappa_l deltaExt *)
   full = firstOrderFull[l, data];
   dFull = scaledError[full["sol2"] /@ aGrid, (kappaL[l] N[phiS[yC[#]], 40] + N[exact[yC[#]], 40]) & /@ aGrid];
   csCheck[label <> " full field on contour = phiR + kappa phiS (101 pts, scaled)", dFull, 10^-12];
   dFullAR = Abs[Re[full["sol2"][aR] - kappaL[l] N[phiS[aR], 40]] - N[exact[aR], 40]]/Max[1, Abs[N[exact[aR]]]];
   imFullAR = Abs[Im[full["sol2"][aR] - kappaL[l] N[phiS[aR], 40]]]/Max[1, Abs[N[exact[aR]]]];
   csCheck[label <> " full field at aR: Re(full - kappa phiS) = phiR (/max(1,|phiR(aR)|))", dFullAR, 10^-15];
   csCheck[label <> " full field at aR: Im(full - kappa phiS) = 0 (/max(1,|phiR(aR)|))", imFullAR, 10^-15];
   (* real line through the tube *)
   line = realLineSolve[l, Function[r, 0], If[data === Automatic, {exact[rMin], exact'[rMin]}, data]];
   dLine = scaledError[line /@ Subdivide[rMin, rMax, 90], N[exact /@ Subdivide[rMin, rMax, 90], 40]];
   csCheck[label <> " real-line solve through the tube vs exact (scaled)", dLine, 10^-12];
   AppendTo[rows, {label, seffMax, dCont, dAR, imAR, dExt, dFull, imFullAR, dLine}];
   res];

Print["\n--- 2a. Injected PhiR = c0 + c1 r cos(theta): l = 0, 1 nonzero, l = 2, 3 identically zero ---"];
resPhys = Association@Table[l -> testCase["PhiR l=" <> ToString[l], l, Automatic, phiRAnalytic[l]], {l, 0, 3}];

Print["\n--- 2b. Generic injected homogeneous field (notebook test: phi(1) = 1, phi'(1) = 0) ---"];
resGen = Association@Table[l -> testCase["generic l=" <> ToString[l], l, {1, 0}, homogeneousWithData[l, 1, {1, 0}]], {l, {2, 3, 6}}];

Print["\n--- 2c. Two-point BVP by linear shooting (values at r = 1 and r = 10) ---"];
bvpRows = {};
Do[
  Module[{exact, srcA, bvp, err, imAR},
   exact = If[l <= 1, phiRAnalytic[l],
     Module[{A, B}, toF = (A z^l + B z^(-l - 1)) /. First@Solve[{A + B == 1, A 10^l + B 10^(-l - 1) == 0}, {A, B}]; Function @@ {{zVar}, toF}]];
   bvp = threeRegionBVP[l, Function[r, 0], seffFunction[l], {exact[rMin], exact[rMax]}];
   err = scaledError[Join[bvp["sol1"] /@ extL, bvp["sol3"] /@ extR], N[exact /@ Join[extL, extR], 40]];
   imAR = Abs[bvp["imAR"][[1]]] + Abs[bvp["imAR"][[2]]];
   csCheck["BVP l=" <> ToString[l] <> If[l <= 1, " (PhiR)", " (phi(1)=1, phi(10)=0, notebook sec. 3)"] <> " exterior error", err, 10^-12];
   csCheck["BVP l=" <> ToString[l] <> " |Im| at aR", imAR, 10^-15];
   AppendTo[bvpRows, {l, err, imAR}]], {l, {0, 1, 2}}];

Print["\n--- 2d. Diagnostic: truncated puncture (Taylor of corrected phiS through order n), generic l = 2 ---"];
Print["    Im of the residual field at aR equals the Im part of the Taylor remainder; the reconstructed exterior phiR must be exact for every n."];
truncRows = {};
Do[
  Module[{res, exact = homogeneousWithData[2, 1, {1, 0}], imRes, dExt},
   res = firstOrderRegular[2, {1, 0}, n];
   imRes = Abs[Im[res["sol2"][aR]]] + Abs[Im[res["sol2"]'[aR]]];
   (* residual = phiR + kappa (phiS - puncture), so its Im part at aR is kappa Im(phiS - puncture)(aR) *)
   Module[{rem = Function @@ {{zVar}, kappaL[2] (singularModeCorrected[2][zVar] - singularModeTaylor[2, n][zVar])}},
    csCheck["n=" <> ToString[n] <> ": Im residual(aR) = kappa Im(phiS - puncture)(aR) (value and derivative)",
     Abs[Im[res["sol2"][aR]] - Im[N[rem[aR], 40]]] + Abs[Im[res["sol2"]'[aR]] - Im[N[rem'[aR], 40]]], 10^-15]];
   dExt = scaledError[res["sol3"] /@ extR, N[exact /@ extR, 40]];
   Print["    n = ", n, ":  |Im residual(aR)| + |Im residual'(aR)| = ", sci[imRes], ",   exterior error (right) = ", sci[dExt]];
   AppendTo[truncRows, {n, imRes, dExt}]], {n, {2, 4, 6, 8, 12}}];
csCheck["truncated puncture: exterior phiR error for all n", Max[truncRows[[All, 3]]], 10^-12];
csCheck["truncated puncture: residual shadow at aR decreases monotonically with n", If[OrderedQ[Reverse[truncRows[[All, 2]]]], 0, 1], 0];

(* ---------------------------------------------------------------- *)
header = {"case", "max|Seff|", "contour_vs_continuation", "aR_val+der_err", "Im_at_aR", "exterior_err", "full_field_err", "full_Im_at_aR", "real_line_err"};
Export[FileNameJoin[{resDir, "stage2_regular_field.csv"}], Prepend[Map[If[NumericQ[#], sci[#], #] &, rows, {2}], header], "CSV"];
Export[FileNameJoin[{resDir, "stage2_bvp_and_truncation.csv"}],
  Join[{{"BVP l", "exterior_err", "Im_at_aR"}}, Map[If[NumericQ[#] && ! IntegerQ[#], sci[#], #] &, bvpRows, {2}],
   {{}, {"puncture order n (l=2 generic)", "Im_residual_at_aR", "exterior_err"}}, Map[If[NumericQ[#] && ! IntegerQ[#], sci[#], #] &, truncRows, {2}]], "CSV"];

Print["\n--- 2e. Plots ---"];
plotCase[res_, exact_, title_] := Plot[{Re[res["sol2"][a]], Im[res["sol2"][a]], Re[N[exact[yC[a]]]], Im[N[exact[yC[a]]]]}, {a, aL, aR},
   PlotStyle -> {{Blue, Thick}, {Red, Thick}, {Black, Dashed}, {Gray, Dashed}},
   PlotLegends -> {"Re \!\(\*SuperscriptBox[\(\[Phi]\), \(R\)]\) contour", "Im \!\(\*SuperscriptBox[\(\[Phi]\), \(R\)]\) contour", "Re continuation", "Im continuation"},
   PlotLabel -> title, AxesLabel -> {"a", ""}, ImageSize -> 420];
errPlot[res_, exact_, title_] := ListLogPlot[{
    {#, Abs[res["sol1"][#] - N[exact[#], 40]] + 10^-30} & /@ Subdivide[rMin, aL, 70],
    {#, Abs[res["sol2"][#] - N[exact[yC[#]], 40]] + 10^-30} & /@ Subdivide[aL, aR, 40],
    {#, Abs[res["sol3"][#] - N[exact[#], 40]] + 10^-30} & /@ Subdivide[aR, rMax, 90]},
   Joined -> True, PlotLegends -> {"region 1 (real)", "contour (vs continuation)", "region 3 (real)"},
   PlotLabel -> title, AxesLabel -> {"r or a", "|error|"}, ImageSize -> 420, PlotRange -> {10^-30, 10^-8}];
Export[FileNameJoin[{resDir, "stage2_phiR_contour.pdf"}],
  Column[{plotCase[resPhys[0], phiRAnalytic[0], "\[ScriptL] = 0, \!\(\*SuperscriptBox[\(\[Phi]\), \(R\)]\) = \!\(\*SubscriptBox[\(c\), \(0\)]\)\!\(\*SqrtBox[\(4 \[Pi]\)]\)"],
    plotCase[resPhys[1], phiRAnalytic[1], "\[ScriptL] = 1, \!\(\*SuperscriptBox[\(\[Phi]\), \(R\)]\) = \!\(\*SubscriptBox[\(c\), \(1\)]\) z \!\(\*SqrtBox[\(4 \[Pi]/3\)]\)"],
    plotCase[resGen[2], homogeneousWithData[2, 1, {1, 0}], "\[ScriptL] = 2, generic homogeneous (notebook test)"]}]];
Export[FileNameJoin[{resDir, "stage2_phiR_errors.pdf"}],
  Column[{errPlot[resPhys[1], phiRAnalytic[1], "\[ScriptL] = 1 error"], errPlot[resGen[2], homogeneousWithData[2, 1, {1, 0}], "\[ScriptL] = 2 generic, error"]}]];
Print["   wrote stage2_phiR_contour.pdf, stage2_phiR_errors.pdf, stage2_regular_field.csv, stage2_bvp_and_truncation.csv"];

Print["\nStage 2 wall time: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage2_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
