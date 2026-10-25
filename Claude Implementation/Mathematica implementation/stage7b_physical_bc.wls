#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 7b: the l'max optimum and shrinking worldtube of Stage 7 (PRL Appendix D, fig_lprime_optimum) with physical boundary
   conditions, psi' = l psi/r at r = 1 and psi' = -(l+1) psi/r at r = 10 (as stage12b), instead of Cauchy data at r = 1.
   Same Lemma-1 decomposition: (c) = response to (truncated - exact source) over [1, 10], (b) = response to the shadow term on the
   tube; the response is now the free-space Green's function G = -(s^2/(2l+1)) r_<^l/r_>^(l+1).  The per-l' tube moments of Stage 7
   ((du, du') at aR from zero data at aL) are converted to the two moments m1 = Int s^(1-l) f, m2 = Int s^(l+2) f, which give
   -m1 r^l/(2l+1) inside and -m2 r^(-l-1)/(2l+1) outside.  Reference: exact total + homogeneous terms meeting the same conditions.
   The Stage 7 (Cauchy) outputs are kept; this writes results/stage7b_* and ../figures/fig_lprime_optimum.pdf.
   Run:  wolframscript -file stage7b_physical_bc.wls                                                                     *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
s7 = Import[FileNameJoin[{dir, "stage7_lprime_optimum.wls"}], "Text"];
s7 = StringTake[s7, {1, StringPosition[s7, "(* ================= main loop"][[1, 1]] - 1}];
s7 = StringReplace[s7, "dir = DirectoryName[$InputFileName];" -> "", 1];
s7 = StringReplace[s7, "Print[\"=== Stage 7" -> "Hold[\"x", 1];
ToExpression[s7];
csResetChecks[]; Off[General::munfl];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
Print["=== Stage 7b: l'max optimum with physical boundary conditions ==="];
GphysV[l_, r_, s_List] := -(s^2/(2 l + 1)) (UnitStep[r - s] s^l/r^(l + 1) + (1 - UnitStep[r - s]) r^l/s^(l + 1));
physResp[l_, nodes_, delta_, grid_] := Table[Total[nodes[[All, 2]] GphysV[l, r, nodes[[All, 1]]] delta], {r, grid}];
physFix[sol_, l_] := Module[{B = (analyticEvalD[sol, 1] - l analyticEval[sol, 1])/(2 l + 1), A = -(analyticEvalD[sol, 10] + (l + 1) analyticEval[sol, 10]/10)/((2 l + 1) 10^(l - 1))},
   Function[r, analyticEval[sol, r] + A r^l + B r^(-l - 1)]];
toPhys[l_, aRv_, {du_, ddu_}, grid_] := Module[{m1 = aRv^-l ((l + 1) du + aRv ddu), m2 = -aRv^(l + 1) (l du - aRv ddu)},
   Table[If[r < r0, -m1 r^l, -m2 r^(-l - 1)]/(2 l + 1), {r, grid}]];
(* check of the moment conversion: a test source on the tube, direct quadrature against G vs conversion of its (du, du') at aR *)
Module[{l = 3, nd = Join[gl[9/2, r0, 30], gl[r0, 11/2, 30]], f, du, direct},
  f = Sin[3 nd[[All, 1]]] + nd[[All, 1]]^2;
  du = {Total[nd[[All, 2]] (nd[[All, 1]]^2/(2 l + 1)) ((11/2)^l nd[[All, 1]]^(-l - 1) - nd[[All, 1]]^l (11/2)^(-l - 1)) f],
        Total[nd[[All, 2]] (nd[[All, 1]]^2/(2 l + 1)) (l (11/2)^(l - 1) nd[[All, 1]]^(-l - 1) + (l + 1) nd[[All, 1]]^l (11/2)^(-l - 2)) f]};
  direct = physResp[l, nd, f, gCom];
  csCheck["moment conversion (du, du') at aR -> free-space response = direct quadrature (relative)", Max[Abs[toPhys[l, 11/2, du, gCom] - direct]]/Max[Abs[direct]], 10^-25]];

rows = {};
Do[
  SetCSParameters[casePars[Q]];
  Do[Module[{c = r0^-(n + 1), nodes, mv, scale, exactSrc, cCurves = <||>, bFixed = <||>},
     nodes = cNodes[12];
     mv = modeVals[Q, n, c, #] & /@ nodes[[All, 1]];
     Print["    ", Q, " n=", n, ": node values cached [", tick[], " s]"];
     Do[
      scale[l] = Max[Abs[physFix[totalSol[Q, l], l] /@ gCom]];
      exactSrc = Table[Module[{s = nodes[[k, 1]]}, pwVal[Sum[pieceAnalyticSource[Q, p, l], {p, {"SS", "SR", "RR"}}], s] - pwVal[{ssPsrc[Q, n, c, l, "real"], ssPsrc[Q, n, c, l, "real"]}, s]], {k, Length[nodes]}];
      Do[cCurves[{l, M}] = physResp[l, nodes, Table[truncSrc[Q, l, M, mv[[k]], nodes[[k, 1]]], {k, Length[nodes]}] - exactSrc, gCom], {M, lpList}];
      Do[bFixed[{l, R}] = Accumulate[bPerLp[Q, n, c, R, l, lpTop, 30]], {R, Rfixed}];
      Do[
       Module[{cc = cCurves[{l, M}], schemes},
        schemes = Join[Table[{ToString[R, InputForm], R, bFixed[{l, R}][[M + 1]]}, {R, Rfixed}],
          {{"r0/l'max", r0/M, Total[bPerLp[Q, n, c, r0/M, l, M, 30]]}}];
        Do[Module[{bh = toPhys[l, r0 + sc[[2]], sc[[3]], gCom]},
           AppendTo[rows, {Q, n, sc[[1]], M, l, Max[Abs[cc + bh]]/scale[l], Max[Abs[bh]]/scale[l], Max[Abs[cc]]/scale[l]}]], {sc, schemes}]],
       {M, lpList}],
      {l, lOut}];
     Print["    ", Q, " n=", n, " done [", tick[], " s]"]], {n, nS}],
  {Q, {"Product", "GradDot"}}];

(* App D numbers: growth per coupled mode of (b) near l' = 30 (n = 3), and the fall of (b) per unit n at l'max = 12 (R = 1/2) *)
growthRows = Flatten[Table[SetCSParameters[casePars[Q]]; Module[{n = 3, c = r0^-4, per, sc},
     per = Table[Module[{pairs = bPerLp[Q, n, c, R, l, 30, 30], s = Max[Abs[physFix[totalSol[Q, l], l] /@ gCom]]}, Table[Max[Abs[toPhys[l, r0 + R, pairs[[lp + 1]], gCom]]]/s, {lp, 0, 30}]], {l, lOut}];
     sc = Max /@ Transpose[per];                                                                  (* per l', max over l *)
     {Q, ToString[R, InputForm], N[(sc[[31]]/sc[[26]])^(1/5), 4], sc[[31]]}], {Q, {"Product", "GradDot"}}, {R, Rfixed}], 1];
nRows = Flatten[Table[SetCSParameters[casePars[Q]]; Table[Module[{c = r0^-(n + 1)},
      {Q, n, Max[Table[Module[{s = Max[Abs[physFix[totalSol[Q, l], l] /@ gCom]]}, Max[Abs[toPhys[l, r0 + 1/2, Total[bPerLp[Q, n, c, 1/2, l, 12, 30]], gCom]]]/s], {l, lOut}]]}], {n, {1, 3, 5, 7}}], {Q, {"Product", "GradDot"}}], 1];
perN = Table[{Q, N[10^Coefficient[Fit[{#[[2]], Log10[#[[3]]]} & /@ Select[nRows, #[[1]] === Q &], {1, x}, x], x], 3]}, {Q, {"Product", "GradDot"}}];
SetCSParameters[casePars["Product"]];

(* ================= tables, figure ================= *)
mxr[Q_, n_, sch_, M_, k_] := Max[Cases[rows, {Q, n, sch, M, __}][[All, k]]];
schemes = {"1/2", "1/4", "r0/l'max"};
md = StringRiffle[Join[{"| case | n | R | quantity | " <> StringRiffle[ToString /@ lpList, " | "] <> " |", "|---|---|---|---|" <> StringRepeat["---|", Length[lpList]]},
    Flatten@Table[Table["| " <> Q <> " | " <> ToString[n] <> " | " <> sch <> " | " <> {"total", "(b)", "(c)"}[[k - 5]] <> " | " <> StringRiffle[sci[mxr[Q, n, sch, #, k]] & /@ lpList, " | "] <> " |",
       {k, 6, 8}], {Q, {"Product", "GradDot"}}, {n, nS}, {sch, schemes}]], "\n"];
optRows = Flatten[Table[Module[{tot = mxr[Q, n, sch, #, 6] & /@ lpList, k},
     k = First[Ordering[tot, 1]]; {Q, n, sch, lpList[[k]], tot[[k]], tot[[-1]], If[OrderedQ[Reverse[tot]], "monotone decreasing", "not monotone"]}], {Q, {"Product", "GradDot"}}, {n, nS}, {sch, schemes}], 2];
md2 = StringRiffle[Join[{"| case | n | R | argmin l'max | min total | total at l'max = 40 | behaviour |", "|---|---|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> #[[3]] <> " | " <> ToString[#[[4]]] <> " | " <> sci[#[[5]]] <> " | " <> sci[#[[6]]] <> " | " <> #[[7]] <> " |") & /@ optRows], "\n"];
md3 = StringRiffle[Join[{"| case | R | growth per l' near l' = 30 (n = 3) | per-mode (b) at l' = 30 |", "|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> #[[2]] <> " | " <> ToString[#[[3]]] <> " | " <> sci[#[[4]]] <> " |") & /@ growthRows], "\n"];
md4 = "(b) at l'max = 12, R = 1/2, falls per unit n by: " <> StringRiffle[(#[[1]] <> " " <> ToString[#[[2]]]) & /@ perN, ", "];
Print["\n", md2, "\n\n", md3, "\n\n", md4];
Export[FileNameJoin[{resDir, "stage7b_summary.md"}],
  "# Stage 7b: l'max optimum and shrinking worldtube, physical boundary conditions\n\nScaled exterior error (max over l = 0..6 and r in [1,3.5] U [6.5,10]).\n\n" <>
   md2 <> "\n\n## growth of (b) per coupled mode\n\n" <> md3 <> "\n\n" <> md4 <> "\n\n## Errors against l'max\n\n" <> md <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage7b_errors.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], sci[#[[6]]], sci[#[[7]]], sci[#[[8]]]} & /@ rows, {"case", "n", "R", "lpmax", "l", "total", "b", "c"}], "CSV"];
$fs = {FontFamily -> "Times", FontSize -> 9};
schStyle = <|"1/2" -> RGBColor[0.1, 0.35, 0.75], "1/4" -> RGBColor[0.85, 0.45, 0.1], "r0/l'max" -> RGBColor[0.15, 0.6, 0.25]|>;
panels = Table[ListLogPlot[Join[Flatten[Table[Table[{M, mxr[Q, n, sch, M, 6]}, {M, lpList}], {n, nS}, {sch, schemes}], 1], Table[Table[{M, mxr[Q, n, "1/2", M, 8]}, {M, lpList}], {n, nS}]],
    Joined -> True, PlotMarkers -> Join[Flatten[Table[{{"\[FilledCircle]", 5}, {"\[FilledSquare]", 5}, {"\[FilledDiamond]", 6}}, {n, nS}], 1], {None, None}],
    PlotStyle -> Join[Flatten[Table[Directive[schStyle[sch], AbsoluteThickness[1.3], If[n == 3, Dashing[{}], Dashed]], {n, nS}, {sch, schemes}], 1],
      Table[Directive[Gray, AbsoluteThickness[0.8], If[n == 3, Dotted, DotDashed]], {n, nS}]],
    Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]'\!\(\*SubscriptBox[\(\[ThinSpace]\), \(max\)]\)", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 212,
    Epilog -> {Inset[Style[If[Q === "Product", "1/\!\(\*SuperscriptBox[\(\[Rho]\), \(2\)]\)", "1/\!\(\*SuperscriptBox[\(\[Rho]\), \(4\)]\)"], $fs], Scaled[{0.97, 0.97}], {Right, Top}]}], {Q, {"Product", "GradDot"}}];
leg = Column[{LineLegend[Values[schStyle], {"R = 1/2", "R = 1/4", "R = \!\(\*SubscriptBox[\(r\), \(0\)]\)/\[ScriptL]'\!\(\*SubscriptBox[\(\[ThinSpace]\), \(max\)]\)"}, LegendLayout -> "Row", LabelStyle -> $fs],
    LineLegend[{Directive[Black], Directive[Black, Dashed], Directive[Gray, Dotted], Directive[Gray, DotDashed]}, {"total, n = 3", "total, n = 5", "(c) only, n = 3", "(c) only, n = 5"}, LegendLayout -> "Row", LabelStyle -> $fs]}, Alignment -> Center];
Export[FileNameJoin[{figDir, "fig_lprime_optimum.pdf"}], Legended[Grid[{panels}, Spacings -> 1], Placed[leg, Below]]];
Export[FileNameJoin[{figDir, "fig_lprime_optimum.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], N[#[[6]], 6], N[#[[7]], 6], N[#[[8]], 6]} & /@ rows, {"case", "n", "R", "lpmax", "l", "total", "b", "c"}], "CSV"];
Print["\nStage 7b wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage7b_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
