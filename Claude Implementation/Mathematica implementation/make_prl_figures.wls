#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 13: PRL figures for the Dzz case, Q = (1/2) d_z^2 (A B) (exact puncture q/rho, PhiR = c0 + c1 z, r0 = 5, tube [4.5, 5.5], amplitude 1).
   Writes ../figures/prl/{prl_intro, prl_source, prl_solution, prl_convergence}.pdf with a CSV beside each.
   prl_convergence reads results/stage12b_routes.csv (routes with physical boundary conditions, from stage12b_physical_bc.wls).
   Run:  wolframscript -file make_prl_figures.wls                                                                       *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
s9 = Import[FileNameJoin[{dir, "stage9_dzz.wls"}], "Text"];
defsA = StringTake[s9, {1, StringPosition[s9, "(* ================= checks of the d_z machinery"][[1, 1]] - 1}];
defsA = StringReplace[StringReplace[defsA, "dir = DirectoryName[$InputFileName];" -> "", 1], "Print[\"=== Stage 9" -> "Hold[\"x", 1];
ToExpression[defsA];                                               (* Stage 6-9 definitions; SetCSParameters[casePars["Product"]] *)
prlDir = FileNameJoin[{ParentDirectory[dir], "figures", "prl"}]; Quiet[CreateDirectory[prlDir]];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
csResetChecks[]; Off[N::meprec, General::munfl];
q = csPar["q"]; c0 = csPar["c0"]; c1 = csPar["c1"];
(* style: 8 pt Times, thin frames, widths 8.6 cm (one column) and 17.8 cm (two columns) *)
cm = 72/2.54; w1 = 8.6 cm; w2 = 17.8 cm;
fs = {FontFamily -> "Times", FontSize -> 8};
frame = {FrameStyle -> Directive[Black, AbsoluteThickness[0.5]], FrameTicksStyle -> Directive[Black, 8], BaseStyle -> fs, Frame -> True, Axes -> False};
colSS = RGBColor[0.12, 0.38, 0.70]; colSR = RGBColor[0.85, 0.40, 0.05]; colRR = RGBColor[0.20, 0.60, 0.25];
colL = <|0 -> RGBColor[0.12, 0.38, 0.70], 2 -> RGBColor[0.85, 0.40, 0.05], 4 -> RGBColor[0.20, 0.60, 0.25]|>;
colRoute = <|"product" -> RGBColor[0.80, 0.15, 0.15], "H0" -> GrayLevel[0.35], "H1" -> RGBColor[0.12, 0.38, 0.70], "H1 + tail" -> RGBColor[0.20, 0.60, 0.25], "(c)" -> Black|>;
tubeShade[{y1_, y2_}] := {GrayLevel[0.9], Rectangle[{aL, y1}, {aR, y2}]};
Print["=== Stage 13: PRL figures -> ", prlDir, " ==="];

(* ---------------- sources, l = 2 ---------------- *)
l = 2;
ssRe = Function @@ {{zVar}, N[ssSourceModeExpr["Dzz", l, zVar, "real"] /. Log[x_^2] :> 2 Log[Abs[x]], 40]};
ssCo = Function @@ {{zVar}, N[ssSourceModeExpr["Dzz", l, zVar, "continued"], 40]};
srPw = pieceAnalyticSource["Dzz", "SR", l];
srRe[r_] := N[If[r < r0, srPw[[1]], srPw[[2]]] /. zVar -> r, 30];
srCo = Function @@ {{zVar}, N[(1/2) dzzExpr[Function[j, If[j < 0, 0, 2 q Sum[KCouplingContinued[j, lp, zVar, "ThetaToH"] phiRAnalytic[lp][zVar], {lp, 0, 1}]]], l], 40]};
cl2 = N[(1/2) dzzDeltaCoeff[pieceAnalyticSource["Product", "SR", #] &, l], 12];
A2 = Module[{u, e = ssSourceModeExpr["Dzz", l, zVar, "real"] /. Log[x_^2] :> 2 Log[Abs[x]]}, N[Coefficient[Normal[Series[(e /. zVar -> r0 + u) /. Abs[u] -> u, {u, 0, 0}]], u, -2], 12]];
csCheck["c_2 = -4 Pi q PhiR(x0) Y_2(0)/r0^2", Abs[cl2 + 4 Pi q (c0 + c1 r0) Sqrt[5/(4 Pi)]/r0^2], 10^-10];
(* SR is bounded at r0: its in/out branches are rational (Coulomb modes times PhiR), with a jump and a delta only *)
csCheck["SR l=2 in/out branches contain no Log and are finite at r0 (bounded: jump + delta only)",
  If[FreeQ[srPw, Log] && And @@ (NumericQ[N[# /. zVar -> r0]] & /@ srPw), 0, 1], 0];
srJumpVal = N[(srPw[[2]] - srPw[[1]]) /. zVar -> r0, 12];
rr = Join[Subdivide[3, 4.99, 400], r0 - 10^-Range[2.5, 6, 0.5], r0 + Reverse[10^-Range[2.5, 6, 0.5]], Subdivide[5.01, 7, 400]];
srcTab = Table[{r, ssRe[r], srRe[r], A2/(r - r0)^2}, {r, rr}];
aa = Subdivide[aL, aR, 200];
coTab = Table[Module[{z = N[yC[a], 30]}, {a, Re[ssCo[z]], Im[ssCo[z]], Re[srCo[z]], Im[srCo[z]]}], {a, aa}];
Module[{yr = {-1.2, 2.2}, pa, pb, srJump},
  pa = Show[
    ListLinePlot[{Select[srcTab[[All, {1, 2}]], #[[1]] < r0 &], Select[srcTab[[All, {1, 2}]], #[[1]] > r0 &], Select[srcTab[[All, {1, 3}]], #[[1]] < r0 &], Select[srcTab[[All, {1, 3}]], #[[1]] > r0 &],
       Select[srcTab[[All, {1, 4}]], #[[1]] < r0 &], Select[srcTab[[All, {1, 4}]], #[[1]] > r0 &]},
     PlotStyle -> {Directive[colSS, AbsoluteThickness[1.2]], Directive[colSS, AbsoluteThickness[1.2]], Directive[colSR, AbsoluteThickness[1.2]], Directive[colSR, AbsoluteThickness[1.2]],
       Directive[colSS, Dashed, AbsoluteThickness[0.7]], Directive[colSS, Dashed, AbsoluteThickness[0.7]]},
     PlotRange -> {{3, 7}, yr}, Evaluate[Sequence @@ frame], FrameLabel -> {"r", "source mode, \[ScriptL] = 2"}, ImageSize -> {w2/2 - 6, 0.62 (w2/2)},
     Prolog -> {tubeShade[yr]},
     Epilog -> {Text[Style["SS", fs, colSS], {4.15, 0.75}], Text[Style["SR", fs, colSR], {3.45, -0.32}],
       Text[Style["A\:2082/\[CapitalDelta]r\.b2", fs, colSS], {5.45, 1.85}, {-1, 0}], Text[Style["(a)", fs], Scaled[{0.04, 0.94}], {-1, 0}]}]];
  pb = ListLinePlot[{coTab[[All, {1, 2}]], coTab[[All, {1, 3}]], coTab[[All, {1, 4}]], coTab[[All, {1, 5}]]},
    PlotStyle -> {Directive[colSS, AbsoluteThickness[1.2]], Directive[colSS, Dashed, AbsoluteThickness[1.0]], Directive[colSR, AbsoluteThickness[1.2]], Directive[colSR, Dashed, AbsoluteThickness[1.0]]},
    PlotRange -> All, Evaluate[Sequence @@ frame], FrameLabel -> {"a  (z = y(a) on the contour)", "continued source, \[ScriptL] = 2"}, ImageSize -> {w2/2 - 6, 0.62 (w2/2)},
    Epilog -> {Text[Style["(b)", fs], Scaled[{0.04, 0.94}], {-1, 0}]}];
  Export[FileNameJoin[{prlDir, "prl_source.pdf"}], Legended[Grid[{{pa, pb}}, Spacings -> 0.5],
     Placed[LineLegend[{Directive[colSS], Directive[colSS, Dashed], Directive[colSR], Directive[colSR, Dashed]}, {"Re SS", "Im SS", "Re SR", "Im SR"}, LabelStyle -> fs, LegendLayout -> {"Row", 1}], Top]]];
  Export[FileNameJoin[{prlDir, "prl_source.csv"}], Join[{{"panel", "x", "SS_or_ReSS", "SR_or_ImSS", "A2_over_Dr2_or_ReSR", "ImSR"}},
     {"a_real_line", N[#[[1]], 8], N[#[[2]], 10], N[#[[3]], 10], N[#[[4]], 10], ""} & /@ srcTab, {"b_contour", N[#[[1]], 8], N[#[[2]], 10], N[#[[3]], 10], N[#[[4]], 10], N[#[[5]], 10]} & /@ coTab,
     {{"delta_coefficient_SR_l2", N[cl2, 12], "", "", "", ""}, {"SR_jump_out_minus_in_at_r0", srJumpVal, "", "", "", ""}, {"A2", N[A2, 12], "", "", "", ""}}], "CSV"]];
Print["    prl_source done [", tick[], " s]"];

(* ---------------- intro figure: the contour in the complex r plane, and the l = 2 source magnitude ---------------- *)
Module[{xr = {3.5, 6.5}, aC = Subdivide[aL, aR, 300], rReal, magR, magC, pa, pb, ymax},
  rReal = Join[Subdivide[3.5, 4.99, 300], r0 - 10^-Range[2.5, 6, 0.5], r0 + Reverse[10^-Range[2.5, 6, 0.5]], Subdivide[5.01, 6.5, 300]];
  magR = {#, ssRe[#] + srRe[#]} & /@ rReal;                                             (* pointwise SS + SR on the real axis (delta terms not shown) *)
  magC = Module[{z = N[yC[#], 30]}, {#, Re[ssCo[z] + srCo[z]]}] & /@ aC;                  (* Re of the continued SS + SR on the contour *)
  joinErr = Table[Module[{z = N[yC[a], 30]}, Abs[Re[ssCo[z] + srCo[z]] - (ssRe[a] + srRe[a])]/Abs[ssRe[a] + srRe[a]]], {a, {aL, aR}}];
  csCheck["prl_intro (b): Re(SS+SR) on the contour joins the real-axis SS+SR at a_L and a_R (relative)", Max[joinErr], 10^-12];
  ymax = 1.15;
  pa = Graphics[{GrayLevel[0.9], Rectangle[{aL, -0.18}, {aR, ymax}],
     Directive[colSS, AbsoluteThickness[1.4]], Line[{{xr[[1]], 0}, {aL, 0}}], Line[{{aR, 0}, {xr[[2]], 0}}],
     Arrowheads[0.045], Arrow[{{3.75, 0}, {3.52, 0}}], Arrow[{{6.25, 0}, {6.48, 0}}],
     Directive[RGBColor[0.80, 0.15, 0.15], AbsoluteThickness[1.4]], Line[{Re[#], Im[#]} & /@ (N[yC[#]] & /@ aC)],
     Black, PointSize[0.018], Point[{r0, 0}], Point[{aL, 0}], Point[{aR, 0}],
     Text[Style[Subscript["r", "0"], fs], {r0, -0.09}], Text[Style[Subscript["a", "L"], fs], {aL, -0.09}], Text[Style[Subscript["a", "R"], fs], {aR, -0.09}],
     Text[Style["to r = 1", fs, colSS], {3.56, 0.12}, {-1, 0}], Text[Style["to r = 10", fs, colSS], {6.44, 0.12}, {1, 0}],
     Text[Style["z = y(a)", fs, RGBColor[0.80, 0.15, 0.15]], {r0 + 0.33, 0.95}], Text[Style["(a)", fs], Scaled[{0.03, 0.9}], {-1, 0}]},
    Frame -> True, FrameStyle -> Directive[Black, AbsoluteThickness[0.5]], FrameTicksStyle -> Directive[Black, 8], BaseStyle -> fs,
    FrameLabel -> {None, "Im z"}, FrameTicks -> {{Automatic, None}, {{Automatic, None}, None}}, PlotRange -> {xr, {-0.18, ymax}}, AspectRatio -> 0.42, ImageSize -> w1];
  pb = Show[Graphics[{GrayLevel[0.9], Rectangle[{aL, -1}, {aR, 3}]}],
    ListLinePlot[{Select[magR, #[[1]] < r0 &], Select[magR, #[[1]] > r0 &], magC},
     PlotStyle -> {Directive[colSS, AbsoluteThickness[1.2]], Directive[colSS, AbsoluteThickness[1.2]], Directive[RGBColor[0.80, 0.15, 0.15], AbsoluteThickness[1.4]]}],
    PlotRange -> {xr, {-1, 3}}, Evaluate[Sequence @@ frame], FrameLabel -> {"a  (r = a on the real axis)", "SS + SR, \[ScriptL] = 2"}, AspectRatio -> 0.55, ImageSize -> w1,
    Epilog -> {Text[Style["(b)", fs], Scaled[{0.03, 0.9}], {-1, 0}], Text[Style["real axis (\[Delta] terms not shown)", fs, colSS], Scaled[{0.03, 0.08}], {-1, 0}],
      Text[Style["source on contour", fs, RGBColor[0.80, 0.15, 0.15]], Scaled[{0.62, 0.12}], {-1, 0}], Text[Style["\[Proportional] \[CapitalDelta]r\:207b\.b2", fs, colSS], Scaled[{0.56, 0.92}], {-1, 0}]}];
  Export[FileNameJoin[{prlDir, "prl_intro.pdf"}], Column[{pa, pb}, Spacings -> 0.2]];
  Export[FileNameJoin[{prlDir, "prl_intro.csv"}], Join[{{"panel", "a", "x", "y_or_magnitude"}},
     {"a_contour", N[#, 8], N[Re[yC[#]], 10], N[Im[yC[#]], 10]} & /@ aC,
     {"b_real_axis_SS_plus_SR", N[#[[1]], 10], N[#[[1]], 10], N[#[[2]], 10]} & /@ magR, {"b_contour_Re_SS_plus_SR", N[#[[1]], 10], "", N[#[[2]], 10]} & /@ magC], "CSV"]];
Print["    prl_intro done [", tick[], " s]"];

(* ---------------- recovered exterior Phi2: contour solves with physical boundary conditions ----------------
   psi' = l psi/r at r = 1 and psi' = -(l+1) psi/r at r = 10 (as stage12b).  Each piece is solved from zero Cauchy data at r = 1
   (three-region contour solve for SS and SR, real line for RR) and alpha r^l is added to meet the condition at r = 10.
   Reference: (1/2) d_z^2 chi + A r^l + B r^(-l-1) satisfying the same conditions.  Point-supported fit: error = alpha [1/rho]_l + beta d_r0 [1/rho]_l. *)
rIn = Subdivide[1, aL, 60]; rOut = Subdivide[aR, 10, 80]; aIn = Subdivide[aL, aR, 100];
solRows = {}; tubeRows = {}; physRows = {};
physFix[sol_, l_] := Module[{B = (analyticEvalD[sol, 1] - l analyticEval[sol, 1])/(2 l + 1), A = -(analyticEvalD[sol, 10] + (l + 1) analyticEval[sol, 10]/10)/((2 l + 1) 10^(l - 1))},
   Function[r, analyticEval[sol, r] + A r^l + B r^(-l - 1)]];
invRho[l_, r_, R0_] := Sqrt[4 Pi/(2 l + 1)] If[r < R0, r^l/R0^(l + 1), R0^l/r^(l + 1)];
pointFit[l_, errs_, grid_] := Module[{X = Table[N[{invRho[l, r, r0], D[invRho[l, r, R0], R0] /. R0 -> r0}, 30], {r, grid}]}, LeastSquares[X, errs]];
Do[Module[{solSS = pieceAnalyticSolution["Dzz", "SS", ll], solSR = pieceAnalyticSolution["Dzz", "SR", ll], solRR = pieceAnalyticSolution["Dzz", "RR", ll], chi, exT, exS, exR, scale,
    rS, rR, lR, fi, fo, srcR, alS, alR, alRR, ext, c1, c2, rrC, grid = Join[rIn, rOut], eS, eR, eT, fitS, fitR, gam, negc},
   chi = solSS + solSR + solRR; exT = physFix[chi, ll]; exS = physFix[solSS, ll]; exR = physFix[solSR, ll];
   scale = Max[Abs[exT /@ gCom]];
   rS = threeRegionSolve[ll, numericFunction[ssSourceModeExpr["Dzz", ll, zVar, "real"], zVar], contourFunction[ssSourceModeExpr["Dzz", ll, zVar, "continued"], zVar], {0, 0}];
   fi = numericFunction[pieceAnalyticSource["Dzz", "SR", ll][[1]], zVar]; fo = numericFunction[pieceAnalyticSource["Dzz", "SR", ll][[2]], zVar];
   srcR[r_?NumericQ] := If[r < r0, fi[r], fo[r]];
   rR = threeRegionSolve[ll, srcR, contourFunction[(1/2) dzzExpr[Function[j, If[j < 0, 0, 2 q Sum[KCouplingContinued[j, lp, zVar, "ThetaToH"] phiRAnalytic[lp][zVar], {lp, 0, 1}]]], ll], zVar], {0, 0}];
   lR = realLineSolve[ll, numericFunction[pieceAnalyticSource["Dzz", "RR", ll][[1]], zVar], {0, 0}];
   ext[res_, r_] := exteriorValue[res, r];
   alS = -(rS["sol3"]'[10] + (ll + 1) rS["sol3"][10]/10)/((2 ll + 1) 10^(ll - 1));
   alR = -(rR["sol3"]'[10] + (ll + 1) rR["sol3"][10]/10)/((2 ll + 1) 10^(ll - 1));
   alRR = -(lR'[10] + (ll + 1) lR[10]/10)/((2 ll + 1) 10^(ll - 1));
   eS = (ext[rS, #] + alS #^ll - exS[#]) & /@ grid; eR = (ext[rR, #] + alR #^ll - exR[#]) & /@ grid;
   eT = (ext[rS, #] + ext[rR, #] + lR[#] + (alS + alR + alRR) #^ll - exT[#]) & /@ grid;
   Do[AppendTo[solRows, {ll, grid[[k]], ext[rS, grid[[k]]] + ext[rR, grid[[k]]] + lR[grid[[k]]] + (alS + alR + alRR) grid[[k]]^ll, exT[grid[[k]]], Abs[eT[[k]]]/scale}], {k, Length[grid]}];
   fitS = pointFit[ll, eS, grid]; fitR = pointFit[ll, eR, grid];
   (* negative control: pointwise SR without its delta misses gam G(r, r0), G the free-space Green's function *)
   gam = (1/2) N[dzzDeltaCoeff[pieceAnalyticSource["Product", "SR", #] &, ll], 40];
   negc = Max[Abs[(gam (-(r0^2/(2 ll + 1)) If[# < r0, #^ll/r0^(ll + 1), r0^ll/#^(ll + 1)])) & /@ gCom]]/scale;
   AppendTo[physRows, {ll, Max[Abs[eS]]/scale, Max[Abs[eR]]/scale, Max[Abs[eT]]/scale, Abs[fitS[[1]]] + Abs[fitS[[2]]], Abs[fitR[[1]]] + Abs[fitR[[2]]], negc}];
   (* inside the tube: psi~(y(a)) = contour solves + RR continued + the alpha z^l correction; RR from zero data = solRR + c1 r^l + c2 r^(-l-1) *)
   {c1, c2} = Module[{v = analyticEval[solRR, 1], d = analyticEvalD[solRR, 1]}, LinearSolve[{{1, 1}, {ll, -(ll + 1)}}, {-v, -d}]];
   Do[Module[{z = N[yC[a], 30], tot}, tot = rS["sol2"][a] + rR["sol2"][a] + N[solRR[[1]] /. zVar -> z, 30] + c1 z^ll + c2 z^(-ll - 1) + (alS + alR + alRR) z^ll;
      AppendTo[tubeRows, {ll, a, Re[tot], Im[tot]}]], {a, aIn}];
   csCheck["l=" <> ToString[ll] <> ": recovered exterior Phi2 = exact, physical boundary conditions (scaled)", Max[Abs[eT]]/scale, 10^-10];
   Print["    l=", ll, ": total ", sci[Max[Abs[eT]]/scale], ", fit SS ", sci[physRows[[-1, 5]]], ", SR ", sci[physRows[[-1, 6]]], ", no-delta control ", sci[negc], " [", tick[], " s]"]], {ll, 0, 6}];
Export[FileNameJoin[{resDir, "stage9b_exact_physical.csv"}], Prepend[{#[[1]], Sequence @@ (sci /@ #[[2 ;; 7]])} & /@ physRows,
   {"l", "SS_err", "SR_err", "total_err", "SS_point_fit_abs_alpha_plus_beta", "SR_point_fit_abs_alpha_plus_beta", "pointwise_SR_no_delta_err"}], "CSV"];
Module[{pa, pb, lsel = {0, 2, 4}, yr, cols7 = Table[Blend[{RGBColor[0.12, 0.38, 0.70], RGBColor[0.20, 0.60, 0.25], RGBColor[0.85, 0.40, 0.05]}, ll/6], {ll, 0, 6}]},
  yr = {Min[Cases[solRows, {ll_ /; MemberQ[lsel, ll], _, _, ex_, _} :> ex]], Max[Cases[solRows, {ll_ /; MemberQ[lsel, ll], _, _, ex_, _} :> ex]]};
  yr = yr + {-0.08, 0.08} (yr[[2]] - yr[[1]]);
  yr = {Min[yr[[1]], Min[Cases[tubeRows, {ll_ /; MemberQ[lsel, ll], _, re_, im_} :> Min[re, im]]]], Max[yr[[2]], Max[Cases[tubeRows, {ll_ /; MemberQ[lsel, ll], _, re_, im_} :> Max[re, im]]]]};
  pa = Show[Graphics[tubeShade[yr]],
    ListLinePlot[Join[Flatten[Table[{Cases[solRows, {ll, r_, _, ex_, _} /; r <= aL :> {r, ex}], Cases[solRows, {ll, r_, _, ex_, _} /; r >= aR :> {r, ex}]}, {ll, lsel}], 1],
       Flatten[Table[{Cases[tubeRows, {ll, a_, re_, _} :> {a, re}], Cases[tubeRows, {ll, a_, _, im_} :> {a, im}]}, {ll, lsel}], 1]],
     PlotStyle -> Join[Flatten[Table[{#, #} &[Directive[colL[ll], AbsoluteThickness[1.0]]], {ll, lsel}]],
       Flatten[Table[{Directive[colL[ll], AbsoluteThickness[1.0]], Directive[colL[ll], Dashed, AbsoluteThickness[0.8]]}, {ll, lsel}]]], PlotRange -> {{1, 10}, yr}, Evaluate[Sequence @@ frame],
     FrameLabel -> {None, "\[CapitalPhi]\:2082,\[ScriptL](r)"}, FrameTicks -> {{Automatic, None}, {{Automatic, None}, None}}, ImageSize -> {w1, 0.55 w1}],
    ListPlot[Table[Cases[solRows, {ll, r_, rec_, _, _} :> {r, rec}][[;; ;; 6]], {ll, lsel}], PlotMarkers -> {"\[EmptyCircle]", 6}, PlotStyle -> Table[colL[ll], {ll, lsel}]],
    Epilog -> {Text[Style["(a)", fs], Scaled[{0.03, 0.92}], {-1, 0}],
      Inset[LineLegend[Table[colL[ll], {ll, lsel}], Table["\[ScriptL] = " <> ToString[ll], {ll, lsel}], LabelStyle -> fs, LegendLayout -> {"Column", 1}], Scaled[{0.03, 0.66}], {Left, Top}],
      Text[Style["lines: exact solution", fs], Scaled[{0.98, 0.62}], {1, 0}], Text[Style["circles: contour solve", fs], Scaled[{0.98, 0.55}], {1, 0}],
      Text[Style["tube: Re \[Psi]\:0303(y(a)) solid, Im dashed", fs], Scaled[{0.98, 0.48}], {1, 0}]},
    PlotRange -> {{1, 10}, yr}, Evaluate[Sequence @@ frame], FrameLabel -> {None, "\[CapitalPhi]\:2082,\[ScriptL]"}, FrameTicks -> {{Automatic, None}, {{Automatic, None}, None}}, AspectRatio -> Full, ImageSize -> {w1, 0.55 w1}, ImagePadding -> {{44, 6}, {6, 6}}];
  pb = ListLogPlot[Flatten[Table[{Cases[solRows, {ll, r_, _, _, e_} /; r <= aL :> {r, e + 10^-30}], Cases[solRows, {ll, r_, _, _, e_} /; r >= aR :> {r, e + 10^-30}]}, {ll, 0, 6}], 1],
    Joined -> True, PlotStyle -> Flatten[Table[{#, #} &[Directive[cols7[[ll + 1]], AbsoluteThickness[0.8]]], {ll, 0, 6}]], PlotRange -> {{1, 10}, {10^-22, 10^-14}}, Evaluate[Sequence @@ frame],
    FrameLabel -> {"a  (r = a outside the tube)", "scaled error"}, ImageSize -> {w1, 0.5 w1}, AspectRatio -> Full, ImagePadding -> {{44, 6}, {32, 4}}, Prolog -> {GrayLevel[0.9], Rectangle[{aL, Log[10.^-22]}, {aR, Log[10.^-14]}]},
    Epilog -> {Text[Style["(b)  \[ScriptL] = 0\[Ellipsis]6", fs], Scaled[{0.03, 0.9}], {-1, 0}]}];
  Export[FileNameJoin[{prlDir, "prl_solution.pdf"}], Column[{pa, pb}, Spacings -> 0.3]];
  Export[FileNameJoin[{prlDir, "prl_solution.csv"}], Join[Prepend[{#[[1]], N[#[[2]], 8], N[#[[3]], 16], N[#[[4]], 16], N[#[[5]], 6]} & /@ solRows, {"l", "r = a (exterior)", "recovered", "exact", "scaled_error"}],
     Prepend[{#[[1]], N[#[[2]], 8], N[#[[3]], 16], N[#[[4]], 16], ""} & /@ tubeRows, {"l", "a (tube, z = y(a))", "Re psi_contour", "Im psi_contour", ""}]], "CSV"]];
Print["    prl_solution done [", tick[], " s]"];

(* ---------------- convergence: Stage 6 puncture, n = 3, 5 ---------------- *)
(* routes with physical boundary conditions (stage12b_physical_bc.wls) *)
numS[x_] := If[StringQ[x], ToExpression[StringReplace[x, "e" -> "*^"]], x];
convRows = {#[[1]], #[[2]], #[[3]], numS[#[[4]]]} & /@ Rest[Import[FileNameJoin[{resDir, "stage12b_routes.csv"}], "CSV"]];
csCheck["H1 = (c) in the physical-BC route table (max relative difference)", Max[Table[Abs[Cases[convRows, {n, "H1", M, e_} :> e][[1]] - Cases[convRows, {n, "(c)", M, e_} :> e][[1]]]/Cases[convRows, {n, "(c)", M, e_} :> e][[1]], {n, {3, 5}}, {M, {4, 8, 12, 16, 20, 30, 40}}]], 10^-6];
Module[{pan},
  pan[n_, lab_, xl_] := ListLogPlot[Table[Cases[convRows, {n, rt, M_, e_} :> {M, e}], {rt, {"product", "H0", "H1", "H1 + tail", "(c)"}}],
    Joined -> True, PlotMarkers -> {{"\[FilledCircle]", 4}, {"\[FilledSquare]", 4}, {"\[FilledDiamond]", 5}, {"\[FilledUpTriangle]", 3}, None},
    PlotStyle -> {Directive[colRoute["product"], AbsoluteThickness[1.0]], Directive[colRoute["H0"], AbsoluteThickness[1.0]], Directive[colRoute["H1"], AbsoluteThickness[1.0]],
      Directive[colRoute["H1 + tail"], AbsoluteThickness[1.0]], Directive[colRoute["(c)"], Dashed, AbsoluteThickness[0.7]]},
    PlotRange -> {{0, 42}, {10^-12, 30}}, Evaluate[Sequence @@ frame], FrameLabel -> {xl, "scaled exterior error"}, ImageSize -> {w1, 0.52 w1},
    Epilog -> {Text[Style[lab <> "  n = " <> ToString[n], fs], Scaled[{0.03, 0.08}], {-1, 0}]}];
  Export[FileNameJoin[{prlDir, "prl_convergence.pdf"}], Column[{pan[3, "(a)", None], pan[5, "(b)", "\[ScriptL]\[Prime]\:2098\:2090\:2093"],
      LineLegend[{colRoute["product"], colRoute["H0"], colRoute["H1"], colRoute["H1 + tail"], Directive[Black, Dashed]}, {"product route (R = \:00bd)", "H0 (no split)", "H1 (k = 0)", "H1 + \[CapitalLambda] tail", "(c) truncation"},
       LabelStyle -> fs, LegendLayout -> {"Row", 3}]}, Spacings -> 0.3, Alignment -> Center]];
  Export[FileNameJoin[{prlDir, "prl_convergence.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], N[#[[4]], 6]} & /@ convRows, {"n", "route", "lpmax", "max_l_scaled_error"}], "CSV"]];
Print["    prl_convergence done [", tick[], " s]"];
csSummary[];
