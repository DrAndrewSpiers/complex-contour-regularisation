#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 8: hybrid split of SR that removes the shadow cross-term.  Stage 6 puncture PhiP = q/rho + c q rho^n, R = 1/2.
     SR = 2Q[PhiP, T_k Phi_res] + 2Q[PhiP, Phi_res - T_k Phi_res],   T_k = 3D Cartesian Taylor polynomial of Phi_res about x0.
   Part 1: closed form (finite Gaunt couplings of [PhiP]_L with the polynomial's modes, l' <= k), continued with K~P (Theta->H)
           times the entire polynomial modes: no shadow cross-term (Im T = 0 on the axis).
   Part 2: product formed on the real line from the real-line residual modes, continued by the Cauchy integral
           s~(z) = (1/(Pi i)) Int_I s(t)/(t - z) dt, I = [4, 6], split at r0, Gauss-Legendre on each side with s at the nodes.
           For z above a piece the cubic Taylor polynomial of s about Re z (from a Chebyshev
           interpolant of s on that side, evaluated on the real interval only) is subtracted and integrated analytically.
           The contour source is then represented by a Chebyshev series in a (NA nodes) for NDSolve.
   H0: Cauchy-continue the whole real SR source (no split).  T_k coefficients from mode sums at r0 (as a self-force code would).
   Exterior error = (c) (real-line l'-truncation, as Stage 7) + Cauchy-step error delta = (contour solve - real-line solve) over the tube.
   Run:  wolframscript -file stage8_hybrid_split.wls                                                                 *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
s7 = Import[FileNameJoin[{dir, "stage7_lprime_optimum.wls"}], "Text"];
s7 = StringTake[s7, {1, StringPosition[s7, "(* ================= main loop"][[1, 1]] - 1}];
s7 = StringReplace[s7, "dir = DirectoryName[$InputFileName];" -> "", 1];
s7 = StringReplace[s7, "Print[\"=== Stage 7" -> "Hold[\"x", 1];
ToExpression[s7];                                           (* Stage 6 + 7 definitions *)
csResetChecks[];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
SetCSParameters[casePars["Product"]]; aL = tubeL[]; aR = tubeR[];        (* R = 1/2 *)
kList = {0, 1, 2}; variants = Join[{"H0"}, ("H1 k=" <> ToString[#]) & /@ kList];
cI = {4, 6}; NC = 60; NA = 24;                                         (* Cauchy interval; Chebyshev degree in t (per side) and per a-piece *)
Print["=== Stage 8: hybrid split; n = ", nS, ", l'max = ", lpList, ", l = ", lOut, ", variants ", variants, " ==="];

(* ---------------- Taylor polynomial T_k of Phi_res about x0 from mode sums at r0 ---------------- *)
Y0[l_] := Sqrt[(2 l + 1)/(4 Pi)];
taylorCoeffs[Q_, n_, c_, M_] := Module[{f0, f1, f2, z = zVar},
   (* real-line residual modes at r0: phi, phi', phi'' (C^(n+1) there, in-branch formulas) *)
   f0 = Sum[Y0[l] evR[resR[n, c, l][[1]], r0], {l, 0, M}];
   f1 = Sum[Y0[l] evR[D[resR[n, c, l][[1]], z], r0], {l, 0, M}];
   f2 = Sum[Y0[l] evR[D[resR[n, c, l][[1]], {z, 2}], r0], {l, 0, M}];
   (* A0 = Phi, A1 = d_z Phi, A2 = d_z^2 Phi, B2 = d_x^2 Phi = (1/r) Phi_r + (1/r^2) Phi_thth at the axis *)
   {f0, f1, f2, Sum[Y0[l] (evR[D[resR[n, c, l][[1]], z], r0]/r0 - l (l + 1)/2 evR[resR[n, c, l][[1]], r0]/r0^2), {l, 0, M}]}];
taylorExact[] := With[{c0 = csPar["c0"], c1 = csPar["c1"], c2 = csPar["c2"]}, N[{c0 + c1 r0 + c2 r0^2, c1 + 2 c2 r0, 2 c2, -c2}, 40]];
(* modes of T_k (polynomial in r, entire): T = A0 + A1 (r mu - r0) + A2/2 (r mu - r0)^2 + B2/2 r^2 (1 - mu^2) *)
tModes[k_, {A0_, A1_, A2_, B2_}] := Module[{mu, cl, poly},
   poly = A0 + If[k >= 1, A1 (zVar mu - r0), 0] + If[k >= 2, A2/2 (zVar mu - r0)^2 + B2/2 zVar^2 (1 - mu^2), 0];
   cl = CoefficientList[Expand[poly], mu];
   Table[Expand[Sum[cl[[m + 1]] muPowerMode[lp, m], {m, 0, Length[cl] - 1}]], {lp, 0, k}]];

(* ---------------- real SR source pieces at a point, from cached mode values (Stage 7 modeVals) ---------------- *)
srTrunc[Q_, l_, M_, mv_, s_] := Module[{a = mv["a"], da = mv["da"], v = mv["v"], dv = mv["dv"]},
   2 Sum[Module[{Ls = Range[Abs[l - lp], l + lp]},
      If[Q === "Product", Sum[gt[l, lp, L] a[[L + 1]], {L, Ls}] v[[lp + 1]],
       Sum[gt[l, lp, L] da[[L + 1]], {L, Ls}] dv[[lp + 1]] + Sum[gt[l, lp, L] angF[l, lp, L] a[[L + 1]], {L, Ls}] v[[lp + 1]]/s^2]], {lp, 0, M}]];
part1Real[Q_, l_, k_, tm_, mv_, s_] := Module[{a = mv["a"], da = mv["da"]},
   2 Sum[Module[{Ls = Range[Abs[l - lp], l + lp], tv = evR[tm[[lp + 1]], s], dtv = evR[D[tm[[lp + 1]], zVar], s]},
      If[Q === "Product", Sum[gt[l, lp, L] a[[L + 1]], {L, Ls}] tv,
       Sum[gt[l, lp, L] da[[L + 1]], {L, Ls}] dtv + Sum[gt[l, lp, L] angF[l, lp, L] a[[L + 1]], {L, Ls}] tv/s^2]], {lp, 0, k}]];
(* continued Part 1 on the a-grid: kv = <|lp -> {K~P, d_r K~P, K~P^theta} values at y(a)|> (cached per l) *)
part1Cont[Q_, l_, k_, tm_, kv_, yA_] := 2 Sum[Module[{tv = evR[tm[[lp + 1]], #] & /@ yA},
     If[Q === "Product", kv[lp][[1]] tv, kv[lp][[2]] (evR[D[tm[[lp + 1]], zVar], #] & /@ yA) + kv[lp][[3]] tv/yA^2]], {lp, 0, k}];

(* ---------------- Chebyshev representation and Cauchy continuation ---------------- *)
chebX[a_, b_, nn_] := N[(a + b)/2 + (b - a)/2 Cos[Pi (Range[0, nn] + 1/2)/(nn + 1)], 50];                  (* first kind: interior nodes *)
(* fixed 40-digit arithmetic: significance tracking is far too pessimistic in these long sums and recurrences *)
chebCoef[f_List] := Block[{$MinPrecision = 40, $MaxPrecision = 40}, Module[{nn = Length[f] - 1, c, ff = SetPrecision[f, 40]},
   c = Table[2/(nn + 1) Sum[ff[[k + 1]] Cos[N[j Pi (k + 1/2)/(nn + 1), 40]], {k, 0, nn}], {j, 0, nn}]; c[[1]] = c[[1]]/2; c]];
chebEval[c_, a_, b_, z_] := Block[{$MinPrecision = 40, $MaxPrecision = 40}, Module[{x = SetPrecision[(2 z - a - b)/(b - a), 40], b1 = 0, b2 = 0, tmp},
   Do[tmp = 2 x b1 - b2 + c[[j + 1]]; b2 = b1; b1 = tmp, {j, Length[c] - 1, 1, -1}];
   x b1 - b2 + c[[1]]]];
pieces = {{cI[[1]], r0}, {r0, cI[[2]]}};
gaussP = gl[#[[1]], #[[2]], 40] & /@ pieces;               (* 40 nodes: exact for the degree-(NC-1) quotient *)
chebDeriv[c_] := Module[{nn = Length[c] - 1, d}, d = ConstantArray[0, nn + 1];
   If[nn >= 1, d[[nn]] = 2 nn c[[nn + 1]]; Do[d[[j]] = If[j + 2 <= nn + 1, d[[j + 2]], 0] + 2 j c[[j + 1]], {j, nn - 1, 1, -1}]; d[[1]] = d[[1]]/2];
   Most[d]];
nG = 80; gaussP = gl[#[[1]], #[[2]], nG] & /@ pieces;
(* Int_a^b Sum_i d_i (t - x)^i/(t - z) dt in closed form ((t - x) = (t - z) + w, w = z - x) *)
polyCauchy[d_, x_, z_, {a_, b_}] := Module[{w = z - x, lg = Log[b - z - I 10^-60] - Log[a - z - I 10^-60]},
   Sum[d[[i + 1]] Sum[Binomial[i, k] w^(i - k) If[k == 0, lg, ((b - z)^k - (a - z)^k)/k], {k, 0, i}], {i, 0, Length[d] - 1}]];
(* sG: s at the Gauss nodes of each piece; dco: per piece the Chebyshev coefficient lists of s and its first three derivatives *)
cauchyAt[sG_, dco_, z_] := Module[{x = Re[z], tot = 0},
   Do[Module[{ab = pieces[[k]], tj = gaussP[[k, All, 1]], wj = gaussP[[k, All, 2]]},
      If[ab[[1]] < x < ab[[2]],             (* subtraction for every z above this piece: exact identity, no switching *)
       Module[{dd = Table[chebEval[dco[[k, m + 1]], ab[[1]], ab[[2]], x] (2/(ab[[2]] - ab[[1]]))^m/m!, {m, 0, 3}], P},
        P = Sum[dd[[m + 1]] (tj - x)^m, {m, 0, 3}];
        tot += Total[wj (sG[[k]] - P)/(tj - z)] + polyCauchy[dd, x, z, ab]],
       tot += Total[wj sG[[k]]/(tj - z)]]], {k, 2}];
   tot/(Pi I)];
derivCoefs[vals_] := Module[{c = chebCoef[vals]}, NestList[chebDeriv, c, 3]];
(* piecewise Chebyshev in a on 40 equal pieces: the composition s~(y(a)) has complex-a singularities about 0.1 from the real segment *)
aBreaks := Subdivide[aL, aR, 40];      (* complex-a singularities (preimages of t = r0, 4, 6) lie ~0.1 from the real a-segment *)
aCheb[nn_] := Flatten[Table[chebX[aBreaks[[k]], aBreaks[[k + 1]], nn], {k, Length[aBreaks] - 1}]];
lineNodes = Join[gl[aL, r0, 40], gl[r0, aR, 40]];
(* contour solve from zero data with the source given by its values at the a-Chebyshev nodes; minus the real-line solve *)
deltaAt[l_, srcVals_, sLine_] := Module[{nb = Length[aBreaks] - 1, cs, f, u, g0, g1, s = lineNodes[[All, 1]], w = lineNodes[[All, 2]]},
   cs = chebCoef /@ Partition[srcVals, Length[srcVals]/nb];
   f[a_?NumericQ] := Module[{k = Min[nb, Max[1, LengthWhile[aBreaks, # <= a &]]]}, chebEval[cs[[k]], aBreaks[[k]], aBreaks[[k + 1]], a]];
   u = solveContourSegment[l, f, {0, 0}];
   g0 = (s^2/(2 l + 1)) (aR^l s^(-l - 1) - s^l aR^(-l - 1)); g1 = (s^2/(2 l + 1)) (l aR^(l - 1) s^(-l - 1) + (l + 1) s^l aR^(-l - 2));
   {Re[u[aR]], Re[u'[aR]]} - {Total[w g0 sLine], Total[w g1 sLine]}];

(* ================= main loop ================= *)
rows = {}; tRows = {};      (* rows: {case, n, variant, M, l, total, delta, c} scaled *)
contourSource[Q_, l_, M_, var_, mvCh_, tCh_, mvG_, tc_, kv_, yA_] := Module[{sCh, sG},
   sCh = Table[srTrunc[Q, l, M, mvCh[[k, i]], tCh[[k, i]]], {k, 2}, {i, Length[tCh[[k]]]}];
   sG = Table[srTrunc[Q, l, M, mvG[[k, i]], gaussP[[k, i, 1]]], {k, 2}, {i, Length[gaussP[[k]]]}];
   If[var === "H0", With[{dco = derivCoefs /@ sCh}, cauchyAt[sG, dco, #] & /@ yA],
    Module[{kk = ToExpression[StringTake[var, -1]], tm, p1Ch, p1G, dco},
     tm = tModes[kk, tc];
     p1Ch = Table[part1Real[Q, l, kk, tm, mvCh[[k, i]], tCh[[k, i]]], {k, 2}, {i, Length[tCh[[k]]]}];
     p1G = Table[part1Real[Q, l, kk, tm, mvG[[k, i]], gaussP[[k, i, 1]]], {k, 2}, {i, Length[gaussP[[k]]]}];
     dco = derivCoefs /@ (sCh - p1Ch);
     part1Cont[Q, l, kk, tm, kv, yA] + (cauchyAt[sG - p1G, dco, #] & /@ yA)]]];
Do[
  SetCSParameters[casePars[Q]];
  Do[Module[{c = r0^-(n + 1), tCh = {chebX[cI[[1]], r0, NC], chebX[r0, cI[[2]], NC]}, aN = aCheb[NA], yA, mvCh, mvG, mvL, nodesC = cNodes[12], mvC},
     yA = N[yC[#], 50] & /@ aN;
     mvCh = Map[modeVals[Q, n, c, #] &, tCh, {2}]; mvL = modeVals[Q, n, c, #] & /@ lineNodes[[All, 1]];
     mvG = Map[modeVals[Q, n, c, #] &, gaussP[[All, All, 1]], {2}];
     mvC = modeVals[Q, n, c, #] & /@ nodesC[[All, 1]];
     Print["    ", Q, " n=", n, ": mode values cached [", tick[], " s]"];
     Do[Module[{tc = taylorCoeffs[Q, n, c, M], te = taylorExact[]}, AppendTo[tRows, {Q, n, M, Sequence @@ Abs[tc - te]}]], {M, lpList}];
     Do[
      Module[{scale = Max[Abs[analyticEval[totalSol[Q, l], #] & /@ gCom]], exactC, kv},
       kv = Association@Table[lp -> {evR[KPc[n, c, l, lp], #] & /@ yA, If[Q === "Product", 0, evR[dKPc[n, c, l, lp], #] & /@ yA],
            If[Q === "Product", 0, evR[KthPc[n, c, l, lp], #] & /@ yA]}, {lp, 0, Max[kList]}];
       exactC = Table[Module[{s = nodesC[[k, 1]]}, pwVal[Sum[pieceAnalyticSource[Q, p, l], {p, {"SS", "SR", "RR"}}], s] -
            pwVal[{ssPsrc[Q, n, c, l, "real"], ssPsrc[Q, n, c, l, "real"]}, s]], {k, Length[nodesC]}];
       Do[
        Module[{eLine, sL = Table[srTrunc[Q, l, M, mvL[[j]], lineNodes[[j, 1]]], {j, Length[lineNodes]}], tc = taylorCoeffs[Q, n, c, M]},
         eLine = cResponse[l, nodesC, Table[truncSrc[Q, l, M, mvC[[k]], nodesC[[k, 1]]], {k, Length[nodesC]}] - exactC, gCom];
         Do[Module[{d, bh},
            d = deltaAt[l, contourSource[Q, l, M, var, mvCh, tCh, mvG, tc, kv, yA], sL];
            bh = Join[ConstantArray[0, Length[gIn]], homOnGrid[l, aR, d, gOut]];
            AppendTo[rows, {Q, n, var, M, l, Max[Abs[eLine + bh]]/scale, Max[Abs[bh]]/scale, Max[Abs[eLine]]/scale}]], {var, variants}]],
        {M, lpList}]],
      {l, lOut}];
     Print["    ", Q, " n=", n, " done [", tick[], " s]"]], {n, nS}],
  {Q, {"Product", "GradDot"}}];

(* ================= checks ================= *)
mxr[Q_, n_, v_, M_, k_] := Max[Cases[rows, {Q, n, v, M, __}][[All, k]]];
Do[csCheck[Q <> " n=" <> ToString[n] <> " " <> v <> ": Cauchy-step error (contour - real line) <= 1e-12 scaled, all l'max",
    Max[Table[mxr[Q, n, v, M, 7], {M, lpList}]], 10^-12];
   csCheck[Q <> " n=" <> ToString[n] <> " " <> v <> ": total follows (c): max over l'max of |total - (c)|/(c) <= 1e-6",
    Max[Table[Abs[mxr[Q, n, v, M, 6] - mxr[Q, n, v, M, 8]]/mxr[Q, n, v, M, 8], {M, lpList}]], 10^-6],
  {Q, {"Product", "GradDot"}}, {n, nS}, {v, variants}];
(* convergence of the Cauchy step: Gauss nodes 80 -> 120 per side, Chebyshev degree 60 -> 90 (t) and 24 -> 36 per a-piece, one setting per case *)
Do[SetCSParameters[casePars[Q]];
  Module[{n = 3, c = r0^-4, M = 12, l = 2, sL, res},
   sL = Table[srTrunc[Q, l, M, modeVals[Q, n, c, s], s], {s, lineNodes[[All, 1]]}];
   res = Table[Block[{nG = pr[[3]], gaussP = gl[#[[1]], #[[2]], pr[[3]]] & /@ pieces},
      Module[{tCh = {chebX[cI[[1]], r0, pr[[1]]], chebX[r0, cI[[2]], pr[[1]]]}, aN = aCheb[pr[[2]]], yA, mvCh, mvG, kv, tc = taylorCoeffs[Q, n, c, M]},
       yA = N[yC[#], 50] & /@ aN; mvCh = Map[modeVals[Q, n, c, #] &, tCh, {2}]; mvG = Map[modeVals[Q, n, c, #] &, gaussP[[All, All, 1]], {2}];
       kv = Association@Table[lp -> {evR[KPc[n, c, l, lp], #] & /@ yA, If[Q === "Product", 0, evR[dKPc[n, c, l, lp], #] & /@ yA],
            If[Q === "Product", 0, evR[KthPc[n, c, l, lp], #] & /@ yA]}, {lp, 0, 2}];
       Table[deltaAt[l, contourSource[Q, l, M, var, mvCh, tCh, mvG, tc, kv, yA], sL], {var, {"H0", "H1 k=2"}}]]], {pr, {{NC, NA, nG}, {3 NC/2, 3 NA/2, 3 nG/2}}}];
   csCheck[Q <> ": Cauchy step converged (Gauss 80 -> 120 per side, Chebyshev 60 -> 90 in t, 24 -> 36 per a-piece; max |delta change|)", Max[Abs[res[[1]] - res[[2]]]], 10^-15]],
  {Q, {"Product", "GradDot"}}];

(* ================= tables and figure ================= *)
s7rows = Import[FileNameJoin[{resDir, "stage7_errors.csv"}], "CSV"];
prod[Q_, n_, M_] := Max[ToExpression /@ Cases[s7rows, {Q, n, "1/2", M, _, tot_, _, _} :> tot]];
md = StringRiffle[Join[{"| case | n | quantity | " <> StringRiffle[ToString /@ lpList, " | "] <> " |", "|---|---|---|" <> StringRepeat["---|", Length[lpList]]},
    Flatten@Table[Join[{"| " <> Q <> " | " <> ToString[n] <> " | Stage 7 product route, R = 1/2 | " <> StringRiffle[sci[prod[Q, n, #]] & /@ lpList, " | "] <> " |",
        "| " <> Q <> " | " <> ToString[n] <> " | (c) | " <> StringRiffle[sci[mxr[Q, n, "H0", #, 8]] & /@ lpList, " | "] <> " |"},
       Table["| " <> Q <> " | " <> ToString[n] <> " | total " <> v <> " | " <> StringRiffle[sci[mxr[Q, n, v, #, 6]] & /@ lpList, " | "] <> " |", {v, variants}],
       Table["| " <> Q <> " | " <> ToString[n] <> " | Cauchy-step error " <> v <> " | " <> StringRiffle[sci[mxr[Q, n, v, #, 7]] & /@ lpList, " | "] <> " |", {v, variants}]],
      {Q, {"Product", "GradDot"}}, {n, nS}]], "\n"];
mdT = StringRiffle[Join[{"| case | n | l'max | err Phi(x0) | err d_z Phi | err d_z^2 Phi | err d_x^2 Phi |", "|---|---|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]]] <> " | " <> StringRiffle[sci /@ #[[4 ;;]], " | "] <> " |") & /@ tRows], "\n"];
Print["\n", md, "\n\nTaylor coefficients from mode sums at r0 (absolute errors vs exact):\n\n", mdT];
Export[FileNameJoin[{resDir, "stage8_summary.md"}],
  "# Stage 8: hybrid split (no shadow cross-term)\n\nScaled exterior errors, max over l = 0..6 (grid and scale as Stage 7), R = 1/2, Cauchy data at r = 1.\n\n" <> md <>
   "\n\n## Taylor coefficients of Phi_res at x0 from mode sums (absolute errors)\n\n" <> mdT <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage8_errors.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], sci[#[[6]]], sci[#[[7]]], sci[#[[8]]]} & /@ rows, {"case", "n", "variant", "lpmax", "l", "total", "cauchy_step", "c"}], "CSV"];
Export[FileNameJoin[{resDir, "stage8_taylor.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]], sci[#[[7]]]} & /@ tRows, {"case", "n", "lpmax", "err_Phi", "err_dz", "err_dzz", "err_dxx"}], "CSV"];
$fs = {FontFamily -> "Times", FontSize -> 9};
vStyle = <|"H0" -> RGBColor[0.15, 0.6, 0.25], "H1 k=0" -> RGBColor[0.55, 0.3, 0.7], "H1 k=1" -> RGBColor[0.85, 0.45, 0.1], "H1 k=2" -> RGBColor[0.1, 0.35, 0.75]|>;
panels = Table[ListLogPlot[Join[Table[Table[{M, prod[Q, n, M]}, {M, lpList}], {n, nS}], Flatten[Table[Table[{M, mxr[Q, n, v, M, 6]}, {M, lpList}], {n, nS}, {v, variants}], 1]],
    Joined -> True, PlotMarkers -> Join[{None, None}, Flatten[Table[{{"\[FilledCircle]", 5}, {"\[EmptySquare]", 6}, {"\[EmptyDiamond]", 7}, {"\[EmptyUpTriangle]", 7}}, {n, nS}], 1]],
    PlotStyle -> Join[Table[Directive[Red, AbsoluteThickness[1.4], If[n == 3, Dashing[{}], Dashed]], {n, nS}],
      Flatten[Table[Directive[vStyle[v], AbsoluteThickness[1], If[n == 3, Dashing[{}], Dashed]], {n, nS}, {v, variants}], 1]],
    Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]'\!\(\*SubscriptBox[\(\[ThinSpace]\), \(max\)]\)", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 212,
    Epilog -> {Inset[Style[If[Q === "Product", "1/\!\(\*SuperscriptBox[\(\[Rho]\), \(2\)]\)", "1/\!\(\*SuperscriptBox[\(\[Rho]\), \(4\)]\)"], $fs], Scaled[{0.97, 0.97}], {Right, Top}]}], {Q, {"Product", "GradDot"}}];
panels2 = Table[ListLogPlot[Table[Table[{M, mxr[Q, 3, v, M, 7] + 10^-30}, {M, lpList}], {v, variants}], Joined -> True, PlotMarkers -> Automatic,
    PlotStyle -> (Directive[vStyle[#], AbsoluteThickness[1]] & /@ variants), Frame -> True, Axes -> False,
    FrameLabel -> {"\[ScriptL]'\!\(\*SubscriptBox[\(\[ThinSpace]\), \(max\)]\)", "Cauchy-step error (n = 3)"}, BaseStyle -> $fs, ImageSize -> 212], {Q, {"Product", "GradDot"}}];
leg = Column[{LineLegend[Join[{Directive[Red, AbsoluteThickness[1.4]]}, Values[vStyle]], Join[{"product route (Stage 7), R = 1/2"}, {"H0: Cauchy, no split", "H1, k = 0", "H1, k = 1", "H1, k = 2"}], LegendLayout -> "Row", LabelStyle -> $fs],
    LineLegend[{Directive[Black], Directive[Black, Dashed]}, {"n = 3", "n = 5"}, LegendLayout -> "Row", LabelStyle -> $fs]}, Alignment -> Center];
Export[FileNameJoin[{figDir, "fig_hybrid.pdf"}], Legended[Grid[{panels, panels2}, Spacings -> {1, 1}], Placed[leg, Below]]];
Export[FileNameJoin[{figDir, "fig_hybrid.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], #[[5]], N[#[[6]], 6], N[#[[7]], 6], N[#[[8]], 6]} & /@ rows, {"case", "n", "variant", "lpmax", "l", "total", "cauchy_step", "c"}], "CSV"];
Print["\nStage 8 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage8_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
