#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 12b: the hybrid-convergence study of Stage 12 (B) (and the route comparison of Stage 9, PRL Fig. 5) with physical
   boundary conditions, psi' = l psi/r at r = 1 and psi' = -(l+1) psi/r at r = 10, instead of Cauchy data at r = 1.
   With these conditions the solution for a source on [1, 10] is the free-space Green's-function integral,
       psi(r) = Int_1^10 G(r, s) S(s) ds,   G(r, s) = -(s^2/(2l+1)) r_<^l / r_>^(l+1),
   so every exterior error is computed as that quadrature of its error source (Lemma 1 + linearity; the same decomposition as
   Stages 9 and 12, with the Cauchy-data Green's function replaced by G):
     (c) = H1 (k = 0): truncated SR + RR pointwise - exact, plus delta(r - r0) of strength (gP1 [= gSR] + gRR - gEx)
     H0            : as (c) but without the SR delta (gRR - gEx)
     product route : (c) plus the shadow term Re S~(x + i0) - pointwise on the tube (l'max <= 20)
   The exact reference is (1/2) d_z^2 chi plus the homogeneous terms that satisfy the same conditions.
   Run:  wolframscript -file stage12b_physical_bc.wls                                                                    *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
s9 = Import[FileNameJoin[{dir, "stage9_dzz.wls"}], "Text"];
defsA = StringTake[s9, {1, StringPosition[s9, "(* ================= checks of the d_z machinery"][[1, 1]] - 1}];
defsA = StringReplace[StringReplace[defsA, "dir = DirectoryName[$InputFileName];" -> "", 1], "Print[\"=== Stage 9" -> "Hold[\"x", 1];
defsB = StringTake[s9, {StringPosition[s9, "(* ================= TEST 2"][[1, 1]], StringPosition[s9, "rows9 = {};"][[1, 1]] - 1}];
defsB = StringReplace[defsB, "Print[\"\\n--- Test 2" -> "Hold[\"x", 1];
ToExpression[defsA]; ToExpression[defsB];
csResetChecks[]; Off[General::munfl];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
$fs = {FontFamily -> "Times", FontSize -> 9};
Print["=== Stage 12b: hybrid convergence with physical boundary conditions ==="];

nB = {1, 3, 5, 7}; jMax = 8; Mm = 40; lB = Range[0, 6]; Mprod = 20; lpRoutes = {4, 8, 12, 16, 20, 30, 40};
gridPlot = Join[Subdivide[1, 9/2, 14], Subdivide[11/2, 10, 18]]; gAll = Union[gCom, gridPlot];
posCom = Flatten[Position[gAll, #] & /@ gCom]; posOut = Flatten[Position[gAll, #] & /@ gOut];
nodesB = cNodes[12, gAll];                                                                  (* covers [1, 10]; grid points are breakpoints *)
Gphys[l_, r_, s_] := -(s^2/(2 l + 1)) If[s < r, s^l/r^(l + 1), r^l/s^(l + 1)];
GphysV[l_, r_, s_List] := -(s^2/(2 l + 1)) (UnitStep[r - s] s^l/r^(l + 1) + (1 - UnitStep[r - s]) r^l/s^(l + 1));
physResp[l_, nodes_, delta_, grid_] := Table[Total[nodes[[All, 2]] GphysV[l, r, nodes[[All, 1]]] delta], {r, grid}];
Module[{l = 3, rr, ss, z},
  csCheck["G solves L G = 0 off the diagonal, unit jump of d_r G, and meets psi' = l psi/r at 1, psi' = -(l+1) psi/r at 10",
   Abs[Simplify[D[Gphys[l, rr, 2], {rr, 2}] + 2/rr D[Gphys[l, rr, 2], rr] - l (l + 1) Gphys[l, rr, 2]/rr^2 /. rr -> 3]] +
    Abs[(D[-(ss^2/(2 l + 1)) ss^l/rr^(l + 1), rr] - D[-(ss^2/(2 l + 1)) rr^l/ss^(l + 1), rr]) /. rr -> ss /. ss -> 7/2] - 1 +
    Abs[(D[-(ss^2/(2 l + 1)) rr^l/ss^(l + 1), rr] - l (-(ss^2/(2 l + 1)) rr^l/ss^(l + 1))/rr) /. {rr -> 1, ss -> 4}] +
    Abs[(D[-(ss^2/(2 l + 1)) ss^l/rr^(l + 1), rr] + (l + 1) (-(ss^2/(2 l + 1)) ss^l/rr^(l + 1))/rr) /. {rr -> 10, ss -> 4}], 10^-30]];
(* exact reference with the same conditions: (1/2) d_z^2 chi + A r^l + B r^(-l-1) *)
exactPhys[l_, grid_] := Module[{chi = Sum[pieceAnalyticSolution["Dzz", p, l], {p, {"SS", "SR", "RR"}}], A, B},
   B = (analyticEvalD[chi, 1] - l analyticEval[chi, 1])/(2 l + 1);
   A = -(analyticEvalD[chi, 10] + (l + 1) analyticEval[chi, 10]/10)/((2 l + 1) 10^(l - 1));
   (analyticEval[chi, #] + A #^l + B #^(-l - 1)) & /@ grid];
nodeFcum[{A_, V_}, Mx_] := Module[{srC, rrSh, tr, z3 = 0 V[[1]]},
   srC = Table[2 prodD[Sum[gt[j, lp, L] A[[L + 1]], {L, Abs[j - lp], j + lp}], V[[lp + 1]]], {j, 0, jMax}, {lp, 0, Mx}];
   tr[j_, l1_, l2_] := With[{g = gt[j, l1, l2]}, If[g == 0, z3, g prodD[V[[l1 + 1]], V[[l2 + 1]]]]];
   rrSh = Table[Sum[tr[j, m, l2], {l2, Max[0, m - j], m}] + Sum[tr[j, l1, m], {l1, Max[0, m - j], m - 1}], {j, 0, jMax}, {m, 0, Mx}];
   {Accumulate /@ srC, Accumulate /@ rrSh}];
LamU[l_, nn_] := 2^nn (2 l + 1)/((2 l + nn + 1) Pochhammer[l + nn/2 + 3/2, -nn]);

(* ================= errors ================= *)
eB = <||>; eH0 = <||>; eProd = <||>; dB = <||>; scaleP = <||>;
Do[Module[{c = r0^-(n + 1), P, AVin, AVout, cin, cout, a0, ep, tubeN, FtR, FtC},
   P = Table[Module[{cum = nodeFcum[nodeAV[n, c, s, Automatic, jMax + Mm, Mm], Mm]},
       Table[(1/2) dzzArr[cum[[1, All, M + 1]] + cum[[2, All, M + 1]], s][[All, 1]], {M, 0, Mm}]], {s, nodesB[[All, 1]]}];
   AVin = nodeAV[n, c, r0, "in", jMax + Mm, Mm, 4]; AVout = nodeAV[n, c, r0, "out", jMax + Mm, Mm, 4];
   cin = nodeFcum[AVin, Mm]; cout = nodeFcum[AVout, Mm];
   a0 = Accumulate[Table[Y0[l] evR[resR[n, c, l][[1]], r0], {l, 0, Mm}]];
   ep = exPairF[n, c];
   If[MemberQ[{3, 5}, n],
    tubeN = Join[gl[aL, r0, 40], gl[r0, aR, 40]];
    FtR = Table[nodeF[nodeAV[n, c, s, Automatic, jMax + Mprod, Mprod], Select[lpRoutes, # <= Mprod &], False], {s, tubeN[[All, 1]]}];
    FtC = Table[nodeF[nodeAVc[n, c, s, jMax + Mprod, Mprod], Select[lpRoutes, # <= Mprod &], False], {s, tubeN[[All, 1]]}]];
   Print["    n=", n, ": nodes cached [", tick[], " s]"];
   Do[Module[{exE = Table[(1/2) dzzExpr[Function[j, If[j < 0, 0, ep[j][[b]]]], l], {b, 2}], exPt, gEx, g0, scale},
      exPt = pwVal[exE, #] & /@ nodesB[[All, 1]]; gEx = (1/2) N[dzzDeltaCoeff[ep, l], 40]; g0 = Gphys[l, #, r0] & /@ gAll;
      scale = Max[Abs[exactPhys[l, gCom]]]; scaleP[{n, l}] = scale;
      Module[{res = Table[Module[{gRR = (1/2) dzzDeltaNum[l, cin[[2, # + 1, M + 1, 2]] &, cout[[2, # + 1, M + 1, 2]] &], tp = tPoly9[0, {a0[[M + 1]], 0, 0, 0, 0, 0}], gP1, gSR, base},
            gP1 = (1/2) dzzDeltaNum[l, p1F[AVin[[1]], tp, 0, r0][[# + 1, 2]] &, p1F[AVout[[1]], tp, 0, r0][[# + 1, 2]] &];
            gSR = (1/2) dzzDeltaNum[l, cin[[1, # + 1, M + 1, 2]] &, cout[[1, # + 1, M + 1, 2]] &];
            base = physResp[l, nodesB, P[[All, M + 1, l + 1]] - exPt, gAll];
            {base + (gP1 + gRR - gEx) g0, base + (gRR - gEx) g0, base + (gSR + gRR - gEx) g0}], {M, 0, Mm}]},
       eB[{n, l}] = res[[All, 1]]/scale; eH0[{n, l}] = res[[All, 2]]/scale;
       csCheck["n=" <> ToString[n] <> " l=" <> ToString[l] <> ": H1 (k = 0) = mode truncation (c), max relative difference over l'max", Max[Abs[res[[All, 1]] - res[[All, 3]]]]/Max[Abs[res[[All, 3]]]], 10^-25];
       (* product route: add the shadow term on the tube, Re(continued source)(x + i0) - pointwise *)
       If[MemberQ[{3, 5}, n], Do[Module[{s = tubeN[[All, 1]], w = tubeN[[All, 2]], fb},
           fb = Table[Re[pt[l, FtC[[i]][M][[1]], s[[i]]]] - pt[l, FtR[[i]][M][[1]], s[[i]]], {i, Length[s]}];
           eProd[{n, l, M}] = (res[[M + 1, 3]] + Table[Total[w GphysV[l, r, s] fb], {r, gAll}])/scale], {M, Select[lpRoutes, # <= Mprod &]}]]];
      dB[{n, l}] = Prepend[Differences[eB[{n, l}]], eB[{n, l}][[1]] - (physResp[l, nodesB, -exPt, gAll] - gEx g0)/scale]], {l, lB}];
   Print["    n=", n, " done [", tick[], " s]"]], {n, nB}];

errM[n_, M_] := Max[Table[Max[Abs[eB[{n, l}][[M + 1, posCom]]]], {l, lB}]];
errH0[n_, M_] := Max[Table[Max[Abs[eH0[{n, l}][[M + 1, posCom]]]], {l, lB}]];
errProd[n_, M_] := Max[Table[Max[Abs[eProd[{n, l, M}][[posCom]]]], {l, lB}]];
slopes = Table[Module[{fitE = Fit[Table[{Log[M], Log[errM[n, M]]}, {M, 20, Mm}], {1, x}, x], lmx, fitD},
    lmx = First@Ordering[Table[-Abs[dB[{n, l}][[Mm + 1, posOut[[1]]]]], {l, lB}], 1] - 1;
    fitD = Fit[Table[{Log[lp], Log[Abs[dB[{n, lmx}][[lp + 1, posOut[[1]]]]]]}, {lp, 25, Mm}], {1, x}, x];
    {n, -Coefficient[fitE, x], -Coefficient[fitD, x], lmx, errM[n, 12], errM[n, 20], errM[n, 40]}], {n, nB}];
Do[csCheck["n=" <> ToString[s[[1]]] <> ": error slope = n + 2 within 0.3 (l'max 20..40)", Abs[s[[2]] - (s[[1]] + 2)], 3/10], {s, slopes}];
Do[csCheck["n=" <> ToString[s[[1]]] <> ": per-mode power = n + 3 within 0.3 (r = 6.5, l' 25..40)", Abs[s[[3]] - (s[[1]] + 3)], 3/10], {s, slopes}];
accel[n_, l_, M_, basis_] := Module[{pw = n + 3, bf, X, Y, coef, tl, wn = If[basis === "Lambda1", 3, 6]},
   bf = Which[basis === "power", Function[lp, {lp^-pw, lp^(-pw - 1)}], basis === "Lambda", Function[lp, {LamU[lp, -pw], LamU[lp, -pw - 2]}], basis === "Lambda1", Function[lp, {LamU[lp, -pw]}]];
   X = N[bf /@ Range[M - wn + 1, M], 30]; Y = dB[{n, l}][[M - wn + 2 ;; M + 1]];
   coef = PseudoInverse[X] . Y;
   tl = If[basis === "power", N[{HurwitzZeta[pw, M + 1], HurwitzZeta[pw + 1, M + 1]}, 30], -N[Sum[bf[lp], {lp, 0, M}], 30]];
   eB[{n, l}][[M + 1]] + tl . coef];
accM[n_, M_, basis_] := Max[Table[Max[Abs[accel[n, l, M, basis][[posCom]]]], {l, lB}]];
accRows = Flatten[Table[{n, M, errM[n, M], accM[n, M, "power"], accM[n, M, "Lambda"], accM[n, M, "Lambda1"]}, {n, nB}, {M, 8, Mm}], 1];
Print["    slopes and acceleration [", tick[], " s]"];

(* ================= Cauchy step of the contour solve, propagated with the physical conditions ================= *)
(* contour - distributional line at aR, (du, du'), from zero data at aL, is the response to the tube-source error;
   with the physical conditions it becomes -m1 r^l/(2l+1) inside and -m2 r^(-l-1)/(2l+1) outside, m1 = aR^-l ((l+1) du + aR du'), m2 = -aR^(l+1) (l du - aR du') *)
cauchyRows = Flatten[Table[Module[{n = 3, c = r0^-4, LA, tc, tCh, aN, yAc, srAt, sCh, sG, sL, AVin, AVout, FIn, FOut, gSR, lineAR, tp, p1At, p1Ch, p1G, src, gP1, d, pred, du, m1, m2, eph},
     LA = jMax + M; tc = taylor9[n, c, M];
     srAt[s_] := pt[l, nodeF[nodeAV[n, c, s, Automatic, LA, M], {M}, False][M][[1]], s];
     AVin = nodeAV[n, c, r0, "in", LA, M]; AVout = nodeAV[n, c, r0, "out", LA, M];
     FIn = nodeF[AVin, {M}, False][M][[1]]; FOut = nodeF[AVout, {M}, False][M][[1]];
     gSR = (1/2) dzzDeltaNum[l, FIn[[# + 1, 2]] &, FOut[[# + 1, 2]] &]; lineAR = gSR {G0[l, aR, r0], dG0[l, aR, r0]};
     tCh = {chebX[cI[[1]], r0, NC], chebX[r0, cI[[2]], NC]}; aN = aCheb[NA]; yAc = N[yC[#], 50] & /@ aN;
     sCh = Map[srAt, tCh, {2}]; sG = Map[srAt, gaussP[[All, All, 1]], {2}]; sL = srAt /@ lineNodes[[All, 1]];
     If[route === "H0", src = cauchyAt[sG, derivCoefs /@ sCh, #] & /@ yAc; gP1 = 0,
      tp = tPoly9[0, tc]; p1At[s_] := pt[l, p1F[nodeAV[n, c, s, Automatic, LA, 0][[1]], tp, 0, s], s];
      p1Ch = Map[p1At, tCh, {2}]; p1G = Map[p1At, gaussP[[All, All, 1]], {2}];
      src = (evR[(1/2) dzzExpr[Function[j, If[j < 0, 0, 2 KPc[n, c, j, 0] tp[[1, 1]]]], l], #] & /@ yAc) + (cauchyAt[sG - p1G, derivCoefs /@ (sCh - p1Ch), #] & /@ yAc);
      gP1 = (1/2) dzzDeltaNum[l, p1F[AVin[[1]], tp, 0, r0][[# + 1, 2]] &, p1F[AVout[[1]], tp, 0, r0][[# + 1, 2]] &]];
     d = deltaAt[l, src, sL] - lineAR; pred = -(gSR - gP1) {G0[l, aR, r0], dG0[l, aR, r0]};
     du = d - pred; m1 = aR^-l ((l + 1) du[[1]] + aR du[[2]]); m2 = -aR^(l + 1) (l du[[1]] - aR du[[2]]);
     eph = Table[If[r < aL, -m1 r^l, -m2 r^(-l - 1)]/(2 l + 1), {r, gCom}];
     {route, n, M, l, Max[Abs[eph]]/scaleP[{n, l}], Max[Abs[eph]]/Max[Abs[If[route === "H0", eH0, eB][{n, l}][[M + 1, posCom]] scaleP[{n, l}]]]}],
    {route, {"H1 k=0", "H0"}}, {M, {8, 20}}, {l, {2, 6}}], 2];
csCheck["Cauchy step (contour - line, physical conditions), scaled, max over cases", Max[cauchyRows[[All, 5]]], 10^-14];
Print["    Cauchy step [", tick[], " s]"];

(* ================= tables, csv ================= *)
md1 = StringRiffle[Join[{"| n | slope (l'max 20-40) | per-mode power (l' 25-40) | error at l'max = 12 | 20 | 40 |", "|---|---|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> ToString[NumberForm[#[[2]], 3]] <> " | " <> ToString[NumberForm[#[[3]], 3]] <> " | " <> StringRiffle[sci /@ #[[5 ;; 7]], " | "] <> " |") & /@ slopes], "\n"];
md2 = StringRiffle[Join[{"| n | l'max | plain | + power (2 terms) | + Lambda (2 terms) | + Lambda (1 term) | gain power | gain Lambda 2 | gain Lambda 1 |", "|---|---|---|---|---|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;; 6]], " | "] <> " | " <> StringRiffle[ToString[NumberForm[#, 3]] & /@ (#[[3]]/#[[4 ;; 6]]), " | "] <> " |") & /@
     Select[accRows, MemberQ[{12, 20, 30}, #[[2]]] &]], "\n"];
md3 = StringRiffle[Join[{"| n | l'max | max error, r <= 3.5 | r in [3.5, 4.5] | r in [5.5, 6.5] | r in [6.5, 10] |", "|---|---|---|---|---|---|"},
    Flatten@Table[Module[{prof = Table[Max[Table[Abs[eB[{n, l}][[M + 1, k]]], {l, lB}]], {k, Length[gAll]}], regMax},
       regMax[{lo_, hi_}] := Max[Table[If[lo <= gAll[[k]] <= hi, prof[[k]], 0], {k, Length[gAll]}]];
       "| " <> ToString[n] <> " | " <> ToString[M] <> " | " <> StringRiffle[sci[regMax[#]] & /@ {{1, 7/2}, {7/2, 9/2}, {11/2, 13/2}, {13/2, 10}}, " | "] <> " |"], {n, {3, 5}}, {M, {12, 20}}]], "\n"];
routeRows = Flatten[Table[Join[If[M <= Mprod, {{n, "product", M, errProd[n, M]}}, {}], {{n, "H0", M, errH0[n, M]}, {n, "H1", M, errM[n, M]}, {n, "(c)", M, errM[n, M]}}], {n, {3, 5}}, {M, lpRoutes}], 2];
routeRows = Join[routeRows, Flatten[Table[If[M >= 12, {{n, "H1 + tail", M, Cases[accRows, {n, M, _, _, _, a_} :> a][[1]]}}, {}], {n, {3, 5}}, {M, 8, Mm}], 2]];
md4 = StringRiffle[Join[{"| n | route | " <> StringRiffle[ToString /@ lpRoutes, " | "] <> " |", "|---|---|" <> StringRepeat["---|", Length[lpRoutes]]},
    Flatten@Table["| " <> ToString[n] <> " | " <> rt <> " | " <> StringRiffle[Function[Mv, If[# === {}, "-", sci[#[[1]]]] &[Cases[routeRows, {n, rt, Mv, e_} :> e]]] /@ lpRoutes, " | "] <> " |",
      {n, {3, 5}}, {rt, {"product", "H0", "H1", "H1 + tail"}}]], "\n"];
md5 = StringRiffle[Join[{"| route | n | l'max | l | Cauchy step, scaled | relative to the route's error |", "|---|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]]] <> " | " <> ToString[#[[4]]] <> " | " <> sci[#[[5]]] <> " | " <> sci[#[[6]]] <> " |") & /@ cauchyRows], "\n"];
Print["\n", md1, "\n\n", md2, "\n\n", md3, "\n\n", md4, "\n\n", md5];
Export[FileNameJoin[{resDir, "stage12b_summary.md"}], "# Stage 12b: hybrid convergence with physical boundary conditions\n\nScaled exterior error (max over l = 0..6 and r in the exterior grid), " <>
   "psi' = l psi/r at r = 1, psi' = -(l+1) psi/r at r = 10.\n\n" <> md1 <> "\n\n## tail acceleration\n\n" <> md2 <> "\n\n## error against r\n\n" <> md3 <>
   "\n\n## routes (PRL Fig. 5)\n\n" <> md4 <> "\n\n## Cauchy step of the contour solve\n\n" <> md5 <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage12b_convergence.csv"}], Prepend[Flatten[Table[{n, M, sci[errM[n, M]]}, {n, nB}, {M, 0, Mm}], 1], {"n", "lpmax", "max_err"}], "CSV"];
Export[FileNameJoin[{resDir, "stage12b_acceleration.csv"}], Prepend[{#[[1]], #[[2]], sci[#[[3]]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]]} & /@ accRows, {"n", "lpmax", "plain", "power_tail", "Lambda_tail", "Lambda1_tail"}], "CSV"];
Export[FileNameJoin[{resDir, "stage12b_permode.csv"}], Prepend[Flatten[Table[{n, l, lp, sci[dB[{n, l}][[lp + 1, posOut[[1]]]]]}, {n, nB}, {l, lB}, {lp, 0, Mm}], 2], {"n", "l", "lp", "d_lp_at_r6.5"}], "CSV"];
Export[FileNameJoin[{resDir, "stage12b_profile.csv"}], Prepend[Flatten[Table[{n, M, N[gAll[[k]], 6], sci[Max[Table[Abs[eB[{n, l}][[M + 1, k]]], {l, lB}]]]}, {n, nB}, {M, {12, 20, 40}}, {k, Length[gAll]}], 2], {"n", "lpmax", "r", "max_l_err"}], "CSV"];
Export[FileNameJoin[{resDir, "stage12b_routes.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], N[#[[4]], 8]} & /@ routeRows, {"n", "route", "lpmax", "max_l_scaled_error"}], "CSV"];
Export[FileNameJoin[{resDir, "stage12b_cauchy_step.csv"}], Prepend[{#[[1]], #[[2]], #[[3]], #[[4]], sci[#[[5]]], sci[#[[6]]]} & /@ cauchyRows, {"route", "n", "lpmax", "l", "scaled", "relative_to_route_error"}], "CSV"];

(* ================= figure (same four panels as Stage 12 (B)) ================= *)
Module[{cols = <|1 -> RGBColor[0.55, 0.3, 0.7], 3 -> RGBColor[0.1, 0.35, 0.75], 5 -> RGBColor[0.15, 0.6, 0.25], 7 -> RGBColor[0.85, 0.45, 0.1]|>, pa, pb, pc, pd, fl},
  fl[s_, n_] := Module[{a = Log[errM[n, 30]] + s Log[30]}, Table[{M, Exp[a - s Log[M]]}, {M, {15, Mm}}]];
  pa = ListLogLogPlot[Join[Table[Table[{M, errM[n, M]}, {M, 2, Mm}], {n, nB}], Table[fl[Select[slopes, #[[1]] == n &][[1, 2]], n], {n, nB}]], Joined -> True,
    PlotStyle -> Join[Table[Directive[cols[n], AbsoluteThickness[1.3]], {n, nB}], Table[Directive[Black, Dashed, AbsoluteThickness[0.6]], {n, nB}]],
    Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]\[Prime]\:2098\:2090\:2093", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
    Epilog -> Table[Inset[Style["n=" <> ToString[s[[1]]] <> ": slope " <> ToString[NumberForm[s[[2]], 3]], $fs, cols[s[[1]]]], Scaled[{0.03, 0.03 + 0.08 (4 - Position[nB, s[[1]]][[1, 1]])}], {Left, Bottom}], {s, slopes}]];
  pb = ListLogLogPlot[Table[Table[{lp, Abs[dB[{s[[1]], s[[4]]}][[lp + 1, posOut[[1]]]]]}, {lp, 2, Mm}], {s, slopes}], Joined -> True,
    PlotStyle -> Table[Directive[cols[n], AbsoluteThickness[1.3]], {n, nB}], Frame -> True, Axes -> False,
    FrameLabel -> {"\[ScriptL]\[Prime]", "|per-mode contribution| at r = 6.5"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
    Epilog -> Table[Inset[Style["n=" <> ToString[s[[1]]] <> ": power " <> ToString[NumberForm[s[[3]], 3]], $fs, cols[s[[1]]]], Scaled[{0.03, 0.03 + 0.08 (4 - Position[nB, s[[1]]][[1, 1]])}], {Left, Bottom}], {s, slopes}]];
  pc = ListLogPlot[Flatten[Table[{Table[{gAll[[k]], Max[Table[Abs[eB[{n, l}][[M + 1, k]]], {l, lB}]] + 10^-30}, {k, Select[Range[Length[gAll]], gAll[[#]] <= 9/2 &]}],
        Table[{gAll[[k]], Max[Table[Abs[eB[{n, l}][[M + 1, k]]], {l, lB}]] + 10^-30}, {k, Select[Range[Length[gAll]], gAll[[#]] >= 11/2 &]}]}, {n, {3, 5}}, {M, {12, 20}}], 2],
    Joined -> True, PlotStyle -> Flatten[Table[{#, #} &[Directive[cols[n], If[M == 12, Dashed, Opacity[1]], AbsoluteThickness[1.2]]], {n, {3, 5}}, {M, {12, 20}}]],
    Frame -> True, Axes -> False, FrameLabel -> {"r", "max\:2097 scaled error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> {{1, 10}, All},
    Prolog -> {GrayLevel[0.9], Rectangle[Scaled[{(aL - 1)/9, 0}], Scaled[{(aR - 1)/9, 1}]]},
    Epilog -> {Black, Inset[Style["tube", $fs], Scaled[{(r0 - 1)/9, 0.5}]], Inset[Style["dashed: \[ScriptL]\[Prime]\:2098\:2090\:2093 = 12, solid: 20", $fs], Scaled[{0.98, 0.03}], {Right, Bottom}]}];
  pd = ListLogLogPlot[Flatten[Table[{{#[[2]], #[[3]]} & /@ Select[accRows, #[[1]] == n &], {#[[2]], #[[5]]} & /@ Select[accRows, #[[1]] == n &], {#[[2]], #[[6]]} & /@ Select[accRows, #[[1]] == n &]}, {n, {3, 5}}], 1],
    Joined -> True, PlotStyle -> Flatten[Table[{Directive[cols[n], AbsoluteThickness[1.3]], Directive[cols[n], Dashed, AbsoluteThickness[1.3]], Directive[cols[n], Dotted, AbsoluteThickness[1.3]]}, {n, {3, 5}}], 1],
    Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]\[Prime]\:2098\:2090\:2093", "scaled exterior error"}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
    Epilog -> {Inset[Style["n = 3 (blue), 5 (green): solid plain;", $fs], Scaled[{0.03, 0.11}], {Left, Bottom}], Inset[Style["dashed + \[CapitalLambda] tail, 2 terms; dotted, 1 term", $fs], Scaled[{0.03, 0.03}], {Left, Bottom}]}];
  Export[FileNameJoin[{figDir, "fig_hybrid_convergence.pdf"}], Grid[{{pa, pb}, {pc, pd}}, Spacings -> {1, 1}]];
  Export[FileNameJoin[{figDir, "fig_hybrid_convergence.csv"}], Prepend[{#[[1]], #[[2]], N[#[[3]], 6], N[#[[4]], 6], N[#[[5]], 6], N[#[[6]], 6]} & /@ accRows, {"n", "lpmax", "plain", "power_tail", "Lambda_tail", "Lambda1_tail"}], "CSV"]];
Print["\nStage 12b wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage12b_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
