#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 14b: real-axis radial source samples needed for the exterior SS field, against the contour's 64/96 nodes (Stage 14 (E)).
   Target: the tube contribution to the exterior SS field, psi_l(r) = (r^l M1 - r^(-l-1) M2)/(2l+1), r in [6.5, 10],
           M_k = FP Int_{aL}^{aR} w_k SS_l ds  (w_1 = s^(1-l), w_2 = s^(l+2)) = Re Int_C w_k SS~_l dz   (Stage 10);
           error = max over l <= Lmax of max_r |psi - psi_ref| / max_r |psi_ref|  (Stage 14 (E) measure).
   Real axis: SS_l = A/Dr^2 + B/Dr + C0 ln|Dr| [+ C1 Dr ln|Dr| + C2 Dr^2 ln|Dr|] + remainder; the listed terms are integrated
   against w_k analytically (finite part), the remainder from the samples.
     R1a : uniform grid r0 +- k h (k = 1..K), endpoint-corrected Gregory rule on each side (value at r0 extrapolated), {A,B,C0}
     R1b : same, with C1 and C2 also subtracted (the best rule one can justify on a uniform grid with the local expansion)
     R2  : split at r0, Gauss panels graded geometrically toward r0 (ratio sigma = 0.2), {A,B,C0}
   Contour (C): Stage 14 shared-node angular quadrature at Gauss panels in a.   Noise: 1e-8 relative on the SS samples.
   Samples: SS_l(r) for all l at once from Legendre Q_j(x) at 40-150 digits (checked against the closed form).
   Run:  wolframscript -file stage14b_radial.wls                                                                       *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]]; Needs["NumericalDifferentialEquationAnalysis`"];
resDir = FileNameJoin[{dir, "results"}]; figDir = FileNameJoin[{ParentDirectory[dir], "figures"}];
SetCSParameters[<|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>];
q = 1; r0 = 5; R = 1/2; aL = r0 - R; aR = r0 + R; Lt = 30; LmaxList = {10, 20, 30}; tols = {10^-8, 10^-12};
csResetChecks[]; Off[General::munfl]; t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
$fs = {FontFamily -> "Times", FontSize -> 9};
gOut = Subdivide[13/2, 10, 14]; gIn = Subdivide[1, 7/2, 12]; gCom = Join[gIn, gOut];
phys = MemberQ[$ScriptCommandLine, "phys"];          (* physical boundary conditions: free-space Green's function, both exterior regions *)
tag = If[phys, "stage14bphys_", "stage14b_"];
Print["=== Stage 14b: radial samples for the exterior SS field, real axis (R1, R2) vs contour ==="];

(* ---------------- SS_l(r), l = 0..Lt, from Legendre Q (real r) ---------------- *)
invZD[p_, z_] := (-1)^p p!/z^(p + 1);
dzArr[Fl_, z_] := Module[{nd = Length[Fl[[1]]] - 1, f}, f[j_, m_] := If[j < 0 || j > Length[Fl] - 1, 0, Fl[[j + 1, m + 1]]];
   Table[Table[aZ[j] (f[j - 1, m + 1] - (j - 1) Sum[Binomial[m, i] f[j - 1, i] invZD[m - i, z], {i, 0, m}]) +
      aZ[j + 1] (f[j + 1, m + 1] + (j + 2) Sum[Binomial[m, i] f[j + 1, i] invZD[m - i, z], {i, 0, m}]), {m, 0, nd - 1}], {j, 0, Length[Fl] - 2}]];
ssQ[r_, prec_] := Module[{rr = SetPrecision[r, prec], x, x1, x2, F},
   x = (rr^2 + r0^2)/(2 rr r0); x1 = (rr^2 - r0^2)/(2 rr^2 r0); x2 = r0/rr^3;
   F = Table[Module[{Q = Re[LegendreQ[j, 0, 3, x]], Qm, Q1, Q2, c = q^2 Sqrt[Pi (2 j + 1)]/r0},
       Qm = If[j == 0, 0, Re[LegendreQ[j - 1, 0, 3, x]]];
       Q1 = If[j == 0, 1/(1 - x^2), j (x Q - Qm)/(x^2 - 1)]; Q2 = (2 x Q1 - j (j + 1) Q)/(1 - x^2);
       c {Q/rr, Q1 x1/rr - Q/rr^2, Q2 x1^2/rr + Q1 x2/rr - 2 Q1 x1/rr^2 + 2 Q/rr^3}], {j, 0, Lt + 2}];
   (1/2) dzArr[dzArr[F, rr], rr][[All, 1]]];
Do[csCheck["SS_l from Legendre Q = closed form (l = " <> ToString[l] <> ", r = 4.7, 5.0025)",
   Max[Table[Abs[ssQ[r, 40][[l + 1]] - (ssSourceModeExpr["Dzz", l, zVar, "real"] /. zVar -> SetPrecision[r, 40] /. Log[x_^2] :> 2 Log[Abs[x]])]/Abs[ssQ[r, 40][[l + 1]]], {r, {47/10, 5 + 1/400}}]], 10^-25], {l, {0, 7, 30}}];

(* ---------------- local expansion: A/D^2 + B/D + (C0 + C1 D + C2 D^2 + C3 D^3) ln|D| + D0 + ... + D4 D^4, fitted at 150 digits ---------------- *)
basis[d_] := {1/d^2, 1/d, Log[Abs[d]], d Log[Abs[d]], d^2 Log[Abs[d]], d^3 Log[Abs[d]], 1, d, d^2, d^3, d^4};
coefFit = Module[{ds = Join[-10^-14 Range[8], 10^-14 Range[8]], vals, X},
   vals = Transpose[ssQ[r0 + #, 150] & /@ ds]; X = basis /@ SetPrecision[ds, 150];
   Table[LeastSquares[X, vals[[l + 1]]][[1 ;; 5]], {l, 0, Lt}]];                        (* {A, B, C0, C1, C2} per l *)
csCheck["local coefficients agree with the Stage 12 series (l = 2: A, B, C0)", Max[Abs[N[coefFit[[3, 1 ;; 3]] - {0.07926654595212022026690057530240366065, -0.04755992757127213216014034518144219639, -0.01585330919042404405338011506048073213}, 30]]], 10^-25];
Print["    samples and local coefficients ready [", tick[], " s]"];

(* ---------------- analytic finite-part integrals of the singular terms against w_k ---------------- *)
nint[f_, rng_] := Quiet[NIntegrate[f, rng, WorkingPrecision -> 40, PrecisionGoal -> 30, AccuracyGoal -> 35, MaxRecursion -> 40], {NIntegrate::izero, NIntegrate::slwcon}];
singInt = Table[Module[{s, ws},
     ws = {s^(1 - l), s^(l + 2)};
     Table[Module[{w0 = w /. s -> r0, w1 = D[w, s] /. s -> r0},
       {nint[Together[(w - w0 - w1 (s - r0))/(s - r0)^2], {s, aL, r0, aR}] - 2 w0/R,           (* FP Int w/D^2, FP Int D^-2 = -2/R *)
        nint[Together[(w - w0)/(s - r0)], {s, aL, r0, aR}],                                (* PV Int w/D *)
        nint[w Log[Abs[s - r0]], {s, aL, r0, aR}], nint[w (s - r0) Log[Abs[s - r0]], {s, aL, r0, aR}], nint[w (s - r0)^2 Log[Abs[s - r0]], {s, aL, r0, aR}]}], {w, ws}]], {l, 0, Lt}];
Print["    singular-term integrals [", tick[], " s]"];
masks = <|"ABC" -> {1, 1, 1, 0, 0}, "ABC+" -> {1, 1, 1, 1, 1}|>;
singVal[l_, mk_, d_] := (mk coefFit[[l + 1]]) . {1/d^2, 1/d, Log[Abs[d]], d Log[Abs[d]], d^2 Log[Abs[d]]};
singMom[l_, mk_] := Table[(mk coefFit[[l + 1]]) . singInt[[l + 1, k]], {k, 2}];
psiOf[l_, {M1_, M2_}] := If[phys, Table[If[r < r0, -M1 r^l, -M2 r^(-l - 1)]/(2 l + 1), {r, gCom}], (#^l M1 - #^(-l - 1) M2)/(2 l + 1) & /@ gOut];
errOf[psis_, ref_, L_] := Max[Table[Max[Abs[N[psis[[l + 1]]] - ref[[l + 1]]]]/Max[Abs[ref[[l + 1]]]], {l, 0, L}]];

(* ---------------- contour reference and contour curve (Stage 14 machinery) ---------------- *)
Srho[z_, rho_List] := With[{d = (z^2 - r0^2 - rho^2)/(2 r0)}, (1/2) q^2 (-2/rho^4 + 8 d^2/rho^6)];
muOf[z_, rho_] := (z^2 + r0^2 - rho^2)/(2 z r0);
legAll[L_, mu_] := Module[{P = {ConstantArray[1., Length[mu]], mu}}, Do[AppendTo[P, ((2 l + 1) mu P[[l + 1]] - l P[[l]])/(l + 1)], {l, 1, L - 1}]; P];
gq[n_] := gq[n] = Transpose[N[GaussianQuadratureWeights[n, 0, 1, 30]]];
nrm = N[Table[Sqrt[Pi (2 l + 1)], {l, 0, Lt}]];
pathNodes[z_, br_, nq_] := Module[{via = If[br === "in", -I/2, I/2], e = If[br === "in", r0 - z, z - r0], t, w}, {t, w} = gq[nq];
   {Join[z + r0 + (via - z - r0) t, via + (e - via) t], Join[(via - z - r0) w, (e - via) w]}];
modesC[z_, br_, nq_, L_] := Module[{rho, w, s, mu}, {rho, w} = pathNodes[z, br, nq]; s = Srho[z, rho]; mu = muOf[z, rho];
   nrm[[;; L + 1]] (legAll[L, mu] . (w s (-rho/(z r0))))];
yA[a_] := a + I (1 - ((a - r0)/R)^2)^6;
aNodes[P_, m_] := Module[{br = Subdivide[aL, aR, P], a}, Flatten[Table[({N[yA[#[[1]]]], N[D[yA[a], a] /. a -> #[[1]]], #[[2]]} & /@ N[GaussianQuadratureWeights[m, br[[k]], br[[k + 1]], 30]]), {k, P}], 1]];
contourPsi[nodes_, noise_: 0, seed_: 1] := Module[{z = nodes[[All, 1]], dzw = nodes[[All, 2]] nodes[[All, 3]], src},
   src = Transpose[modesC[#, "out", 128, Lt] & /@ z];
   If[noise != 0, SeedRandom[seed]; src = src (1 + noise RandomVariate[NormalDistribution[], Dimensions[src]])];
   Table[psiOf[l, {Re[Total[z^(1 - l) src[[l + 1]] dzw]], Re[Total[z^(l + 2) src[[l + 1]] dzw]]}], {l, 0, Lt}]];
psiRef = contourPsi[aNodes[80, 30]];
Print["    contour reference [", tick[], " s]"];

(* ---------------- R1: uniform grid, Gregory rule on each side ---------------- *)
gam = {1/12, 1/24, 19/720, 3/160, 863/60480, 275/24192};
gregoryW[K_] := Module[{w = ConstantArray[1, K + 1]}, w[[1]] = 1/2; w[[K + 1]] = 1/2;
   Do[Do[w[[i + 1]] -= gam[[j]] (-1)^j (-1)^(j - i) Binomial[j, i]; w[[K + 1 - i]] -= gam[[j]] (-1)^i Binomial[j, i], {i, 0, j}], {j, Length[gam]}]; w];
csCheck["Gregory weights (6 differences) integrate x^m, m <= 7, exactly on [0, 40]", Max[Table[Abs[gregoryW[40] . Table[If[k == 0 && m == 0, 1, k^m], {k, 0, 40}] - 40^(m + 1)/(m + 1)]/40^(m + 1), {m, 0, 7}]], 0];
Kmax = 3200; hmin = R/Kmax;
gridL = Table[r0 - k hmin, {k, Kmax}]; gridR = Table[r0 + k hmin, {k, Kmax}];
ssGridL = ssQ[#, 40] & /@ gridL; ssGridR = ssQ[#, 40] & /@ gridR;
Print["    uniform-grid samples (", 2 Kmax, " points) [", tick[], " s]"];
r1Psi[K_, mk_, noise_: 0, seed_: 1] := Module[{st = Kmax/K, h = R/K, W = gregoryW[K], xi},
   If[noise != 0, SeedRandom[seed]; xi = RandomVariate[NormalDistribution[], {2, K, Lt + 1}]];
   Table[Module[{mom = {0, 0}},
      Do[Module[{pts = If[side == 1, gridL, gridR][[st ;; Kmax ;; st]], ss = If[side == 1, ssGridL, ssGridR][[st ;; Kmax ;; st, l + 1]], rem, f0, sgn = If[side == 1, -1, 1]},
         If[noise != 0, ss = ss (1 + noise xi[[side, All, l + 1]])];
         rem = ss - (singVal[l, mk, # - r0] & /@ pts);
         f0 = InterpolatingPolynomial[Transpose[{Range[6], rem[[1 ;; 6]]}], 0];                          (* value at r0, extrapolated *)
         mom += Table[h (W . (Prepend[rem, f0] Prepend[pts, r0]^wexp)), {wexp, {1 - l, l + 2}}]], {side, 2}];
      psiOf[l, mom + singMom[l, mk]]], {l, 0, Lt}]];
r1Rows = Flatten[Table[{v, K, 2 K, Table[errOf[r1Psi[K, masks[v]], psiRef, L], {L, LmaxList}]}, {v, {"ABC", "ABC+"}}, {K, {25, 50, 100, 200, 400, 800, 1600, 3200}}], 1];
Print["    R1 [", tick[], " s]"];

(* ---------------- R2: graded Gauss panels toward r0 ---------------- *)
sig = 1/5;
r2Nodes[nl_, m_] := Module[{br = Join[R sig^Range[0, nl], {0}], g},
   Flatten[Table[Module[{gw = N[GaussianQuadratureWeights[m, br[[k + 1]], br[[k]], 60], 60]}, gw], {k, nl + 1}], 1]];        (* {Dr, weight} on (0, R) *)
r2Psi[nl_, m_, mk_, noise_: 0, seed_: 1] := Module[{nd = r2Nodes[nl, m], xi},
   If[noise != 0, SeedRandom[seed]; xi = RandomVariate[NormalDistribution[], {2, Length[nd], Lt + 1}]];
   Module[{ss = Table[ssQ[r0 + sg #, 80] & /@ nd[[All, 1]], {sg, {-1, 1}}]},
     Table[Module[{mom = {0, 0}},
        Do[Module[{pts = r0 + sg nd[[All, 1]], s1 = ss[[If[sg == -1, 1, 2], All, l + 1]], rem},
           If[noise != 0, s1 = s1 (1 + noise xi[[If[sg == -1, 1, 2], All, l + 1]])];
           rem = s1 - (singVal[l, mk, # - r0] & /@ pts);
           mom += Table[Total[nd[[All, 2]] rem pts^wexp], {wexp, {1 - l, l + 2}}]], {sg, {-1, 1}}];
        psiOf[l, mom + singMom[l, mk]]], {l, 0, Lt}]]];
r2Rows = Flatten[Table[{"R2", nl, m, 2 m (nl + 1), Table[errOf[r2Psi[nl, m, masks["ABC"]], psiRef, L], {L, LmaxList}]}, {nl, {2, 4, 6, 8, 10, 12, 14}}, {m, {4, 6, 8, 12, 16}}], 1];
csCheck["R2 (finest) = contour reference (exterior, l <= 30, relative): FP on the real axis = contour", Min[r2Rows[[All, 5, 3]]], 10^-12];
Print["    R2 [", tick[], " s]"];

(* ---------------- contour curve ---------------- *)
cRows = Flatten[Table[{"contour", P, m, P m, Table[errOf[contourPsi[aNodes[P, m]], psiRef, L], {L, LmaxList}]}, {P, {1, 2, 3, 4, 6, 8}}, {m, {6, 8, 12, 16, 20, 24}}], 1];
Print["    contour scan [", tick[], " s]"];

(* ---------------- points needed, and noise ---------------- *)
best[rowsF_, k_, tol_] := Module[{ok = Select[rowsF, #[[-1, k]] <= tol &]}, If[ok === {}, Missing["not reached"], First[SortBy[ok, #[[-2]] &]]]];
r1Fmt = {#[[1]], #[[2]], #[[3]], #[[4]]} & /@ r1Rows;                                             (* {variant, K, points, errs} *)
rowsBy = <|"R1a" -> Select[r1Fmt, #[[1]] === "ABC" &], "R1b" -> Select[r1Fmt, #[[1]] === "ABC+" &],
   "R2" -> ({#[[1]], {#[[2]], #[[3]]}, #[[4]], #[[5]]} & /@ r2Rows), "contour" -> ({#[[1]], {#[[2]], #[[3]]}, #[[4]], #[[5]]} & /@ cRows)|>;
need = Flatten[Table[Module[{b = best[rowsBy[mth], k, tol]}, {mth, LmaxList[[k]], tol, If[MissingQ[b], "not reached", b[[3]]], If[MissingQ[b], Min[rowsBy[mth][[All, 4, k]]], b[[4, k]]], If[MissingQ[b], None, b[[2]]]}],
     {mth, {"contour", "R1a", "R1b", "R2"}}, {k, Length[LmaxList]}, {tol, tols}], 2];
noiseRows = Flatten[Table[Module[{b = best[rowsBy[mth], 3, tol], e},
      If[MissingQ[b], {mth, tol, "-", "-"},
       e = Max[Table[errOf[Switch[mth, "contour", contourPsi[aNodes @@ b[[2]], 10^-8, sd], "R1a", r1Psi[b[[2]], masks["ABC"], 10^-8, sd], "R1b", r1Psi[b[[2]], masks["ABC+"], 10^-8, sd],
            "R2", r2Psi[b[[2, 1]], b[[2, 2]], masks["ABC"], 10^-8, sd]], psiRef, Max[LmaxList]], {sd, 3}]];
       {mth, tol, b[[3]], e}]], {mth, {"contour", "R1a", "R1b", "R2"}}, {tol, tols}], 1];
Print["    noise [", tick[], " s]"];
csCheck["contour reaches 1e-12 at Lmax = 30 within the scan", If[NumberQ[Cases[need, {"contour", 30, 10^-12, p_, __} :> p][[1]]], 0, 1], 0];

(* ---------------- tables, csv, figure ---------------- *)
md1 = StringRiffle[Join[{"| method | Lmax | target | radial samples | error reached | configuration |", "|---|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> sci[N[#[[3]]]] <> " | " <> ToString[#[[4]]] <> " | " <> sci[#[[5]]] <> " | " <> ToString[#[[6]]] <> " |") & /@ need], "\n"];
md2 = StringRiffle[Join[{"| method | configuration chosen for (Lmax = 30) | radial samples | error with 1e-8 relative noise on SS (max of 3 seeds) |", "|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> sci[N[#[[2]]]] <> " | " <> ToString[#[[3]]] <> " | " <> If[NumberQ[#[[4]]], sci[#[[4]]], "-"] <> " |") & /@ noiseRows], "\n"];
Print["\n", md1, "\n\n", md2];
Export[FileNameJoin[{resDir, tag <> "summary.md"}], "# Stage 14b: radial samples for the exterior SS field\n\n" <> md1 <> "\n\n## 1e-8 relative noise on the SS samples\n\n" <> md2 <> "\n", "Text"];
Export[FileNameJoin[{resDir, tag <> "points_vs_error.csv"}], Prepend[Flatten[Table[{mth, ToString[#[[2]]], #[[3]], Sequence @@ (sci /@ #[[4]])} & /@ rowsBy[mth], {mth, Keys[rowsBy]}], 1],
   {"method", "configuration", "radial_samples", "err_L10", "err_L20", "err_L30"}], "CSV"];
Export[FileNameJoin[{resDir, tag <> "needed.csv"}], Prepend[{#[[1]], #[[2]], sci[N[#[[3]]]], #[[4]], sci[#[[5]]], ToString[#[[6]]]} & /@ need, {"method", "Lmax", "target", "radial_samples", "error", "configuration"}], "CSV"];
Export[FileNameJoin[{resDir, tag <> "noise.csv"}], Prepend[{#[[1]], sci[N[#[[2]]]], #[[3]], If[NumberQ[#[[4]]], sci[#[[4]]], "-"]} & /@ noiseRows, {"method", "target", "radial_samples", "err_with_1e-8_noise"}], "CSV"];
Module[{front, cols = <|"contour" -> RGBColor[0.15, 0.6, 0.25], "R1a" -> RGBColor[0.85, 0.45, 0.1], "R1b" -> RGBColor[0.75, 0.1, 0.45], "R2" -> RGBColor[0.1, 0.35, 0.75]|>, p},
  front[mth_] := Module[{pts = SortBy[{#[[3]], #[[4, 3]]} & /@ rowsBy[mth], First], out = {}, cur = Infinity}, Do[If[pt[[2]] < cur, AppendTo[out, pt]; cur = pt[[2]]], {pt, pts}]; out];
  p = ListLogLogPlot[front /@ Keys[cols], Joined -> True, PlotMarkers -> {"\[FilledCircle]", 5}, PlotStyle -> (Directive[#, AbsoluteThickness[1.3]] & /@ Values[cols]),
    Frame -> True, Axes -> False, FrameLabel -> {"radial source samples", "exterior SS error (l \[LessEqual] 30)"}, BaseStyle -> $fs, ImageSize -> 300, PlotRange -> {All, {10^-15, 1}},
    PlotLegends -> Placed[LineLegend[Values[cols], {"contour (Gauss panels in a)", "R1a: uniform grid, Gregory, {A,B,C\:2080} subtracted", "R1b: same, {A,B,C\:2080,C\:2081,C\:2082}", "R2: graded Gauss panels, {A,B,C\:2080}"}, LabelStyle -> $fs], Below],
    Epilog -> {Dashed, Gray, Line[{{1, 10^-8}, {10^4, 10^-8}}], Line[{{1, 10^-12}, {10^4, 10^-12}}]}];
  Export[FileNameJoin[{figDir, If[phys, "fig_radial.pdf", "fig_radial_cauchy.pdf"]}], p];
  Export[FileNameJoin[{figDir, If[phys, "fig_radial.csv", "fig_radial_cauchy.csv"]}], Prepend[Flatten[Table[{mth, #[[1]], N[#[[2]], 6]} & /@ front[mth], {mth, Keys[cols]}], 1], {"method", "radial_samples", "err_L30_pareto"}], "CSV"]];
Print["\nStage 14b wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, tag <> "checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
