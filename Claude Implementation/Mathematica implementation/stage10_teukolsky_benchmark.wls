#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 10: contour vs the mode-level Hadamard finite part of Musaeus, Pound & Upton (arXiv:2609.10521), Sec. V.D.
   Radial operator  [r^2 d^2 - 2 r d - (l-1)(l+2)] psi = S(r),  Green's function eq. (140) (regular at 0 and infinity).
   S = A/Dr^4 + B/Dr^3 + C/Dr^2 + D/Dr + E/(1 + (r/r0)^2) + F Log|Dr| r0^3/(r0^3 + r^3),  Dr = r - r0.
   P1  : eq. (141) with the boundary-term formulas (A2)-(A5), cutoff a.
   P1' : independent Hadamard finite part by Taylor subtraction (checks the transcription of (A2)-(A5)).
   P2  : the l = 2 closed form, eq. (142) (pure singular terms).
   P3  : contour, Re Int_C G(r,z) S~(z) dz over the tube + real-line integral outside; P3b: contour ODE solve.
   Run:  wolframscript -file stage10_teukolsky_benchmark.wls                                                      *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]]; Needs["NumericalDifferentialEquationAnalysis`"];
resDir = FileNameJoin[{dir, "results"}]; figDir = FileNameJoin[{ParentDirectory[dir], "figures"}];
csResetChecks[];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
wp = 40; pg = 25;
r0 = 5; lList = Range[2, 6]; aList = {1/4, 1/2, 1}; rObs = {1, 2, 3, 7, 9, 12, 20};
tubes = {{1/2, 1}, {1/2, 1/4}, {1/2, 2}, {1/4, 1}, {1, 1}};            (* contour {R, amplitude} *)
sets = <|"A" -> {1, 0, 0, 0, 0, 0}, "B" -> {0, 1, 0, 0, 0, 0}, "C" -> {0, 0, 1, 0, 0, 0}, "D" -> {0, 0, 0, 1, 0, 0},
   "mixed" -> {7/10, -13/10, 21/10, 2/5, 3/2, -4/5}|>;
Print["=== Stage 10: Teukolsky-type benchmark; l = ", lList, ", cutoffs a = ", aList, ", tubes {R, amp} = ", tubes, " ==="];

(* ---------------- Green's function, eq. (140) ---------------- *)
gOutB[l_, r_, s_] := -(1/((2 l + 1) s)) (s/r)^(l - 1);           (* s < r *)
gInB[l_, r_, s_] := -(1/((2 l + 1) s)) (r/s)^(l + 2);             (* s > r *)
Lop[l_, f_, r_] := r^2 D[f, {r, 2}] - 2 r D[f, r] - (l - 1) (l + 2) f;
Do[Module[{r, s},
   csCheck["l=" <> ToString[l] <> ": eq. (140) branches solve the homogeneous equation; continuity; r'^2 [d_r G] = 1",
    Abs[Simplify[Lop[l, gOutB[l, r, s], r]]] + Abs[Simplify[Lop[l, gInB[l, r, s], r]]] + Abs[Simplify[(gOutB[l, r, s] - gInB[l, r, s]) /. r -> s]] +
     Abs[Simplify[s^2 (D[gOutB[l, r, s], r] - D[gInB[l, r, s], r]) /. r -> s] - 1], 0]], {l, lList}];
(* r-derivative of order rd of G, as a function of s, valid on one side of r *)
gD[l_, r_, s_, rd_, side_] := Module[{rr}, D[If[side === "out", gOutB[l, rr, s], gInB[l, rr, s]], {rr, rd}] /. rr -> r];
gFull[l_, r_, s_?NumericQ, rd_] := If[s < r, gD[l, r, s, rd, "out"], gD[l, r, s, rd, "in"]];

(* ---------------- source ---------------- *)
sing[cf_, d_] := cf[[1]]/d^4 + cf[[2]]/d^3 + cf[[3]]/d^2 + cf[[4]]/d;
rem[cf_, s_, lg_] := cf[[5]]/(1 + (s/r0)^2) + cf[[6]] lg r0^3/(r0^3 + s^3);
sReal[cf_, s_] := sing[cf, s - r0] + rem[cf, s, Log[Abs[s - r0]]];
sCont[cf_, z_] := sing[cf, z - r0] + rem[cf, z, Log[z - r0]];
nint[f_, rng_] := Quiet[NIntegrate[f, rng, WorkingPrecision -> wp, PrecisionGoal -> pg, AccuracyGoal -> Infinity, MaxRecursion -> 40], NIntegrate::izero];   (* izero: exactly vanishing log integrals *)
(* real-line integral of G S over (0, lo) and (hi, infinity), with r as a breakpoint *)
outer[l_, r_, cf_, rd_, lo_, hi_] := outer[l, r, cf, rd, lo, hi] = Module[{s},
   nint[gFull[l, r, s, rd] sReal[cf, s], Join[{s, 0}, If[r < lo, {r}, {}], {lo}]] + nint[gFull[l, r, s, rd] sReal[cf, s], Join[{s, hi}, If[r > hi, {r}, {}], {Infinity}]]];

(* ---------------- P1: eq. (141) with (A2)-(A5) ---------------- *)
fpPaper[g_, s_, cf_, a_] := Module[{gk, li, A = cf[[1]], B = cf[[2]], C = cf[[3]], Dd = cf[[4]]},
   gk[k_, x_] := D[g, {s, k}] /. s -> x;
   li[k_] := If[D[g, {s, k}] === 0, 0, nint[Log[Abs[s - r0]/a] D[g, {s, k}], {s, r0 - a, r0, r0 + a}]];
   -Dd li[1] +
    (-(C/a) (gk[0, r0 + a] + gk[0, r0 - a]) - C li[2]) +
    (-(B/(2 a^2)) (gk[0, r0 + a] - gk[0, r0 - a]) - (B/(2 a)) (gk[1, r0 + a] + gk[1, r0 - a]) - (1/2) B li[3]) +
    (-(A/(3 a^3)) (gk[0, r0 + a] + gk[0, r0 - a]) - (A/(6 a^2)) (gk[1, r0 + a] - gk[1, r0 - a]) - (A/(6 a)) (gk[2, r0 + a] + gk[2, r0 - a]) - (1/6) A li[4])];
(* P1': Hadamard finite part by subtracting the Taylor polynomial of g at r0 (exact rational algebra, no cancellation) *)
fpTaylor[g_, s_, cf_, a_] := Module[{x, gx, tk},
   gx = g /. s -> r0 + x;
   tk[k_] := SeriesCoefficient[gx, {x, 0, k}];
   Sum[If[cf[[5 - n]] == 0, 0, cf[[5 - n]] (nint[Together[(gx - Sum[tk[k] x^k, {k, 0, n - 1}])/x^n], {x, -a, 0, a}] +
        Sum[If[EvenQ[k - n], tk[k] 2 a^(k - n + 1)/(k - n + 1), 0], {k, 0, n - 1}])], {n, 1, 4}]];
p1[l_, r_, cf_, a_, rd_: 0, fp_: fpPaper] := Module[{s, g},
   g = gD[l, r, s, rd, If[r > r0, "out", "in"]];
   fp[g, s, cf, a] + If[cf[[5]] == 0 && cf[[6]] == 0, 0, nint[g rem[cf, s, Log[Abs[s - r0]]], {s, r0 - a, r0, r0 + a}]] + outer[l, r, cf, rd, r0 - a, r0 + a]];

(* ---------------- P2: eq. (142), l = 2, pure singular terms ---------------- *)
p2[cf_, r_] := With[{A = cf[[1]], B = cf[[2]], C = cf[[3]], Dd = cf[[4]], L1 = Log[Abs[1 - r0/r]], L2 = Log[Abs[(r0 - r)/r0]]},
   Dd/(60 r r0^5) (12 r^5 L1 + r r0 (12 r^3 + 6 r^2 r0 + 4 r r0^2 + 3 r0^3) - 12 r0^5 L2) -
    C/(12 r0^6) (12 r^4 L1 + r0 (12 r^3 + 6 r^2 r0 + 4 r r0^2 + 3 r0^3)) -
    B/(4 r0^7) (r0 (-12 r^4 + 6 r^3 r0 + 2 r^2 r0^2 + r r0^3 + r0^4)/(r - r0) - 12 r^4 L1) -
    A/(12 r0^8) (84 r^4 L1 + r0 (84 r^5 - 126 r^4 r0 + 28 r^3 r0^2 + 7 r^2 r0^3 + 2 r r0^4 + 3 r0^5)/(r - r0)^2)];

(* ---------------- P3: contour integral ---------------- *)
yT[a_, R_, amp_] := a + I amp (1 - ((a - r0)/R)^2)^6;
contourNodes[R_, amp_] := contourNodes[R, amp] = Module[{br = Subdivide[r0 - R, r0 + R, 40], a},
   Flatten[Table[Module[{gw = GaussianQuadratureWeights[30, br[[k]], br[[k + 1]], 60]},
      {N[yT[#[[1]], R, amp], 60], N[D[yT[a, R, amp], a] /. a -> #[[1]], 60], #[[2]]} & /@ gw], {k, Length[br] - 1}], 1]];
p3[l_, r_, cf_, {R_, amp_}, rd_: 0] := Module[{z, g, f},
   g = gD[l, r, z, rd, If[r > r0, "out", "in"]]; f = Function @@ {{z}, g sCont[cf, z]};
   Re[Total[(f[#[[1]]] #[[2]] #[[3]]) & /@ contourNodes[R, amp]]] + outer[l, r, cf, rd, r0 - R, r0 + R]];

(* ================= main comparison ================= *)
rows = {}; p2rows = {}; aRows = {}; tRows = {};
Do[Do[Module[{cf = sets[key], v1, v1T, v1a, v3, sc},
    v1a = Table[Table[p1[l, r, cf, a], {r, rObs}], {a, aList}]; v1 = v1a[[2]]; sc = Max[Abs[v1]];          (* a = 1/2 is the reference *)
    v1T = Table[p1[l, r, cf, 1/2, 0, fpTaylor], {r, rObs}];
    v3 = Table[Table[p3[l, r, cf, tb], {r, rObs}], {tb, tubes}];
    AppendTo[rows, {key, l, sc, Max[Abs[v1T - v1]]/sc, Max[Abs[# - v1] & /@ v1a]/sc, Max[Abs[# - v1] & /@ v3]/sc}];
    Do[AppendTo[aRows, {key, l, aList[[i]], Max[Abs[v1a[[i]] - v1]]/sc}], {i, Length[aList]}];
    Do[AppendTo[tRows, {key, l, tubes[[i]], Max[Abs[v3[[i]] - v1]]/sc}], {i, Length[tubes]}];
    csCheck[key <> " l=" <> ToString[l] <> ": (A2)-(A5) = Taylor-subtracted Hadamard finite part (scaled)", Max[Abs[v1T - v1]]/sc, 10^-20];
    csCheck[key <> " l=" <> ToString[l] <> ": P1 independent of the cutoff a in {1/4, 1/2, 1} (scaled)", Max[Abs[# - v1] & /@ v1a]/sc, 10^-20];
    csCheck[key <> " l=" <> ToString[l] <> ": contour P3 = P1 for all {R, amp} (scaled)", Max[Abs[# - v1] & /@ v3]/sc, 10^-20];
    If[l == 2 && cf[[5]] == 0 && cf[[6]] == 0, Module[{v2 = Table[N[p2[cf, r], wp], {r, rObs}]},
       AppendTo[p2rows, {key, Table[{rObs[[i]], N[v1[[i]], 12], N[v2[[i]], 12]}, {i, Length[rObs]}], Max[Abs[v2 - v1]]/sc}];
       csCheck[key <> " l=2: eq. (142) = P1 (scaled)", Max[Abs[v2 - v1]]/sc, 10^-20]]];
    Print["    ", key, " l=", l, ": scale ", sci[sc], ", P1'-P1 ", sci[rows[[-1, 4]]], ", a-dependence ", sci[rows[[-1, 5]]], ", P3-P1 ", sci[rows[[-1, 6]]], "  [", tick[], " s]"]],
   {key, Keys[sets]}], {l, lList}];

(* ================= P3b: contour ODE solve (mixed source, l = 3, R = 1/2, amp = 1) ================= *)
Quiet[Module[{l = 3, cf = sets["mixed"], R = 1/2, amp = 1, aL, aR, k, psi, a, rr, s1, s2, s3, d1, v7, v9, ref7, ref9, yp, ypp},
  aL = r0 - R; aR = r0 + R; k = (l - 1) (l + 2);
  d1 = {p1[l, 1, cf, 1/2, 0], p1[l, 1, cf, 1/2, 1]};
  s1 = NDSolveValue[{rr^2 psi''[rr] - 2 rr psi'[rr] - k psi[rr] == sReal[cf, rr], psi[1] == d1[[1]], psi'[1] == d1[[2]]}, psi, {rr, 1, aL},
    WorkingPrecision -> 30, PrecisionGoal -> 20, AccuracyGoal -> 35, MaxSteps -> 10^6, InterpolationOrder -> All];
  yp[x_] := D[yT[a, R, amp], a] /. a -> x; ypp[x_] := D[yT[a, R, amp], {a, 2}] /. a -> x;
  s2 = NDSolveValue[{yT[a, R, amp]^2 (psi''[a]/yp[a]^2 - psi'[a] ypp[a]/yp[a]^3) - 2 yT[a, R, amp] psi'[a]/yp[a] - k psi[a] == sCont[cf, yT[a, R, amp]],
     psi[aL] == s1[aL], psi'[aL] == s1'[aL]}, psi, {a, aL, aR}, WorkingPrecision -> 30, PrecisionGoal -> 20, AccuracyGoal -> 35, MaxSteps -> 10^6, InterpolationOrder -> All];
  s3 = NDSolveValue[{rr^2 psi''[rr] - 2 rr psi'[rr] - k psi[rr] == sReal[cf, rr], psi[aR] == Re[s2[aR]], psi'[aR] == Re[s2'[aR]]}, psi, {rr, aR, 9},
    WorkingPrecision -> 30, PrecisionGoal -> 20, AccuracyGoal -> 35, MaxSteps -> 10^6, InterpolationOrder -> All];
  {ref7, ref9} = {p1[l, 7, cf, 1/2], p1[l, 9, cf, 1/2]}; {v7, v9} = {s3[7], s3[9]};
  odeRow = {Abs[v7 - ref7]/Abs[ref7], Abs[v9 - ref9]/Abs[ref9], Abs[Im[s2[aR]]], Abs[s1[3] - p1[l, 3, cf, 1/2]]/Abs[p1[l, 3, cf, 1/2]]};
  csCheck["P3b mixed l=3: contour ODE solve = P1 at r = 7, 9 (relative)", Max[odeRow[[1 ;; 2]]], 10^-14];
  Print["    P3b (contour ODE): relative error at r = 7, 9: ", sci /@ odeRow[[1 ;; 2]], "; Im psi(aR) = ", sci[odeRow[[3]]], "  [", tick[], " s]"]], NDSolveValue::precw];      (* precw: Cauchy data carry ~30 digits *)

(* ================= tables ================= *)
md1 = StringRiffle[Join[{"| source | l | scale max abs(psi) | P1' (Taylor FP) - P1 | max over a of P1(a) - P1(1/2) | max over {R, amp} of P3 - P1 |", "|---|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;;]], " | "] <> " |") & /@ rows], "\n"];
md2 = StringRiffle[Join[{"| source | r | P1 | eq. (142) |", "|---|---|---|---|"},
    Flatten[Table[("| " <> pr[[1]] <> " | " <> ToString[#[[1]]] <> " | " <> ToString[NumberForm[#[[2]], 12]] <> " | " <> ToString[NumberForm[#[[3]], 12]] <> " |") & /@ pr[[2]], {pr, p2rows}]]], "\n"];
md3 = "P3b (contour ODE solve, mixed source, l = 3, R = 1/2, amp = 1): relative error " <> sci[odeRow[[1]]] <> " (r = 7), " <> sci[odeRow[[2]]] <> " (r = 9); Im psi(aR) = " <> sci[odeRow[[3]]] <> ".";
Print["\n", md1, "\n\n", md2, "\n\n", md3];
Export[FileNameJoin[{resDir, "stage10_summary.md"}], "# Stage 10: contour vs mode-level Hadamard finite part (arXiv:2609.10521, Sec. V.D)\n\n" <>
   "Differences are scaled by max over r in " <> ToString[rObs] <> " of abs(P1).\n\n" <> md1 <> "\n\n## eq. (142), l = 2\n\n" <> md2 <> "\n\n" <> md3 <> "\n", "Text"];
Export[FileNameJoin[{resDir, "stage10_comparison.csv"}], Prepend[{#[[1]], #[[2]], Sequence @@ (sci /@ #[[3 ;;]])} & /@ rows, {"source", "l", "scale", "P1T_minus_P1", "max_a_dependence", "max_contour_minus_P1"}], "CSV"];
Export[FileNameJoin[{resDir, "stage10_cutoff_tube.csv"}], Join[{{"source", "l", "a", "diff"}}, {#[[1]], #[[2]], ToString[#[[3]], InputForm], sci[#[[4]]]} & /@ aRows,
    {{"source", "l", "R_amp", "diff"}}, {#[[1]], #[[2]], ToString[#[[3]], InputForm], sci[#[[4]]]} & /@ tRows], "CSV"];
Export[FileNameJoin[{resDir, "stage10_eq142.csv"}], Prepend[Flatten[Table[Prepend[#, pr[[1]]] & /@ pr[[2]], {pr, p2rows}], 1], {"source", "r", "P1", "eq142"}], "CSV"];
Print["\nStage 10 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, "stage10_checks.csv"}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
