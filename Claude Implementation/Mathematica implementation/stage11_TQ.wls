#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 11: the T_Q term of Musaeus, Pound & Upton (arXiv:2609.10521, eqs. 22-24, 52, 105-112) in the d_z^2 toy.
   Exact puncture q/rho, PhiR = c0 + c1 z, SR = q d_z^2 (PhiR/rho)  (= the Stage 9 "Dzz" SR).
   (a) toy T_Q := canonical SR - 4D Hadamard finite part (ball-regularised) of the pointwise SR, measured with test functions.
   (b) the Stage 9 mode-level delta coefficients c_l split into the delta content of the 4D finite part and of T_Q;
       the mode-level finite part of the pointwise modes (eq. 105) is the slab-regularised distribution.
   (c) the paper's prescription end to end (mode-level FP via (A2)-(A3) with the Coulomb-mode Green's function, + T_Q),
       against the exact (1/2) d_z^2 chi and the contour.   (d) errors without T_Q and with pointwise modes.
   Run:  wolframscript -file stage11_TQ.wls                                                                           *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]]; Needs["NumericalDifferentialEquationAnalysis`"];
resDir = FileNameJoin[{dir, "results"}];
csResetChecks[];
t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
SetCSParameters[<|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>];
q = 1; r0 = 5; c0 = 1; c1 = 1/10; Phi0 = c0 + c1 r0; wp = 40; pg = 20;
phys = MemberQ[$ScriptCommandLine, "phys"];        (* physical boundary conditions psi' = l psi/r at 1, -(l+1) psi/r at 10, for (c), (d) *)
physFix[sol_, l_] := Module[{B = (analyticEvalD[sol, 1] - l analyticEval[sol, 1])/(2 l + 1), A = -(analyticEvalD[sol, 10] + (l + 1) analyticEval[sol, 10]/10)/((2 l + 1) 10^(l - 1))},
   Function[r, analyticEval[sol, r] + A r^l + B r^(-l - 1)]];
lList = Range[0, 6]; gIn = Subdivide[1, 7/2, 12]; gOut = Subdivide[13/2, 10, 14]; gCom = Join[gIn, gOut];
Y0[l_] := Sqrt[(2 l + 1)/(4 Pi)];
cL[l_] := -4 Pi q Phi0 Y0[l]/r0^2;                                     (* Stage 9 mode-level delta coefficient of SR *)
nint[f_, rng_] := Quiet[NIntegrate[f, rng, WorkingPrecision -> wp, PrecisionGoal -> pg, AccuracyGoal -> Infinity, MaxRecursion -> 40], NIntegrate::izero];
Print["=== Stage 11: T_Q in the d_z^2 toy (exact puncture, l = ", lList, ") ==="];
Off[N::meprec];        (* exact zeros *)
Do[csCheck["l=" <> ToString[l] <> ": c_l = (1/2) delta coefficient of d_z^2 [2 q PhiR/rho]_l (Stage 9) = -4 Pi q PhiR(x0) Y_l(0)/r0^2",
   Abs[N[(1/2) dzzDeltaCoeff[pieceAnalyticSource["Product", "SR", #] &, l] - cL[l], 40]], 10^-35], {l, lList}];
On[N::meprec];

(* ================= (a), (b): test-function pairings ================= *)
Print["\n--- (a), (b): pairings with test functions phi = g(r) P_lt(cos theta) ---"];
r1 = 26/5; w = 4/5; gT[r_] := Exp[-(r - r1)^2/w^2];
(* the regulators: angular weight of the 1/rho^3 term of SR, (3 mu^2 - 1), under ball (rho^B) and slab (|Delta r|^B') regularisation *)
Module[{Bp, ang, mu},
  ang = Simplify[2 Pi Integrate[Abs[mu]^Bp (3 mu^2 - 1), {mu, -1, 1}, Assumptions -> Bp > -1]];
  csCheck["ball (rho/rho0)^B finite part has no pole: Int (3 mu^2 - 1) dOmega = 0", Abs[2 Pi Integrate[3 mu^2 - 1, {mu, -1, 1}]], 0];
  csCheck["slab |Dr|^B': Int |mu|^B' (3 mu^2 - 1) dOmega = 8 Pi B'/((B'+1)(B'+3))", Abs[Simplify[ang - 8 Pi Bp/((Bp + 1) (Bp + 3))]], 0];
  (* F(B, B') ~ phi(x0) q PhiR(x0) eps^(B+B')/(B+B') ang(B'): FP_B' then B -> 0 gives 0 (ball); B = 0 then FP_B' gives ang/B' -> 8 Pi/3 *)
  csCheck["slab - ball = lim_{B'->0} ang/B' = 8 Pi/3 (= (2/3) 4 Pi)", Abs[Limit[ang/Bp, Bp -> 0] - 8 Pi/3], 0]];
abRows = {};
Do[Module[{phiMode, dzzPhi, canon, slab, pv, srHalf, fpw, f3, phi3, rr, mu, rho, ref, kTQ, kPV, tq, inner},
   phiMode = Sqrt[4 Pi/(2 lt + 1)] gT[zVar];
   srHalf[j_] := pieceAnalyticSource["Product", "SR", j]/2;                                   (* [q PhiR/rho]_j in/out *)
   dzzPhi[j_] := dzzExpr[Function[i, If[i == lt, phiMode, 0]], j];
   (* canonical: <d_z^2 (q PhiR/rho), phi> = <q PhiR/rho, d_z^2 phi>, mode by mode (integrable) *)
   canon = Sum[If[dzzPhi[j] === 0, 0, Module[{e = srHalf[j] dzzPhi[j] zVar^2},
       nint[e[[1]], {zVar, 1/100, r0}] + nint[e[[2]], {zVar, r0, 14}]]], {j, Max[0, lt - 2], lt + 2}];
   (* slab = mode-level pairing of the pointwise modes (what eq. 105 gives) *)
   fpw = pieceAnalyticSource["Dzz", "SR", lt];
   slab = nint[fpw[[1]] phiMode zVar^2, {zVar, 1/100, r0}] + nint[fpw[[2]] phiMode zVar^2, {zVar, r0, 14}];
   (* ball-regularised principal value (= Riesz finite part), in coordinates (rho, mu) centred at x0 *)
   rr[rh_, m_] := Sqrt[r0^2 + rh^2 + 2 r0 rh m];
   f3[rh_, m_] := q ((c0 + c1 (r0 + rh m)) (3 m^2 - 1)/rh^3 - 2 c1 m/rh^2);
   phi3[rh_, m_] := gT[rr[rh, m]] LegendreP[lt, (r0 + rh m)/rr[rh, m]];
   (* the mu-integral first: the 1/rho^3 term has zero spherical mean, so rho^2 Int f phi dmu = O(rho) and the rho-integral is ordinary *)
   inner[rh_?NumericQ] := NIntegrate[rh^2 f3[rh, mu] phi3[rh, mu], {mu, -1, 1}, WorkingPrecision -> 45, PrecisionGoal -> 25, MaxRecursion -> 30];
   pv = 2 Pi NIntegrate[inner[rho], {rho, 0, 1/5, 1, 5, 51/5, 20}, WorkingPrecision -> 30, PrecisionGoal -> 15, MaxRecursion -> 30];
   ref = cL[lt] Sqrt[4 Pi/(2 lt + 1)] gT[r0] r0^2;                   (* sum_l c_l phi_l(r0) r0^2 = -4 Pi q PhiR(x0) phi(x0) *)
   tq = canon - pv; kTQ = tq/ref; kPV = (pv - slab)/ref;
   AppendTo[abRows, {lt, N[canon, 20], N[pv, 20], N[slab, 20], N[tq/(-(4 Pi/3) q Phi0 gT[r0]), 20], N[kPV, 20], N[kTQ, 20], N[(slab - canon)/(4 Pi q Phi0 gT[r0]), 20]}];
   csCheck["lt=" <> ToString[lt] <> ": canonical - ball PV = -(4 Pi/3) q PhiR(x0) phi(x0)  (relative)", Abs[tq/(-(4 Pi/3) q Phi0 gT[r0]) - 1], 10^-10];
   csCheck["lt=" <> ToString[lt] <> ": slab (mode-level, pointwise modes) - canonical = +4 Pi q PhiR(x0) phi(x0)  (relative)", Abs[(slab - canon)/(4 Pi q Phi0 gT[r0]) - 1], 10^-15];
   csCheck["lt=" <> ToString[lt] <> ": delta content of the 4D finite part = (2/3) c_l", Abs[kPV - 2/3], 10^-10];
   csCheck["lt=" <> ToString[lt] <> ": delta content of T_Q = (1/3) c_l", Abs[kTQ - 1/3], 10^-10];
   Print["    lt=", lt, ": kappa_PV = ", N[kPV, 12], ", kappa_TQ = ", N[kTQ, 12], "  [", tick[], " s]"]], {lt, lList}];

(* ================= (c), (d): prescriptions end to end ================= *)
Print["\n--- (c), (d): end to end, Cauchy data at r = 1, exterior r in [6.5, 10] ---"];
(* SS pointwise mode near r0: A/Dr^2 + B/Dr + O(log); coefficients from the Laurent series on each side *)
ssAB[l_] := ssAB[l] = Module[{e = ssSourceModeExpr["Dzz", l, zVar, "real"], u, sp, sm, cp, cm},
   sp = Normal[Series[(e /. zVar -> r0 + u) /. Log[x_^2] :> 2 Log[Abs[x]] /. Abs[u] -> u, {u, 0, 0}]];
   sm = Normal[Series[(e /. zVar -> r0 - u) /. Log[x_^2] :> 2 Log[Abs[x]] /. Abs[-u] -> u /. Abs[u] -> u, {u, 0, 0}]];
   cp = {Coefficient[sp, u, -2], Coefficient[sp, u, -1]}; cm = {Coefficient[sm, u, -2], -Coefficient[sm, u, -1]};
   {N[cp, 50], Max[Abs[N[cp - cm, 50]]], FreeQ[cp, u]}];
(* finite part over [r0-a, r0+a] of wt(s) (A/Dr^2 + B/Dr), eqs. (A2)-(A3) with wt in place of G *)
fpAB[wt_, s_, {A_, B_}, a_] := Module[{li},
   li[k_] := nint[Log[Abs[s - r0]/a] D[wt, {s, k}], {s, r0 - a, r0, r0 + a}];
   -(A/a) ((wt /. s -> r0 + a) + (wt /. s -> r0 - a)) - A li[2] - B li[1]];
yT[a_] := a + I (1 - ((a - r0)/(1/2))^2)^6;
cNodes = Module[{br = Subdivide[r0 - 1/2, r0 + 1/2, 40], a}, Flatten[Table[({N[yT[#[[1]]], 60], N[D[yT[a], a] /. a -> #[[1]], 60], #[[2]]} & /@
       GaussianQuadratureWeights[30, br[[k]], br[[k + 1]], 60]), {k, Length[br] - 1}], 1]];
srContExpr[l_] := (1/2) dzzExpr[Function[j, If[j < 0, 0, 2 q Sum[KCouplingContinued[j, lp, zVar, "ThetaToH"] phiRAnalytic[lp][zVar], {lp, 0, 1}]]], l];
cRows = {};
Do[Module[{chi, scale, data, hom, wts, s, ss, AB, srP, rr, mom, momC, psi, ex, vars, errs, kap, AB2, tail},
   chi = Sum[pieceAnalyticSolution["Dzz", p, l], {p, {"SS", "SR", "RR"}}];
   scale = If[phys, Max[Abs[physFix[chi, l] /@ gCom]], Max[Abs[analyticEval[chi, #] & /@ gCom]]];
   data = {analyticEval[chi, 1], analyticEvalD[chi, 1]};
   hom = Module[{A, B}, (A #^l + B #^(-l - 1)) & /. First@Solve[{A + B == data[[1]], l A - (l + 1) B == data[[2]]}, {A, B}]];
   wts = {s^(1 - l), s^(l + 2)};                                      (* G0(r, s) = (r^l s^(1-l) - r^(-l-1) s^(l+2))/(2l+1), s < r *)
   ss = ssSourceModeExpr["Dzz", l, s, "real"] /. Log[x_^2] :> 2 Log[Abs[x]];
   AB2 = ssAB[l]; AB = AB2[[1]];
   csCheck["l=" <> ToString[l] <> ": SS pointwise mode = A/Dr^2 + B/Dr + O(log) with the same A, B on both sides, no log/Dr terms", AB2[[2]] + If[AB2[[3]], 0, 1], 10^-40];
   srP = pieceAnalyticSource["Dzz", "SR", l] /. zVar -> s; rr = pieceAnalyticSource["Dzz", "RR", l][[1]] /. zVar -> s;
   (* moments of the pointwise SS (mode-level FP, cutoff a), SR (pointwise), RR over [1, 6] *)
   mom[a_] := Table[Module[{w8 = wts[[k]]},
       nint[w8 ss, {s, 1, r0 - a}] + nint[w8 ss, {s, r0 + a, 6}] + nint[w8 (ss - AB[[1]]/(s - r0)^2 - AB[[2]]/(s - r0)), {s, r0 - a, r0, r0 + a}] + fpAB[w8, s, AB, a] +
        nint[w8 srP[[1]], {s, 1, r0}] + nint[w8 srP[[2]], {s, r0, 6}] + nint[w8 rr, {s, 1, 6}]], {k, 2}];
   (* contour: Re Int_C wt (SS~ + SR~) dz over the tube + real-line pieces *)
   momC = Table[Module[{w8 = wts[[k]], fC = Function @@ {{s}, wts[[k]] ((ssSourceModeExpr["Dzz", l, zVar, "continued"] + srContExpr[l]) /. zVar -> s)}},
       Re[Total[(fC[#[[1]]] #[[2]] #[[3]]) & /@ cNodes]] + nint[w8 (ss + srP[[1]]), {s, 1, r0 - 1/2}] + nint[w8 (ss + srP[[2]]), {s, r0 + 1/2, 6}] + nint[w8 rr, {s, 1, 6}]], {k, 2}];
   (* exterior solution from the moments over [1, 6], the smooth source on [6, r], and a delta(r - r0) of strength gam *)
   tail = Association@Table[r -> Table[nint[wts[[k]] (ss + srP[[2]] + rr), {s, 6, r}], {k, 2}], {r, gOut}];
   psi[m_, gam_, r_] := Module[{mm = m + gam (wts /. s -> r0) + tail[r]},
      hom[r] + (r^l mm[[1]] - r^(-l - 1) mm[[2]])/(2 l + 1)];
   ex = analyticEval[chi, #] & /@ gOut;
   If[phys, Module[{smooth, t610, lo},
      (* free-space Green's function: psi(r) = -(r^(-l-1) Int_1^r s^(l+2) S + r^l Int_r^10 s^(1-l) S)/(2l+1); moments over [1, 6] as above *)
      smooth[k_, a_, b_] := If[b <= a, 0, nint[wts[[k]] (ss + If[a >= r0, srP[[2]], srP[[1]]] + rr), {s, a, b}]];
      t610 = Table[smooth[k, 6, 10], {k, 2}];
      psi[m_, gam_, r_] := Module[{mm = m + gam (wts /. s -> r0)},
         -(If[r > 6, r^(-l - 1) (mm[[2]] + smooth[2, 6, r]) + r^l smooth[1, r, 10],
             r^(-l - 1) smooth[2, 1, r] + r^l (mm[[1]] - smooth[1, 1, r] + t610[[1]])])/(2 l + 1)];
      ex = physFix[chi, l] /@ gCom]];
   kap = Select[abRows, #[[1]] == l &][[1]];                          (* measured kappa_PV, kappa_TQ for this l *)
   vars = {{"contour", momC, 0}, {"Pf modes + T_Q (paper, eq. 24)", mom[1/2], (kap[[6]] + kap[[7]]) cL[l]}, {"Pf modes + T_Q, a = 1/4", mom[1/4], (kap[[6]] + kap[[7]]) cL[l]},
     {"Pf modes, no T_Q", mom[1/2], kap[[6]] cL[l]}, {"pointwise modes (eq. 105) + T_Q", mom[1/2], kap[[7]] cL[l]}, {"pointwise modes, no delta", mom[1/2], 0}};
   errs = Table[Max[Abs[(psi[v[[2]], v[[3]], #] & /@ If[phys, gCom, gOut]) - ex]]/scale, {v, vars}];
   AppendTo[cRows, Prepend[errs, l]];
   csCheck["l=" <> ToString[l] <> ": contour = exact (1/2) d_z^2 chi (scaled)", errs[[1]], 10^-15];
   csCheck["l=" <> ToString[l] <> ": paper prescription (mode FP of SS via (A2)-(A3) + Pf-mode delta + T_Q) = exact (scaled)", errs[[2]], 10^-9];
   csCheck["l=" <> ToString[l] <> ": paper prescription independent of the cutoff a (scaled)", errs[[3]], 10^-9];
   csCheck["l=" <> ToString[l] <> ": omitting T_Q misses (1/3) c_l G (relative to the response)", Abs[errs[[4]]/(Abs[cL[l]] If[phys, Max[Abs[(r0^2/(2 l + 1)) If[# < r0, #^l/r0^(l + 1), r0^l/#^(l + 1)]] & /@ gCom],
       Max[Abs[(r0^(1 - l) #^l - r0^(l + 2) #^(-l - 1))/(2 l + 1)] & /@ gOut]]/(3 scale)) - 1], 10^-8];
   Print["    l=", l, ": ", StringRiffle[sci /@ errs, ", "], "  [", tick[], " s]"]], {l, lList}];
varNames = {"contour", "Pf modes + T_Q (paper, eq. 24)", "same, a = 1/4", "Pf modes, no T_Q", "pointwise modes (eq. 105) + T_Q", "pointwise modes, no delta"};

(* ================= tables ================= *)
md1 = StringRiffle[Join[{"| l_test | canonical | ball PV (4D Pf) | slab (mode-level, pointwise) | T_Q / (-(4 Pi/3) q PhiR phi(x0)) | kappa_PV | kappa_TQ | (slab - canonical)/(4 Pi q PhiR phi(x0)) |", "|---|---|---|---|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> StringRiffle[ToString[NumberForm[#, 12]] & /@ #[[2 ;;]], " | "] <> " |") & /@ abRows], "\n"];
md2 = StringRiffle[Join[{"| l | c_l (Stage 9) | 4D Pf part (2/3) | T_Q part (1/3) |", "|---|---|---|---|"},
    ("| " <> ToString[#[[1]]] <> " | " <> ToString[NumberForm[N[cL[#[[1]]]], 8]] <> " | " <> ToString[NumberForm[N[#[[6]] cL[#[[1]]]], 8]] <> " | " <> ToString[NumberForm[N[#[[7]] cL[#[[1]]]], 8]] <> " |") & /@ abRows], "\n"];
md3 = StringRiffle[Join[{"| l | " <> StringRiffle[varNames, " | "] <> " |", "|---|" <> StringRepeat["---|", Length[varNames]]},
    ("| " <> ToString[#[[1]]] <> " | " <> StringRiffle[sci /@ #[[2 ;;]], " | "] <> " |") & /@ cRows], "\n"];
Print["\n(a) pairings:\n", md1, "\n\n(b) decomposition of c_l:\n", md2, "\n\n(c), (d) scaled exterior errors against (1/2) d_z^2 chi:\n", md3];
Export[FileNameJoin[{resDir, If[phys, "stage11phys_summary.md", "stage11_summary.md"]}], "# Stage 11: T_Q in the d_z^2 toy (exact puncture)\n\n## (a) test-function pairings\n\n" <> md1 <>
   "\n\n## (b) mode-level delta coefficients\n\n" <> md2 <> "\n\n## (c), (d) scaled exterior errors\n\n" <> md3 <> "\n", "Text"];
Export[FileNameJoin[{resDir, If[phys, "stage11phys_pairings.csv", "stage11_pairings.csv"]}], Prepend[{#[[1]], Sequence @@ (N[#, 16] & /@ #[[2 ;;]])} & /@ abRows, {"l_test", "canonical", "ball_PV", "slab", "TQ_ratio", "kappa_PV", "kappa_TQ", "slab_minus_canonical_ratio"}], "CSV"];
Export[FileNameJoin[{resDir, If[phys, "stage11phys_prescriptions.csv", "stage11_prescriptions.csv"]}], Prepend[{#[[1]], Sequence @@ (sci /@ #[[2 ;;]])} & /@ cRows, Prepend[varNames, "l"]], "CSV"];
Print["\nStage 11 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, If[phys, "stage11phys_checks.csv", "stage11_checks.csv"]}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
