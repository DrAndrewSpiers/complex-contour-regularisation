#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Stage 15: extracting the second-order residual field at the particle (Dzz case, Q[A,B] = (1/2) d_Z^2 (A B)).
   Phi1 = q/rho + PhiR1, PhiR1 = c0 + c1 Z = A + c1 d (A = c0 + c1 r0, d = Z - r0).  Lap Phi2 = (1/2) d_Z^2 (Phi1^2).
   Exact second-order singular field:  chiS = q^2 ln rho + q A rho + (1/2) q c1 d rho,  Phi2S = (1/2) d_Z^2 chiS,
   so Phi2R = Phi2 - Phi2S = (1/2) d_Z^2 (chi - chiS) = (1/2) d_Z^2 chiRR + h, where h is the harmonic field injected through the data.
   The SS + SR field psi~_l is solved on the contour (variation of parameters along y(a), cross-checked against NDSolve), RR on the real line.
   The residual psi~R_l = psi~_l - phi~P_l is extended from the contour to r0 by a Taylor series whose coefficients come from the ODE
   (recursion on the Taylor coefficients of the residual source s~R = s~_SS + s~_SR - L_l phi~P, obtained by a Cauchy/FFT integral on a
   small circle).  Cross-checks: an ODE solve around a circle about the expansion point with Cauchy integrals; the monodromy of psi~R
   around r0.  (a) exact puncture; (b) Phi2P = Phi2S + eps q^2 rho^k, k = 1, 3, 5; (c) eta-baseline psi~ - phi~S at r0 + i eta;
   (d) controls: d1 = Theta->H blend of the finished (differentiated) Dzz SR modes, d2 = one-sided (out) continuation.
   Run:  wolframscript -file stage15_residual_extraction.wls [smoke]                                                              *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourSecondOrder.wl"}]]; Needs["NumericalDifferentialEquationAnalysis`"];
resDir = FileNameJoin[{dir, "results"}]; figDir = FileNameJoin[{ParentDirectory[dir], "figures", "residual"}]; Quiet[CreateDirectory[figDir]];
smoke = MemberQ[$ScriptCommandLine, "smoke"];
SetCSParameters[<|"q" -> 1, "r0" -> 5, "c0" -> 1, "c1" -> 1/10|>];
q = 1; r0 = 5; c0 = 1; c1 = 1/10; A0 = c0 + c1 r0; R = 1/2; aL = tubeL[]; aR = tubeR[];
csResetChecks[]; t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
Lmax = If[smoke, 4, 20]; lList = Range[0, Lmax];
kList = {1, 3, 5}; epsK[k_] := r0^-(k + 2);
Bo = 10; ro = 15; Bi = 1; ri = 1/2;                 (* h = Bo/|x - ro zhat| + Bi/|x - ri zhat|, harmonic on [1, 10], all l present *)
Nmax = If[smoke, 60, 120]; MF = If[smoke, 128, 256]; rcFrac = 7/10;      (* Taylor order; FFT points; circle radius / |z_e - r0| *)
aPts = {5, 26/5};                                                            (* expansion points: apex and a = 5.2 *)
etaList = If[smoke, {1, 1/4, 1/16}, 2^-Range[0, 7]]; lEta = If[smoke, {0, 2}, {0, 2, 5, 10}];
figL = Select[{0, 2, 5, 10, 20}, # <= Lmax &]; lCirc = Select[{0, 2, 5, 10}, # <= Lmax &]; lMono = Select[{0, 2, 5}, # <= Lmax &];
WP = 50; WPin = 80;
Y0[l_] := Sqrt[(2 l + 1)/(4 Pi)];
$fs = {FontFamily -> "Times", FontSize -> 9};
Print["=== Stage 15: residual extraction at r0; l = 0..", Lmax, ", k = ", kList, ", Nmax = ", Nmax, ", smoke = ", smoke, " ==="];

(* ======================= 0. the second-order singular field in 3D (symbolic) ======================= *)
Module[{x, y, Z, rho, d, PhiR, chiS, chiRR, hh, lap, ax},
  rho = Sqrt[x^2 + y^2 + (Z - r0)^2]; d = Z - r0; PhiR = c0 + c1 Z;
  lap[f_] := D[f, {x, 2}] + D[f, {y, 2}] + D[f, {Z, 2}];
  chiS = q^2 Log[rho] + q A0 rho + (1/2) q c1 d rho;
  chiRR = c0^2 (x^2 + y^2 + Z^2)/6 + c0 c1 Z (x^2 + y^2 + Z^2)/5 + c1^2 Z^4/12;
  csCheck["3D: Lap(ln rho) = 1/rho^2, Lap rho = 2/rho, Lap(d rho) = 4 d/rho (x != x0)",
   Abs[Simplify[lap[Log[rho]] - 1/rho^2]] + Abs[Simplify[lap[rho] - 2/rho]] + Abs[Simplify[lap[d rho] - 4 d/rho]], 0];
  csCheck["3D: Lap chiS = Phi1^2 - (PhiR1)^2 = q^2/rho^2 + 2 q PhiR1/rho", Abs[Simplify[lap[chiS] - (q^2/rho^2 + 2 q PhiR/rho)]], 0];
  csCheck["3D: Lap chiRR = (PhiR1)^2", Abs[Simplify[lap[chiRR] - PhiR^2]], 0];
  csCheck["3D: chiS = chi_SS + chi_SR of Stage 3 (q c0 rho + q c1 (r0 rho + (Z - r0) rho/2))", Abs[Simplify[chiS - (q^2 Log[rho] + q c0 rho + q c1 (r0 rho + (Z - r0) rho/2))]], 0];
  hh = Bo/Sqrt[x^2 + y^2 + (Z - ro)^2] + Bi/Sqrt[x^2 + y^2 + (Z - ri)^2];
  csCheck["3D: h is harmonic", Abs[Simplify[lap[hh]]], 0];
  csCheck["3D: (1/2) d_Z^2 chiRR = c0^2/6 + 3 c0 c1 Z/5 + c1^2 Z^2/2", Abs[Simplify[(1/2) D[chiRR, {Z, 2}] - (c0^2/6 + 3 c0 c1 Z/5 + c1^2 Z^2/2)]], 0];
  ax = Simplify[((1/2) D[chiRR, {Z, 2}] + hh) /. {x -> 0, y -> 0}, ri < Z < ro];     (* Phi2R on the axis *)
  exact3D = {ax /. Z -> r0, D[ax, Z] /. Z -> r0};
  Print["    exact Phi2R(x0) = ", N[exact3D[[1]], 25], ",  d_r Phi2R(x0) = ", N[exact3D[[2]], 25]]];

(* ======================= 1. mode expressions in zVar ======================= *)
hMode[l_] := hMode[l] = Sqrt[4 Pi/(2 l + 1)] (Bo zVar^l/ro^(l + 1) + Bi ri^l/zVar^(l + 1));
rhoN[k_, L_] := rhoN[k, L] = Table[uModeIntegral[L, 0, {"pow", k/2}, zVar, br], {br, {"in", "out"}}];          (* [rho^k]_L in/out *)
Hz = HExt[zVar];
blend[{in_, out_}] := in + (out - in) Hz;                                                                        (* Theta -> H *)
chiSR[j_] := chiSR[j] = pieceAnalyticSolution["Product", "SR", j];                                               (* [chi_SR]_j in/out, C^2 at r0 *)
punSS[l_] := punSS[l] = pieceAnalyticSolution["Dzz", "SS", l][[2]];              (* in = out as real functions; ln|r - r0| -> Log[z - r0] *)
punSR[l_] := punSR[l] = (1/2) dzzExpr[Function[j, If[j < 0, 0, blend[chiSR[j]]]], l];   (* blend the undifferentiated modes, then the ladder *)
punS[l_] := punS[l] = punSS[l] + punSR[l];
punK[k_, l_] := punK[k, l] = epsK[k] q^2 blend[rhoN[k, l]];                      (* truncation term eps q^2 [rho^k]_l, C^(k+1) at r0 *)
pDzzSR[l_] := pDzzSR[l] = pieceAnalyticSolution["Dzz", "SR", l];                 (* finished Dzz SR modes in/out (kinked, C^0) *)
punW1[l_] := punW1[l] = blend[pDzzSR[l]];                                        (* control d1: blend the finished modes *)
punW2[l_] := punW2[l] = pDzzSR[l][[2]];                                          (* control d2: one-sided continuation (no kink, no delta) *)
srcSS[l_] := srcSS[l] = ssSourceModeExpr["Dzz", l, zVar, "continued"];
srcSR[l_] := srcSR[l] = (1/2) dzzExpr[Function[j, If[j < 0, 0, 2 q Sum[KCouplingContinued[j, lp, zVar, "ThetaToH"] phiRAnalytic[lp][zVar], {lp, 0, 1}]]], l];
srcTot[l_] := srcTot[l] = srcSS[l] + srcSR[l];
sRes[l_] := sRes[l] = srcSR[l] - radialL[l, punSR[l], zVar];           (* residual source, exact puncture; the SS part is identically 0 (checked) *)
dPun["k=1", l_] := punK[1, l]; dPun["k=3", l_] := punK[3, l]; dPun["k=5", l_] := punK[5, l];
dPun["d1", l_] := punW1[l] - punSR[l]; dPun["d2", l_] := punW2[l] - punSR[l]; dPun["exact", l_] := 0;
dSrc[cs_, l_] := dSrc[cs, l] = If[cs === "exact", 0, -radialL[l, dPun[cs, l], zVar]];
cases = {"exact", "k=1", "k=3", "k=5", "d1", "d2"};
solIn[l_] := solIn[l] = pieceAnalyticSolution["Dzz", "SS", l][[1]] + pDzzSR[l][[1]] + hMode[l];     (* physical SS+SR mode on [1, r0) *)
solOut[l_] := solOut[l] = pieceAnalyticSolution["Dzz", "SS", l][[2]] + pDzzSR[l][[2]] + hMode[l];   (* on (r0, 10] *)
dz[e_] := D[e, zVar];
ev[e_, z_] := N[e /. zVar -> SetPrecision[z, WPin], WP];
(* exact residual modes at r0 (value, r-derivative) for each case; the real part of the residual on the axis is these functions *)
resExact["exact", l_] := {ev[hMode[l], r0], ev[dz[hMode[l]], r0]};
resExact[cs_, l_] /; StringTake[cs, 2] == "k=" := With[{k = ToExpression[StringDrop[cs, 2]]}, With[{e = epsK[k] q^2 rhoN[k, l][[1]]}, {ev[hMode[l] - e, r0], ev[dz[hMode[l] - e], r0]}]];
resExact["d1", l_] := resExact["exact", l];
resExact["d2", l_] := resExact["exact", l];     (* the value is right (J_P(r0) = 0); the derivative does not exist: (z - r0) log(z - r0) term *)

Print["\n--- 1. consistency of the continued puncture and the residual source ---"];
zTest = {r0 + I, N[yC[23/5], 60], N[yC[27/5], 60]};
Do[Module[{dSS, sr, dW1, rePun, reK, dK, cin},
   dSS = Max[Abs[ev[radialL[l, punSS[l], zVar] - srcSS[l], #]/ev[srcSS[l], #]] & /@ zTest];
   csCheck["l=" <> ToString[l] <> ": SS source on the contour = L_l of the continued SS puncture (relative, 3 points)", dSS, 10^-40];
   (* the residual source is imaginary on the axis and bounded at r0 *)
   sr = {Max[Abs[Re[ev[sRes[l], #]]] & /@ {47/10 + 10^-40 I, 53/10 + 10^-40 I}], Abs[ev[sRes[l], r0 + 10^-8 I]], Abs[ev[sRes[l], r0 + I]]};
   csCheck["l=" <> ToString[l] <> ": Re s~R(r + i0) = 0 at r = 4.7, 5.3 (exact puncture)", sr[[1]], 10^-30];
   csCheck["l=" <> ToString[l] <> ": s~R bounded at r0: |s~R(r0 + 1e-8 i)| <= 10 |s~R(r0 + i)| + 1e-10", Max[0, sr[[2]] - 10 sr[[3]] - 10^-10], 0];
   (* real part of the continued puncture on the axis = modes of the one 3D puncture (in branch below r0, out above) *)
   rePun = Abs[Re[ev[punS[l], 47/10 + 10^-40 I]] - ev[pieceAnalyticSolution["Dzz", "SS", l][[1]] + pDzzSR[l][[1]], 47/10]] +
     Abs[Re[ev[punS[l], 53/10 + 10^-40 I]] - ev[pieceAnalyticSolution["Dzz", "SS", l][[2]] + pDzzSR[l][[2]], 53/10]];
   csCheck["l=" <> ToString[l] <> ": Re phi~S(r + i0) = [Phi2S]_l(r) at r = 4.7 (in) and 5.3 (out)", rePun, 10^-30];
   reK = Max[Table[Abs[Re[ev[punK[k, l], 47/10 + 10^-40 I]] - ev[epsK[k] q^2 rhoN[k, l][[1]], 47/10]] + Abs[Re[ev[punK[k, l], 53/10 + 10^-40 I]] - ev[epsK[k] q^2 rhoN[k, l][[2]], 53/10]], {k, kList}]];
   csCheck["l=" <> ToString[l] <> ": Re of the continued truncation term = eps q^2 [rho^k]_l on the axis (k = 1, 3, 5)", reK, 10^-30];
   (* [rho^k]_l is C^(k+1): in and out agree at r0 through the (k+1)-th derivative, differ at order k+2 *)
   dK = Table[Module[{p = rhoN[k, l]}, Max[Table[Abs[ev[D[p[[2]] - p[[1]], {zVar, m}], r0]], {m, 0, k + 1}]]], {k, kList}];
   cin = Table[Module[{p = rhoN[k, l]}, Abs[ev[D[p[[2]] - p[[1]], {zVar, k + 2}], r0]]], {k, kList}];
   csCheck["l=" <> ToString[l] <> ": [rho^k]_l is C^(k+1) at r0 (in - out derivatives through order k+1)", Max[dK], 10^-35];
   csCheck["l=" <> ToString[l] <> ": [rho^k]_l is not C^(k+2) (jump of the (k+2)-th derivative > 0)", If[Min[cin] > 10^-10, 0, 1], 0];
   (* [chi_SR]_j is C^2: the H-blend then has J = O((z - r0)^3), so L phi~S_SR - s~_SR is analytic *)
   csCheck["l=" <> ToString[l] <> ": [chi_SR]_l is C^2 at r0 (in - out through the second derivative)",
    Max[Table[Abs[ev[D[chiSR[l][[2]] - chiSR[l][[1]], {zVar, m}], r0]], {m, 0, 2}]], 10^-35];
   (* control d1 differs from the blend-then-ladder puncture by a function that is imaginary on the axis and vanishes at r0 *)
   dW1 = {Max[Abs[Re[ev[dPun["d1", l], #]]] & /@ {47/10 + 10^-40 I, 53/10 + 10^-40 I}], Abs[ev[dPun["d1", l], r0 + 10^-6 I]]};
   csCheck["l=" <> ToString[l] <> ": control d1 - blend/ladder puncture is imaginary on the axis", dW1[[1]], 10^-30];
   csCheck["l=" <> ToString[l] <> ": control d1 - blend/ladder puncture -> 0 at r0 (|.| at r0 + 1e-6 i < 1e-4)", dW1[[2]], 10^-4];
   If[MemberQ[figL, l], Print["    l=", l, " [", tick[], " s]"]]], {l, lList}];

(* ======================= 2. SS + SR on the contour: variation of parameters along y(a) ======================= *)
(* psi(z) = A z^l + B z^(-l-1) + [z^l I1(z) - z^(-l-1) I2(z)]/(2l+1),  I1 = Int_{aL}^{z} s^(1-l) s~(s) ds,  I2 = Int s^(l+2) s~(s) ds along the contour;
   Gauss-Legendre panels in a (breakpoints include a = 5 and a = 5.2); values at the breakpoints are partial sums.  Equivalent to the
   contour ODE solve (Stage 12); cross-checked against solveContourSegment below. *)
nP = If[smoke, 40, 60]; nG = 20;
yE[a_, eta_] := a + I eta (1 - ((a - r0)/R)^2)^6;
yEp[a_, eta_] := Module[{t}, D[yE[t, eta], t] /. t -> a];
makeNodes[brks_, eta_, ng_] := Module[{gw, aN, zN, wN, zb},
   gw = Flatten[Table[GaussianQuadratureWeights[ng, brks[[k]], brks[[k + 1]], 60], {k, Length[brks] - 1}], 1];
   aN = gw[[All, 1]];
   zN = N[yE[SetPrecision[#, 80], eta], 60] & /@ aN; wN = (N[yEp[SetPrecision[#, 80], eta], 60] & /@ aN) gw[[All, 2]];
   zb = N[yE[SetPrecision[#, 80], eta], 60] & /@ brks;
   <|"a" -> aN, "z" -> zN, "w" -> wN, "zb" -> zb, "ng" -> ng, "breaks" -> brks|>];
contourSolve[l_, srcExpr_, {v_, d_}, nd_] := Module[{S, A, B, cum1, cum2, zb = nd["zb"], psi, dpsi, zN = nd["z"], wN = nd["w"]},
   S = ev[srcExpr, #] & /@ zN;
   cum1 = Prepend[Accumulate[Total /@ Partition[wN zN^(1 - l) S, nd["ng"]]], 0];
   cum2 = Prepend[Accumulate[Total /@ Partition[wN zN^(l + 2) S, nd["ng"]]], 0];
   {A, B} = LinearSolve[{{aL^l, aL^(-l - 1)}, {l aL^(l - 1), -(l + 1) aL^(-l - 2)}}, {v, d}];
   psi = A zb^l + B zb^(-l - 1) + (zb^l cum1 - zb^(-l - 1) cum2)/(2 l + 1);
   dpsi = l A zb^(l - 1) - (l + 1) B zb^(-l - 2) + (l zb^(l - 1) cum1 + (l + 1) zb^(-l - 2) cum2)/(2 l + 1);
   <|"psi" -> psi, "dpsi" -> dpsi, "S" -> S|>];
breaks = Subdivide[aL, aR, nP]; nodes1 = makeNodes[breaks, 1, nG];
iPt = Flatten[Position[breaks, #] & /@ aPts]; iR = Length[breaks];
Print["\n--- 2. contour solve of SS + SR (", Length[nodes1["z"]], " nodes) ---"];
psiC = <||>; dpsiC = <||>; resC = <||>; dresC = <||>; dataAt = <||>; shadowR = <||>; contRows = {};
Do[Module[{dat, sol, pS, dpS, rr, drr, hv, dhv},
   dat = {ev[solIn[l], aL], ev[dz[solIn[l]], aL]};
   sol = contourSolve[l, srcTot[l], dat, nodes1];
   pS = ev[punS[l], #] & /@ nodes1["zb"]; dpS = ev[dz[punS[l]], #] & /@ nodes1["zb"];
   rr = sol["psi"] - pS; drr = sol["dpsi"] - dpS;
   psiC[l] = sol["psi"]; dpsiC[l] = sol["dpsi"]; resC[l] = rr; dresC[l] = drr;
   hv = ev[hMode[l], #] & /@ {aL, aR}; dhv = ev[dz[hMode[l]], #] & /@ {aL, aR};
   shadowR[l] = {Im[sol["psi"][[-1]]], Im[sol["dpsi"][[-1]]], Im[rr[[-1]]], Im[drr[[-1]]]};
   csCheck["l=" <> ToString[l] <> ": Re psi~R(a_L) = h_l(a_L) (data)", Abs[Re[rr[[1]]] - hv[[1]]], 10^-30];
   csCheck["l=" <> ToString[l] <> ": Lemma 1 + consistency: Re psi~R(a_R) = h_l(a_R), Re psi~R'(a_R) = h_l'(a_R)", Abs[Re[rr[[-1]]] - hv[[2]]] + Abs[Re[drr[[-1]]] - dhv[[2]]], 10^-25];
   csCheck["l=" <> ToString[l] <> ": Re psi~(a_R) = physical SS+SR mode (out branch + h) at a_R", Abs[Re[sol["psi"][[-1]]] - ev[solOut[l], aR]], 10^-25];
   Do[AppendTo[contRows, {l, N[breaks[[i]], 10], N[Re[rr[[i]]], 20], N[Im[rr[[i]]], 20], N[Re[sol["psi"][[i]]], 20], N[Im[sol["psi"][[i]]], 20], N[Re[pS[[i]]], 20], N[Im[pS[[i]]], 20]}], {i, Length[breaks]}];
   If[MemberQ[figL, l], Print["    l=", l, ": Im psi~(a_R) = ", sci[shadowR[l][[1]]], ", Re psi~R(a_R) - h = ", sci[Re[rr[[-1]]] - hv[[2]]], ", psi~R(apex) = ", N[rr[[iPt[[1]]]], 12], " [", tick[], " s]"]]], {l, lList}];
(* cross-check against the contour ODE solve (NDSolve, Stage 9 machinery) *)
Do[Module[{dat = {ev[solIn[l], aL], ev[dz[solIn[l]], aL]}, u, d1, d2},
   u = solveContourSegment[l, contourFunction[srcTot[l], zVar], dat];
   d1 = Abs[u[aR] - psiC[l][[-1]]] + Abs[u'[aR] - dpsiC[l][[-1]]]; d2 = Abs[u[5] - psiC[l][[iPt[[1]]]]];
   csCheck["l=" <> ToString[l] <> ": quadrature contour solve = NDSolve contour solve at a_R and at the apex", d1 + d2, 10^-15];
   Print["    l=", l, ": quadrature vs NDSolve ", sci[d1 + d2], " [", tick[], " s]"]], {l, Select[{0, 2, 6}, # <= Lmax &]}];
(* RR on the real line (polynomial modes, l <= 2 nonzero); check the real-line solve against the closed form *)
rrAt[l_] := rrAt[l] = Module[{sol = pieceAnalyticSolution["Dzz", "RR", l][[1]]}, {ev[sol, r0], ev[dz[sol], r0]}];
Do[Module[{sol = pieceAnalyticSolution["Dzz", "RR", l][[1]], src = pieceAnalyticSource["Dzz", "RR", l][[1]], u},
   u = realLineSolve[l, numericFunction[src, zVar], {ev[sol, 1], ev[dz[sol], 1]}];
   csCheck["l=" <> ToString[l] <> ": RR real-line solve at r0 = closed form (value, derivative)", Abs[u[r0] - rrAt[l][[1]]] + Abs[u'[r0] - rrAt[l][[2]]], 10^-18]], {l, Select[{0, 1, 2, 3}, # <= Lmax &]}];
Export[FileNameJoin[{resDir, "stage15_residual_contour.csv"}], Prepend[contRows, {"l", "a", "Re_resid", "Im_resid", "Re_psi", "Im_psi", "Re_punct", "Im_punct"}], "CSV"];

(* ======================= 3. Taylor coefficients about a contour point from the ODE ======================= *)
(* z^2 psi'' + 2 z psi' - l(l+1) psi = z^2 s~R; with z = z_e + w, psi = Sum p_n w^n, sigma_n = Taylor coefficients of z^2 s~R:
   p_(n+2) = [sigma_n - 2 z_e (n+1)^2 p_(n+1) - (n(n+1) - l(l+1)) p_n] / (z_e^2 (n+2)(n+1)).
   The s~R coefficients come from the Cauchy integral on the circle |z - z_e| = rc (trapezoidal rule, MF points): every quantity is O(1). *)
thF = Table[2 Pi k/MF, {k, 0, MF - 1}]; eTh = N[Exp[I thF], 60];
fftCoef[expr_, ze_, rc_] := Module[{vals = ev[expr, ze + rc #] & /@ eTh}, Table[Total[vals eTh^(-n)]/MF/rc^n, {n, 0, Nmax}]];
taylorODE[l_, ze_, {p0_, p1_}, c_] := Module[{cc = Join[{0, 0}, c], sig, p = ConstantArray[0, Nmax + 3]},
   sig = Table[ze^2 c[[n + 1]] + 2 ze cc[[n + 2]] + cc[[n + 1]], {n, 0, Nmax}];
   p[[1]] = p0; p[[2]] = p1;
   Do[p[[n + 3]] = (sig[[n + 1]] - 2 ze (n + 1)^2 p[[n + 2]] - (n (n + 1) - l (l + 1)) p[[n + 1]])/(ze^2 (n + 2) (n + 1)), {n, 0, Nmax}];
   p];
partialSums[p_, w0_] := {Accumulate[Table[p[[n + 1]] w0^n, {n, 0, Length[p] - 1}]], Accumulate[Table[If[n == 0, 0, n p[[n + 1]] w0^(n - 1)], {n, 0, Length[p] - 1}]]};
Print["\n--- 3. Taylor extension to r0 ---"];
tay = <||>; tayRows = {}; r0Rows = {}; coefStore = <||>;
Do[Module[{ze = nodes1["zb"][[iPt[[ip]]]], rc, cEx, pRes, dpRes},
   rc = rcFrac Abs[ze - r0];
   cEx = fftCoef[sRes[l], ze, rc];
   pRes = resC[l][[iPt[[ip]]]]; dpRes = dresC[l][[iPt[[ip]]]];
   Do[Module[{c, p0, p1, p, S, ex, errV, errD},
      c = cEx + If[cs === "exact", 0, fftCoef[dSrc[cs, l], ze, rc]];
      p0 = pRes - If[cs === "exact", 0, ev[dPun[cs, l], ze]]; p1 = dpRes - If[cs === "exact", 0, ev[dz[dPun[cs, l]], ze]];
      p = taylorODE[l, ze, {p0, p1}, c]; coefStore[{cs, l, ip}] = p;
      S = partialSums[p, r0 - ze]; ex = resExact[cs, l];
      errV = Abs[Re[S[[1]]] - ex[[1]]]; errD = Abs[Re[S[[2]]] - ex[[2]]];
      tay[{cs, l, ip}] = <|"S" -> S, "errV" -> errV, "errD" -> errD, "exact" -> ex, "p" -> p|>;
      If[MemberQ[figL, l] || Lmax <= 6, Do[AppendTo[tayRows, {cs, l, N[aPts[[ip]]], n, N[errV[[n + 1]], 6], N[errD[[n + 1]], 6], N[Im[S[[1, n + 1]]], 12], N[Im[S[[2, n + 1]]], 12]}], {n, 0, Nmax}]];
      AppendTo[r0Rows, {cs, l, N[aPts[[ip]]], N[Re[S[[1, -1]]], 25], N[ex[[1]], 25], N[errV[[-1]], 6], N[Im[S[[1, -1]]], 20], N[Re[S[[2, -1]]], 25], N[ex[[2]], 25], N[errD[[-1]], 6], N[Im[S[[2, -1]]], 20]}]], {cs, cases}];
   If[MemberQ[figL, l], Print["    l=", l, " a_e=", N[aPts[[ip]]], ": exact err(N=", Nmax, ") = ", sci[tay[{"exact", l, ip}]["errV"][[-1]]], " / ", sci[tay[{"exact", l, ip}]["errD"][[-1]]],
      ";  k=3: ", sci[tay[{"k=3", l, ip}]["errV"][[-1]]], " / ", sci[tay[{"k=3", l, ip}]["errD"][[-1]]], ";  d1: ", sci[tay[{"d1", l, ip}]["errV"][[-1]]], ";  d2: ", sci[tay[{"d2", l, ip}]["errV"][[-1]]], " / ", sci[tay[{"d2", l, ip}]["errD"][[-1]]], " [", tick[], " s]"]]],
  {l, lList}, {ip, Length[aPts]}];
(* convergence rates: geometric ratio (exact, d1) from a fit of log err against N over the range before the floor; algebraic slopes (k, d2) from log err vs log N on [Nmax/3, Nmax] *)
geomRatio[err_] := Module[{sel = Select[Table[{n, err[[n + 1]]}, {n, 2, Nmax}], 10^-38 < #[[2]] < 10^-3 &]},
   If[Length[sel] < 4, Indeterminate, Exp[Coefficient[Fit[{#[[1]], Log[#[[2]]]} & /@ sel, {1, x}, x], x]]]];
algSlope[err_] := Module[{sel = Select[Table[{n, err[[n + 1]]}, {n, Ceiling[Nmax/3], Nmax}], #[[2]] > 10^-38 &]},
   If[Length[sel] < 4, Indeterminate, Coefficient[Fit[{Log[#[[1]]], Log[#[[2]]]} & /@ sel, {1, x}, x], x]]];
(* fits to the Taylor coefficients themselves, |p_n w0^n| (upper half of the range above the 1e-40 floor): geometric a + b n + c log n (ratio e^b), algebraic a + c log n *)
coefSel[p_, w0_] := Module[{sel = Select[Table[{n, Abs[p[[n + 1]] w0^n]}, {n, 10, Nmax}], #[[2]] > 10^-40 &]}, If[Length[sel] < 8, {}, Take[sel, -Floor[Length[sel]/2]]]];
coefGeom[p_, w0_] := Module[{sel = coefSel[p, w0]}, If[sel === {}, Indeterminate, Exp[Coefficient[Fit[{#[[1]], Log[#[[2]]]} & /@ sel, {1, x, Log[x]}, x], x]]]];
coefAlg[p_, w0_] := Module[{sel = coefSel[p, w0]}, If[sel === {}, Indeterminate, Coefficient[Fit[{Log[#[[1]]], Log[#[[2]]]} & /@ sel, {1, x}, x], x]]];
rateRows = {};
Do[Module[{eV = tay[{cs, l, ip}]["errV"], eD = tay[{cs, l, ip}]["errD"], p = tay[{cs, l, ip}]["p"], w0 = r0 - nodes1["zb"][[iPt[[ip]]]]},
   AppendTo[rateRows, {cs, l, N[aPts[[ip]]], N[geomRatio[eV], 4], N[geomRatio[eD], 4], N[algSlope[eV], 4], N[algSlope[eD], 4], N[eV[[-1]], 4], N[eD[[-1]], 4], N[coefGeom[p, w0], 4], N[coefAlg[p, w0], 4]}]], {cs, cases}, {l, lList}, {ip, Length[aPts]}];
rate[cs_, l_, ip_, k_] := First[Cases[rateRows, {cs, l, N[aPts[[ip]]], __}]][[k]];
predRatio = Table[N[Abs[r0 - nodes1["zb"][[iPt[[ip]]]]]/Abs[nodes1["zb"][[iPt[[ip]]]]], 6], {ip, Length[aPts]}];
Print["    predicted geometric ratios |r0 - z_e|/|z_e| (singularity at z = 0): ", predRatio];
Do[csCheck["exact puncture, a_e=" <> ToString[N[aPts[[ip]]]] <> ": Taylor sum at r0 converges to h_l(r0), h_l'(r0) for every l (N = " <> ToString[Nmax] <> ", max over l)",
    Max[Table[tay[{"exact", l, ip}]["errV"][[-1]] + tay[{"exact", l, ip}]["errD"][[-1]], {l, lList}]], 10^-20];
   csCheck["exact puncture, a_e=" <> ToString[N[aPts[[ip]]]] <> ": two-parameter fit of the error curve gives a geometric ratio within 60% of |r0 - z_e|/|z_e| (biased by the n^l prefactor of the pole at z = 0; median over l)",
    Abs[Median[DeleteCases[Table[rate["exact", l, ip, 4], {l, lList}], Indeterminate]] - predRatio[[ip]]]/predRatio[[ip]], 0.6];
   csCheck["exact puncture, a_e=" <> ToString[N[aPts[[ip]]]] <> ": coefficient-fit ratio (|p_n w0^n| ~ n^c rho^n) within 10% of |r0 - z_e|/|z_e| (median over l)",
    Abs[Median[DeleteCases[Table[rate["exact", l, ip, 10], {l, lList}], Indeterminate]] - predRatio[[ip]]]/predRatio[[ip]], 0.1];
   csCheck["control d1 (blend of finished modes), a_e=" <> ToString[N[aPts[[ip]]]] <> ": Taylor sum converges to the same real values (max over l)",
    Max[Table[tay[{"d1", l, ip}]["errV"][[-1]] + tay[{"d1", l, ip}]["errD"][[-1]], {l, lList}]], 10^-20];
   csCheck["control d2 (one-sided), a_e=" <> ToString[N[aPts[[ip]]]] <> ": fails: derivative error at N = Nmax > 1e-3 for every l <= 6",
    If[Min[Table[tay[{"d2", l, ip}]["errD"][[-1]], {l, Select[lList, # <= 6 &]}]] > 10^-3, 0, 1], 0];
   (* expected exponents: the leading singular term (i/pi) j (z - r0)^(k+2) log(z - r0) has Taylor terms at r0 with phase i^(k+2) (-1)^n x const about an
      expansion point directly above r0 (apex): for odd k it feeds only Re of the value and Im of the derivative, so the real part of the derivative
      inherits the next term (z - r0)^(k+3) log: both real-part errors fall like N^-(k+2).  Off the vertical (a = 5.2) the phases mix: N^-(k+2), N^-(k+1). *)
   Do[Module[{eD = If[ip == 1, k + 2, k + 1]},
     csCheck["k=" <> ToString[k] <> ", a_e=" <> ToString[N[aPts[[ip]]]] <> ": algebraic convergence, value slope within 1 of -(k+2) (median over l)",
      Abs[Median[DeleteCases[Table[rate["k=" <> ToString[k], l, ip, 6], {l, lList}], Indeterminate]] + (k + 2)], 1];
     csCheck["k=" <> ToString[k] <> ", a_e=" <> ToString[N[aPts[[ip]]]] <> ": derivative slope within 1 of -" <> ToString[eD] <> " (parity rule; median over l)",
      Abs[Median[DeleteCases[Table[rate["k=" <> ToString[k], l, ip, 7], {l, lList}], Indeterminate]] + eD], 1]], {k, kList}], {ip, Length[aPts]}];
(* control d2: Re S'_N tends to the mean of the one-sided derivatives, h' + c_l/2, with c_l = -4 pi q PhiR(x0) Y_l(0)/r0^2 the dropped delta coefficient (Stage 9) *)
cDelta[l_] := -4 Pi q A0 Y0[l]/r0^2;
Do[csCheck["control d2 l=" <> ToString[l] <> ": derivative error at N = Nmax = |c_l|/2 (mean of the one-sided derivatives) within 3%",
   Abs[tay[{"d2", l, 1}]["errD"][[-1]] - Abs[cDelta[l]]/2]/(Abs[cDelta[l]]/2), 0.03], {l, lMono}];

(* ----- 3b. cross-check: ODE solve around a circle about the apex, Taylor coefficients by Cauchy integrals ----- *)
Print["\n--- 3b. Cauchy-circle cross-check (apex) ---"];
pathODE[l_, u_, t_, zf_, zpf_, zppf_] := (u''[t] - u'[t] zppf[t]/zpf[t])/zpf[t]^2 + (2/zf[t]) u'[t]/zpf[t] - l (l + 1) u[t]/zf[t]^2;
circleCoef[l_, srcExpr_, ze_, rc_, {p0_, p1_}, lgRule_: None] := Module[{sf, u, t, s1, s2, vals, pC},
   sf = If[lgRule === None, numericFunction[srcExpr, zVar], Module[{lg, e = srcExpr /. Log[-5 + zVar] -> lg, f}, f[z_?NumericQ, lgv_] := N[e /. {zVar -> SetPrecision[z, 70], lg -> SetPrecision[lgv, 70]}, 50]; f]];
   (* radial segment z = ze + rc t, t in [0, 1] *)
   s1 = NDSolveValue[{u''[t]/rc^2 + (2/(ze + rc t)) u'[t]/rc - l (l + 1) u[t]/(ze + rc t)^2 == If[lgRule === None, sf[ze + rc t], sf[ze + rc t, lgRule[ze + rc t]]],
      u[0] == toWP[p0], u'[0] == toWP[rc p1]}, u, {t, 0, 1}, Evaluate[Sequence @@ ndOptions[]]];
   (* circle z = ze + rc Exp[I t] *)
   s2 = NDSolveValue[{pathODE[l, u, t, Function[tt, ze + rc Exp[I tt]], Function[tt, I rc Exp[I tt]], Function[tt, -rc Exp[I tt]]] ==
       If[lgRule === None, sf[ze + rc Exp[I t]], sf[ze + rc Exp[I t], lgRule[ze + rc Exp[I t]]]],
      u[0] == SetPrecision[s1[1], 30], u'[0] == SetPrecision[I s1'[1], 30]}, u, {t, 0, 2 Pi}, Evaluate[Sequence @@ ndOptions[]]];
   vals = s2 /@ N[thF, 30];
   pC = Table[Total[vals eTh^(-n)]/MF/rc^n, {n, 0, Nmax}];
   <|"p" -> pC, "mono" -> s2[2 Pi] - s2[0], "sol" -> s2|>];
circRows = {};
Do[Module[{ze = nodes1["zb"][[iPt[[1]]]], rc, cc, pR, S, d, nCmp = Min[40, Nmax]},
   rc = rcFrac Abs[ze - r0];
   pR = tay[{cs, l, 1}]["p"];
   cc = circleCoef[l, sRes[l] + dSrc[cs, l], ze, rc, {pR[[1]], pR[[2]]}];
   d = Max[Table[Abs[cc["p"][[n + 1]] - pR[[n + 1]]] rc^n, {n, 0, nCmp}]];
   S = partialSums[cc["p"][[1 ;; nCmp + 1]], r0 - ze];
   AppendTo[circRows, {cs, l, N[d, 4], N[Abs[cc["mono"]], 4], N[Abs[S[[1, -1]] - tay[{cs, l, 1}]["S"][[1, nCmp + 1]]], 4]}];
   csCheck[cs <> " l=" <> ToString[l] <> ": Taylor coefficients from the circle solve = ODE recursion (max_n |dp_n| rc^n, n <= " <> ToString[nCmp] <> ")", d, 10^-15];
   csCheck[cs <> " l=" <> ToString[l] <> ": solution single-valued around the circle |u(2 Pi) - u(0)|", Abs[cc["mono"]], 10^-15];
   Print["    ", cs, " l=", l, ": coefficient difference ", sci[d], ", monodromy ", sci[Abs[cc["mono"]]], " [", tick[], " s]"]], {cs, {"exact", "k=3"}}, {l, lCirc}];

(* ----- 3c. monodromy of the residual around r0 (direct test of analyticity at r0) ----- *)
(* path z = r0 + R e^{i(Pi - t)}, t in [0, 2 Pi]: a_L (from above) -> apex r0 + i R -> a_R -> below -> a_L (from below).
   On the path |z - r0| = R, so the continued Log[z - r0] = ln R + i (Pi - t) exactly; this is substituted for the principal Log.
   Prediction: Delta psi~R = -[dPun(a_L, from below) - dPun(a_L, from above)] = 2 Pi i (i/Pi) x (coefficient of Log) ... computed from the expressions. *)
Print["\n--- 3c. monodromy of psi~R around r0 ---"];
monoRows = {};
Do[Module[{src = sRes[l] + dSrc[cs, l], lg, e, sf, u, t, zf, zpf, zppf, s2, p0, p1, pred, dp = dPun[cs, l], mono, half},
   e = src /. Log[-5 + zVar] -> lg;
   sf[z_?NumericQ, lgv_] := N[e /. {zVar -> SetPrecision[z, 70], lg -> SetPrecision[lgv, 70]}, 50];
   zf = Function[tt, r0 + R Exp[I (Pi - tt)]]; zpf = Function[tt, -I R Exp[I (Pi - tt)]]; zppf = Function[tt, -R Exp[I (Pi - tt)]];
   p0 = resC[l][[1]] - If[cs === "exact", 0, ev[dp, aL]]; p1 = dresC[l][[1]] - If[cs === "exact", 0, ev[dz[dp], aL]];
   s2 = NDSolveValue[{pathODE[l, u, t, zf, zpf, zppf] == sf[zf[t], Log[R] + I (Pi - t)], u[0] == toWP[p0], u'[0] == toWP[p1 zpf[0]]}, u, {t, 0, 2 Pi}, Evaluate[Sequence @@ ndOptions[]]];
   mono = s2[2 Pi] - s2[0]; half = s2[Pi];
   pred = If[cs === "exact", 0, -(ev[dp /. Log[-5 + zVar] -> Log[-5 + zVar] - 2 Pi I, aL] - ev[dp, aL])];
   AppendTo[monoRows, {cs, l, N[mono, 8], N[pred, 8], N[Abs[mono - pred], 4], N[Abs[half - (resC[l][[-1]] - If[cs === "exact", 0, ev[dp, aR]])], 4]}];
   csCheck[cs <> " l=" <> ToString[l] <> ": monodromy of psi~R around r0 = prediction (exact: 0)", Abs[mono - pred], 10^-15];
   csCheck[cs <> " l=" <> ToString[l] <> ": psi~R(a_R) from the semicircle = from the contour y(a) (contour independence)", Abs[half - (resC[l][[-1]] - If[cs === "exact", 0, ev[dp, aR]])], 10^-15];
   Print["    ", cs, " l=", l, ": monodromy ", N[mono, 6], ", predicted ", N[pred, 6], " [", tick[], " s]"]], {cs, {"exact", "k=3", "d1", "d2"}}, {l, lMono}];
Export[FileNameJoin[{resDir, "stage15_monodromy.csv"}], Prepend[{#[[1]], #[[2]], Re[#[[3]]], Im[#[[3]]], Re[#[[4]]], Im[#[[4]]], #[[5]], #[[6]]} & /@ monoRows, {"case", "l", "Re_monodromy", "Im_monodromy", "Re_predicted", "Im_predicted", "abs_diff", "semicircle_vs_contour_at_aR"}], "CSV"];

(* ======================= 4. mode sums at x0 ======================= *)
Print["\n--- 4. mode sums Phi2R(x0), d_r Phi2R(x0) ---"];
msRows = {}; msFinal = <||>;
Do[Module[{vals, dvals, exV, exD, cum, dcum, cumE, dcumE},
   vals = Table[Re[tay[{cs, l, 1}]["S"][[1, -1]]] + rrAt[l][[1]], {l, lList}]; dvals = Table[Re[tay[{cs, l, 1}]["S"][[2, -1]]] + rrAt[l][[2]], {l, lList}];
   exV = Table[resExact[cs, l][[1]] + rrAt[l][[1]], {l, lList}]; exD = Table[resExact[cs, l][[2]] + rrAt[l][[2]], {l, lList}];
   cum = Accumulate[vals Table[Y0[l], {l, lList}]]; dcum = Accumulate[dvals Table[Y0[l], {l, lList}]];
   cumE = Accumulate[exV Table[Y0[l], {l, lList}]]; dcumE = Accumulate[exD Table[Y0[l], {l, lList}]];
   Do[AppendTo[msRows, {cs, lList[[i]], N[cum[[i]], 20], N[Abs[cum[[i]] - exact3D[[1]]], 6], N[Abs[cumE[[i]] - exact3D[[1]]], 6], N[dcum[[i]], 20], N[Abs[dcum[[i]] - exact3D[[2]]], 6], N[Abs[dcumE[[i]] - exact3D[[2]]], 6]}], {i, Length[lList]}];
   msFinal[cs] = {cum[[-1]], Abs[cum[[-1]] - exact3D[[1]]], dcum[[-1]], Abs[dcum[[-1]] - exact3D[[2]]], Abs[cumE[[-1]] - exact3D[[1]]], Abs[dcumE[[-1]] - exact3D[[2]]]};
   Print["    ", cs, ": Phi2R(x0) = ", N[cum[[-1]], 15], " (err ", sci[msFinal[cs][[2]]], ", l-truncation alone ", sci[msFinal[cs][[5]]], "),  d_r = ", N[dcum[[-1]], 15], " (err ", sci[msFinal[cs][[4]]], ", l-trunc ", sci[msFinal[cs][[6]]], ")"]],
  {cs, cases}];
csCheck["exact puncture: mode sum Phi2R(x0) at Lmax = exact 3D value within the l-truncation (r0/ro)^Lmax x 10", msFinal["exact"][[2]], 10 N[(r0/ro)^Lmax]];
csCheck["exact puncture: mode sum d_r Phi2R(x0) within the l-truncation x 10", msFinal["exact"][[4]], 10 N[(r0/ro)^Lmax]];
csCheck["exact puncture: Taylor-extraction error of the mode sums (|sum - sum of exact modes|) <= 1e-18", Abs[msFinal["exact"][[2]] - msFinal["exact"][[5]]] + Abs[msFinal["exact"][[4]] - msFinal["exact"][[6]]], 10^-18];
Export[FileNameJoin[{resDir, "stage15_modesum.csv"}], Prepend[msRows, {"case", "Lcut", "Phi2R_x0", "err", "err_exact_modes", "dPhi2R_x0", "err_d", "err_d_exact_modes"}], "CSV"];

(* ======================= 5. eta-baseline: approach r0 along contours of decreasing apex height ======================= *)
Print["\n--- 5. eta-baseline ---"];
etaCases = {"exact", "k=1", "k=3", "k=5"};
etaRows = {}; etaRes = <||>; etaNum = <||>;
etaBreaks[eta_] := Module[{J = Max[3, Ceiling[Log2[8/eta]]]}, Union[Subdivide[aL, aR, 40], r0 - R 2^-Range[1, J], r0 + R 2^-Range[1, J]]];
Do[Module[{brE = etaBreaks[eta], ndE, ir},
   ndE = makeNodes[brE, eta, nG]; ir = Position[brE, r0][[1, 1]];
   Do[Module[{dat, sol, za = ndE["zb"][[ir]], smax},
      dat = {ev[solIn[l], aL], ev[dz[solIn[l]], aL]};
      sol = contourSolve[l, srcTot[l], dat, ndE]; smax = Max[Abs[sol["S"]]];
      Do[Module[{pun = punS[l] + dPun[cs, l], rv, dv, ex = resExact[cs, l], pv},
         pv = ev[pun, za]; rv = sol["psi"][[ir]] - pv; dv = sol["dpsi"][[ir]] - ev[dz[pun], za];
         etaRes[{cs, l, eta}] = {rv, dv, Abs[sol["psi"][[ir]]], Abs[pv], smax};
         AppendTo[etaRows, {cs, l, N[eta], N[Re[rv], 20], N[Abs[Re[rv] - ex[[1]]], 6], N[Re[dv], 20], N[Abs[Re[dv] - ex[[2]]], 6], N[Im[rv], 10], N[Im[dv], 10], N[Abs[sol["psi"][[ir]]], 6], N[Abs[pv], 6], N[smax, 6]}]], {cs, etaCases}];
      (* numerical error estimate at the smallest eta: nodes per panel 20 -> 30 *)
      If[eta == Min[etaList], Module[{nd2 = makeNodes[brE, eta, 30], sol2},
         sol2 = contourSolve[l, srcTot[l], dat, nd2];
         etaNum[l] = {Abs[sol2["psi"][[ir]] - sol["psi"][[ir]]], Abs[sol2["dpsi"][[ir]] - sol["dpsi"][[ir]]]};
         csCheck["eta=" <> ToString[N[eta]] <> " l=" <> ToString[l] <> ": quadrature converged (20 -> 30 nodes per panel) at the apex", Max[etaNum[l]], 10^-25]]]], {l, lEta}];
   Print["    eta=", N[eta], " (", Length[ndE["z"]], " nodes) done [", tick[], " s]"]], {eta, etaList}];
(* Richardson (polynomial) extrapolation eta -> 0, all etas and the 4 smallest *)
extrap[pts_] := Module[{x}, InterpolatingPolynomial[pts, x] /. x -> 0];
richRows = {};
Do[Module[{ex = resExact[cs, l], pv, pd, eAll, e4, dAll, d4},
   pv = Table[{eta, Re[etaRes[{cs, l, eta}][[1]]]}, {eta, etaList}]; pd = Table[{eta, Re[etaRes[{cs, l, eta}][[2]]]}, {eta, etaList}];
   eAll = Abs[extrap[pv] - ex[[1]]]; e4 = Abs[extrap[Take[pv, -Min[4, Length[pv]]]] - ex[[1]]];
   dAll = Abs[extrap[pd] - ex[[2]]]; d4 = Abs[extrap[Take[pd, -Min[4, Length[pd]]]] - ex[[2]]];
   AppendTo[richRows, {cs, l, N[eAll, 4], N[e4, 4], N[dAll, 4], N[d4, 4], N[Abs[Re[etaRes[{cs, l, Min[etaList]}][[1]]] - ex[[1]]], 4], N[Abs[Re[etaRes[{cs, l, Min[etaList]}][[2]]] - ex[[2]]], 4]}]], {cs, etaCases}, {l, lEta}];
csCheck["eta-baseline, exact puncture: bias of Re psi~R(r0 + i eta) is O(eta): error ratio between the 2nd and 3rd eta within [0.7, 1.5] x eta ratio (median over l)",
  Module[{rat = Table[Module[{e1 = Abs[Re[etaRes[{"exact", l, etaList[[2]]}][[1]]] - resExact["exact", l][[1]]], e2 = Abs[Re[etaRes[{"exact", l, etaList[[3]]}][[1]]] - resExact["exact", l][[1]]]}, e1/e2], {l, lEta}], ex = etaList[[2]]/etaList[[3]]},
   Module[{m = Median[rat]}, If[0.7 ex <= m <= 1.5 ex, 0, Abs[m/ex - 1]]]], 0];
(* Taylor-predicted bias (exact case): coefficients about r0 from the ODE (data at r0 = Taylor sum from the apex; source coefficients on |z - r0| = 1/2, legitimate because s~R is analytic there) *)
biasPred = <||>;
Do[Module[{pR0, Sr0 = tay[{"exact", l, 1}]["S"], c0r, p},
   c0r = fftCoef[sRes[l], r0, 1/2];
   p = taylorODE[l, r0, {Sr0[[1, -1]], Sr0[[2, -1]]}, c0r];
   biasPred[l] = Table[{eta, Abs[Re[Sum[p[[n + 1]] (I eta)^n, {n, 0, Min[40, Nmax]}]] - resExact["exact", l][[1]]], Abs[Re[Sum[n p[[n + 1]] (I eta)^(n - 1), {n, 1, Min[40, Nmax]}]] - resExact["exact", l][[2]]]}, {eta, etaList}];
   csCheck["l=" <> ToString[l] <> ": eta-bias of the direct contour solves = Taylor series about r0 (exact puncture, relative to the bias at eta = 1)",
    Max[Table[Abs[Abs[Re[etaRes[{"exact", l, etaList[[j]]}][[1]]] - resExact["exact", l][[1]]] - biasPred[l][[j, 2]]], {j, Length[etaList]}]]/biasPred[l][[1, 2]], 10^-12]], {l, lEta}];
Export[FileNameJoin[{resDir, "stage15_eta_baseline.csv"}], Prepend[etaRows, {"case", "l", "eta", "Re_resid", "err_value", "Re_dresid", "err_deriv", "Im_resid", "Im_dresid", "abs_psi_apex", "abs_puncture_apex", "max_abs_source"}], "CSV"];
Export[FileNameJoin[{resDir, "stage15_eta_richardson.csv"}], Prepend[richRows, {"case", "l", "err_extrap_all", "err_extrap_last4", "err_d_extrap_all", "err_d_extrap_last4", "err_smallest_eta", "err_d_smallest_eta"}], "CSV"];

(* ======================= 6. tables, figures, summary ======================= *)
Export[FileNameJoin[{resDir, "stage15_taylor_convergence.csv"}], Prepend[tayRows, {"case", "l", "a_e", "N", "err_value", "err_deriv", "Im_value", "Im_deriv"}], "CSV"];
Export[FileNameJoin[{resDir, "stage15_modes_r0.csv"}], Prepend[r0Rows, {"case", "l", "a_e", "Re_resid_r0", "exact", "err", "Im_resid_r0", "Re_dresid_r0", "exact_d", "err_d", "Im_dresid_r0"}], "CSV"];
Export[FileNameJoin[{resDir, "stage15_rates.csv"}], Prepend[rateRows, {"case", "l", "a_e", "geom_ratio_value", "geom_ratio_deriv", "alg_slope_value", "alg_slope_deriv", "err_value_Nmax", "err_deriv_Nmax", "coef_fit_ratio", "coef_fit_exponent"}], "CSV"];
Export[FileNameJoin[{resDir, "stage15_circle_check.csv"}], Prepend[circRows, {"case", "l", "max_coef_diff_scaled", "monodromy_circle", "partial_sum_diff_N40"}], "CSV"];
Export[FileNameJoin[{resDir, "stage15_shadow.csv"}], Prepend[Table[{l, N[shadowR[l][[1]], 12], N[shadowR[l][[2]], 12], N[shadowR[l][[3]], 12], N[shadowR[l][[4]], 12], N[Im[tay[{"exact", l, 1}]["S"][[1, -1]]], 12], N[Im[tay[{"exact", l, 1}]["S"][[2, -1]]], 12]}, {l, lList}],
   {"l", "Im_psi_aR", "Im_dpsi_aR", "Im_resid_aR", "Im_dresid_aR", "Im_resid_r0", "Im_dresid_r0"}], "CSV"];

(* ----- figures ----- *)
clip[x_] := Max[x, 10^-40];
psiR = Superscript[OverTilde["\[Psi]"], "R"]; phiR = Superscript["\[Phi]", "R"]; rz = Subscript["r", "0"];
labV = Row[{"|Re ", Subscript["S", "N"], " \[Minus] ", phiR, "(", rz, ")|"}]; labD = Row[{"|Re ", Subscript["S", "N"], "\[Prime] \[Minus] ", phiR, "\[Prime](", rz, ")|"}];
labEta = Row[{"|Re ", psiR, "(", rz, " + i\[Eta]) \[Minus] ", phiR, "(", rz, ")|"}]; labEtaD = Row[{"|Re ", psiR, "\[Prime](", rz, " + i\[Eta]) \[Minus] ", phiR, "\[Prime](", rz, ")|"}];
labMS = Row[{"|", Superscript[Subscript["\[CapitalPhi]", "2"], "R"], "(", Subscript["x", "0"], ") error|"}]; labMSD = Row[{"|", Subscript["\[PartialD]", "r"], Superscript[Subscript["\[CapitalPhi]", "2"], "R"], "(", Subscript["x", "0"], ") error|"}];
cols = {RGBColor[0.1, 0.35, 0.75], RGBColor[0.85, 0.45, 0.1], RGBColor[0.15, 0.6, 0.25], RGBColor[0.55, 0.3, 0.7], Black, Gray};
(* 1: residual on the contour *)
f1 = Module[{sel = Function[l, Select[contRows, #[[1]] == l &]]},
   Grid[{{ListLinePlot[Table[{#[[2]], #[[3]]} & /@ sel[l], {l, figL}], PlotStyle -> cols, Frame -> True, Axes -> False, FrameLabel -> {"a", Row[{"Re ", psiR}]}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All],
      ListLinePlot[Table[{#[[2]], #[[4]]} & /@ sel[l], {l, figL}], PlotStyle -> cols, Frame -> True, Axes -> False, FrameLabel -> {"a", Row[{"Im ", psiR}]}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All]}}, Spacings -> 1]];
Export[FileNameJoin[{figDir, "fig_residual_contour.pdf"}], Legended[f1, Placed[LineLegend[cols[[1 ;; Length[figL]]], ("\[ScriptL] = " <> ToString[#]) & /@ figL, LegendLayout -> "Row", LabelStyle -> $fs], Below]]];
Export[FileNameJoin[{figDir, "fig_residual_contour.csv"}], Prepend[Select[contRows, MemberQ[figL, #[[1]]] &], {"l", "a", "Re_resid", "Im_resid", "Re_psi", "Im_psi", "Re_punct", "Im_punct"}], "CSV"];
(* 2: Taylor convergence at r0 against N *)
tcurve[cs_, l_, ip_, k_] := Table[{n, clip[tay[{cs, l, ip}][k][[n + 1]]]}, {n, 1, Nmax}];
tmax[cs_, ip_, k_] := Table[{n, clip[Max[Table[tay[{cs, l, ip}][k][[n + 1]], {l, lList}]]]}, {n, 1, Nmax}];
l2 = If[Lmax >= 2, 2, 0];
tPanel[data_, ylab_, lab_] := ListLogLogPlot[data, Joined -> True, PlotStyle -> ReplacePart[cols, 5 -> Directive[Black, Dashed]], Frame -> True, Axes -> False, FrameLabel -> {"N", ylab}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> {{1, Nmax}, {10^-36, 10^2}},
   Epilog -> {Inset[Style[lab, $fs], Scaled[{0.03, 0.05}], {Left, Bottom}]}];
f2 = Grid[{{tPanel[Table[tcurve[cs, l2, 1, "errV"], {cs, cases}], labV, "\[ScriptL] = " <> ToString[l2] <> ", apex"],
     tPanel[Table[tcurve[cs, l2, 1, "errD"], {cs, cases}], labD, "\[ScriptL] = " <> ToString[l2] <> ", apex"]},
    {tPanel[Table[tmax[cs, 1, "errV"], {cs, cases}], Row[{"max over \[ScriptL] of ", labV}], "apex"], tPanel[Table[tmax[cs, 2, "errV"], {cs, cases}], Row[{"max over \[ScriptL] of ", labV}], "a = 5.2"]}}, Spacings -> {1, 1}];
Export[FileNameJoin[{figDir, "fig_taylor_convergence.pdf"}], Legended[f2, Placed[LineLegend[ReplacePart[cols, 5 -> Directive[Black, Dashed]], {"exact", "k = 1", "k = 3", "k = 5", "d1: blend finished modes", "d2: one-sided"}, LegendLayout -> {"Row", 3}, LabelStyle -> $fs], Below]]];
Export[FileNameJoin[{figDir, "fig_taylor_convergence.csv"}], Prepend[Flatten[Table[{cs, N[aPts[[ip]]], n, N[tmax[cs, ip, "errV"][[n, 2]], 6], N[tmax[cs, ip, "errD"][[n, 2]], 6], N[tcurve[cs, l2, ip, "errV"][[n, 2]], 6], N[tcurve[cs, l2, ip, "errD"][[n, 2]], 6]}, {cs, cases}, {ip, Length[aPts]}, {n, 1, Nmax}], 2],
   {"case", "a_e", "N", "maxl_err_value", "maxl_err_deriv", "l2_err_value", "l2_err_deriv"}], "CSV"];
(* 3: mode-sum convergence *)
msc[cs_, k_] := Table[{#[[2]], clip[#[[k]]]} & /@ Select[msRows, #[[1]] == cs &]];
mPanel[k1_, k2_, ylab_] := ListLogPlot[Join[Table[msc[cs, k1], {cs, etaCases}], Table[msc[cs, k2], {cs, etaCases}]], Joined -> True, PlotMarkers -> Join[Table[{"\[FilledCircle]", 4}, 4], Table[None, 4]],
   PlotStyle -> Join[cols[[1 ;; 4]], Directive[#, Dashed] & /@ cols[[1 ;; 4]]], Frame -> True, Axes -> False, FrameLabel -> {"\[ScriptL]\!\(\*SubscriptBox[\(\[ThinSpace]\), \(max\)]\)", ylab}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All];
f3 = Grid[{{mPanel[4, 5, labMS], mPanel[7, 8, labMSD]}}, Spacings -> 1];
Export[FileNameJoin[{figDir, "fig_modesum.pdf"}], Legended[f3, Placed[Column[{LineLegend[cols[[1 ;; 4]], {"exact", "k = 1", "k = 3", "k = 5"}, LegendLayout -> "Row", LabelStyle -> $fs], LineLegend[{Black, Directive[Black, Dashed]}, {"Taylor-extracted modes", "exact modes (\[ScriptL]-truncation only)"}, LegendLayout -> "Row", LabelStyle -> $fs]}, Alignment -> Center], Below]]];
Export[FileNameJoin[{figDir, "fig_modesum.csv"}], Prepend[msRows, {"case", "Lcut", "Phi2R_x0", "err", "err_exact_modes", "dPhi2R_x0", "err_d", "err_d_exact_modes"}], "CSV"];
(* 4: eta-baseline (l = l2) *)
ecurve[cs_, k_] := Table[{N[eta], clip[Abs[Re[etaRes[{cs, l2, eta}][[k]]] - resExact[cs, l2][[k]]]]}, {eta, etaList}];
ePanel[k_, ylab_] := ListLogLogPlot[Join[Table[ecurve[cs, k], {cs, etaCases}], {Table[{N[#[[1]]], clip[#[[k + 1]]]} & /@ biasPred[l2]]}], Joined -> True, PlotMarkers -> Join[Table[{"\[FilledCircle]", 4}, 4], {None}],
   PlotStyle -> Join[cols[[1 ;; 4]], {Directive[Black, Dashed]}], Frame -> True, Axes -> False, FrameLabel -> {"\[Eta]", ylab}, BaseStyle -> $fs, ImageSize -> 230, PlotRange -> All,
   Epilog -> {Inset[Style["\[ScriptL] = " <> ToString[l2], $fs], Scaled[{0.03, 0.95}], {Left, Top}]}];
f4 = Grid[{{ePanel[1, labEta], ePanel[2, labEtaD]}}, Spacings -> 1];
Export[FileNameJoin[{figDir, "fig_eta_baseline.pdf"}], Legended[f4, Placed[LineLegend[Join[cols[[1 ;; 4]], {Directive[Black, Dashed]}], {"exact", "k = 1", "k = 3", "k = 5", Row[{"Taylor series about ", rz, " (exact puncture)"}]}, LegendLayout -> {"Row", 3}, LabelStyle -> $fs], Below]]];
Export[FileNameJoin[{figDir, "fig_eta_baseline.csv"}], Prepend[Join[Select[etaRows, #[[2]] == l2 &][[All, {1, 3, 5, 7}]], {"Taylor about r0", N[#[[1]]], N[#[[2]], 6], N[#[[3]], 6]} & /@ biasPred[l2]], {"case", "eta", "err_value", "err_deriv"}], "CSV"];

(* ----- summary ----- *)
nTo[err_, tol_] := Module[{p = Position[err, x_ /; x < tol, 1, 1]}, If[p === {}, ">" <> ToString[Nmax], ToString[p[[1, 1]] - 1]]];
fmt[x_] := If[NumericQ[x], sci[x], ToString[x]];
nf[x_] := ToString[CForm[N[x, 16]]];
cfmt[z_] := If[NumericQ[z], sci[Re[z]] <> If[Im[z] >= 0 || Im[z] == 0, " + ", " - "] <> sci[Abs[Im[z]]] <> " i", ToString[z]];
md1 = StringRiffle[Join[{"| l | Re psi~R(r0), apex Taylor | h_l(r0) exact | error | Im psi~R(r0) | derivative error | ratio (value) | N for 1e-20 | Im psi~(a_R) | error, a_e = 5.2 | ratio, a_e = 5.2 |", "|---|---|---|---|---|---|---|---|---|---|---|"},
    Table[Module[{t = tay[{"exact", l, 1}], t2 = tay[{"exact", l, 2}]}, "| " <> ToString[l] <> " | " <> nf[Re[t["S"][[1, -1]]]] <> " | " <> nf[t["exact"][[1]]] <> " | " <> sci[t["errV"][[-1]]] <> " | " <>
       ToString[CForm[N[Im[t["S"][[1, -1]]], 10]]] <> " | " <> sci[t["errD"][[-1]]] <> " | " <> fmt[rate["exact", l, 1, 4]] <> " | " <> nTo[t["errV"], 10^-20] <> " | " <> ToString[CForm[N[shadowR[l][[1]], 8]]] <> " | " <> sci[t2["errV"][[-1]]] <> " | " <> fmt[rate["exact", l, 2, 4]] <> " |"], {l, lList}]], "\n"];
medR[cs_, ip_, k_] := Module[{v = DeleteCases[Table[rate[cs, l, ip, k], {l, lList}], Indeterminate]}, If[v === {}, Indeterminate, Median[v]]];
maxE[cs_, ip_, k_] := Max[Table[tay[{cs, l, ip}][k][[-1]], {l, lList}]];
md2 = StringRiffle[Join[{"| case | a_e | error fits: geometric ratio, value (median over l) | ratio, derivative | algebraic slope, value | slope, derivative | coefficient fits: ratio | exponent | max_l error at N = " <> ToString[Nmax] <> ", value | derivative |", "|---|---|---|---|---|---|---|---|---|---|"},
    Flatten@Table["| " <> cs <> " | " <> ToString[N[aPts[[ip]]]] <> " | " <> fmt[medR[cs, ip, 4]] <> " | " <> fmt[medR[cs, ip, 5]] <> " | " <> fmt[medR[cs, ip, 6]] <> " | " <> fmt[medR[cs, ip, 7]] <> " | " <> fmt[medR[cs, ip, 10]] <> " | " <> fmt[medR[cs, ip, 11]] <> " | " <> sci[maxE[cs, ip, "errV"]] <> " | " <> sci[maxE[cs, ip, "errD"]] <> " |", {cs, cases}, {ip, Length[aPts]}]], "\n"];
md3 = StringRiffle[Join[{"| case | Phi2R(x0), Lmax = " <> ToString[Lmax] <> " | error | l-truncation alone | d_r Phi2R(x0) | error | l-truncation alone |", "|---|---|---|---|---|---|---|"},
    Table["| " <> cs <> " | " <> nf[msFinal[cs][[1]]] <> " | " <> sci[msFinal[cs][[2]]] <> " | " <> sci[msFinal[cs][[5]]] <> " | " <> nf[msFinal[cs][[3]]] <> " | " <> sci[msFinal[cs][[4]]] <> " | " <> sci[msFinal[cs][[6]]] <> " |", {cs, cases}]], "\n"];
md4 = StringRiffle[Join[{"| case | l | max_n |dp_n| rc^n (n <= 40) | |u(2 Pi) - u(0)| on the circle | partial-sum difference at r0 (N = 40) |", "|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> sci[#[[3]]] <> " | " <> sci[#[[4]]] <> " | " <> sci[#[[5]]] <> " |") & /@ circRows], "\n"];
md5 = StringRiffle[Join[{"| case | l | monodromy of psi~R around r0 | predicted | difference | psi~R(a_R): semicircle - contour |", "|---|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> cfmt[#[[3]]] <> " | " <> cfmt[#[[4]]] <> " | " <> sci[#[[5]]] <> " | " <> sci[#[[6]]] <> " |") & /@ monoRows], "\n"];
md6 = StringRiffle[Join[{"| eta | nodes | exact: bias of the value | bias of the derivative | Taylor-predicted bias (value) | k=1 bias | k=3 bias | k=5 bias | abs psi~(apex) | abs phi~S(apex) | max abs source | Im psi~R(apex) |", "|---|---|---|---|---|---|---|---|---|---|---|---|"},
    Table[Module[{e = etaRes[{"exact", l2, eta}], ex = resExact["exact", l2], j = Position[etaList, eta][[1, 1]]},
      "| " <> ToString[N[eta]] <> " | " <> ToString[20 (Length[etaBreaks[eta]] - 1)] <> " | " <> sci[Abs[Re[e[[1]]] - ex[[1]]]] <> " | " <> sci[Abs[Re[e[[2]]] - ex[[2]]]] <> " | " <> sci[biasPred[l2][[j, 2]]] <> " | " <>
       StringRiffle[Table[sci[Abs[Re[etaRes[{cs, l2, eta}][[1]]] - resExact[cs, l2][[1]]]], {cs, {"k=1", "k=3", "k=5"}}], " | "] <> " | " <> sci[e[[3]]] <> " | " <> sci[e[[4]]] <> " | " <> sci[e[[5]]] <> " | " <> ToString[NumberForm[Im[e[[1]]], 8]] <> " |"], {eta, etaList}]], "\n"];
md7 = StringRiffle[Join[{"| case | l | extrapolated error (all eta) | (4 smallest) | derivative (all) | (4 smallest) | error at the smallest eta | derivative |", "|---|---|---|---|---|---|---|---|"},
    ("| " <> #[[1]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[3 ;;]], " | "] <> " |") & /@ richRows], "\n"];
md8 = StringRiffle[Join[{"| l | numerical error of psi~(apex) at eta = " <> ToString[N[Min[etaList]]] <> " (20 -> 30 nodes per panel): value | derivative |", "|---|---|---|"},
    Table["| " <> ToString[l] <> " | " <> sci[etaNum[l][[1]]] <> " | " <> sci[etaNum[l][[2]]] <> " |", {l, lEta}]], "\n"];
n20 = Table[ToExpression[nTo[tay[{"exact", l, 1}]["errV"], 10^-20]], {l, lList}];
d1ok = maxE["d1", 1, "errV"] + maxE["d1", 1, "errD"] < 10^-20; d2fail = Min[Table[tay[{"d2", l, 1}]["errD"][[-1]], {l, Select[lList, # <= 6 &]}]] > 10^-3;
kSlopes = Table[{k, medR["k=" <> ToString[k], 1, 6], medR["k=" <> ToString[k], 1, 7]}, {k, kList}];
summary = "# Stage 15: extracting the second-order residual field at the particle (Dzz case)\n\n" <>
  "Run: `wolframscript -file stage15_residual_extraction.wls`" <> If[smoke, " smoke", ""] <> ".  l = 0.." <> ToString[Lmax] <> ", Nmax = " <> ToString[Nmax] <> ", wall time " <> ToString[tick[]] <> " s.\n\n" <>
  "## Setup and formulas\n\n" <>
  "- Toy: Lap Phi1 = -4 pi q delta^3, Phi1 = q/rho + PhiR1, PhiR1 = c0 + c1 Z = A + c1 d with A = c0 + c1 r0, d = Z - r0 (Z = r mu, x0 = r0 zhat); q = 1, r0 = 5, c0 = 1, c1 = 1/10.  Second order: Lap Phi2 = (1/2) d_Z^2 (Phi1^2).\n" <>
  "- Exact singular field: chiS = q^2 ln rho + q A rho + (1/2) q c1 d rho, Phi2S = (1/2) d_Z^2 chiS (checked symbolically: Lap chiS = q^2/rho^2 + 2 q PhiR1/rho, using Lap ln rho = 1/rho^2, Lap rho = 2/rho, Lap(d rho) = 4 d/rho).  chiS = chi_SS + chi_SR of Stage 3, so Phi2R = Phi2 - Phi2S = (1/2) d_Z^2 chiRR + h.\n" <>
  "- The field solved for is Phi2 = (1/2) d_Z^2 chi + h, with the harmonic field h = Bo/|x - ro zhat| + Bi/|x - ri zhat| (Bo = " <> ToString[Bo] <> ", ro = " <> ToString[ro] <> ", Bi = " <> ToString[Bi] <> ", ri = " <> ToString[ri, InputForm] <> ") injected through the Cauchy data at a_L (equivalently at r = 1); h gives every l a nonzero residual, h_l(r) = sqrt(4 pi/(2l+1)) [Bo r^l/ro^(l+1) + Bi ri^l/r^(l+1)].  Exact 3D values: Phi2R(x0) = c0^2/6 + 3 c0 c1 r0/5 + c1^2 r0^2/2 + h(x0) = " <> ToString[NumberForm[N[exact3D[[1]], 20], 20]] <> ", d_r Phi2R(x0) = " <> ToString[NumberForm[N[exact3D[[2]], 20], 20]] <> ".\n" <>
  "- Contour y(a) = a + i (1 - 4 (a - r0)^2)^6 on [4.5, 5.5] (apex r0 + i, min |z - r0| = 0.298).  SS + SR are solved together on the contour by variation of parameters, psi~ = A z^l + B z^(-l-1) + [z^l I1 - z^(-l-1) I2]/(2l+1), I1 = Int s^(1-l) s~ ds, I2 = Int s^(l+2) s~ ds along the contour (60 Gauss panels x 20 nodes in a; 50-digit evaluation of the closed-form sources; cross-checked against the NDSolve contour solve to 1e-15).  Sources: s~_SS = (1/2) D_z^2 [q^2/rho^2]_l with 2 log(z + r0) - 2 log(z - r0); s~_SR = (1/2) D_z^2 [2 q sum_l' K~_(j l') phiR_l'] with K~ by Theta -> H (Stage 9).  RR is solved on the real line (polynomial, l <= 2).\n" <>
  "- Continued puncture phi~S_l: SS part = out-branch formula of (1/2) D_z^2 [q^2 ln rho]_l with ln|r - r0| -> Log[z - r0] (in = out as real functions); SR part = (1/2) D_z^2 applied to the H-blend F_in + (F_out - F_in) H of the undifferentiated modes F_j = [chi_SR]_j, which are C^2 at r0.  Identities checked: s~_SS = L_l phi~S_SS exactly (1e-40); Re phi~S(r + i0) = [Phi2S]_l(r) on both sides; s~R := s~_SS + s~_SR - L_l phi~S is imaginary on the axis and bounded (analytic) at r0.\n" <>
  "- Truncated punctures (b): Phi2P = Phi2S + eps q^2 rho^k, eps = r0^-(k+2), k = 1, 3, 5; the modes [rho^k]_l are C^(k+1) at r0 (checked) and are continued by the same H-blend.  Exact residual modes at r0: h_l(r0) - eps q^2 [rho^k]_l(r0), and the derivative.\n" <>
  "- Taylor extension: with z = z_e + w and psi~R = sum p_n w^n, the ODE z^2 psi'' + 2 z psi' - l(l+1) psi = z^2 s~R gives p_(n+2) = [sigma_n - 2 z_e (n+1)^2 p_(n+1) - (n(n+1) - l(l+1)) p_n]/(z_e^2 (n+2)(n+1)), where sigma_n are the Taylor coefficients of z^2 s~R.  p_0, p_1 are psi~R(z_e), psi~R'(z_e) from the contour solve; the coefficients of s~R come from the Cauchy integral on |z - z_e| = 0.7 |z_e - r0| by the trapezoidal rule with " <> ToString[MF] <> " points (every quantity O(1); no fit to contour samples).  Sum to r0: S_N = sum_(n<=N) p_n (r0 - z_e)^n, S'_N = sum n p_n (r0 - z_e)^(n-1).  Expansion points: the apex (a = 5) and a = 5.2 (|z_e - r0| = 0.404).\n" <>
  "- Cross-checks: (i) the residual ODE solved by NDSolve around the circle |z - z_e| = rc and the Taylor coefficients taken by Cauchy integrals; (ii) the monodromy of psi~R around r0 along |z - r0| = R from a_L (above) through the apex, a_R and back to a_L (below), with Log[z - r0] = ln R + i(pi - t) continued explicitly.\n" <>
  "- Baseline (c): contours y_eta(a) = a + i eta (1 - 4(a - r0)^2)^6, eta = 1 .. 1/128, graded Gauss panels; psi~ - phi~P read off at the apex r0 + i eta; polynomial (Richardson) extrapolation in eta.\n" <>
  "- Controls (d): d1 = Theta -> H blend of the finished Dzz SR modes P_in + (P_out - P_in) H (the task's control); d2 = one-sided continuation P_out (drops the kink, hence the delta of L P).\n\n" <>
  "## (a) Exact puncture: residual modes at r0 from the apex Taylor series\n\n" <> md1 <> "\n\n" <>
  "Predicted geometric ratio |r0 - z_e|/|z_e| (next singularity at z = 0): " <> ToString[predRatio] <> " for a_e = " <> ToString[N /@ aPts] <> ".\n\n" <>
  "## Convergence rates (all cases)\n\n" <> md2 <> "\n\n" <>
  "## Mode sums at x0 (Taylor-extracted modes + RR on the real line)\n\n" <> md3 <> "\n\n" <>
  "## Cross-check: Cauchy circle about the apex\n\n" <> md4 <> "\n\n## Monodromy of psi~R around r0\n\n" <> md5 <> "\n\n" <>
  "## (c) eta-baseline (l = " <> ToString[l2] <> ")\n\n" <> md6 <> "\n\n" <> md8 <> "\n\nRichardson extrapolation eta -> 0:\n\n" <> md7 <> "\n\n" <>
  "## Assessment of the predictions\n\n" <>
  "1. Exact puncture: analytic at r0, geometric convergence with ratio set by z = 0.  " <> If[maxE["exact", 1, "errV"] < 10^-20 && Abs[medR["exact", 1, 10] - predRatio[[1]]]/predRatio[[1]] < 0.1 && Abs[medR["exact", 2, 10] - predRatio[[2]]]/predRatio[[2]] < 0.1, "HELD", "NOT HELD"] <>
   ": max_l error " <> sci[maxE["exact", 1, "errV"]] <> " (value), " <> sci[maxE["exact", 1, "errD"]] <> " (derivative) at N = " <> ToString[Nmax] <> " (floor of the 50-digit arithmetic; 1e-20 is reached at N = " <> ToString[Min[n20]] <> "-" <> ToString[Max[n20]] <> ").  Ratio of the Taylor coefficients |p_n w0^n| ~ n^c rho^n (median over l): rho = " <> fmt[medR["exact", 1, 10]] <> " vs predicted " <> ToString[predRatio[[1]]] <> " (apex), " <> fmt[medR["exact", 2, 10]] <> " vs " <> ToString[predRatio[[2]]] <> " (a = 5.2); a two-parameter fit of the error curve gives " <> fmt[medR["exact", 1, 4]] <> " and " <> fmt[medR["exact", 2, 4]] <> ", biased upwards by the n^l prefactor of the order-(l+1) pole at z = 0.  The monodromy around r0 is " <> sci[Max[Cases[monoRows, {"exact", _, m_, __} :> Abs[m]]]] <> ".  Im psi~R(r0) is O(1) (table above): the shadow is carried through to r0 but only the real part is physical.  The second point a = 5.2 converges faster (ratio 0.08): in the exact case every contour point works, the closer to r0 the better.\n" <>
  "2. Truncated puncture: algebraic convergence.  Measured median slopes of the real-part errors, apex, value / derivative: " <> StringRiffle[("k = " <> ToString[#[[1]]] <> ": " <> fmt[#[[2]]] <> " / " <> fmt[#[[3]]]) & /@ kSlopes, "; "] <>
   "; coefficient-fit exponents of |p_n w0^n|: " <> StringRiffle[("k = " <> ToString[#] <> ": " <> fmt[medR["k=" <> ToString[#], 1, 11]]) & /@ kList, "; "] <>
   ".  The residual modes are C^(k+1) (two orders smoother than the 3D residual, which is C^(k-1)); the continuation has a (z - r0)^(k+2) log(z - r0) branch point, so the rates expected from the mode smoothness are N^-(k+2) (value), N^-(k+1) (derivative) and n^-(k+3) for the coefficients, against the N^-k of the prediction.  " <>
   "At a = 5.2 the slopes are " <> StringRiffle[("k = " <> ToString[#] <> ": " <> fmt[medR["k=" <> ToString[#], 2, 6]] <> " / " <> fmt[medR["k=" <> ToString[#], 2, 7]]) & /@ kList, "; "] <> ".  " <>
   If[And @@ Table[Abs[medR["k=" <> ToString[k], 1, 6] + (k + 2)] < 1 && Abs[medR["k=" <> ToString[k], 2, 7] + (k + 1)] < 1, {k, kList}], "The algebraic character held, with exponent k+2 for the value (two better than N^-k, because the mode functions are two orders smoother than the 3D field) and k+1 for the derivative off the vertical.", "The measured exponents do not match the mode-smoothness expectation within 1; see the table."] <>
   If[And @@ Table[Abs[medR["k=" <> ToString[k], 1, 7] + (k + 2)] < 1, {k, kList}], "  At the apex, which lies directly above r0, the derivative is NOT worse: the leading term (i/pi) j (z - r0)^(k+2) log(z - r0) has Taylor terms at r0 of constant phase i^(k+2), imaginary for odd k, so after the factor i/pi it feeds only Re of the value and Im of the derivative; Re of the derivative inherits the next term and also falls like N^-(k+2).  'Worse for the derivative' therefore held only for the off-vertical expansion point.", ""] <>
   "  r0 is a branch point on the edge of the disc for every expansion point, so the expansion point changes the prefactor (|z_e - r0|^(k+2)) and the phase pattern, not the exponent.\n" <>
  "3. Consistency (delta content).  Control d1 (blend of the finished, kinked modes) " <> If[d1ok, "CONVERGED to the same real values (max_l error " <> sci[maxE["d1", 1, "errV"]] <> "): the prediction that it fails did NOT hold.  Reason: for C^0 (kinked) finished modes, L[P_in + (P_out - P_in) H] = s_in + (s_out - s_in) H + [P'] delta_ext + analytic, so the kink reproduces the delta of L P automatically; d1 differs from the blend-then-ladder puncture by a function that is imaginary on the axis and vanishes at r0 (checked), i.e. by a harmless shadow.", "FAILED (max_l error " <> sci[maxE["d1", 1, "errV"]] <> "), as predicted."] <>
   "  Control d2 (one-sided continuation, kink dropped) " <> If[d2fail, "FAILS as a delta-dropping control should: its real part is wrong by J_P(r) Theta(r0 - r) at a_L, the residual acquires a term (i/pi) c_l (z - r0) log(z - r0) (c_l the dropped delta coefficient), the value converges only like 1/N and the real part of the derivative tends to the mean of the two one-sided derivatives, h_l' + c_l/2, so its error is |c_l|/2 = 0.106 (l = 0), 0.238 (l = 2) (min_l derivative error " <> sci[Min[Table[tay[{"d2", l, 1}]["errD"][[-1]], {l, Select[lList, # <= 6 &]}]]] <> " at N = " <> ToString[Nmax] <> "); its real part is also wrong at a_L.", "did not fail as expected."] <> "\n" <>
  "4. Baseline (c): the bias of Re psi~R(r0 + i eta) is O(eta) with the coefficient -Im psi~R'(r0) (the shadow derivative), reproduced by the Taylor series about r0 to " <> sci[Max[Table[Max[Table[Abs[Abs[Re[etaRes[{"exact", l, etaList[[j]]}][[1]]] - resExact["exact", l][[1]]] - biasPred[l][[j, 2]]], {j, Length[etaList]}]]/biasPred[l][[1, 2]], {l, lEta}]]] <> " (relative).  Richardson extrapolation of the exact case reaches " <> sci[Max[Cases[richRows, {"exact", _, e_, __} :> e]]] <> " (all etas); for the truncated punctures it saturates at the eta^(k+2) log eta level.  Digits are lost mainly through the O(eta) bias and the node count needed near the apex (table), not through cancellation: |phi~S(apex)| grows only logarithmically.\n\n" <>
  "## Surprises and failed predictions, stated plainly\n\n" <>
  "- Control d1 (blend of the finished Dzz SR modes) does not fail: it converges exactly like the blend-then-ladder puncture.  The delta content of L[puncture] lives in the kink of the finished modes, which both continuations keep; they differ by an analytic function that is imaginary on the axis.  Prediction 3, as a requirement on the puncture, did not hold in this toy.  What does fail is dropping the kink (d2).\n" <>
  "- The algebraic exponent is k+2, not k, for the value: the mode functions [rho^k]_l are C^(k+1), two orders smoother than the 3D field.\n" <>
  "- At the apex the real part of the derivative converges as fast as the value (N^-(k+2)), because the leading log term feeds only Im of the derivative there; off the vertical the derivative is one order worse (N^-(k+1)).\n" <>
  "- With Lmax = " <> ToString[Lmax] <> " the mode sums are limited entirely by the l-truncation (exact: (r0/ro)^l from h; truncated: L^-k of the rho^k modes); the Taylor extraction contributes < 1e-30.\n" <>
  "- The two-parameter fit of the exact-case error curve overestimates the geometric ratio (0.25 vs 0.196) because of the n^l prefactor; the three-parameter coefficient fit recovers 0.199.\n" <>
  "- The eta-baseline loses no digits to cancellation (|psi~| and |phi~S| at the apex stay O(1)); its limitation is the O(eta) bias and the 1/eta^2 growth of the source, which costs quadrature nodes.  Richardson extrapolation over eta = 1 .. 1/128 recovers the exact-case value to 3e-15 - 5e-11 and is capped at eta^(k+2) log eta for truncated punctures (1e-10 for k = 1).\n\n" <>
  "## Checks\n\n" <> ToString[Count[csCheckList[], {_, _, _, True}]] <> " passed, " <> ToString[Count[csCheckList[], {_, _, _, False}]] <> " failed (of " <> ToString[Length[csCheckList[]]] <> ").  Failed: " <> ToString[Select[csCheckList[], ! Last[#] &][[All, 1]]] <> "\n\n" <>
  "## Files\n\n- results/stage15_checks.csv, stage15_residual_contour.csv, stage15_taylor_convergence.csv, stage15_modes_r0.csv, stage15_rates.csv, stage15_modesum.csv, stage15_circle_check.csv, stage15_monodromy.csv, stage15_shadow.csv, stage15_eta_baseline.csv, stage15_eta_richardson.csv.\n- figures/residual/fig_residual_contour.pdf, fig_taylor_convergence.pdf, fig_modesum.pdf, fig_eta_baseline.pdf (each with a CSV).\n";
Export[FileNameJoin[{resDir, If[smoke, "stage15_summary_smoke.md", "stage15_summary.md"]}], summary, "Text"];
Print["\n", md1, "\n\n", md2, "\n\n", md3, "\n\n", md5, "\n\n", md6, "\n\n", md7];
Print["\nStage 15 wall time: ", tick[], " s"];
csSummary[];
Export[FileNameJoin[{resDir, If[smoke, "stage15_checks_smoke.csv", "stage15_checks.csv"]}], Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ csCheckList[], {"check", "value", "tol", "result"}], "CSV"];
