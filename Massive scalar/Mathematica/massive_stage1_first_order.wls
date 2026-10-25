#!/usr/bin/env wolframscript
(* ::Package:: *)
(* Massive scalar, Stage 1: first order.  (Lap - m^2) Phi1 = -4 Pi q delta^3,  Phi1 = q Exp[-m rho]/rho + PhiR.
   (a) homogeneous solutions u_l, v_l (modified spherical Bessel), Wronskian, Green's function, the massive Coulomb mode three ways;
   (b) the continued singular mode (T3): Theta -> H blend of u_l(r<) v_l(r>): (i) Re = mode of the 3D puncture on both sides,
       (ii) L psi - deltaExt = (I/Pi) E0 with E0 real-analytic on the axis (violates (ii)); the corrected mode psi - (I/Pi) P, L P = E0,
       satisfies both; continuity along the contour; DSolve attempt (the massless App. C route);
   (c) first-order residual solves on the contour with S_eff = kappa (deltaExt - L phiP) for the Theta->H, corrected and truncated punctures
       PhiP_n (series of Exp[-m rho]/rho through rho^n): shadow at aR and exterior exactness;
   (d) smoothness of the residual Phi1 - PhiP_n: first non-smooth term rho^(n+2) (odd power) -> C^(n+1) in 3D, modes C^(n+3);
   (e) m = 0 reproduces the massless package (closed-form phiS_0 of the massless Stage 1; shadows and exterior errors of its Stage 2).
   Run:  wolframscript -file massive_stage1_first_order.wls [smoke]                                                           *)

$HistoryLength = 0;
dir = DirectoryName[$InputFileName];
Get[FileNameJoin[{dir, "ContourMassive.wl"}]];
resDir = FileNameJoin[{ParentDirectory[dir], "results"}]; figDir = FileNameJoin[{ParentDirectory[dir], "figures"}];
smoke = MemberQ[$ScriptCommandLine, "smoke"];
csResetChecks[]; t0 = AbsoluteTime[]; tick[] := Round[AbsoluteTime[] - t0, 0.1];
mList = {1/5, 1, 0}; lList = If[smoke, {0, 2}, Range[0, 6]]; nList = {1, 3, 5};
r0 = 5; q = 1; rMin = 1; rMax = 10; aL = 9/2; aR = 11/2;
Print["=== Massive Stage 1: first order; m = ", mList, ", l = ", lList, " ==="];
Off[N::meprec]; Off[General::munfl];
rows = {}; shadowRows = {}; smoothRows = {}; dsolveRows = {};
mStr[m_] := ToString[m, InputForm];
(* P: real particular solution of L P = E0 (E0 analytic, removable 0/0 at r0: Taylor fallback there), from zero data at rMin *)
solveP[l_] := Module[{E0 = singularResidualE[l], ser, z, E0s, u, a, P, Pc},
   ser = Normal[Series[E0[z], {z, r0, 8}]];
   E0s[x_?NumericQ] := If[Abs[x - r0] < 10^-4, N[ser /. z -> SetPrecision[x, 60], 40], N[E0[SetPrecision[x, 60]], 40]];
   P = NDSolveValue[{realL[l, u, a] == E0s[a], u[rMin] == 0, u'[rMin] == 0}, u, {a, rMin, rMax}, WorkingPrecision -> 30, PrecisionGoal -> 22, AccuracyGoal -> 30, MaxSteps -> 10^6, InterpolationOrder -> All];
   Pc = solveContourSegment[l, numericFunction[E0[yC[a]], a], {P[aL], P'[aL]}];
   {P, Pc, E0}];

Do[SetCMParameters[<|"m" -> m|>];
 Print["\n--- m = ", m, " ---"];
 (* (a) homogeneous solutions *)
 Do[csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": Wronskian u v' - u' v = -(2l+1)/r^2 (r = 1, 3, 10)",
    Max[Table[Abs[N[uFun[l, r] vFunD[l, r] - uFunD[l, r] vFun[l, r] + (2 l + 1)/r^2, 30]], {r, {1, 3, 10}}]], 10^-25], {l, {0, 1, 2, 6, 12}}];
 Do[Module[{z}, csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": L_l u_l = 0 and L_l v_l = 0 at r = 2, 7",
     Max[Table[Abs[N[(radialL[l, #, z]) /. z -> r, 30]], {r, {2, 7}}] & /@ {uFun[l, z], vFun[l, z]}], 10^-25]], {l, {0, 3}}];
 Module[{l = 2, s = 4, z, G}, G[r_] := greenPhys[l, r, s];
  csCheck["m=" <> mStr[m] <> ": Green's function: L G = 0 off s, unit jump of d_r G at s, Robin conditions at 1 and 10 (l = 2, s = 4)",
   Abs[N[radialL[l, -(s^2/(2 l + 1)) uFun[l, s] vFun[l, z], z] /. z -> 7, 30]] + Abs[N[radialL[l, -(s^2/(2 l + 1)) uFun[l, z] vFun[l, s], z] /. z -> 2, 30]] +
    Abs[N[-(s^2/(2 l + 1)) (uFun[l, s] vFunD[l, s] - uFunD[l, s] vFun[l, s]) - 1, 30]] +
    Abs[N[robinResidual[l, G[1], -(s^2/(2 l + 1)) uFunD[l, 1] vFun[l, s], 1], 30]] + Abs[N[robinResidual[l, G[10], -(s^2/(2 l + 1)) uFun[l, s] vFunD[l, 10], 10], 30]], 10^-25]];
 Do[Module[{num, cf, sub, r = rr},
    num = Sqrt[Pi (2 l + 1)] NIntegrate[Exp[-m Sqrt[r^2 + r0^2 - 2 r r0 mu]]/Sqrt[r^2 + r0^2 - 2 r r0 mu] LegendreP[l, mu], {mu, -1, 1}, WorkingPrecision -> 40, PrecisionGoal -> 28];
    cf = N[Sqrt[4 Pi/(2 l + 1)] coulombMassive[l, r], 40];
    sub = If[m == 0, cf, evR[logReal[expRhoModePair[-1, m, l][[If[r < r0, 1, 2]]]], r]];      (* the exponential form is not defined at m = 0 *)
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> " r=" <> ToString[r] <> ": [Exp[-m rho]/rho]_l: angular quadrature = Bessel form = rho-substitution form", Abs[num - cf] + Abs[sub - cf], 10^-25]],
  {l, If[smoke, {0, 2}, {0, 1, 2, 5}]}, {rr, {3, 7}}];
 csCheck["m=" <> mStr[m] <> ": PhiR(x0) = c0 + c1 r0 and PhiR is a homogeneous solution (modes l = 0, 1)",
  Abs[N[phiR3D[r0, 1] - phiRAtParticle[], 30]] + Max[Table[Abs[N[radialL[l, phiRMode[l][z], z] /. z -> 3, 30]], {l, 0, 1}]], 10^-25];

 (* (b) continued singular mode *)
 Do[Module[{f = singularModeThetaH[l], E0 = singularResidualE[l], z, pts, reErr, resAxis, cont, P, Pc, corr, lres, jump, u, a},
    pts = Join[Subdivide[1, 9/2, 8], Subdivide[11/2, 10, 8]];
    reErr = Max[Abs[N[Re[f[#]] - coulombModeUnit[l, #], 40]] & /@ pts];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (i) Re of the Theta->H singular mode = G_l on both sides of r0 (17 points)", reErr, 10^-30];
    (* (ii): L psi - deltaExt = (I/Pi) E0, E0 real on the axis and analytic at r0 *)
    resAxis = Table[N[(radialL[l, f[z], z] - deltaExt[z]) /. z -> r, 40], {r, pts}];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (ii) Theta->H: Re(L psi - deltaExt) = 0 on the axis", Max[Abs[Re[resAxis]]], 10^-30];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": (ii) Theta->H: L psi - deltaExt = (I/Pi) E0 on the contour (5 points)",
     Max[Table[Abs[N[(radialL[l, f[z], z] - deltaExt[z] - I/Pi E0[z]) /. z -> N[yC[a], 50], 40]], {a, {46/10, 48/10, 5, 52/10, 54/10}}]], 10^-30];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": E0(r0) = (3/2) J''(r0) + 2/r0 = -1/r0 for every l and m (analytic at r0)",
     Abs[N[E0[r0 + 10^-12] + 1/r0, 30]], 10^-10];
    (* corrected mode: psi - (I/Pi) P with L P = E0 on the real line (P real), continued along the contour by the ODE.
       L phiS = deltaExt then holds identically; what is tested is that the continuation of P along the contour returns to the real
       solution at aR (no 2 Pi I jump: E0 is analytic) and that the contour ODE residual of the interpolant is small *)
    {P, Pc} = solveP[l][[1 ;; 2]];
    lres = Table[With[{aa = SetPrecision[a, 30]}, Abs[N[contourL[l, Pc, aa] - E0[yC[aa]], 25]]], {a, {46/10, 48/10, 5, 52/10, 54/10}}];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": corrected mode psi - (I/Pi) P: contour ODE residual of P~ (interpolant, 5 points, relative to max|E0|)", Max[lres]/Max[Table[Abs[N[E0[yC[a]], 20]], {a, {46/10, 48/10, 5, 52/10, 54/10}}]], 10^-8];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": corrected mode: P~ continued along the contour returns to the real P at aR (value and derivative, relative; Im = 0)", (Abs[Pc[aR] - P[aR]] + Abs[Pc'[aR] - P'[aR]])/(Abs[P[aR]] + Abs[P'[aR]]), 10^-12];
    corr[a_?NumericQ] := f[yC[a]] - I/Pi Pc[a];
    (* continuity of the Theta -> H mode along the contour: NDSolve from aL with deltaExt source reproduces f(y(a)) *)
    cont = solveContourSegment[l, numericFunction[deltaExt[yC[a]] + I/Pi E0[yC[a]], a], {f[aL], Module[{t}, D[f[t], t] /. t -> aL]}];
    csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": Theta->H mode is continuous along the contour (ODE from aL reproduces it at aR)", Abs[cont[aR] - f[aR]], 10^-15];
    jump = Max[Abs[Im[resAxis]]];
    AppendTo[rows, {m, l, reErr, Max[Abs[Re[resAxis]]], jump, N[E0[r0 + 10^-12], 6], Max[lres], Abs[Im[f[aR]]], Abs[Im[f[aR] - I/Pi Pc[aR]]]}];
    Print["    l=", l, ": Theta->H shadow source max|Im(L psi - deltaExt)| on the axis = ", sci[jump], "; Im psi(aR) = ", sci[Im[f[aR]]], "; corrected Im(aR) = ", sci[Im[f[aR] - I/Pi Pc[aR]]], "  [", tick[], " s]"]],
  {l, lList}];

 (* DSolve attempt (massless App. C): does DSolve return a particular solution, and what homogeneous offset D_l = Re psi_p(r+i0) - G_l does it carry? *)
 If[m != 0, Do[Module[{psi, r, sol, t, pp, z, Dv, al, be, fitRes, pts = {3/2, 2, 3, 4, 6, 7, 8, 9}},
    sol = TimeConstrained[Quiet[DSolveValue[psi''[r] + (2/r) psi'[r] - (l (l + 1)/r^2 + m^2) psi[r] == I/(Pi (r - r0)), psi[r], r]], 120, $Aborted];
    If[sol === $Aborted || ! FreeQ[sol, DSolveValue] || ! FreeQ[sol, Integrate] || ! FreeQ[sol, Inactive],
     t = If[sol === $Aborted, "timed out (120 s)", "unevaluated"]; AppendTo[dsolveRows, {m, l, t, "-", "-", "-"}],
     pp = Function @@ {{z}, (sol /. C[_] -> 0) /. r -> z};
     Dv = Table[N[Re[pp[SetPrecision[x, 80] + I 10^-40]], 40] - N[coulombModeUnit[l, x], 40], {x, pts}];      (* from above *)
     {al, be} = LinearSolve[N[{{uFun[l, 2], vFun[l, 2]}, {uFun[l, 3], vFun[l, 3]}}, 40], Dv[[{2, 3}]]];
     fitRes = Max[Abs[Dv - Table[N[al uFun[l, x] + be vFun[l, x], 40], {x, pts}]]];
     t = "closed form (" <> StringTake[ToString[InputForm[sol /. C[_] -> 0]], UpTo[120]] <> "...)";
     AppendTo[dsolveRows, {m, l, t, sci[al], sci[be], sci[fitRes]}];
     csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": DSolve particular solution: Re psi_p(r+i0) - G_l is a homogeneous solution alpha u_l + beta v_l (8 points, both sides)", fitRes, 10^-25]];
    Print["    DSolve l=", l, ": ", t, "  [", tick[], " s]"]], {l, {0, 2}}]];

 (* (c) first-order residual solves with S_eff = kappa (deltaExt - L phiP): PhiR injected through Cauchy data at rMin.
    The reference singular field at the tube ends is the corrected mode (conditions (i) and (ii)); the residual on the contour is
    phiR~ + kappa (phiS~ - phiP~), and at aR the puncture is added back.  Exterior = Re of the data at aR. *)
 Do[Module[{phiA, d, f, E0, P, Pc, t, exF, exD, pts, names},
    pts = Join[Subdivide[1, 9/2, 8], Subdivide[11/2, 10, 8]];
    phiA = toFunction[uFun[l, zVar]/uFun[l, r0], zVar];       (* test regular field for every l: the homogeneous solution regular at 0, normalised to 1 at r0 *)
    d = {phiA[rMin], phiA'[rMin]};
    f = singularModeThetaH[l]; {P, Pc, E0} = solveP[l];
    exF = Function @@ {{a}, f[yC[a]] - I/Pi Pc[a]}; exD[aa_] := (D[f[t], t] /. t -> aa) - I/Pi Pc'[aa];      (* corrected mode and its z-derivative at the ends (y' = 1 there) *)
    names = Join[{"ThetaH", "corrected"}, ("n=" <> ToString[#]) & /@ nList];
    (* ThetaH / corrected: S_eff = kappa (deltaExt - L phiP), reference singular field = corrected mode; the residual on the contour
       is compared with the entire phiR~ (its imaginary part is the residual's shadow).  Truncated: S_eff = kappa L(f - phiP_n) with f
       the Theta->H exact mode, so that the residual equals the closed form phiR~ + kappa (f - phiP_n)~ (massless Stage 6 check). *)
    Do[Module[{pf, pdL, pdR, sL, sR, pL, pR, s1, s2, s3, vR, dR, seff, shadow, err, pfun, a, aGrid = Subdivide[aL, aR, 24], ref, cmp},
       Which[v === "ThetaH", pf = Function @@ {{a}, f[yC[a]]}; pdL = D[f[t], t] /. t -> aL; pdR = D[f[t], t] /. t -> aR;
         seff = numericFunction[kappaL[l] (deltaExt[yC[a]] - contourL[l, pf, a]), a]; sL = kappaL[l] {exF[aL], exD[aL]}; sR = kappaL[l] {exF[aR], exD[aR]};
         ref = Table[phiA[N[yC[aa], 40]] - kappaL[l] I/Pi Pc[aa], {aa, aGrid}],                 (* the residual's shadow is the continued -(I/Pi) kappa P *)
        v === "corrected", pf = exF; pdL = exD[aL]; pdR = exD[aR]; seff = Function[x, 0]; sL = kappaL[l] {exF[aL], exD[aL]}; sR = kappaL[l] {exF[aR], exD[aR]};
         ref = Table[phiA[N[yC[aa], 40]], {aa, aGrid}],
        True, Module[{nn = ToExpression[StringTake[v, {3, -1}]]}, pfun = toFunction[punctureModeCont[punctureCoeffs[nn], l]/kappaL[l], zVar];
          pf = Function @@ {{a}, pfun[yC[a]]}; pdL = D[pfun[t], t] /. t -> aL; pdR = D[pfun[t], t] /. t -> aR;
          seff = numericFunction[kappaL[l] (contourL[l, Function @@ {{a}, f[yC[a]]}, a] - contourL[l, pf, a]), a];
          sL = kappaL[l] {f[aL], D[f[t], t] /. t -> aL}; sR = kappaL[l] {f[aR], D[f[t], t] /. t -> aR};
          ref = Table[With[{zz = N[yC[aa], 40]}, phiA[zz] + kappaL[l] (f[zz] - pfun[zz])], {aa, aGrid}]]];
       pL = kappaL[l] {pf[aL], pdL}; pR = kappaL[l] {pf[aR], pdR};
       s1 = solveRealSegment[l, Function[r, 0], {rMin, aL}, d];
       s2 = solveContourSegment[l, seff, {s1[aL] + sL[[1]] - pL[[1]], s1'[aL] + sL[[2]] - pL[[2]]}];
       vR = s2[aR] + pR[[1]] - sR[[1]]; dR = s2'[aR] + pR[[2]] - sR[[2]];
       s3 = solveRealSegment[l, Function[r, 0], {aR, rMax}, {Re[vR], Re[dR]}];
       err = Max[Abs[Table[If[r <= aL, s1[r], s3[r]] - phiA[r], {r, pts}]]]/Max[Abs[phiA /@ pts], 10^-30];
       cmp = Max[Abs[Table[s2[aGrid[[k]]] - ref[[k]], {k, Length[aGrid]}]]]/Max[Abs[ref]];
       shadow = Max[Abs[Table[Im[s2[aGrid[[k]]] - phiA[N[yC[aGrid[[k]]], 40]]], {k, Length[aGrid]}]]]/Max[Abs[ref]];
       AppendTo[shadowRows, {m, l, v, err, cmp, shadow}];
       csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> " puncture " <> v <> ": exterior regular field exact (scaled)", err, 10^-14];
       csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> " puncture " <> v <> ": contour residual = closed form " <> Switch[v, "ThetaH", "phiR~ - (I/Pi) kappa P~", "corrected", "phiR~", _, "phiR~ + kappa (f - phiP_n)~"] <> " (scaled, 25 points)", cmp, 10^-14];
       If[v === "ThetaH", csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": Theta->H puncture violates (ii): the contour residual carries an O(1) shadow (> 1e-3 scaled)", If[shadow > 10^-3, 0, 1], 0]];
       If[v === "corrected", csCheck["m=" <> mStr[m] <> " l=" <> ToString[l] <> ": corrected puncture: contour residual has no shadow", shadow, 10^-14]]],
      {v, names}];
    Print["    l=", l, ": residual solves done [", tick[], " s]"]], {l, If[smoke, {0, 2}, {0, 1, 2, 4, 6}]}];

 (* (d) smoothness of the residual modes [Phi1 - PhiP_n]_l' at r0: derivative jumps (out - in)^(k)(r0), k = 0..n+5 (m = 0: the residual is PhiR, nothing to test) *)
 If[m != 0, Do[Module[{e = residualExactPair[punctureCoeffs[n], lp], jumps, first, lead},
    jumps = Table[Abs[N[(D[e[[2]], {zVar, k}] - D[e[[1]], {zVar, k}]) /. zVar -> r0, 40]], {k, 0, n + 5}];
    first = Position[jumps, x_ /; x > 10^-20, 1, 1];
    AppendTo[smoothRows, {m, n, lp, Sequence @@ (sci /@ jumps)}];
    csCheck["m=" <> mStr[m] <> " n=" <> ToString[n] <> " l'=" <> ToString[lp] <> ": residual modes are C^(n+3) at r0 (derivatives 0..n+3 continuous, jump at order n+4)",
     If[Max[jumps[[1 ;; n + 4]]] < 10^-20 && jumps[[n + 5]] > 10^-20, 0, 1], 0]], {n, nList}, {lp, If[smoke, {2}, {0, 2, 5}]}]];
 Print["    smoothness done [", tick[], " s]"],
 {m, mList}];

(* (e) m = 0: the massless package's closed-form phiS_0 (Stage 1 of the massless test) and its Stage-2 numbers *)
Print["\n--- m = 0: closed-form phiS_0 of the massless Stage 1 ---"];
SetCMParameters[<|"m" -> 0|>];
Module[{l = 0, z, psiML, f, E0, P, Pc, u, a, pts, d1, d2},
  (* massless: psi_p - D_0 with psi_p = (I/(2 Pi)) (75 + z^2)/z + (5 I/Pi) Log[5 - z] - (25 I/Pi) Log[z - 5]/z (canonical logs), D_0 = 5 + 25/z *)
  psiML = Function[z, (I/(2 Pi)) (75 + z^2)/z + (5 I/Pi) (Log[z - 5] - I Pi) - (25 I/Pi) Log[z - 5]/z - 5 - 25/z];
  f = singularModeThetaH[0]; E0 = singularResidualE[0];
  csCheck["m=0 l=0: Theta->H mode = -r0 + (r0 - r0^2/z) H(z) (massless App. C form)", Max[Table[Abs[N[f[zz] - (-5 + (5 - 25/zz) HExt[zz]), 30]], {zz, {3, 7, 51/10 + I/2}}]], 10^-25];
  csCheck["m=0 l=0: E0 = -1/z, i.e. L[Theta->H mode] = deltaExt - (I/Pi)/z (massless App. C (ii))", Max[Table[Abs[N[E0[zz] + 1/zz, 30]], {zz, {3, 7, 51/10 + I/2}}]], 10^-25];
  (* corrected mode = Theta->H + (I/(2 Pi)) z up to I x (real homogeneous); compare with the massless closed form: real parts equal, Im difference homogeneous *)
  pts = {2, 3, 7, 8};
  d1 = Table[N[Re[psiML[r]] - Re[f[r] + I/(2 Pi) r], 30], {r, pts}];
  d2 = Table[N[Im[psiML[r]] - Im[f[r] + I/(2 Pi) r], 30], {r, pts}];
  csCheck["m=0 l=0: corrected mode has the real part of the massless closed-form phiS_0", Max[Abs[d1]], 10^-25];
  csCheck["m=0 l=0: Im difference to the massless phiS_0 is a homogeneous solution A + B/r", Module[{A, B}, {A, B} = {A, B} /. First@Solve[{A + B/2 == d2[[1]], A + B/3 == d2[[2]]}, {A, B}]; Max[Abs[d2 - (A + B/pts)]]], 10^-25];
  csCheck["m=0 l=0: L[Theta->H + (I/(2 Pi)) z] = deltaExt (the massless fix)", Max[Table[Abs[N[(radialL[0, f[z] + I/(2 Pi) z, z] - deltaExt[z]) /. z -> zz, 30]], {zz, {3, 51/10 + I/2}}]], 10^-25]];
csCheck["m=0: punctureCoeffs[n] reduces to q/rho (all mass terms vanish)", If[punctureCoeffs[5] === {{-1, 1}}, 0, 1], 0];

(* ================= tables ================= *)
md1 = StringRiffle[Join[{"| m | l | (i) max abs(Re psi - G_l) | (ii) max abs(Re(L psi - deltaExt)) | max abs(Im(L psi - deltaExt)) on axis | E0(r0) | corrected: max abs(L phiS - deltaExt) on contour | Im psi_ThetaH(aR) | Im psi_corr(aR) |", "|---|---|---|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> StringRiffle[sci /@ #[[{3, 4, 5}]], " | "] <> " | " <> ToString[#[[6]]] <> " | " <> StringRiffle[sci /@ #[[{7, 8, 9}]], " | "] <> " |") & /@ rows], "\n"];
md2 = StringRiffle[Join[{"| m | l | puncture | exterior phiR error (scaled) | contour residual vs closed form (scaled) | residual shadow on the contour (scaled) |", "|---|---|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> #[[3]] <> " | " <> sci[#[[4]]] <> " | " <> sci[#[[5]]] <> " | " <> sci[#[[6]]] <> " |") & /@ shadowRows], "\n"];
md3 = StringRiffle[Join[{"| m | n | l' | jumps of (out - in)^(k)(r0), k = 0 ... n+5 |", "|---|---|---|---|"},
    ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> ToString[#[[3]]] <> " | " <> StringRiffle[#[[4 ;;]], ", "] <> " |") & /@ smoothRows], "\n"];
md4 = StringRiffle[Join[{"| m | l | DSolve outcome | alpha (u_l) | beta (v_l) | fit residual |", "|---|---|---|---|---|---|"}, ("| " <> mStr[#[[1]]] <> " | " <> ToString[#[[2]]] <> " | " <> #[[3]] <> " | " <> ToString[#[[4]]] <> " | " <> ToString[#[[5]]] <> " | " <> ToString[#[[6]]] <> " |") & /@ dsolveRows], "
"];
Print["\n", md1, "\n\n", md2, "\n\n", md3, "\n\n", md4];
Export[FileNameJoin[{resDir, "massive_stage1_summary.md"}], "# Massive Stage 1: first order\n\n## Continued singular mode (T3 (i), (ii))\n\n" <> md1 <> "\n\n## Residual solves on the contour\n\n" <> md2 <>
   "\n\n## Smoothness of the residual modes at r0\n\n" <> md3 <> "\n\n## DSolve\n\n" <> md4 <> "\n", "Text"];
Export[FileNameJoin[{resDir, "massive_stage1_singular_mode.csv"}], Prepend[{mStr[#[[1]]], #[[2]], sci[#[[3]]], sci[#[[4]]], sci[#[[5]]], #[[6]], sci[#[[7]]], sci[#[[8]]], sci[#[[9]]]} & /@ rows, {"m", "l", "re_mismatch", "re_residual_axis", "im_residual_axis", "E0_r0", "corrected_residual_contour", "Im_psi_thetaH_aR", "Im_psi_corrected_aR"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage1_residual_solves.csv"}], Prepend[{mStr[#[[1]]], #[[2]], #[[3]], sci[#[[4]]], sci[#[[5]]], sci[#[[6]]]} & /@ shadowRows, {"m", "l", "puncture", "exterior_error", "contour_vs_closed_form", "contour_shadow"}], "CSV"];
Export[FileNameJoin[{resDir, "massive_stage1_smoothness.csv"}], Prepend[MapAt[mStr, #, 1] & /@ smoothRows, Join[{"m", "n", "lp"}, ("jump_k" <> ToString[#]) & /@ Range[0, 10]][[1 ;; Length[smoothRows[[1]]]]]], "CSV"];
Print["\nStage 1 wall time: ", tick[], " s"];
csSummary[];
exportChecks[FileNameJoin[{resDir, "massive_stage1_checks.csv"}]];
