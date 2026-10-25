(* ::Package:: *)

(* ::Title:: *)
(*ContourMassive` -- complex-contour regularisation, massive scalar toy model*)

(* ::Text:: *)
(*Adapted from ContourSecondOrder.wl (massless test, Claude Implementation/Mathematica implementation/).  Nothing is imported*)
(*from that folder; everything needed is copied here and modified for the screened operator (Laplacian - m^2).*)
(**)
(*Model.  (Lap - m^2) Phi1 = -4 Pi q delta^3(x - x0),  Phi1 = q Exp[-m rho]/rho + PhiR,  (Lap - m^2) PhiR = 0 (regular at 0).*)
(*        (Lap - m^2) Phi2 = (1/2) d_Z^2 (Phi1^2),  Z = r mu the Cartesian axis coordinate.  Exact: Phi2 = (1/2) d_Z^2 chi, (Lap - m^2) chi = Phi1^2.*)
(*Modes.  Phi = Sum phi_l(r) Y_l0,  [f]_l = Sqrt[Pi (2l+1)] Int f P_l dmu.  Radial operator L_l psi = psi'' + 2 psi'/r - (l(l+1)/r^2 + m^2) psi.*)
(*        Homogeneous solutions, normalised to the massless powers:  u_l(r) = (2l+1)!! i_l(m r)/m^l -> r^l,  v_l(r) = k_l(m r) m^(l+1)/(2l-1)!! -> r^(-l-1),*)
(*        with Wronskian u v' - u' v = -(2l+1)/r^2.  Green's function (regular at 0, decaying at infinity): G_l(r,s) = -(s^2/(2l+1)) u_l(r<) v_l(r>).*)
(*        [q Exp[-m rho]/rho]_L = q Sqrt[4 Pi/(2L+1)] u_L(r<) v_L(r>)  (the massive Coulomb mode; its derivative jump at r0 is m-independent).*)
(*        Regular field: PhiR = c0 u_0(r)/u_0(r0) + c1 r0 (u_1(r)/u_1(r0)) mu, so that PhiR(x0) = c0 + c1 r0 for every m and PhiR -> c0 + c1 r mu as m -> 0.*)
(*Contour.  z = y(a) = a + I amp (1 - ((a-r0)/radius)^2)^power on [aL, aR]; H(z) = 1 + (I/Pi) Log[z - r0]; deltaExt = I/(Pi (z - r0)).*)
(*          Real-axis values are boundary values from above.  Complex radial variable: zVar (symbolic) or z; Cartesian axis coordinate: Z or r mu, never z.*)
(*Logs.   Closed-form modes carry the symbol LOG for Log|r - r0|.  Substitute LOG -> Log[Abs[zVar - r0]] (real axis), LOG -> Log[zVar - r0] (out branch,*)
(*        continued) or LOG -> Log[zVar - r0] - I Pi (in branch, canonical: same analytic function in the upper half plane, from-above values on the axis).*)
(*        A Theta -> H blend in + (out - in) H of two closed-form branches must use logOut on both, so that out - in is real-analytic on the axis.*)
(*m = 0.  Every function falls back to the massless expressions, so m = 0 reproduces the massless test exactly.*)

Needs["NumericalDifferentialEquationAnalysis`"];
BeginPackage["ContourMassive`", {"NumericalDifferentialEquationAnalysis`"}];

(* ---------------- parameters and utilities ---------------- *)
$CMParameters::usage = "$CMParameters: association of model and numerical parameters (q, r0, m, rMin, rMax, radius, amp, power, c0, c1, wp, ...).";
SetCMParameters::usage = "SetCMParameters[assoc] merges assoc into the defaults, stores it and clears caches.";
csPar::usage = "csPar[key] returns a parameter.";
clearCMCaches::usage = "clearCMCaches[] clears memoised results.";
toWP::usage = "toWP[x] converts x to WorkingPrecision csPar[\"wp\"].";
sci::usage = "sci[x] formats a number in short scientific notation.";
csCheck::usage = "csCheck[name, value, tol] records and prints PASS/FAIL for |value| <= tol.";
csCheckList::usage = "csCheckList[] returns {name, value, tol, pass} for every recorded check.";
csResetChecks::usage = "csResetChecks[] clears the recorded checks.";
csSummary::usage = "csSummary[] prints the PASS/FAIL summary; True if all passed.";
exportChecks::usage = "exportChecks[file] writes the check list as CSV.";
scaledError::usage = "scaledError[num, ref] = Max|num - ref|/Max|ref| (absolute if ref is identically 0).";
zVar::usage = "zVar: the formal complex radial variable in symbolic expressions.";
LOG::usage = "LOG: placeholder for Log|r - r0| in closed-form modes (see the package header).";
logReal::usage = "logReal[expr] substitutes LOG -> Log[(zVar - r0)^2]/2 (= Log|r - r0| on the real axis, differentiable).";
logOut::usage = "logOut[expr] substitutes LOG -> Log[zVar - r0] (continued out branch).";
logIn::usage = "logIn[expr] substitutes LOG -> Log[zVar - r0] - I Pi (continued in branch, canonical: from-above values on r < r0).  NB: the Theta -> H blend in + (out - in) H needs logOut on BOTH branches (so that out - in is real-analytic on the axis); logOut[in] differs from logIn[in] by the pure shadow I Pi A.";
evR::usage = "evR[expr, r] evaluates an expression in zVar at r (70-digit input, 50-digit result).";
numericFunction::usage = "numericFunction[expr, var] returns a black-box numeric function of var at extra precision.";
contourFunction::usage = "contourFunction[expr, z] returns f[a] = expr at z = y(a), y(a) computed numerically first.";
dfact::usage = "dfact[n] = n!! for odd n >= -1.";
toFunction::usage = "toFunction[expr, var] = Function[var, expr].";

(* ---------------- contour, operators, homogeneous solutions ---------------- *)
tubeL::usage = "tubeL[] = r0 - radius."; tubeR::usage = "tubeR[] = r0 + radius.";
yC::usage = "yC[a]: the contour z = y(a)."; yCp::usage = "yCp[a] = y'(a)."; yCpp::usage = "yCpp[a] = y''(a).";
radialL::usage = "radialL[l, expr, z]: L_l (massive) applied to expr in z.";
realL::usage = "realL[l, u, r]: L_l u at r for NDSolve.";
contourL::usage = "contourL[l, u, a]: L_l in the contour parameter a.";
HExt::usage = "HExt[z] = 1 + (I/Pi) Log[z - r0]."; deltaExt::usage = "deltaExt[z] = I/(Pi (z - r0)).";
kappaL::usage = "kappaL[l] = -(q/r0^2) Sqrt[4 Pi (2l+1)].";
uFun::usage = "uFun[l, z] = (2l+1)!! i_l(m z)/m^l (-> z^l at m = 0): the homogeneous solution regular at 0.";
vFun::usage = "vFun[l, z] = k_l(m z) m^(l+1)/(2l-1)!! (-> z^(-l-1)): the decaying homogeneous solution.";
uFunD::usage = "uFunD[l, z] = d/dz uFun."; vFunD::usage = "vFunD[l, z] = d/dz vFun.";
greenPhys::usage = "greenPhys[l, r, s] = -(s^2/(2l+1)) u_l(r<) v_l(r>): Green's function with the physical conditions.";
coulombMassive::usage = "coulombMassive[L, r] = u_L(r<) v_L(r>) for real r (massive analogue of r<^L/r>^(L+1)).";
coulombMassivePair::usage = "coulombMassivePair[L] = {in, out} expressions in zVar for u_L(r<) v_L(r>).";
coulombMassiveCont::usage = "coulombMassiveCont[L, z] = in + (out - in) H(z): the Theta -> H continuation.";
coulombModeUnit::usage = "coulombModeUnit[l, r] = G_l(r, r0) = -(r0^2/(2l+1)) coulombMassive[l, r].";
singularModeThetaH::usage = "singularModeThetaH[l] = Function[z, G_l continued by Theta -> H] (correct real part; L psi - deltaExt is imaginary on the axis).";
singularResidualE::usage = "singularResidualE[l] = Function[z, E0(z)]: L[ThetaH mode] - deltaExt = (I/Pi) E0, E0 analytic and real on the axis.";

(* ---------------- regular field ---------------- *)
phiRMode::usage = "phiRMode[l] = Function[z, phiR_l(z)]: modes of the injected regular field (only l = 0, 1 nonzero).";
phiRAtParticle::usage = "phiRAtParticle[] = c0 + c1 r0.";
phiR3D::usage = "phiR3D[r, mu] = PhiR at (r, mu).";

(* ---------------- closed-form modes ---------------- *)
rhoModePair::usage = "rhoModePair[k, L] = {in, out} expressions in zVar for [rho^k]_L (k integer, k >= -2; LOG for the log at k = -2).";
expRhoModePair::usage = "expRhoModePair[p, c, L] = {in, out} expressions for [Exp[-c rho] rho^p]_L, p >= -2 (LOG carries the log when p = -2).";
eEnt::usage = "eEnt[w] = Ei(w) - EulerGamma - Log[w]: the entire part of the exponential integral (eEnt = eEntK[0]).";
eEntK::usage = "eEntK[k][w]: k-th derivative of eEnt (entire); derivatives of eEnt are returned in this form.";
ssUndiffPair::usage = "ssUndiffPair[coeffs, j] = {in, out} modes of (PhiP)^2 (coeffs = \"exact\" or a puncture coefficient list), LOG forms.";
srUndiffPair::usage = "srUndiffPair[coeffs, j, lpList] = {in, out} modes of 2 PhiP Phi_res with Phi_res modes from residualExactPair (coeffs = \"exact\": 2 q Exp[-m rho] PhiR/rho via Bessel couplings).";
rrUndiffPair::usage = "rrUndiffPair[coeffs, j, lpList] = {in, out} modes of (Phi_res)^2.";
chiGreen::usage = "chiGreen[j, srcPair, grid] = {chi, chi', chi'', chi'''} at the grid points for L_j chi = src (physical conditions), by Green's-function quadrature (NIntegrate, DoubleExponential, breakpoint at r0).";
exactReference::usage = "exactReference[lmax, grid, pieces] = <|l -> {values, derivatives}|> of Phi2_l = (1/2) d_Z^2 chi + A u_l + B v_l (physical conditions) on grid (which must contain rMin and rMax), chi from the undifferentiated pieces (\"SS\", \"SR\", \"RR\" with the exact puncture).";
muPowerMode::usage = "muPowerMode[l, m] = Sqrt[Pi(2l+1)] Int mu^m P_l dmu.";
gaunt::usage = "gaunt[l1, l2, l3] = Int Y_l1 Y_l2 Y_l3 dOmega (m = 0).";
Y0::usage = "Y0[l] = Y_l0(theta = 0) = Sqrt[(2l+1)/(4 Pi)].";
punctureCoeffs::usage = "punctureCoeffs[n] = list of {k, p_k} with PhiP_n = q Sum p_k rho^k: the series of Exp[-m rho]/rho through rho^n (n >= -1; n = -1 is q/rho).";
punctureModePair::usage = "punctureModePair[coeffs, L] = {in, out} modes of q Sum p_k rho^k (coeffs = {{k, p_k}, ...}).";
punctureModeCont::usage = "punctureModeCont[coeffs, L] = Theta -> H continuation of punctureModePair (each kinked term blended).";
exactPunctureModePair::usage = "exactPunctureModePair[L] = {in, out} modes of q Exp[-m rho]/rho (massive Coulomb mode).";
residualExactPair::usage = "residualExactPair[coeffs, lp] = {in, out} modes of Phi_res = PhiR + q (Exp[-m rho]/rho - Sum p_k rho^k).";
residualExactCont::usage = "residualExactCont[coeffs, lp] = Theta -> H continuation of residualExactPair.";

(* ---------------- d_Z ladder ---------------- *)
aZ::usage = "aZ[l] = l/Sqrt[(2l-1)(2l+1)].";
dzExpr::usage = "dzExpr[F, l] = (d_Z F)_l for mode expressions F[j] in zVar.";
dzzExpr::usage = "dzzExpr[F, l] = (d_Z^2 F)_l.";
dzzDeltaCoeff::usage = "dzzDeltaCoeff[Fpair, l]: coefficient of delta(r - r0) in the distributional (d_Z^2 F)_l for continuous modes Fpair[j] = {in, out}; pass the raw LOG forms (the log part, identical on both sides, is removed internally).";
prodD::usage = "prodD[a, b]: Leibniz product of derivative lists.";
dzArr::usage = "dzArr[Fl, z]: d_Z on lists of derivative lists Fl[[j+1]] = {F_j, F_j', ...} (valid to one order less).";
dzzArr::usage = "dzzArr[Fl, z] = dzArr[dzArr[Fl, z], z].";
dzzDeltaNum::usage = "dzzDeltaNum[l, dIn, dOut]: delta coefficient from one-sided first derivatives dIn[j], dOut[j].";

(* ---------------- 3D sources and complex-radius angular quadrature ---------------- *)
ssH3::usage = "ssH3[coeffs][z, mu, rho] = {h, d_r h, d_r^2 h} (fixed mu) for h = (PhiP)^2 with PhiP = q Sum p_k rho^k; coeffs = \"exact\" gives q^2 Exp[-2 m rho]/rho^2.";
srH3::usage = "srH3[coeffs, T][z, mu, rho] = {h, d_r h, d_r^2 h} for h = 2 PhiP T where T[z, mu] is a 3D function given as {T, d_r T, d_r^2 T}.";
angularNodes::usage = "angularNodes[nq] sets the Gauss rule (nq nodes per segment) used by the angular quadrature.";
branchModes::usage = "branchModes[z, br, h3, jm] = list over j <= jm of {M, M', M''}: the angular quadrature of h3 at complex z along the rho path of branch br (\"in\": via -I/2 to r0 - z; \"out\": via I/2 to z - r0).";
residueModes::usage = "residueModes[z, h3, jm]: residues at rho = 0 of the three rho-integrands (small-circle trapezoidal rule), per j.";
blendedModes::usage = "blendedModes[z, h3, jm] = {F~, F~', F~''} per j: Fin* + (Fout - Fin*) H with the H', H'' terms, where Fin* = (in-branch quadrature) + I Pi A puts the log of the in branch on the same branch as the out branch (A = residue at rho = 0); Fin* differs from the canonical in branch by the pure shadow I Pi A.";
quadModeDzz::usage = "quadModeDzz[l, z, h3] = (1/2) D_Z^2 of the blended quadrature modes: the continued source mode at z.";
quadModePlain::usage = "quadModePlain[l, z, h3, br] = (1/2) D_Z^2 of one branch alone (no blend).";

(* ---------------- Cauchy continuation of a real-axis function (from the massless Stage 8) ---------------- *)
chebX::usage = "chebX[a, b, n]: Chebyshev (first kind) nodes on [a, b].";
chebCoef::usage = "chebCoef[vals]: Chebyshev coefficients from values at chebX nodes (40-digit arithmetic).";
chebEval::usage = "chebEval[c, a, b, x]: Clenshaw evaluation.";
chebDeriv::usage = "chebDeriv[c]: coefficients of the derivative (on [-1, 1]).";
derivCoefs::usage = "derivCoefs[vals]: Chebyshev coefficients of s and its first three derivatives.";
cauchySetup::usage = "cauchySetup[{a, r0, b}, nG]: sets the Cauchy interval pieces and Gauss nodes; returns the association used by cauchyAt.";
cauchyAt::usage = "cauchyAt[cs, sG, dco, z]: (1/(Pi I)) Int s(t)/(t - z) dt with the cubic Taylor subtraction about Re z.";

(* ---------------- solvers ---------------- *)
ndOptions::usage = "ndOptions[]: NDSolve options.";
solveRealSegment::usage = "solveRealSegment[l, src, {r1, r2}, {v, d}]: L_l u = src on a real segment with Cauchy data at r1.";
solveContourSegment::usage = "solveContourSegment[l, srcA, {v, d}]: L_l u = srcA(a) along the contour with data at aL.";
threeRegionSolve::usage = "threeRegionSolve[l, srcReal, srcContour, {v1, d1}]: real [rMin, aL], contour, real [aR, rMax] (data Re u, Re u' at aR).";
threeRegionPhys::usage = "threeRegionPhys[l, srcReal, srcContour]: three-region solve with the physical conditions (zero data at rMin plus alpha u_l fixed by the Robin condition at rMax).";
realLinePhys::usage = "realLinePhys[l, srcIn, srcOut]: real-line solve split at r0 with the physical conditions.";
exteriorValue::usage = "exteriorValue[res, r]: value of a three-region solution at an exterior point.";
robinResidual::usage = "robinResidual[l, f, df, r]: the Robin-condition residual at r = rMin (u_l) or rMax (v_l).";
physFix::usage = "physFix[l, vals1, der1, vals10, der10]: {A, B} such that adding A u_l + B v_l meets the physical conditions.";
physResp::usage = "physResp[l, nodes, delta, grid]: Int G_l(r, s) delta(s) ds on the grid from quadrature nodes {s, w}.";
glNodes::usage = "glNodes[a, b, n]: Gauss-Legendre nodes and weights on [a, b] at 50 digits.";
panelNodes::usage = "panelNodes[breaks, n]: Gauss panels between consecutive breakpoints.";

Begin["`Private`"];

(* ================= parameters ================= *)
$CMDefaults = <|"q" -> 1, "r0" -> 5, "m" -> 1, "rMin" -> 1, "rMax" -> 10, "radius" -> 1/2, "amp" -> 1, "power" -> 6,
   "c0" -> 1, "c1" -> 1/10, "wp" -> 30, "precisionGoal" -> 20, "accuracyGoal" -> 35, "extraPrecision" -> 40|>;
$CMParameters = $CMDefaults;
$cache = <||>;
clearCMCaches[] := ($cache = <||>;);
SetCMParameters[new_Association] := ($CMParameters = Join[$CMDefaults, new]; clearCMCaches[]; $CMParameters);
csPar[k_String] := $CMParameters[k];
SetAttributes[memo, HoldRest];
memo[key_, expr_] := With[{k = {key, $CMParameters}}, If[KeyExistsQ[$cache, k], $cache[k], $cache[k] = expr]];

toWP[x_] := If[Precision[x] === Infinity, N[x, csPar["wp"]], SetPrecision[x, csPar["wp"]]];
dfact[n_Integer] := If[n <= 0, 1, n dfact[n - 2]];

sci[x_] := Module[{y = If[Precision[x] === Infinity, Quiet@N[Abs[x]], Abs[x]], e},
   Which[! NumericQ[y], ToString[x], Precision[x] === Infinity && x == 0, "0", y == 0 && MachineNumberQ[y], "0.(machine)",
    y == 0, "<1e" <> ToString[-Round[Accuracy[y]]], y < 10^-300, "<1e-300",
    True, e = Floor[Log10[y]];
    (If[Re[N[x]] < 0 && Im[N[x]] == 0, "-", ""]) <> ToString[NumberForm[SetPrecision[y/10^e, MachinePrecision], {3, 2}]] <> "e" <> ToString[e]]];
scaledError[num_List, ref_List] := With[{mx = Max[Abs[ref]]}, If[mx == 0, Max[Abs[num - ref]], Max[Abs[num - ref]]/mx]];

$checks = {};
csResetChecks[] := ($checks = {};);
csCheckList[] := $checks;
csCheck[name_String, value_, tol_] := Module[{v = Abs[value], pass},
   pass = TrueQ[NumericQ[v] && If[Precision[v] =!= Infinity && v == 0, 10^-Accuracy[v] <= tol || tol == 0 && Accuracy[v] > 300, v <= tol]];
   AppendTo[$checks, {name, v, tol, pass}];
   Print[If[pass, "PASS  ", "FAIL  "], name, ":  ", sci[v], "   (tol ", sci[tol], ")"]; pass];
csSummary[] := Module[{np = Count[$checks, {_, _, _, True}], nf = Count[$checks, {_, _, _, False}]},
   Print["------------------------------------------------------------"];
   Print["SUMMARY: ", np, " passed, ", nf, " failed (of ", Length[$checks], ")"];
   If[nf > 0, Print["FAILED: ", Column[Select[$checks, ! Last[#] &][[All, 1]]]]]; nf == 0];
exportChecks[file_String] := Export[file, Prepend[{#[[1]], sci[#[[2]]], sci[#[[3]]], If[#[[4]], "PASS", "FAIL"]} & /@ $checks, {"check", "value", "tol", "result"}], "CSV"];

toFunction[expr_, var_Symbol] := Function @@ {{var}, expr};
evR[e_, r_] := N[e /. zVar -> SetPrecision[r, 70], 50];
numericFunction[expr_, var_Symbol] := With[{prec = csPar["wp"] + csPar["extraPrecision"]}, Module[{f}, f[x_?NumericQ] := N[expr /. var -> SetPrecision[x, prec], prec]; f]];
contourFunction[expr_, var_Symbol] := With[{prec = csPar["wp"] + csPar["extraPrecision"]}, Module[{f}, f[x_?NumericQ] := N[expr /. var -> N[yC[SetPrecision[x, prec]], prec], prec]; f]];
logReal[e_] := e /. LOG -> Log[(zVar - csPar["r0"])^2]/2;        (* = Log|r - r0| for real r, and differentiable *)
logOut[e_] := e /. LOG -> Log[zVar - csPar["r0"]];
logIn[e_] := e /. LOG -> Log[zVar - csPar["r0"]] - I Pi;

(* ================= contour and operators ================= *)
tubeL[] := csPar["r0"] - csPar["radius"]; tubeR[] := csPar["r0"] + csPar["radius"];
yC[a_] := a + I csPar["amp"] (1 - ((a - csPar["r0"])/csPar["radius"])^2)^csPar["power"];
yCp[a_] := Module[{t}, D[yC[t], t] /. t -> a]; yCpp[a_] := Module[{t}, D[yC[t], {t, 2}] /. t -> a];
radialL[l_, expr_, z_] := D[expr, {z, 2}] + 2/z D[expr, z] - (l (l + 1)/z^2 + csPar["m"]^2) expr;
realL[l_, u_, r_] := u''[r] + (2/r) u'[r] - (l (l + 1)/r^2 + csPar["m"]^2) u[r];
contourL[l_, u_, a_] := (1/yCp[a]^2) u''[a] - (yCpp[a]/yCp[a]^3) u'[a] + (2/yC[a]) (1/yCp[a]) u'[a] - (l (l + 1)/yC[a]^2 + csPar["m"]^2) u[a];
HExt[z_] := 1 + I/Pi Log[z - csPar["r0"]]; deltaExt[z_] := I/(Pi (z - csPar["r0"]));
kappaL[l_] := -(csPar["q"]/csPar["r0"]^2) Sqrt[4 Pi (2 l + 1)];

(* homogeneous solutions; i_l(x) = Sqrt[Pi/(2x)] I_(l+1/2)(x), k_l(x) = Sqrt[2/(Pi x)] K_(l+1/2)(x) *)
uFun[l_, z_] := With[{m = csPar["m"]}, If[m == 0, z^l, dfact[2 l + 1] Sqrt[Pi/(2 m z)] BesselI[l + 1/2, m z]/m^l]];
vFun[l_, z_] := With[{m = csPar["m"]}, If[m == 0, z^(-l - 1), Sqrt[2/(Pi m z)] BesselK[l + 1/2, m z] m^(l + 1)/dfact[2 l - 1]]];
uFunD[l_, z_] := Module[{t}, D[uFun[l, t], t] /. t -> z]; vFunD[l_, z_] := Module[{t}, D[vFun[l, t], t] /. t -> z];
greenPhys[l_, r_, s_] := -(s^2/(2 l + 1)) If[s < r, uFun[l, s] vFun[l, r], uFun[l, r] vFun[l, s]];
coulombMassive[L_, r_] := With[{r0 = csPar["r0"]}, If[r < r0, uFun[L, r] vFun[L, r0], uFun[L, r0] vFun[L, r]]];
coulombMassivePair[L_] := With[{r0 = csPar["r0"]}, {uFun[L, zVar] vFun[L, r0], uFun[L, r0] vFun[L, zVar]}];
coulombMassiveCont[L_, z_] := With[{p = coulombMassivePair[L] /. zVar -> z}, p[[1]] + (p[[2]] - p[[1]]) HExt[z]];
coulombModeUnit[l_, r_] := -(csPar["r0"]^2/(2 l + 1)) coulombMassive[l, r];
singularModeThetaH[l_] := toFunction[-(csPar["r0"]^2/(2 l + 1)) coulombMassiveCont[l, zVar], zVar];
(* L[F + J H] - deltaExt with F, J homogeneous, J(r0) = 0, J'(r0) = 1 (unit derivative jump): = (I/Pi) E0, E0 = (2J' - 1)/(z - r0) - J/(z - r0)^2 + 2J/(z (z - r0)) *)
singularResidualE[l_] := Module[{p = coulombMassivePair[l], J, r0 = csPar["r0"]},
   J = -(r0^2/(2 l + 1)) (p[[2]] - p[[1]]);
   toFunction[(2 D[J, zVar] - 1)/(zVar - r0) - J/(zVar - r0)^2 + 2 J/(zVar (zVar - r0)), zVar]];

(* ================= regular field ================= *)
phiRMode[l_Integer] := With[{r0 = csPar["r0"], c0 = csPar["c0"], c1 = csPar["c1"]},
   Which[l == 0, toFunction[c0 Sqrt[4 Pi] uFun[0, zVar]/uFun[0, r0], zVar],
    l == 1, toFunction[c1 r0 Sqrt[4 Pi/3] uFun[1, zVar]/uFun[1, r0], zVar], True, toFunction[0, zVar]]];
phiRAtParticle[] := csPar["c0"] + csPar["c1"] csPar["r0"];
phiR3D[r_, mu_] := With[{r0 = csPar["r0"]}, csPar["c0"] uFun[0, r]/uFun[0, r0] + csPar["c1"] r0 (uFun[1, r]/uFun[1, r0]) mu];

(* ================= closed-form modes ================= *)
Y0[l_] := Sqrt[(2 l + 1)/(4 Pi)];
gaunt[l1_, l2_, l3_] := If[Abs[l1 - l2] <= l3 <= l1 + l2 && EvenQ[l1 + l2 + l3], Sqrt[(2 l1 + 1) (2 l2 + 1) (2 l3 + 1)/(4 Pi)] ThreeJSymbol[{l1, 0}, {l2, 0}, {l3, 0}]^2, 0];
muPowerMode[l_, mm_] := Sqrt[Pi (2 l + 1)] Integrate[mu^mm LegendreP[l, mu], {mu, -1, 1}];
(* eEnt = eEntK[0]; eEntK[k] is its k-th derivative (entire).  Near 0 the Taylor series is used; elsewhere Ei - EulerGamma - Log (with the
   value on the negative real axis made continuous from above) and D[(Exp[w] - 1)/w, {w, k - 1}]. *)
eEnt = eEntK[0];
eEntK /: Derivative[n_Integer][eEntK[k_]] := eEntK[k + n];
eEntClosed[k_] := eEntClosed[k] = Module[{w}, Function @@ {{w}, D[(Exp[w] - 1)/w, {w, k - 1}]}];
eEntK[k_Integer][w_?NumericQ] := Which[w == 0, If[k == 0, 0, 1/k],
   Abs[w] < 1/2, With[{pr = If[Precision[w] === Infinity, 60, Precision[w] + 5]}, Sum[N[w, pr]^(n - k)/(n (n - k)!), {n, Max[1, k], 70}]],
   k == 0, ExpIntegralEi[w] - EulerGamma - Log[w] + If[TrueQ[Im[w] == 0 && Re[w] < 0], I Pi, 0],
   True, eEntClosed[k][w]];
(* Antiderivative of rho^n Exp[-c rho] (n integer >= -1), as a function of rho; the log at n = -1 is returned as the symbol LOG plus Log[c] *)
expAnti[n_, c_, rho_] := Which[n == -1, EulerGamma + Log[c] + LOG + eEnt[-c rho], n >= 0, -(n!/c^(n + 1)) Exp[-c rho] Sum[(c rho)^i/i!, {i, 0, n}]];
expAntiOuter[n_, c_, rho_] := Which[n == -1, EulerGamma + Log[c] + Log[rho] + eEnt[-c rho], n >= 0, -(n!/c^(n + 1)) Exp[-c rho] Sum[(c rho)^i/i!, {i, 0, n}]];     (* at rho = z + r0, Re > 0: analytic (ExpIntegralEi itself sits on its cut there) *)
(* generic: [f(rho)]_L with f = u^kk (u = rho^2) via the u-substitution; uAnti gives the antiderivative in w = rho of 2 rho^(2j+1) f(rho^2) *)
powAnti[kk_, j_, w_] := If[j + kk == -1, 2 LOG, w^(2 (j + kk + 1))/(j + kk + 1)];
powAntiOuter[kk_, j_, w_] := If[j + kk == -1, 2 Log[w], w^(2 (j + kk + 1))/(j + kk + 1)];
modeBySubst[L_, antiIn_, antiOut_] := Module[{mu, u, a, b, coeffs, r0 = csPar["r0"], s, t},
   a = zVar^2 + r0^2; b = 2 zVar r0;
   coeffs = CoefficientList[Expand[LegendreP[L, mu] /. mu -> (a - u)/b], u];
   t = zVar + r0;
   Table[Sqrt[Pi (2 L + 1)] (1/b) Sum[coeffs[[j + 1]] (antiOut[j, t] - antiIn[j, If[br == 1, r0 - zVar, zVar - r0]]), {j, 0, Length[coeffs] - 1}], {br, 2}]];
rhoModePair[k_Integer, L_Integer] := memo[{"rhoMode", k, L}, Expand /@ modeBySubst[L, powAnti[k/2, #1, #2] &, powAntiOuter[k/2, #1, #2] &]];
expRhoModePair[p_Integer, c_, L_Integer] := memo[{"expRhoMode", p, c, L},
   modeBySubst[L, 2 expAnti[2 #1 + p + 1, c, #2] &, 2 expAntiOuter[2 #1 + p + 1, c, #2] &]];
(* puncture: PhiP_n = q Sum_{k=-1}^{n} (-m)^(k+1)/(k+1)! rho^k; at m = 0 only q/rho survives *)
punctureCoeffs[n_Integer] := With[{m = csPar["m"]}, Select[Table[{k, If[k == -1, 1, (-m)^(k + 1)/(k + 1)!]}, {k, -1, n}], #[[2]] != 0 &]];
punctureModePair[coeffs_List, L_Integer] := csPar["q"] Sum[c[[2]] rhoModePair[c[[1]], L], {c, coeffs}];
contPair[{in_, out_}, z_] := (in /. zVar -> z) + ((out - in) /. zVar -> z) HExt[z];
punctureModeCont[coeffs_List, L_Integer] := contPair[punctureModePair[coeffs, L], zVar];
exactPunctureModePair[L_Integer] := csPar["q"] Sqrt[4 Pi/(2 L + 1)] coulombMassivePair[L];
residualExactPair[coeffs_List, lp_Integer] := ConstantArray[phiRMode[lp][zVar], 2] + exactPunctureModePair[lp] - punctureModePair[coeffs, lp];
residualExactCont[coeffs_List, lp_Integer] := contPair[residualExactPair[coeffs, lp], zVar];

(* ================= d_Z ladder ================= *)
aZ[l_] := If[l <= 0, 0, l/Sqrt[(2 l - 1) (2 l + 1)]];
dzExpr[F_, l_] := If[l < 0, 0, If[l >= 1, aZ[l] (D[F[l - 1], zVar] - (l - 1) F[l - 1]/zVar), 0] + aZ[l + 1] (D[F[l + 1], zVar] + (l + 2) F[l + 1]/zVar)];
dzzExpr[F_, l_] := dzExpr[Function[j, dzExpr[F, j]], l];
(* The log term C(r) Log|r - r0| of a mode is the same analytic C on both sides and carries no delta, so LOG (or an explicit Log[zVar - r0])
   is removed before the branches are differenced; only the kinks of the regular parts contribute. *)
dzzDeltaCoeff[Fpair_, l_] := Module[{J, r0 = csPar["r0"], F},
   F[i_] := If[i < 0, {0, 0}, Fpair[i] /. {LOG -> 0, Log[zVar - r0] -> 0}];
   J[j_] := If[j < 0, 0, Module[{e = {dzExpr[Function[i, F[i][[1]]], j], dzExpr[Function[i, F[i][[2]]], j]}},
      Re[N[(e[[2]] - e[[1]]) /. zVar -> r0, 60]]]];
   If[l >= 1, aZ[l] J[l - 1], 0] + aZ[l + 1] J[l + 1]];
prodD[a_, b_] := Table[Sum[Binomial[mm, i] a[[i + 1]] b[[mm - i + 1]], {i, 0, mm}], {mm, 0, Length[a] - 1}];
invZD[p_, z_] := (-1)^p p!/z^(p + 1);
dzArr[Fl_, z_] := Module[{nd = Length[Fl[[1]]] - 1, f},
   f[j_, mm_] := If[j < 0 || j > Length[Fl] - 1, 0, Fl[[j + 1, mm + 1]]];
   Table[Table[aZ[j] (f[j - 1, mm + 1] - (j - 1) Sum[Binomial[mm, i] f[j - 1, i] invZD[mm - i, z], {i, 0, mm}]) +
      aZ[j + 1] (f[j + 1, mm + 1] + (j + 2) Sum[Binomial[mm, i] f[j + 1, i] invZD[mm - i, z], {i, 0, mm}]), {mm, 0, nd - 1}], {j, 0, Length[Fl] - 2}]];
dzzArr[Fl_, z_] := dzArr[dzArr[Fl, z], z];
dzzDeltaNum[l_, dIn_, dOut_] := Module[{J},
   J[j_] := If[j < 0, 0, If[j >= 1, aZ[j] (dOut[j - 1] - dIn[j - 1]), 0] + aZ[j + 1] (dOut[j + 1] - dIn[j + 1])];
   If[l >= 1, aZ[l] J[l - 1], 0] + aZ[l + 1] J[l + 1]];

(* ================= 3D functions and their r-derivatives at fixed mu ================= *)
(* for f(rho): d_r f = f'(rho) (r - r0 mu)/rho; d_r^2 f = f'' (r - r0 mu)^2/rho^2 + f' (1/rho - (r - r0 mu)^2/rho^3) *)
rhoDerivs[f_, z_, mu_, rho_] := With[{w = z - csPar["r0"] mu}, {f[[1]], f[[2]] w/rho, f[[3]] w^2/rho^2 + f[[2]] (1/rho - w^2/rho^3)}];
polyRho[coeffs_, rho_] := {Sum[c[[2]] rho^c[[1]], {c, coeffs}], Sum[c[[2]] c[[1]] rho^(c[[1]] - 1), {c, coeffs}], Sum[c[[2]] c[[1]] (c[[1]] - 1) rho^(c[[1]] - 2), {c, coeffs}]};
ssH3["exact"] := Function[{z, mu, rho}, With[{m = csPar["m"], q = csPar["q"]},
    rhoDerivs[q^2 {Exp[-2 m rho]/rho^2, -Exp[-2 m rho] (2 m/rho^2 + 2/rho^3), Exp[-2 m rho] (4 m^2/rho^2 + 8 m/rho^3 + 6/rho^4)}, z, mu, rho]]];
ssH3[coeffs_List] := Function[{z, mu, rho}, With[{q = csPar["q"]}, Module[{p = polyRho[coeffs, rho]},
     rhoDerivs[q^2 {p[[1]]^2, 2 p[[1]] p[[2]], 2 p[[2]]^2 + 2 p[[1]] p[[3]]}, z, mu, rho]]]];
pH3[coeffs_, z_, mu_, rho_] := If[coeffs === "exact", With[{m = csPar["m"]}, rhoDerivs[csPar["q"] {Exp[-m rho]/rho, -Exp[-m rho] (m/rho + 1/rho^2), Exp[-m rho] (m^2/rho + 2 m/rho^2 + 2/rho^3)}, z, mu, rho]],
   rhoDerivs[csPar["q"] polyRho[coeffs, rho], z, mu, rho]];
srH3[coeffs_, T_] := Function[{z, mu, rho}, Module[{P = pH3[coeffs, z, mu, rho], t = T[z, mu]},
    2 {P[[1]] t[[1]], P[[2]] t[[1]] + P[[1]] t[[2]], P[[3]] t[[1]] + 2 P[[2]] t[[2]] + P[[1]] t[[3]]}]];

(* ================= complex-radius angular quadrature ================= *)
nq = 80; {tq, wq} = Transpose[GaussianQuadratureWeights[80, 0, 1, 40]];
angularNodes[n_Integer] := (nq = n; {tq, wq} = Transpose[GaussianQuadratureWeights[n, 0, 1, 40]];);
muOf[z_, rho_] := (z^2 + csPar["r0"]^2 - rho^2)/(2 z csPar["r0"]);
pathPts[z_, br_] := With[{via = If[br === "in", -I/2, I/2], e = If[br === "in", csPar["r0"] - z, z - csPar["r0"]], r0 = csPar["r0"]},
   Join[Transpose[{z + r0 + (via - z - r0) tq, (via - z - r0) wq}], Transpose[{via + (e - via) tq, (e - via) wq}]]];
projModes[z_, pts_, h3_, jm_] := Module[{rho = pts[[All, 1]], w = pts[[All, 2]], mu, vals, leg, r0 = csPar["r0"]},
   mu = muOf[z, rho];
   vals = (# (-rho/(z r0)) w) & /@ h3[z, mu, rho];
   leg = Table[LegendreP[j, mu], {j, 0, jm}];
   Table[Sqrt[Pi (2 j + 1)] (leg[[j + 1]] . #) & /@ vals, {j, 0, jm}]];
branchModes[z_, br_, h3_, jm_] := projModes[z, pathPts[z, br], h3, jm];
nRes = 64;
residueModes[z_, h3_, jm_] := Module[{rho = (1/2) Exp[2 Pi I (Range[nRes] - 1/2)/nRes], w},
   w = rho/nRes;                                       (* (1/(2 Pi I)) Int g drho = (1/N) Sum g(rho_k) rho_k *)
   projModes[z, Transpose[{rho, w}], h3, jm]];
blendedModes[z_, h3_, jm_] := Module[{pin = branchModes[z, "in", h3, jm], pout = branchModes[z, "out", h3, jm], A = residueModes[z, h3, jm],
    H0 = HExt[z], H1 = deltaExt[z], H2 = -I/(Pi (z - csPar["r0"])^2)},
   Table[Module[{a = pin[[j + 1]] + I Pi A[[j + 1]], d}, d = pout[[j + 1]] - a;
     {a[[1]] + d[[1]] H0, a[[2]] + d[[2]] H0 + d[[1]] H1, a[[3]] + d[[3]] H0 + 2 d[[2]] H1 + d[[1]] H2}], {j, 0, jm}]];
quadModeDzz[l_, z_, h3_] := (1/2) dzzArr[blendedModes[z, h3, l + 2], z][[l + 1, 1]];
quadModePlain[l_, z_, h3_, br_] := (1/2) dzzArr[branchModes[z, br, h3, l + 2], z][[l + 1, 1]];

(* ================= Cauchy continuation (massless Stage 8) ================= *)
chebX[a_, b_, nn_] := N[(a + b)/2 + (b - a)/2 Cos[Pi (Range[0, nn] + 1/2)/(nn + 1)], 50];
chebCoef[f_List] := Block[{$MinPrecision = 40, $MaxPrecision = 40}, Module[{nn = Length[f] - 1, c, ff = SetPrecision[f, 40]},
   c = Table[2/(nn + 1) Sum[ff[[k + 1]] Cos[N[j Pi (k + 1/2)/(nn + 1), 40]], {k, 0, nn}], {j, 0, nn}]; c[[1]] = c[[1]]/2; c]];
chebEval[c_, a_, b_, z_] := Block[{$MinPrecision = 40, $MaxPrecision = 40}, Module[{x = SetPrecision[(2 z - a - b)/(b - a), 40], b1 = 0, b2 = 0, tmp},
   Do[tmp = 2 x b1 - b2 + c[[j + 1]]; b2 = b1; b1 = tmp, {j, Length[c] - 1, 1, -1}]; x b1 - b2 + c[[1]]]];
chebDeriv[c_] := Module[{nn = Length[c] - 1, d}, d = ConstantArray[0, nn + 1];
   If[nn >= 1, d[[nn]] = 2 nn c[[nn + 1]]; Do[d[[j]] = If[j + 2 <= nn + 1, d[[j + 2]], 0] + 2 j c[[j + 1]], {j, nn - 1, 1, -1}]; d[[1]] = d[[1]]/2]; Most[d]];
derivCoefs[vals_] := Module[{c = chebCoef[vals]}, NestList[chebDeriv, c, 3]];
glNodes[a_, b_, n_] := GaussianQuadratureWeights[n, a, b, 50];
panelNodes[breaks_List, n_] := Join @@ Table[glNodes[breaks[[k]], breaks[[k + 1]], n], {k, Length[breaks] - 1}];
cauchySetup[{a_, r0_, b_}, nG_] := Module[{pieces = {{a, r0}, {r0, b}}}, <|"pieces" -> pieces, "gauss" -> (glNodes[#[[1]], #[[2]], nG] & /@ pieces)|>];
polyCauchy[d_, x_, z_, {a_, b_}] := Module[{w = z - x, lg = Log[b - z - I 10^-60] - Log[a - z - I 10^-60]},
   Sum[d[[i + 1]] Sum[Binomial[i, k] w^(i - k) If[k == 0, lg, ((b - z)^k - (a - z)^k)/k], {k, 0, i}], {i, 0, Length[d] - 1}]];
cauchyAt[cs_, sG_, dco_, z_] := Module[{x = Re[z], tot = 0, pieces = cs["pieces"], gaussP = cs["gauss"]},
   Do[Module[{ab = pieces[[k]], tj = gaussP[[k, All, 1]], wj = gaussP[[k, All, 2]]},
      If[ab[[1]] < x < ab[[2]],
       Module[{dd = Table[chebEval[dco[[k, mm + 1]], ab[[1]], ab[[2]], x] (2/(ab[[2]] - ab[[1]]))^mm/mm!, {mm, 0, 3}], P},
        P = Sum[dd[[mm + 1]] (tj - x)^mm, {mm, 0, 3}];
        tot += Total[wj (sG[[k]] - P)/(tj - z)] + polyCauchy[dd, x, z, ab]],
       tot += Total[wj sG[[k]]/(tj - z)]]], {k, 2}];
   tot/(Pi I)];

(* ================= solvers ================= *)
ndOptions[] := {WorkingPrecision -> csPar["wp"], PrecisionGoal -> csPar["precisionGoal"], AccuracyGoal -> csPar["accuracyGoal"], MaxSteps -> 10^6, InterpolationOrder -> All};
solveRealSegment[l_, src_, {r1_, r2_}, {v_, d_}] := Module[{u, r},
   NDSolveValue[{realL[l, u, r] == src[r], u[r1] == toWP[v], u'[r1] == toWP[d]}, u, {r, Min[r1, r2], Max[r1, r2]}, Evaluate[Sequence @@ ndOptions[]]]];
solveContourSegment[l_, srcA_, {v_, d_}] := Module[{u, a},
   NDSolveValue[{contourL[l, u, a] == srcA[a], u[tubeL[]] == toWP[v], u'[tubeL[]] == toWP[d]}, u, {a, tubeL[], tubeR[]}, Evaluate[Sequence @@ ndOptions[]]]];
threeRegionSolve[l_, srcReal_, srcContour_, {v1_, d1_}] := Module[{s1, s2, s3, vR, dR, aL = tubeL[], aR = tubeR[]},
   s1 = solveRealSegment[l, srcReal, {csPar["rMin"], aL}, {v1, d1}];
   s2 = solveContourSegment[l, srcContour, {s1[aL], s1'[aL]}];
   vR = s2[aR]; dR = s2'[aR];
   s3 = solveRealSegment[l, srcReal, {aR, csPar["rMax"]}, {Re[vR], Re[dR]}];
   <|"sol1" -> s1, "sol2" -> s2, "sol3" -> s3, "valR" -> vR, "derR" -> dR, "imAR" -> {Im[vR], Im[dR]}|>];
exteriorValue[res_Association, r_] := If[r <= tubeL[], res["sol1"][r], res["sol3"][r]];
robinResidual[l_, f_, df_, r_] := If[r == csPar["rMin"], df - f uFunD[l, r]/uFun[l, r], df - f vFunD[l, r]/vFun[l, r]];
(* zero data at rMin satisfy the Robin condition there; add alpha u_l to meet the one at rMax *)
alphaFromRobin[l_, f10_, df10_] := Module[{rM = csPar["rMax"], lv}, lv = N[vFunD[l, rM]/vFun[l, rM], 50]; (lv f10 - df10)/(N[uFunD[l, rM], 50] - lv N[uFun[l, rM], 50])];
threeRegionPhys[l_, srcReal_, srcContour_] := Module[{res = threeRegionSolve[l, srcReal, srcContour, {0, 0}], al, rM = csPar["rMax"], uf},
   al = alphaFromRobin[l, res["sol3"][rM], res["sol3"]'[rM]];
   uf = Function[r, Evaluate[al uFun[l, r]]];
   Join[res, <|"alpha" -> al, "ext" -> Function[r, exteriorValue[res, r] + al uFun[l, r]], "contour" -> Function[a, res["sol2"][a] + al uFun[l, yC[a]]]|>]];
realLinePhys[l_, srcIn_, srcOut_] := Module[{s1, s2, al, rM = csPar["rMax"], r0 = csPar["r0"]},
   s1 = solveRealSegment[l, srcIn, {csPar["rMin"], r0}, {0, 0}];
   s2 = solveRealSegment[l, srcOut, {r0, rM}, {s1[r0], s1'[r0]}];
   al = alphaFromRobin[l, s2[rM], s2'[rM]];
   <|"f" -> Function[r, If[r <= r0, s1[r], s2[r]] + al uFun[l, r]], "d" -> Function[r, If[r <= r0, s1'[r], s2'[r]] + al uFunD[l, r]], "alpha" -> al|>];
physFix[l_, f1_, df1_, f10_, df10_] := Module[{A, B, r1 = csPar["rMin"], r10 = csPar["rMax"], sol},
   sol = First@Solve[{df1 + A uFunD[l, r1] + B vFunD[l, r1] == (f1 + A uFun[l, r1] + B vFun[l, r1]) uFunD[l, r1]/uFun[l, r1],
       df10 + A uFunD[l, r10] + B vFunD[l, r10] == (f10 + A uFun[l, r10] + B vFun[l, r10]) vFunD[l, r10]/vFun[l, r10]}, {A, B}];
   {A, B} /. sol];
physResp[l_, nodes_, delta_, grid_] := Module[{s = nodes[[All, 1]], w = nodes[[All, 2]], us, vs},
   us = uFun[l, #] & /@ s; vs = vFun[l, #] & /@ s;
   Table[-(1/(2 l + 1)) (vFun[l, r] Total[w s^2 us delta UnitStep[r - s]] + uFun[l, r] Total[w s^2 vs delta (1 - UnitStep[r - s])]), {r, grid}]];

(* ================= undifferentiated second-order source modes ================= *)
ssUndiffPair["exact", j_Integer] := memo[{"ssU", j}, csPar["q"]^2 If[csPar["m"] == 0, rhoModePair[-2, j], expRhoModePair[-2, 2 csPar["m"], j]]];      (* m = 0: the massless [1/rho^2]_j *)
ssUndiffPair[coeffs_List, j_Integer] := memo[{"ssUP", coeffs, j}, csPar["q"]^2 Sum[c1[[2]] c2[[2]] rhoModePair[c1[[1]] + c2[[1]], j], {c1, coeffs}, {c2, coeffs}]];
srUndiffPair["exact", j_Integer, lpList_: {0, 1}] := memo[{"srU", j, lpList},
   2 csPar["q"] Sum[Sum[gaunt[j, lp, L] Sqrt[4 Pi/(2 L + 1)] coulombMassivePair[L], {L, Abs[j - lp], j + lp}] phiRMode[lp][zVar], {lp, lpList}]];
srUndiffPair[coeffs_List, j_Integer, lpList_] := memo[{"srUP", coeffs, j, lpList},
   2 Sum[Sum[gaunt[j, lp, L] punctureModePair[coeffs, L], {L, Abs[j - lp], j + lp}] residualExactPair[coeffs, lp], {lp, lpList}]];
rrUndiffPair["exact", j_Integer, lpList_: {0, 1}] := memo[{"rrU", j, lpList},
   ConstantArray[Sum[gaunt[j, l1, l2] phiRMode[l1][zVar] phiRMode[l2][zVar], {l1, lpList}, {l2, lpList}], 2]];
rrUndiffPair[coeffs_List, j_Integer, lpList_] := memo[{"rrUP", coeffs, j, lpList},
   Sum[If[gaunt[j, l1, l2] == 0, 0, gaunt[j, l1, l2] residualExactPair[coeffs, l1] residualExactPair[coeffs, l2]], {l1, lpList}, {l2, lpList}]];

(* ================= exact reference by Green's-function quadrature ================= *)
(* chi_j(r) = -(1/(2j+1)) [v_j(r) I1(r) + u_j(r) I2(r)],  I1 = Int_rMin^r s^2 u_j S, I2 = Int_r^rMax s^2 v_j S; cumulative over the sorted grid *)
(* the double-exponential rule handles the logarithmic end-point singularity at r0 but reports slow convergence there; its accuracy is
   established by the cross-checks of Stage 2 (independent rule, 3D closed form, NDSolve), so the messages are silenced *)
nint[f_, a_, b_] := If[b <= a, 0, Module[{s, r0 = csPar["r0"]},
    Quiet[NIntegrate[f[s], Evaluate[{s, Sequence @@ Union[{a, b}, Select[{r0}, a < # < b &]]}], Method -> "DoubleExponential", WorkingPrecision -> 40, PrecisionGoal -> 28, MaxRecursion -> 14], {NIntegrate::eincr, NIntegrate::slwcon, NIntegrate::ncvi, NIntegrate::ncvb}]]];
chiGreen[j_, {in_, out_}, grid_] := Module[{r0 = csPar["r0"], fIn, fOut, S, dS, g = Sort[grid], I1, I2, u, v, du, dv, chi, dchi, m = csPar["m"], rMin = csPar["rMin"], rMax = csPar["rMax"]},
   If[in === 0 && out === 0, Return[ConstantArray[{0, 0, 0, 0}, Length[g]]]];
   fIn = numericFunction[logReal[in], zVar]; fOut = numericFunction[logReal[out], zVar];
   S[r_?NumericQ] := If[r < r0, fIn[r], fOut[r]];
   dS = Module[{dIn = numericFunction[D[logReal[in], zVar], zVar], dOut = numericFunction[D[logReal[out], zVar], zVar]}, Function[r, If[r < r0, dIn[r], dOut[r]]]];
   u = Function[r, Evaluate[N[uFun[j, r], 50]]]; v = Function[r, Evaluate[N[vFun[j, r], 50]]];
   I1 = Accumulate[Prepend[Table[nint[Function[s, s^2 u[s] S[s]], g[[k]], g[[k + 1]]], {k, Length[g] - 1}], nint[Function[s, s^2 u[s] S[s]], rMin, g[[1]]]]];
   I2 = Reverse[Accumulate[Reverse[Append[Table[nint[Function[s, s^2 v[s] S[s]], g[[k]], g[[k + 1]]], {k, Length[g] - 1}], nint[Function[s, s^2 v[s] S[s]], g[[-1]], rMax]]]]];
   du = Function[r, Evaluate[N[uFunD[j, r], 50]]]; dv = Function[r, Evaluate[N[vFunD[j, r], 50]]];
   chi = Table[-(1/(2 j + 1)) (v[g[[k]]] I1[[k]] + u[g[[k]]] I2[[k]]), {k, Length[g]}];
   dchi = Table[-(1/(2 j + 1)) (dv[g[[k]]] I1[[k]] + du[g[[k]]] I2[[k]]), {k, Length[g]}];
   Table[With[{r = g[[k]], c = chi[[k]], dc = dchi[[k]]}, Module[{d2, d3},
      d2 = S[r] - 2 dc/r + (j (j + 1)/r^2 + m^2) c;
      d3 = dS[r] - 2 d2/r + 2 dc/r^2 + (j (j + 1)/r^2 + m^2) dc - 2 j (j + 1) c/r^3;
      {c, dc, d2, d3}]], {k, Length[g]}]];
exactReference[lmax_, grid_, pieces_: {"SS", "SR", "RR"}] := Module[{g = Sort[grid], chis, res = <||>, p1, p10, AB, rMin = csPar["rMin"], rMax = csPar["rMax"]},
   chis = Table[chiGreen[j, Sum[Switch[p, "SS", ssUndiffPair["exact", j], "SR", srUndiffPair["exact", j], "RR", rrUndiffPair["exact", j]], {p, pieces}], g], {j, 0, lmax + 2}];
   p1 = Position[g, rMin][[1, 1]]; p10 = Position[g, rMax][[1, 1]];
   Do[Module[{vals = Table[(1/2) dzzArr[chis[[All, k]], g[[k]]][[l + 1]], {k, Length[g]}]},
     AB = physFix[l, vals[[p1, 1]], vals[[p1, 2]], vals[[p10, 1]], vals[[p10, 2]]];
     (* the reference is real; the closed-form sources carry inexact zero imaginary parts, which are dropped here *)
     res[l] = <|"grid" -> g, "val" -> Table[Re[vals[[k, 1]] + AB[[1]] uFun[l, g[[k]]] + AB[[2]] vFun[l, g[[k]]]], {k, Length[g]}],
       "der" -> Table[Re[vals[[k, 2]] + AB[[1]] uFunD[l, g[[k]]] + AB[[2]] vFunD[l, g[[k]]]], {k, Length[g]}], "AB" -> Re[AB]|>], {l, 0, lmax}];
   res];

End[];
EndPackage[];
