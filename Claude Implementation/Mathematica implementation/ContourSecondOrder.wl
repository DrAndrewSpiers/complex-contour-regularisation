(* ::Package:: *)

(* ::Title:: *)
(*ContourSecondOrder` -- complex-contour S/R split, first and second order (toy model)*)

(* ::Text:: *)
(*Flat 3D Poisson toy model of the complex-contour regularisation project (A. Spiers, UCD).*)
(*Spec: "Claude Implementation/singular-field-consistency.tex".*)
(**)
(*Conventions*)
(*  Particle charge q at x0 = r0 zhat, rho = |x - x0|.  Phi = Sum_l phi_l(r) Y_l0(theta), Y_l0 orthonormal.*)
(*  Mode of a 3D function: [f]_l(r) = Sqrt[Pi (2l+1)] Integrate[f P_l(mu), {mu,-1,1}].*)
(*  Radial operator L_l psi = psi'' + 2 psi'/r - l(l+1) psi/r^2.  [Laplacian f]_l = L_l [f]_l.*)
(*  First order: Laplacian Phi1 = -4 Pi q delta^3, mode source kappa_l delta(r-r0), kappa_l = -(q/r0^2) Sqrt[4 Pi (2l+1)].*)
(*  Unit-strength Coulomb mode G_l = -(r0^2/(2l+1)) r_<^l/r_>^(l+1);  kappa_l G_l = [q/rho]_l.*)
(*  Extended delta: H(z) = 1 + (I/Pi) Log[z-r0],  deltaExt(z) = I/(Pi (z-r0)).  Real-axis values are boundary values from above.*)
(*  Contour z = y(a) = a + I amp (1-((a-r0)/radius)^2)^power on [aL,aR] = [r0-radius, r0+radius], real outside.*)
(*  Second order: Laplacian Phi2 = Q[Phi1,Phi1];  Q = "Product" (Q[A,B] = A B, the 1/rho^2 case) and, since Stage 5,*)
(*  Q = "GradDot" (Grad A . Grad B, the 1/rho^4 case).  For "GradDot" the regular-field inputs are pairs {phi, dphi/dr}.*)
(*  Complex radial variable: z.  Cartesian z-coordinate is written r mu (never z) to avoid a clash.*)

BeginPackage["ContourSecondOrder`"];

(* ---------------- parameters and utilities ---------------- *)
$CSParameters::usage = "$CSParameters is the association of model and numerical parameters.";
SetCSParameters::usage = "SetCSParameters[assoc] merges assoc into the defaults, stores it in $CSParameters and clears caches.";
csPar::usage = "csPar[key] returns a parameter (q, r0, rMin, rMax, radius, amp, power, c0, c1, wp, precisionGoal, accuracyGoal).";
clearCSCaches::usage = "clearCSCaches[] clears all memoised results.";
toWP::usage = "toWP[x] converts x to WorkingPrecision csPar[\"wp\"].";
sci::usage = "sci[x] formats a number as a short string in scientific notation.";
csCheck::usage = "csCheck[name, value, tol] records and prints PASS/FAIL for |value| <= tol.";
scaledError::usage = "scaledError[num, ref] = Max|num - ref| / Max|ref| (absolute error if ref is identically zero).";
csCheckList::usage = "csCheckList[] returns the list of recorded checks {name, value, tol, pass}.";
csResetChecks::usage = "csResetChecks[] clears the recorded checks.";
csSummary::usage = "csSummary[] prints a PASS/FAIL summary of the recorded checks and returns True if all passed.";
zVar::usage = "zVar is the formal complex radial variable used in symbolic expressions returned by the package.";

(* ---------------- contour, operators, distributions ---------------- *)
tubeL::usage = "tubeL[] = r0 - radius (left end of the worldtube).";
tubeR::usage = "tubeR[] = r0 + radius (right end of the worldtube).";
yC::usage = "yC[a] is the contour z = y(a).";
yCp::usage = "yCp[a] = y'(a).";
yCpp::usage = "yCpp[a] = y''(a).";
radialL::usage = "radialL[l, expr, z] applies the radial operator L_l to the expression expr in the variable z.";
realL::usage = "realL[l, u, r] is L_l acting on the function u at r (for NDSolve).";
contourL::usage = "contourL[l, u, a] is L_l written in the contour parameter a for u(a) = psi(y(a)) (the notebook's odeOpComplex).";
HExt::usage = "HExt[z] = 1 + (I/Pi) Log[z - r0].";
deltaExt::usage = "deltaExt[z] = I/(Pi (z - r0)).";
coulombRatio::usage = "coulombRatio[L, r] = r_<^L/r_>^(L+1) for real r.";
coulombModeUnit::usage = "coulombModeUnit[l, r] = G_l(r), unit-strength Coulomb mode.";
kappaL::usage = "kappaL[l] = -(q/r0^2) Sqrt[4 Pi (2l+1)].";
coulombRatioContinued::usage = "coulombRatioContinued[L, z, method] continues r_<^L/r_>^(L+1) off the real axis. method = \"ThetaToH\" (Theta -> H) or \"CorrectedSingular\" (uses the corrected first-order singular mode of Stage 1).";

(* ---------------- branch handling ---------------- *)
canonicalizeLogs::usage = "canonicalizeLogs[expr, z] rewrites every Log[c1 z + c0] (real c1, c0) as Log[|c1|] + Log[z - z*] - I Pi Boole[c1 < 0], the form that is analytic in the upper half plane and whose principal-branch value on the real axis is the boundary value from above.";
fromAbove::usage = "fromAbove[f, r] evaluates the function f at r + I eps (eps = 10^-40) with 80-digit arithmetic: the boundary value from above.";

(* ---------------- Stage 1: first-order singular mode ---------------- *)
dsolveSingularMode::usage = "dsolveSingularMode[l] returns Function[z, psi_p] where psi_p = DSolveValue[L_l psi == deltaExt] /. C[_] -> 0, unmodified.";
singularModeCanonical::usage = "singularModeCanonical[l] is dsolveSingularMode[l] with logs canonicalised (same analytic function in the upper half plane).";
singularModeD::usage = "singularModeD[l] = {alpha, beta}: D_l(r) = Re psi_p(r+i0) - G_l(r) = alpha r^l + beta r^(-l-1), fitted exactly at r = 2, 3.";
singularModeDValue::usage = "singularModeDValue[l, r] = Re psi_p(r+i0) - G_l(r), exact, for real r (from the canonical form).";
singularModeCorrected::usage = "singularModeCorrected[l] returns Function[z, phiS_l] with phiS_l = psi_p - D_l (unit strength).";
singularModeTaylor::usage = "singularModeTaylor[l, n] returns the truncated puncture: series of phiS_l about r0 through order n.";

(* ---------------- solvers ---------------- *)
ndOptions::usage = "ndOptions[] returns the NDSolve options used throughout.";
solveRealSegment::usage = "solveRealSegment[l, src, {r1, r2}, {v, d}] solves L_l u = src(r) on the real segment from r1 to r2 with u(r1) = v, u'(r1) = d.";
solveContourSegment::usage = "solveContourSegment[l, srcA, {v, d}] solves L_l u = srcA(a) along the contour a in [aL, aR] with u(aL) = v, u'(aL) = d.";
threeRegionSolve::usage = "threeRegionSolve[l, srcReal, srcContour, {v1, d1}] solves on [rMin,aL] (real), the contour, and [aR,rMax] (real, data Re u(aR), Re u'(aR)). Returns an association with sol1, sol2, sol3, valR, derR, imAR.";
threeRegionBVP::usage = "threeRegionBVP[l, srcReal, srcContour, {vLeft, vRight}] solves the two-point problem u(rMin) = vLeft, u(rMax) = vRight by linear shooting on u'(rMin).";
realLineSolve::usage = "realLineSolve[l, src, {v, d}] solves L_l u = src on [rMin, rMax] on the real line (Cauchy data at rMin).";
realLineBVP::usage = "realLineBVP[l, src, {vLeft, vRight}] solves the two-point problem on the real line by linear shooting.";
exteriorValue::usage = "exteriorValue[res, r] evaluates a three-region solution at an exterior point r.";
numericFunction::usage = "numericFunction[expr, var] returns a black-box function f[x_?NumericQ] evaluating expr at extra precision.";
contourFunction::usage = "contourFunction[expr, z] returns a black-box function f[a_?NumericQ] that evaluates expr at z = y(a): y(a) is computed numerically first (building expr[y(a)] symbolically triggers polynomial expansion and catastrophic cancellation).";

(* ---------------- Stage 2: first-order regular field ---------------- *)
phiRAnalytic::usage = "phiRAnalytic[l] = Function[z, phiR_l(z)] for PhiR = c0 + c1 (r mu): only l = 0, 1 nonzero.";
homogeneousWithData::usage = "homogeneousWithData[l, r1, {v, d}] = Function[z, A z^l + B z^(-l-1)] with the given data at r1.";
seffFunction::usage = "seffFunction[l, puncture] returns the effective source kappa_l (deltaExt - L_l[puncture]) on the contour, computed as in the notebook.";
firstOrderRegular::usage = "firstOrderRegular[l, data, puncture] runs the notebook's three-region regular-field solve. data = Automatic uses phiRAnalytic. puncture = \"exact\" (corrected phiS), n (Taylor of phiS through order n) or \"ThetaToH\" (G_l with Theta -> H; negative control).";
firstOrderFull::usage = "firstOrderFull[l, data] solves the full field (source kappa_l deltaExt) on the contour starting from regular + singular data at aL.";

(* ---------------- Stage 3: second-order sources and couplings ---------------- *)
gaunt::usage = "gaunt[l1, l2, l3] = Integrate[Y_l1 Y_l2 Y_l3 dOmega] (m = 0).";
KCoupling::usage = "KCoupling[l, lp, r] = Integrate[Y_l Y_lp / rho dOmega] (real r).";
KCouplingContinued::usage = "KCouplingContinued[l, lp, z, method] is the continuation of KCoupling (see coulombRatioContinued).";
WLegendre::usage = "WLegendre[l, x] = W_(l-1)(x) = Sum[P_(k-1)(x) P_(l-k)(x)/k, {k,1,l}].";
legendreQNeumann::usage = "legendreQNeumann[l, x, logTerm] = P_l(x) logTerm/2 - W_(l-1)(x), with logTerm = Log[(x+1)/(x-1)].";
xOfZ::usage = "xOfZ[z] = (z^2 + r0^2)/(2 z r0).";
ssSourceModeExpr::usage = "ssSourceModeExpr[Q, l, z, branch] is the SS source mode; branch = \"real\" (uses Log|r - r0|) or \"continued\" (uses Log[z - r0]).";
ssSourceReal::usage = "ssSourceReal[Q, l] returns a numeric black-box function of real r.";
ssSourceContour::usage = "ssSourceContour[Q, l] returns a function of the contour parameter a.";
srSourceReal::usage = "srSourceReal[Q, l, phiRfuns] with phiRfuns = <|lp -> function of r|> returns the SR source mode on the real line.";
srSourceContour::usage = "srSourceContour[Q, l, phiRtilde, method] with phiRtilde = <|lp -> function of a|> returns the SR source on the contour.";
rrSourceReal::usage = "rrSourceReal[Q, l, phiRfuns] returns the RR source mode (Gaunt coupling) on the real line.";
bilinearModes::usage = "bilinearModes[Q, l, aModes, bModes] couples two mode sets (associations l -> value) into mode l of Q[A,B].";

(* ---------------- analytic modes for checks ---------------- *)
uModeIntegral::usage = "uModeIntegral[l, m, fType, r, branch] = Sqrt[Pi(2l+1)] Integrate[mu^m P_l(mu) f(r^2+r0^2-2 r r0 mu), {mu,-1,1}] in closed form, branch \"in\" (r<r0) or \"out\" (r>r0). fType: \"sqrt\" (u^(1/2)), \"halfLog\" (Log[u]/2), \"invSqrt\" (u^(-1/2)), \"inv\" (1/u), \"inv2\", \"invSqrt3\", or {\"pow\", k} (u^k = rho^(2k)).";
muPowerMode::usage = "muPowerMode[l, m] = Sqrt[Pi(2l+1)] Integrate[mu^m P_l(mu), {mu,-1,1}].";
pieceAnalyticSolution::usage = "pieceAnalyticSolution[Q, piece, l] = {inExpr, outExpr} in zVar: modes of the analytic particular solution (SS: q^2 Log rho; SR: q c0 rho + q c1 (r0 rho + (r mu - r0) rho/2); RR: c0^2 r^2/6 + c0 c1 r^3 mu/5 + c1^2 r^4 mu^4/12).";
pieceAnalyticSource::usage = "pieceAnalyticSource[Q, piece, l] = {inExpr, outExpr} in zVar: modes of the source of each piece by direct (closed-form) angular integration.";
analyticEval::usage = "analyticEval[{in, out}, r] evaluates a piecewise analytic mode at real r (exact input, 30 digits).";
analyticEvalD::usage = "analyticEvalD[{in, out}, r] evaluates the r-derivative.";
analytic3D::usage = "analytic3D[piece, r, theta] (Q = Product) or analytic3D[Q, piece, r, theta] is the 3D analytic particular solution for piece = \"SS\", \"SR\", \"RR\".";
legendreQNeumannDx::usage = "legendreQNeumannDx[l, x, logTerm, inv] = Q_l'(x) = P_l'(x) logTerm/2 - P_l(x) inv - W_(l-1)'(x), with logTerm = Log[(x+1)/(x-1)] and inv = 1/(x^2-1).";
invX2m1::usage = "invX2m1[z] = 1/(x(z)^2 - 1) = 4 z^2 r0^2/((z-r0)^2 (z+r0)^2), meromorphic.";
dCoulombRatio::usage = "dCoulombRatio[L, r] = d/dr (r_<^L/r_>^(L+1)) for real r != r0 (jump at r0).";
dCoulombRatioContinued::usage = "dCoulombRatioContinued[L, z, method]: \"ThetaToH\" (piecewise derivative with Theta -> H) or \"CorrectedSingular\" (d/dz of the corrected singular mode).";
KdCoupling::usage = "KdCoupling[l, lp, r] = Integrate[d_r(1/rho) Y_l Y_lp dOmega] = d/dr K_{l lp}(r).";
KthetaCoupling::usage = "KthetaCoupling[l, lp, r] = Integrate[d_theta(1/rho) d_theta(Y_lp) Y_l dOmega] = Sum_L <l lp L> [1/rho]_L (L(L+1) + lp(lp+1) - l(l+1))/2.";
KdCouplingContinued::usage = "KdCouplingContinued[l, lp, z, method] continues KdCoupling off the axis.";
KthetaCouplingContinued::usage = "KthetaCouplingContinued[l, lp, z, method] continues KthetaCoupling off the axis.";
aZ::usage = "aZ[l] = l/Sqrt[(2l-1)(2l+1)] (0 for l <= 0): coefficient of the d/dz mode coupling.";
dzExpr::usage = "dzExpr[F, l] = (d_z F)_l for mode expressions F[j] in zVar (one branch): a_l (F'_(l-1) - (l-1) F_(l-1)/z) + a_(l+1) (F'_(l+1) + (l+2) F_(l+1)/z).";
dzzExpr::usage = "dzzExpr[F, l] = (d_z^2 F)_l (pointwise, one branch).";
dzzDeltaCoeff::usage = "dzzDeltaCoeff[Fpair, l]: coefficient of delta(r - r0) in the distributional (d_z^2 F)_l for continuous modes Fpair[j] = {in, out}.";
phiRAtParticle::usage = "phiRAtParticle[] = PhiR(x0) = c0 + c1 r0 + c2 r0^2."; 

Begin["`Private`"];

(* ================= parameters ================= *)
$CSDefaults = <|"q" -> 1, "r0" -> 5, "rMin" -> 1, "rMax" -> 10, "radius" -> 1/2, "amp" -> 1, "power" -> 6,
   "c0" -> 1, "c1" -> 1/10, "c2" -> 0, "wp" -> 30, "precisionGoal" -> 20, "accuracyGoal" -> 35, "extraPrecision" -> 40|>;
(* accuracyGoal: with 22 the l >= 7 second-order modes, which are ~(r/r0)^l small at r = 1, lost relative accuracy in
   region 1 (scaled exterior error 5e-10 at l = 11, 2e-9 at l = 12).  Infinity (pure relative control) breaks NDSolve when
   a component is exactly zero (e.g. phi'(1) = 0); 40 fails the same way; 35 works for every case tested. *)
$CSParameters = $CSDefaults;
$cache = <||>;
clearCSCaches[] := ($cache = <||>;);
SetCSParameters[new_Association] := ($CSParameters = Join[$CSDefaults, new]; clearCSCaches[]; $CSParameters);
csPar[k_String] := $CSParameters[k];

SetAttributes[memo, HoldRest];
memo[key_, expr_] := With[{k = {key, $CSParameters}},
   If[KeyExistsQ[$cache, k], $cache[k], $cache[k] = expr]];

toWP[x_] := If[Precision[x] === Infinity, N[x, csPar["wp"]], SetPrecision[x, csPar["wp"]]];

sci[x_] := Module[{y = If[Precision[x] === Infinity, Quiet@N[Abs[x]], Abs[x]], e},
   Which[
    ! NumericQ[y], ToString[x],
    Precision[x] === Infinity && x == 0, "0",
    y == 0 && MachineNumberQ[y], "0.(machine)",
    y == 0, "<1e" <> ToString[-Round[Accuracy[y]]],
    y < 10^-300, "<1e-300",
    True,
    e = Floor[Log10[y]];
    (If[Re[N[x]] < 0 && Im[N[x]] == 0, "-", ""]) <>
     ToString[NumberForm[SetPrecision[y/10^e, MachinePrecision], {3, 2}]] <> "e" <> ToString[e]]];

scaledError[num_List, ref_List] := With[{m = Max[Abs[ref]]}, If[m == 0, Max[Abs[num - ref]], Max[Abs[num - ref]]/m]];

$checks = {};
csResetChecks[] := (scaledError[num_List, ref_List] := With[{m = Max[Abs[ref]]}, If[m == 0, Max[Abs[num - ref]], Max[Abs[num - ref]]/m]];

$checks = {};);
csCheckList[] := $checks;
csCheck[name_String, value_, tol_] := Module[{v = Abs[value], pass},
   (* an inexact zero only certifies |value| < 10^-Accuracy: compare that bound with the tolerance *)
   pass = TrueQ[NumericQ[v] && If[Precision[v] =!= Infinity && v == 0, 10^-Accuracy[v] <= tol || tol == 0 && Accuracy[v] > 300, v <= tol]];
   AppendTo[$checks, {name, v, tol, pass}];
   Print[If[pass, "PASS  ", "FAIL  "], name, ":  ", sci[v], "   (tol ", sci[tol], ")"];
   pass];
csSummary[] := Module[{np = Count[$checks, {_, _, _, True}], nf = Count[$checks, {_, _, _, False}]},
   Print["------------------------------------------------------------"];
   Print["SUMMARY: ", np, " passed, ", nf, " failed (of ", Length[$checks], ")"];
   If[nf > 0, Print["FAILED: ", Column[Select[$checks, ! Last[#] &][[All, 1]]]]];
   nf == 0];

toFunction[expr_, var_Symbol] := Function @@ {{var}, expr};

numericFunction[expr_, var_Symbol] := With[{prec = csPar["wp"] + csPar["extraPrecision"]},
   Module[{f},
    f[x_?NumericQ] := N[expr /. var -> SetPrecision[x, prec], prec];
    f]];

contourFunction[expr_, var_Symbol] := With[{prec = csPar["wp"] + csPar["extraPrecision"]},
   Module[{f},
    f[x_?NumericQ] := N[expr /. var -> N[yC[SetPrecision[x, prec]], prec], prec];
    f]];

(* ================= contour and operators ================= *)
tubeL[] := csPar["r0"] - csPar["radius"];
tubeR[] := csPar["r0"] + csPar["radius"];
yC[a_] := a + I csPar["amp"] (1 - ((a - csPar["r0"])/csPar["radius"])^2)^csPar["power"];
yCp[a_] := Module[{t}, D[yC[t], t] /. t -> a];
yCpp[a_] := Module[{t}, D[yC[t], {t, 2}] /. t -> a];

radialL[l_, expr_, z_] := D[expr, {z, 2}] + 2/z D[expr, z] - l (l + 1)/z^2 expr;
realL[l_, u_, r_] := u''[r] + (2/r) u'[r] - (l (l + 1)/r^2) u[r];
contourL[l_, u_, a_] := (1/yCp[a]^2) u''[a] - (yCpp[a]/yCp[a]^3) u'[a] +
    (2/yC[a]) (1/yCp[a]) u'[a] - (l (l + 1)/yC[a]^2) u[a];

HExt[z_] := 1 + I/Pi Log[z - csPar["r0"]];
deltaExt[z_] := I/(Pi (z - csPar["r0"]));
coulombRatio[L_, r_] := Piecewise[{{r^L/csPar["r0"]^(L + 1), r < csPar["r0"]}}, csPar["r0"]^L/r^(L + 1)];
coulombModeUnit[l_, r_] := -(csPar["r0"]^2/(2 l + 1)) coulombRatio[l, r];
kappaL[l_] := -(csPar["q"]/csPar["r0"]^2) Sqrt[4 Pi (2 l + 1)];

coulombRatioContinued[L_, z_, "ThetaToH"] := With[{r0 = csPar["r0"]},
   z^L/r0^(L + 1) + (r0^L/z^(L + 1) - z^L/r0^(L + 1)) HExt[z]];
coulombRatioContinued[L_, z_, "CorrectedSingular"] :=
  -((2 L + 1)/csPar["r0"]^2) singularModeCorrected[L][z];

(* ================= branch handling ================= *)
canonicalizeLogs[expr_, z_Symbol] := expr /. Log[arg_] /; ! FreeQ[arg, z] :> canonLog[arg, z];
canonLog[arg_, z_] := Module[{p = Expand[arg], c1, c0, zs},
   If[PolynomialQ[p, z] && Exponent[p, z] == 1,
    c1 = Coefficient[p, z, 1]; c0 = Coefficient[p, z, 0]; zs = -c0/c1;
    If[TrueQ[Element[c1, Reals]] && TrueQ[Element[zs, Reals]],
     Log[Abs[c1]] + Log[z - zs] + If[TrueQ[c1 < 0], -I Pi, 0],
     Log[arg]],
    Log[arg]]];

fromAbove[f_, r_?NumericQ, eps_: 10^-40, prec_: 80] := N[f[SetPrecision[r, prec] + I eps], prec - 30];

(* ================= Stage 1 ================= *)
dsolveSingularMode[l_Integer?NonNegative] := memo[{"psiP", l},
   Module[{psi, r, sol},
    sol = DSolveValue[psi''[r] + (2/r) psi'[r] - (l (l + 1)/r^2) psi[r] == I/(Pi (r - csPar["r0"])), psi[r], r] /. C[_] -> 0;
    toFunction[sol /. r -> zVar, zVar]]];

singularModeCanonical[l_Integer?NonNegative] := memo[{"psiPcanon", l},
   toFunction[canonicalizeLogs[dsolveSingularMode[l][zVar], zVar], zVar]];

singularModeDValue[l_, r_?NumericQ] :=
  Simplify[ComplexExpand[Re[singularModeCanonical[l][r]]] - coulombModeUnit[l, r]];

singularModeD[l_Integer?NonNegative] := memo[{"Dfit", l},
   Module[{r1 = 2, r2 = 3, d1, d2, al, be},
    d1 = singularModeDValue[l, r1]; d2 = singularModeDValue[l, r2];
    {al, be} = {al, be} /. First@Solve[{al r1^l + be r1^(-l - 1) == d1, al r2^l + be r2^(-l - 1) == d2}, {al, be}];
    Simplify[{al, be}]]];

singularModeCorrected[l_Integer?NonNegative] := memo[{"phiS", l},
   With[{ab = singularModeD[l]},
    toFunction[singularModeCanonical[l][zVar] - ab[[1]] zVar^l - ab[[2]] zVar^(-l - 1), zVar]]];

singularModeTaylor[l_Integer?NonNegative, n_Integer?Positive] := memo[{"phiSTaylor", l, n},
   toFunction[Normal[Series[singularModeCorrected[l][zVar], {zVar, csPar["r0"], n}]], zVar]];

(* ================= solvers ================= *)
(* InterpolationOrder -> All: dense output of full method order.  With the default (low-order Hermite) dense output,
   solutions that the integrator crosses in very few steps (e.g. the polynomial-like SR modes for r < r0) showed
   1e-12 interpolation noise between steps although the step values were accurate to 1e-21. *)
ndOptions[] := {WorkingPrecision -> csPar["wp"], PrecisionGoal -> csPar["precisionGoal"],
   AccuracyGoal -> csPar["accuracyGoal"], MaxSteps -> 10^6, InterpolationOrder -> All};

solveRealSegment[l_, src_, {r1_, r2_}, {v_, d_}] := Module[{u, r},
   NDSolveValue[{realL[l, u, r] == src[r], u[r1] == toWP[v], u'[r1] == toWP[d]}, u,
    {r, Min[r1, r2], Max[r1, r2]}, Evaluate[Sequence @@ ndOptions[]]]];

solveContourSegment[l_, srcA_, {v_, d_}] := Module[{u, a},
   NDSolveValue[{contourL[l, u, a] == srcA[a], u[tubeL[]] == toWP[v], u'[tubeL[]] == toWP[d]}, u,
    {a, tubeL[], tubeR[]}, Evaluate[Sequence @@ ndOptions[]]]];

threeRegionSolve[l_, srcReal_, srcContour_, {v1_, d1_}] := Module[{s1, s2, s3, vR, dR, aL = tubeL[], aR = tubeR[]},
   (* y'(aL) = y'(aR) = 1 and y'' = 0 there, so d/da = d/dz at the matching points *)
   s1 = solveRealSegment[l, srcReal, {csPar["rMin"], aL}, {v1, d1}];
   s2 = solveContourSegment[l, srcContour, {s1[aL], s1'[aL]}];
   vR = s2[aR]; dR = s2'[aR];
   s3 = solveRealSegment[l, srcReal, {aR, csPar["rMax"]}, {Re[vR], Re[dR]}];
   <|"sol1" -> s1, "sol2" -> s2, "sol3" -> s3, "valR" -> vR, "derR" -> dR, "imAR" -> {Im[vR], Im[dR]}|>];

threeRegionBVP[l_, srcReal_, srcContour_, {vLeft_, vRight_}] := Module[{f0, f1, p, res, rMax = csPar["rMax"]},
   f0 = threeRegionSolve[l, srcReal, srcContour, {vLeft, 0}]["sol3"][rMax];
   f1 = threeRegionSolve[l, srcReal, srcContour, {vLeft, 1}]["sol3"][rMax];
   p = (toWP[vRight] - f0)/(f1 - f0);
   res = threeRegionSolve[l, srcReal, srcContour, {vLeft, p}];
   Join[res, <|"shootDerivative" -> p, "bvpResidual" -> res["sol3"][rMax] - vRight|>]];

realLineSolve[l_, src_, {v_, d_}] := solveRealSegment[l, src, {csPar["rMin"], csPar["rMax"]}, {v, d}];

realLineBVP[l_, src_, {vLeft_, vRight_}] := Module[{f0, f1, p, rMax = csPar["rMax"]},
   f0 = realLineSolve[l, src, {vLeft, 0}][rMax];
   f1 = realLineSolve[l, src, {vLeft, 1}][rMax];
   p = (toWP[vRight] - f0)/(f1 - f0);
   realLineSolve[l, src, {vLeft, p}]];

exteriorValue[res_Association, r_] := If[r <= tubeL[], res["sol1"][r], res["sol3"][r]];

(* ================= Stage 2 ================= *)
phiRAnalytic[l_Integer] := Which[
   l == 0, Function @@ {{zVar}, csPar["c0"] Sqrt[4 Pi]},
   l == 1, Function @@ {{zVar}, csPar["c1"] Sqrt[4 Pi/3] zVar},
   l == 2, Function @@ {{zVar}, csPar["c2"] Sqrt[4 Pi/5] zVar^2},
   True, Function @@ {{zVar}, 0}];

homogeneousWithData[l_, r1_, {v_, d_}] := Module[{A, B, sol},
   sol = First@Solve[{A r1^l + B r1^(-l - 1) == v, l A r1^(l - 1) - (l + 1) B r1^(-l - 2) == d}, {A, B}];
   toFunction[(A zVar^l + B zVar^(-l - 1)) /. sol, zVar]];

punctureFunction[l_, "exact"] := singularModeCorrected[l];
punctureFunction[l_, n_Integer] := singularModeTaylor[l, n];
(* G_l with Theta -> H: correct real part, but violates (C2) (note section 6.2); used only as a negative control *)
punctureFunction[l_, "ThetaToH"] := toFunction[-(csPar["r0"]^2/(2 l + 1)) coulombRatioContinued[l, zVar, "ThetaToH"], zVar];

(* The notebook's Seff[a] = complexSource[y[a]] - odeOpComplex[localS, a], times kappa_l *)
seffFunction[l_, punct_: "exact"] := memo[{"seff", l, punct},
   Module[{a, ps = punctureFunction[l, punct], expr},
    expr = kappaL[l] (deltaExt[yC[a]] - contourL[l, Function @@ {{a}, ps[yC[a]]}, a]);
    numericFunction[expr, a]]];

(* puncture values at the matching points: boundary values from above (canonical logs) *)
punctureData[l_, punct_, r_] := Module[{ps = punctureFunction[l, punct], z},
   kappaL[l] {ps[r], D[ps[z], z] /. z -> r}];

firstOrderRegular[l_Integer, data_: Automatic, punct_: "exact"] := memo[{"phiR1", l, data, punct},
   Module[{phiA = phiRAnalytic[l], d, res, aL = tubeL[], aR = tubeR[], pL, pR, sL, sR, vR, dR, s3},
    d = If[data === Automatic, {phiA[csPar["rMin"]], phiA'[csPar["rMin"]]}, data];
    (* Region 2 data: valL + exactValL - taylorValL (notebook), exact = phiS, taylor = puncture *)
    pL = punctureData[l, punct, aL]; sL = punctureData[l, "exact", aL];
    pR = punctureData[l, punct, aR]; sR = punctureData[l, "exact", aR];
    res = Module[{s1, s2},
      s1 = solveRealSegment[l, Function[r, 0], {csPar["rMin"], aL}, d];
      s2 = solveContourSegment[l, seffFunction[l, punct], {s1[aL] + sL[[1]] - pL[[1]], s1'[aL] + sL[[2]] - pL[[2]]}];
      vR = s2[aR] + pR[[1]] - sR[[1]]; dR = s2'[aR] + pR[[2]] - sR[[2]];
      s3 = solveRealSegment[l, Function[r, 0], {aR, csPar["rMax"]}, {Re[vR], Re[dR]}];
      <|"sol1" -> s1, "sol2" -> s2, "sol3" -> s3, "valR" -> vR, "derR" -> dR, "imAR" -> {Im[vR], Im[dR]}|>];
    res]];

firstOrderFull[l_Integer, data_: Automatic] := memo[{"phiFull1", l, data},
   Module[{phiA = phiRAnalytic[l], d, s1, sL, src, s2, aL = tubeL[], a},
    d = If[data === Automatic, {phiA[csPar["rMin"]], phiA'[csPar["rMin"]]}, data];
    s1 = solveRealSegment[l, Function[r, 0], {csPar["rMin"], aL}, d];
    sL = punctureData[l, "exact", aL];
    src = numericFunction[kappaL[l] deltaExt[yC[a]], a];
    s2 = solveContourSegment[l, src, {s1[aL] + sL[[1]], s1'[aL] + sL[[2]]}];
    <|"sol1" -> s1, "sol2" -> s2|>]];

(* ================= Stage 3: couplings and sources ================= *)
gaunt[l1_, l2_, l3_] := If[Abs[l1 - l2] <= l3 <= l1 + l2 && EvenQ[l1 + l2 + l3],
   Sqrt[(2 l1 + 1) (2 l2 + 1) (2 l3 + 1)/(4 Pi)] ThreeJSymbol[{l1, 0}, {l2, 0}, {l3, 0}]^2, 0];

KCoupling[l_, lp_, r_] := Sum[gaunt[l, lp, L] Sqrt[4 Pi/(2 L + 1)] coulombRatio[L, r], {L, Abs[l - lp], l + lp}];
KCouplingContinued[l_, lp_, z_, method_] :=
  Sum[gaunt[l, lp, L] Sqrt[4 Pi/(2 L + 1)] coulombRatioContinued[L, z, method], {L, Abs[l - lp], l + lp}];

WLegendre[l_, x_] := Sum[LegendreP[k - 1, x] LegendreP[l - k, x]/k, {k, 1, l}];
legendreQNeumann[l_, x_, logTerm_] := (1/2) LegendreP[l, x] logTerm - WLegendre[l, x];
xOfZ[z_] := (z^2 + csPar["r0"]^2)/(2 z csPar["r0"]);

notImplemented[Q_] := (Print["Q = ", Q, " is not implemented (1/rho^4 case deferred)."]; $Failed);

(* --- SS --- *)
ssSourceModeExpr["Product", l_, z_, "real"] := With[{r0 = csPar["r0"]},
   csPar["q"]^2 Sqrt[Pi (2 l + 1)]/(z r0) legendreQNeumann[l, xOfZ[z], 2 Log[z + r0] - 2 Log[Abs[z - r0]]]];
ssSourceModeExpr["Product", l_, z_, "continued"] := With[{r0 = csPar["r0"]},
   csPar["q"]^2 Sqrt[Pi (2 l + 1)]/(z r0) legendreQNeumann[l, xOfZ[z], 2 Log[z + r0] - 2 Log[z - r0]]];
ssSourceModeExpr[Q_, ___] := notImplemented[Q];

ssSourceReal[Q_, l_] := Module[{r}, numericFunction[ssSourceModeExpr[Q, l, r, "real"], r]];
ssSourceContour[Q_, l_] := Module[{z}, contourFunction[ssSourceModeExpr[Q, l, z, "continued"], z]];

(* --- generic bilinear coupling of two mode sets --- *)
bilinearModes["Product", l_, aModes_Association, bModes_Association] :=
  Sum[gaunt[l, l1, l2] aModes[l1] bModes[l2], {l1, Keys[aModes]}, {l2, Keys[bModes]}];
bilinearModes[Q_, ___] := notImplemented[Q];

(* --- SR --- *)
srSourceReal["Product", l_, phiRfuns_Association] := Module[{f, q = csPar["q"]},
   f[r_?NumericQ] := 2 q Sum[KCoupling[l, lp, r] phiRfuns[lp][r], {lp, Keys[phiRfuns]}];
   f];
srSourceContour["Product", l_, phiRtilde_Association, method_] := Module[{f, q = csPar["q"], a, kexpr},
   kexpr = AssociationMap[Function[lp, Module[{z}, contourFunction[KCouplingContinued[l, lp, z, method], z]]], Keys[phiRtilde]];
   f[x_?NumericQ] := 2 q Sum[kexpr[lp][x] phiRtilde[lp][x], {lp, Keys[phiRtilde]}];
   f];
srSourceReal[Q_, ___] := notImplemented[Q];
srSourceContour[Q_, ___] := notImplemented[Q];

(* --- RR --- *)
rrSourceReal["Product", l_, phiRfuns_Association] := Module[{f},
   f[r_?NumericQ] := bilinearModes["Product", l, Map[#[r] &, phiRfuns], Map[#[r] &, phiRfuns]];
   f];
rrSourceReal[Q_, ___] := notImplemented[Q];

(* ================= analytic modes ================= *)
(* Antiderivative of u^j f(u), evaluated at u = w^2 with w > 0 *)
uAnti["sqrt", j_, w_] := w^(2 j + 3)/(j + 3/2);
uAnti["invSqrt", j_, w_] := w^(2 j + 1)/(j + 1/2);
uAnti["halfLog", j_, w_] := (1/2) (w^(2 j + 2) 2 Log[w]/(j + 1) - w^(2 j + 2)/(j + 1)^2);
uAnti["inv", j_, w_] := If[j == 0, 2 Log[w], w^(2 j)/j];
uAnti["inv2", j_, w_] := Which[j == 0, -1/w^2, j == 1, 2 Log[w], True, w^(2 j - 2)/(j - 1)];
uAnti["invSqrt3", j_, w_] := w^(2 j - 1)/(j - 1/2);
uAnti[{"pow", k_}, j_, w_] := If[j + k == -1, 2 Log[w], w^(2 (j + k + 1))/(j + k + 1)];   (* general f(u) = u^k, i.e. rho^(2k) (Stage 6) *)

uModeIntegral[l_, m_, fType_, r_, branch_] := Module[{mu, u, a, b, s, t, coeffs, res, r0 = csPar["r0"]},
   a = r^2 + r0^2; b = 2 r r0;
   coeffs = CoefficientList[Expand[(mu^m LegendreP[l, mu]) /. mu -> (a - u)/b], u];
   res = (1/b) Sum[coeffs[[j + 1]] (uAnti[fType, j, t] - uAnti[fType, j, s]), {j, 0, Length[coeffs] - 1}];
   res = Sqrt[Pi (2 l + 1)] res /. {t -> r + r0, s -> If[branch === "in", r0 - r, r - r0]};
   Simplify[res, If[branch === "in", 0 < r < r0, r > r0]]];

muPowerMode[l_, m_] := Sqrt[Pi (2 l + 1)] Integrate[mu^m LegendreP[l, mu], {mu, -1, 1}];

pieceAnalyticSolution["Product", "SS", l_] := memo[{"solSS", l},
   Table[csPar["q"]^2 uModeIntegral[l, 0, "halfLog", zVar, br], {br, {"in", "out"}}]];
pieceAnalyticSolution["Product", "SR", l_] := memo[{"solSR", l},
   With[{q = csPar["q"], c0 = csPar["c0"], c1 = csPar["c1"], r0 = csPar["r0"]},
    Table[Simplify[q (c0 + c1 r0/2) uModeIntegral[l, 0, "sqrt", zVar, br] + (q c1/2) zVar uModeIntegral[l, 1, "sqrt", zVar, br]],
     {br, {"in", "out"}}]]];
pieceAnalyticSolution["Product", "RR", l_] := memo[{"solRR", l},
   With[{c0 = csPar["c0"], c1 = csPar["c1"]},
    ConstantArray[Simplify[c0^2 zVar^2/6 muPowerMode[l, 0] + c0 c1 zVar^3/5 muPowerMode[l, 1] + c1^2 zVar^4/12 muPowerMode[l, 4]], 2]]];
pieceAnalyticSolution[Q_, ___] := notImplemented[Q];

pieceAnalyticSource["Product", "SS", l_] := memo[{"srcSS", l},
   Table[csPar["q"]^2 uModeIntegral[l, 0, "inv", zVar, br], {br, {"in", "out"}}]];
pieceAnalyticSource["Product", "SR", l_] := memo[{"srcSR", l},
   With[{q = csPar["q"], c0 = csPar["c0"], c1 = csPar["c1"]},
    Table[Simplify[2 q (c0 uModeIntegral[l, 0, "invSqrt", zVar, br] + c1 zVar uModeIntegral[l, 1, "invSqrt", zVar, br])],
     {br, {"in", "out"}}]]];
pieceAnalyticSource["Product", "RR", l_] := memo[{"srcRR", l},
   With[{c0 = csPar["c0"], c1 = csPar["c1"]},
    ConstantArray[Simplify[c0^2 muPowerMode[l, 0] + 2 c0 c1 zVar muPowerMode[l, 1] + c1^2 zVar^2 muPowerMode[l, 2]], 2]]];
pieceAnalyticSource[Q_, ___] := notImplemented[Q];

(* Quiet N::meprec: raised only when a mode is an exact zero that N cannot certify (it then returns 0 with accuracy ~80) *)
analyticEval[{in_, out_}, r_] := Quiet[N[If[r < csPar["r0"], in, out] /. zVar -> r, 30], N::meprec];
analyticEvalD[{in_, out_}, r_] := N[D[If[r < csPar["r0"], in, out], zVar] /. zVar -> r, 30];

analytic3D[piece_, r_, th_] := Module[{q = csPar["q"], c0 = csPar["c0"], c1 = csPar["c1"], r0 = csPar["r0"], rho, zc},
   rho = Sqrt[r^2 + r0^2 - 2 r r0 Cos[th]]; zc = r Cos[th];
   Switch[piece,
    "SS", q^2 Log[rho],
    "SR", q c0 rho + q c1 (r0 rho + (zc - r0) rho/2),
    "RR", c0^2 r^2/6 + c0 c1 zc r^2/5 + c1^2 zc^4/12]];

(* ================= Q = "GradDot" (Grad A . Grad B): the 1/rho^4 case (Stage 5) ================= *)
phiRAtParticle[] := csPar["c0"] + csPar["c1"] csPar["r0"] + csPar["c2"] csPar["r0"]^2;

(* SS: [1/rho^4]_l = -Sqrt[Pi(2l+1)] Q_l'(x)/(2 r^2 r0^2) (derivative of Neumann's formula) *)
legendreQNeumannDx[l_, x_, logTerm_, inv_] := Module[{t},
   ((1/2) D[LegendreP[l, t], t] logTerm - LegendreP[l, t] inv - D[WLegendre[l, t], t]) /. t -> x];
invX2m1[z_] := With[{r0 = csPar["r0"]}, 4 z^2 r0^2/((z - r0)^2 (z + r0)^2)];
ssSourceModeExpr["GradDot", l_, z_, "real"] := With[{r0 = csPar["r0"]},
   -csPar["q"]^2 Sqrt[Pi (2 l + 1)]/(2 z^2 r0^2) legendreQNeumannDx[l, xOfZ[z], 2 Log[z + r0] - 2 Log[Abs[z - r0]], invX2m1[z]]];
ssSourceModeExpr["GradDot", l_, z_, "continued"] := With[{r0 = csPar["r0"]},
   -csPar["q"]^2 Sqrt[Pi (2 l + 1)]/(2 z^2 r0^2) legendreQNeumannDx[l, xOfZ[z], 2 Log[z + r0] - 2 Log[z - r0], invX2m1[z]]];

(* SR couplings: d_r K and the angular coupling K^theta *)
dCoulombRatio[L_, r_] := Piecewise[{{L r^(L - 1)/csPar["r0"]^(L + 1), r < csPar["r0"]}}, -(L + 1) csPar["r0"]^L/r^(L + 2)];
dCoulombRatioContinued[L_, z_, "ThetaToH"] := With[{r0 = csPar["r0"]},
   L z^(L - 1)/r0^(L + 1) + (-(L + 1) r0^L/z^(L + 2) - L z^(L - 1)/r0^(L + 1)) HExt[z]];
dCoulombRatioContinued[L_, z_, "CorrectedSingular"] := Module[{t},
   -((2 L + 1)/csPar["r0"]^2) D[singularModeCorrected[L][t], t] /. t -> z];
angFactor[l_, lp_, L_] := (L (L + 1) + lp (lp + 1) - l (l + 1))/2;
KdCoupling[l_, lp_, r_] := Sum[gaunt[l, lp, L] Sqrt[4 Pi/(2 L + 1)] dCoulombRatio[L, r], {L, Abs[l - lp], l + lp}];
KthetaCoupling[l_, lp_, r_] := Sum[gaunt[l, lp, L] Sqrt[4 Pi/(2 L + 1)] angFactor[l, lp, L] coulombRatio[L, r], {L, Abs[l - lp], l + lp}];
KdCouplingContinued[l_, lp_, z_, method_] :=
  Sum[gaunt[l, lp, L] Sqrt[4 Pi/(2 L + 1)] dCoulombRatioContinued[L, z, method], {L, Abs[l - lp], l + lp}];
KthetaCouplingContinued[l_, lp_, z_, method_] :=
  Sum[gaunt[l, lp, L] Sqrt[4 Pi/(2 L + 1)] angFactor[l, lp, L] coulombRatioContinued[L, z, method], {L, Abs[l - lp], l + lp}];

(* SR: [2 q Grad(1/rho).Grad PhiR]_l = 2q Sum_lp [ (d_r K) d_r phiR_lp + K^theta phiR_lp / r^2 ];  phiRfuns = <|lp -> {f, df}|> *)
srSourceReal["GradDot", l_, phiRfuns_Association] := Module[{f, q = csPar["q"]},
   f[r_?NumericQ] := 2 q Sum[KdCoupling[l, lp, r] phiRfuns[lp][[2]][r] + KthetaCoupling[l, lp, r] phiRfuns[lp][[1]][r]/r^2, {lp, Keys[phiRfuns]}];
   f];
(* on the contour d_r phi~ = (d/da phi~)/y'(a); phiRtilde = <|lp -> function of a|> (InterpolatingFunction) *)
srSourceContour["GradDot", l_, phiRtilde_Association, method_] := Module[{f, q = csPar["q"], kd, kt, ypF},
   kd = AssociationMap[Function[lp, Module[{z}, contourFunction[KdCouplingContinued[l, lp, z, method], z]]], Keys[phiRtilde]];
   kt = AssociationMap[Function[lp, Module[{z}, contourFunction[KthetaCouplingContinued[l, lp, z, method]/z^2, z]]], Keys[phiRtilde]];
   ypF = Module[{a}, numericFunction[yCp[a], a]];
   f[x_?NumericQ] := 2 q Sum[kd[lp][x] phiRtilde[lp]'[x]/ypF[x] + kt[lp][x] phiRtilde[lp][x], {lp, Keys[phiRtilde]}];
   f];

(* RR and generic coupling: [Grad A.Grad B]_l = Sum <l l1 l2> (a'_l1 b'_l2 + a_l1 b_l2 (l1(l1+1)+l2(l2+1)-l(l+1))/(2 r^2)) *)
bilinearModes["GradDot", l_, aModes_Association, bModes_Association, r_] :=
  Sum[gaunt[l, l1, l2] (aModes[l1][[2]] bModes[l2][[2]] + angFactor[l, l1, l2] aModes[l1][[1]] bModes[l2][[1]]/r^2), {l1, Keys[aModes]}, {l2, Keys[bModes]}];
rrSourceReal["GradDot", l_, phiRfuns_Association] := Module[{f},
   f[r_?NumericQ] := With[{m = Map[{#[[1]][r], #[[2]][r]} &, phiRfuns]}, bilinearModes["GradDot", l, m, m, r]];
   f];

(* analytic particular solutions and closed-form source modes *)
pieceAnalyticSolution["GradDot", "SS", l_] := memo[{"solSS4", l},
   Table[(csPar["q"]^2/2) uModeIntegral[l, 0, "inv", zVar, br], {br, {"in", "out"}}]];
pieceAnalyticSource["GradDot", "SS", l_] := memo[{"srcSS4", l},
   Table[csPar["q"]^2 uModeIntegral[l, 0, "inv2", zVar, br], {br, {"in", "out"}}]];
(* SR: q (PhiR(x) - PhiR(x0))/rho;  PhiR - PhiR(x0) = c0 - PhiR0 + c1 r mu + c2 r^2 (3 mu^2 - 1)/2 *)
pieceAnalyticSolution["GradDot", "SR", l_] := memo[{"solSR4", l},
   With[{q = csPar["q"], c0 = csPar["c0"], c1 = csPar["c1"], c2 = csPar["c2"], p0 = phiRAtParticle[]},
    Table[Simplify[q ((c0 - p0 - c2 zVar^2/2) uModeIntegral[l, 0, "invSqrt", zVar, br] + c1 zVar uModeIntegral[l, 1, "invSqrt", zVar, br] +
         (3/2) c2 zVar^2 uModeIntegral[l, 2, "invSqrt", zVar, br])], {br, {"in", "out"}}]]];
(* SR source: 2q Grad(1/rho).Grad PhiR = -2q [ -c2 (x^2+y^2) + (z - r0)(c1 + 2 c2 z) ]/rho^3 *)
pieceAnalyticSource["GradDot", "SR", l_] := memo[{"srcSR4", l},
   With[{q = csPar["q"], c1 = csPar["c1"], c2 = csPar["c2"], r0 = csPar["r0"]},
    Module[{mu, cl},
     cl = CoefficientList[-c2 zVar^2 (1 - mu^2) + (zVar mu - r0) (c1 + 2 c2 zVar mu), mu];
     Table[Simplify[-2 q Sum[cl[[k + 1]] uModeIntegral[l, k, "invSqrt3", zVar, br], {k, 0, Length[cl] - 1}]], {br, {"in", "out"}}]]]];
pieceAnalyticSolution["GradDot", "RR", l_] := memo[{"solRR4", l},
   With[{c0 = csPar["c0"], c1 = csPar["c1"], c2 = csPar["c2"]},
    Module[{mu, cl}, cl = CoefficientList[Expand[(c0 + c1 zVar mu + c2 zVar^2 (3 mu^2 - 1)/2)^2/2], mu];
     ConstantArray[Simplify[Sum[cl[[k + 1]] muPowerMode[l, k], {k, 0, Length[cl] - 1}]], 2]]]];
pieceAnalyticSource["GradDot", "RR", l_] := memo[{"srcRR4", l},
   With[{c1 = csPar["c1"], c2 = csPar["c2"]},
    Module[{mu, cl}, cl = CoefficientList[Expand[c2^2 zVar^2 (1 - mu^2) + (c1 + 2 c2 zVar mu)^2], mu];
     ConstantArray[Simplify[Sum[cl[[k + 1]] muPowerMode[l, k], {k, 0, Length[cl] - 1}]], 2]]]];
analytic3D["GradDot", piece_, r_, th_] := Module[{q = csPar["q"], c0 = csPar["c0"], c1 = csPar["c1"], c2 = csPar["c2"], r0 = csPar["r0"], rho, phiR},
   rho = Sqrt[r^2 + r0^2 - 2 r r0 Cos[th]]; phiR = c0 + c1 r Cos[th] + c2 r^2 LegendreP[2, Cos[th]];
   Switch[piece, "SS", q^2/(2 rho^2), "SR", q (phiR - phiRAtParticle[])/rho, "RR", phiR^2/2]];

(* ================= Q = "Dzz": Q[A,B] = (1/2) d_z^2 (A B)  (Stage 9) ================= *)
(* d_z (f Y_l) = a_(l+1) (f' - l f/r) Y_(l+1) + a_l (f' + (l+1) f/r) Y_(l-1),  a_l = l/Sqrt[(2l-1)(2l+1)] *)
aZ[l_] := If[l <= 0, 0, l/Sqrt[(2 l - 1) (2 l + 1)]];
dzExpr[F_, l_] := If[l < 0, 0, If[l >= 1, aZ[l] (D[F[l - 1], zVar] - (l - 1) F[l - 1]/zVar), 0] + aZ[l + 1] (D[F[l + 1], zVar] + (l + 2) F[l + 1]/zVar)];
dzzExpr[F_, l_] := dzExpr[Function[j, dzExpr[F, j]], l];
(* delta(r - r0) coefficient: (d_z F)_j has a jump J_j = a_j [F'_(j-1)] + a_(j+1) [F'_(j+1)]; d_z of that jump gives a_l J_(l-1) + a_(l+1) J_(l+1) *)
dzzDeltaCoeff[Fpair_, l_] := Module[{J},
   J[j_] := If[j < 0, 0, Module[{e = {dzExpr[Function[i, If[i < 0, 0, Fpair[i][[1]]]], j], dzExpr[Function[i, If[i < 0, 0, Fpair[i][[2]]]], j]}},
      (e[[2]] - e[[1]]) /. zVar -> csPar["r0"]]];
   Simplify[If[l >= 1, aZ[l] J[l - 1], 0] + aZ[l + 1] J[l + 1]]];
(* all Dzz pieces are (1/2) d_z^2 of the Product pieces (Dzz SR = d_z^2(PhiP PhiR) = (1/2) d_z^2 (2 PhiP PhiR)); pointwise forms *)
pieceAnalyticSolution["Dzz", piece_, l_] := memo[{"solDzz", piece, l},
   Table[With[{br = br}, (1/2) dzzExpr[Function[j, If[j < 0, 0, pieceAnalyticSolution["Product", piece, j][[br]]]], l]], {br, 2}]];
pieceAnalyticSource["Dzz", piece_, l_] := memo[{"srcDzz", piece, l},
   Table[With[{br = br}, (1/2) dzzExpr[Function[j, If[j < 0, 0, pieceAnalyticSource["Product", piece, j][[br]]]], l]], {br, 2}]];
(* real branch: Log|r - r0| -> Log[(r - r0)^2]/2 so that it can be differentiated *)
ssSourceModeExpr["Dzz", l_, z_, br_] := (1/2) dzzExpr[Function[j, If[j < 0, 0, ssSourceModeExpr["Product", j, zVar, br] /. Log[Abs[x_]] :> Log[x^2]/2]], l] /. zVar -> z;

End[];
EndPackage[];
