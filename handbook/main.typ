#import "template.typ": *

#show: handbook.with(
  title: "From Optimal Transport to Trajectory Inference",
  subtitle: "Entropic transport, Schrödinger bridges and mean-field Langevin dynamics: a handbook for readers who know OT",
  author: "Eren Can",
  date: "October 2026",
  abstract: [
    Trajectory inference (TI) asks how individuals in a population move over time when the population can only be observed through _snapshots_, as in single-cell RNA sequencing, where measuring a cell destroys it. Mathematically, TI is a problem of selecting a law on path space from a sequence of its marginals, and optimal transport (OT) enters at exactly that point: OT is a principled way of choosing couplings. This handbook develops the theory of TI from the perspective of a reader who already knows OT and Wasserstein gradient flows. We go from Kantorovich's problem to entropic OT and Sinkhorn's algorithm, from there to the static and dynamic Schrödinger problems, and arrive at the min-entropy estimator of Lavenant et al., a multi-marginal Schrödinger problem that is consistent for drift-diffusion processes. We then explain how Chizat et al. compute it with mean-field Langevin dynamics, a Wasserstein gradient flow on a product of probability spaces whose global convergence rests on _flat_ rather than _displacement_ convexity, and how unbalanced transport (Wasserstein–Fisher–Rao / Hellinger–Kantorovich) accounts for cell division and death. Throughout, "Connection" boxes make explicit the bridges between these topics, and the last section collects open problems where the OT and gradient-flow viewpoints are likely to be decisive.
  ],
)

= Why trajectory inference is a transport problem <sec:intro>

== The data: snapshots of a population

Let $cal(X) subset bb(R)^d$ be a state space (in biology: gene-expression space, typically after dimensionality reduction) and $Omega = C([0,1], cal(X))$ the space of continuous paths. A population of particles evolves according to an unknown stochastic process $(X_t)_(t in [0,1])$ with law $P in cal(P)(Omega)$. Writing $e_t : omega |-> omega(t)$ for the evaluation map, the time marginals are $P_t = (e_t)_(\#) P$ and the two-time couplings are $P_(s,t) = (e_s, e_t)_(\#) P$.

#definition(title: "Trajectory inference")[
  Given times $0 <= t_1 < dots.c < t_T <= 1$ and empirical measures
  $ hat(mu)_i = 1/n_i sum_(j=1)^(n_i) delta_(X_(t_i, j)), quad i = 1, dots, T, $
  where all samples are independent (*each particle is observed at most once*), estimate the law $P in cal(P)(Omega)$.
]

The final clause is what separates TI from classical statistics of stochastic processes, where one observes whole (discretised) trajectories. In single-cell RNA sequencing the measurement is destructive, so the data carry information about marginals only @schiebinger2019optimal @weinreb2018fundamental.

== The missing information is a coupling

Marginals never determine the law of a process. Both of the following have $P_t = cal(N)(0,1)$ for every $t$:
(a) the static process $X_t equiv X_0 tilde cal(N)(0,1)$;
(b) the stationary Ornstein–Uhlenbeck process $dif X_t = -X_t dif t + sqrt(2) dif B_t$.
Their couplings $P_(s,t)$ are, respectively, the identity coupling and a Gaussian coupling with correlation $e^(-|t-s|)$. More generally, for a Markov process the law is determined by its marginals _together with_ its transition kernels, so the unknown in TI is exactly the sequence of couplings $(P_(t_i, t_(i+1)))_i$, or rather a choice among the many couplings compatible with the data.

#connection(title: "TI is coupling selection")[
  Choosing a coupling $pi in Pi(mu, nu)$ between two given marginals by a variational principle is precisely what optimal transport does. TI therefore "is" an OT problem as soon as one fixes a principle for picking couplings. The whole theory below can be read as a sequence of increasingly faithful answers to the question _which_ variational principle is the right one for noisy, branching populations, and _how_ to compute its solution.
]

== A first answer: Waddington-OT

The first OT-based method, Waddington-OT @schiebinger2019optimal, couples each pair of consecutive snapshots by (unbalanced, entropically regularised) OT, $gamma_i in argmin_(gamma in Pi(hat(mu)_i, hat(mu)_(i+1))) integral c dif gamma + epsilon H(gamma)$, and composes the resulting transition kernels into trajectories. It was strikingly successful on stem-cell reprogramming data. However, it was a heuristic: why entropic regularisation, why this cost, how to choose $epsilon$, and what is estimated when samples are scarce?

#cite(<lavenant2021towards>, form: "prose") answered these questions by identifying the right object: the *min-entropy estimator*, a law on paths minimising relative entropy with respect to Brownian motion while fitting the snapshots. It is consistent for drift-diffusion processes, and it reduces to Waddington-OT in a limiting regime (@sec:ti). #cite(<chizat2022trajectory>, form: "prose") then showed how to compute it with a grid-free, provably convergent particle method (@sec:mfl). The route there passes through a chain of classical ideas from OT:

#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: (left, left, left),
    stroke: none,
    inset: (x: 6pt, y: 4pt),
    table.hline(stroke: 0.8pt),
    table.header([*Step*], [*Mathematical object*], [*Role in TI*]),
    table.hline(stroke: 0.5pt),
    [@sec:ot], [Kantorovich OT, Benamou–Brenier, Otto calculus], [zero-noise coupling principle; geometry for gradient flows],
    [@sec:eot], [entropic OT, Sinkhorn, Schrödinger potentials], [noise-aware couplings; algorithms; first variations],
    [@sec:sb], [static and dynamic Schrödinger problem], [the "most likely" evolution of a Brownian cloud],
    [@sec:ti], [multi-marginal Schrödinger problem with soft constraints], [the min-entropy estimator],
    [@sec:mfl], [mean-field Langevin dynamics], [a convergent gradient-flow solver],
    [@sec:uot], [unbalanced OT (WFR / HK)], [cell division and death],
    table.hline(stroke: 0.8pt),
  ),
  caption: [Roadmap of the handbook.],
) <tab:roadmap>

#figure(
  image("figures/paths.png", width: 100%),
  caption: [Trajectory inference on a simulated 10-dimensional bifurcating process (from the companion notebook). _Left:_ hidden ground-truth paths (grey) and the observed snapshots, coloured by time. _Middle:_ marginals and trajectories reconstructed by the min-entropy estimator, computed with mean-field Langevin dynamics. _Right:_ first coordinate over time, ground truth vs. inferred sample paths.],
) <fig:paths>

= Optimal transport essentials <sec:ot>

This section fixes notation and recalls the facts used later; see @villani2009optimal @santambrogio2015optimal @peyre2019computational for details. Throughout, $c(x,y) = 1/2 norm(x-y)^2$ unless stated otherwise.

== Static OT and duality

For $mu, nu in cal(P)_2(bb(R)^d)$, the Kantorovich problem and its dual read
$ OT(mu, nu) = min_(pi in Pi(mu, nu)) integral c dif pi = sup_(phi(x) + psi(y) <= c(x,y)) integral phi dif mu + integral psi dif nu. $ <eq:ot>
At the optimum one may take $psi = phi^c$, where $phi^c (y) = inf_x c(x,y) - phi(x)$ is the $c$-transform. For the quadratic cost and $mu$ absolutely continuous, Brenier's theorem gives a unique optimal plan $pi = ("id", T)_(\#) mu$ with $T = nabla u$, $u = 1/2 norm(dot)^2 - phi$ convex. In particular $T(x) = x - nabla phi(x)$, an identity that reappears in entropic form in @sec:eot.

== Dynamic OT and the geometry of $cal(P)_2$

Benamou and Brenier @benamou2000computational showed that
$ W_2^2(mu, nu) = min_((rho, v)) integral_0^1 integral norm(v_t)^2 dif rho_t dif t quad "s.t." quad partial_t rho_t + nabla dot (rho_t v_t) = 0, quad rho_0 = mu, thick rho_1 = nu, $ <eq:bb>
with $W_2^2 = 2 OT$. The minimiser is the _displacement interpolation_ $rho_t = ((1-t) "id" + t T)_(\#) mu$ @mccann1997convexity: particles move on straight lines at constant speed and never cross.

Otto's formal Riemannian calculus @otto2001geometry @ambrosio2008gradient reads @eq:bb as a geodesic problem: tangent vectors at $rho$ are velocity fields $v = nabla phi$, with metric $norm(v)_rho^2 = integral norm(v)^2 dif rho$. The Wasserstein gradient of a functional $cal(F)$ is $nabla_(W_2) cal(F)(rho) = nabla (delta cal(F)) / (delta rho)$, where $(delta cal(F)) / (delta rho)$ is the first variation, and the *Wasserstein gradient flow* is
$ partial_t rho_t = nabla dot (rho_t nabla (delta cal(F))/(delta rho)[rho_t]). $ <eq:wgf>
Its time discretisation is the JKO scheme @jordan1998variational, $rho_(k+1) in argmin_rho W_2^2(rho, rho_k)/(2h) + cal(F)(rho)$.

#example(title: "Fokker–Planck as a gradient flow")[
  For $cal(F)(rho) = integral Psi dif rho + tau/2 integral rho log rho$, @eq:wgf is the Fokker–Planck equation $partial_t rho = nabla dot (rho nabla Psi) + tau/2 Delta rho$, i.e. the law of the Langevin diffusion $dif X_t = -nabla Psi(X_t) dif t + sqrt(tau) dif B_t$ @jordan1998variational. If $Psi$ is $lambda$-convex, $cal(F)$ is $lambda$-displacement convex @mccann1997convexity and the flow contracts in $W_2$ at rate $lambda$, the starting point of the evolution-variational-inequality (EVI) theory of @ambrosio2008gradient.
]

#connection(title: "OT already is a model of trajectories")[
  Displacement interpolation is a trajectory model: it reconstructs a deterministic, cyclically monotone flow between two snapshots. It is the *zero-noise* answer to TI. Its defects for biology are just as instructive: trajectories never cross, nothing is random, the interpolation between three or more snapshots is a broken geodesic with no memory of the dynamics, and the estimation of $W_2$ (and of its maps) from samples suffers from the curse of dimensionality, at rate $n^(-1\/d)$. Each defect is cured by adding noise, which is where entropy enters.
]

= Entropic OT and Sinkhorn's algorithm <sec:eot>

== Definition, dual and the Schrödinger system

#definition(title: "Entropic optimal transport")[
  For $epsilon > 0$,
  $ OT_epsilon (mu, nu) = min_(pi in Pi(mu, nu)) integral c dif pi + epsilon H(pi | mu times.o nu), quad H(pi | rho) = integral log (dif pi)/(dif rho) dif pi. $ <eq:eot>
]

Strict convexity gives a unique optimiser. Lagrangian duality yields an unconstrained concave dual @peyre2019computational @nutz2021introduction:
$ OT_epsilon (mu, nu) = sup_(phi, psi) integral phi dif mu + integral psi dif nu - epsilon integral e^((phi(x) + psi(y) - c(x,y))\/epsilon) dif mu(x) dif nu(y) + epsilon, $ <eq:eot-dual>
and the primal–dual relation $dif pi^star = e^((phi plus.o psi - c)\/epsilon) dif (mu times.o nu)$. The optimality conditions are the *Schrödinger system*
$ phi(x) = -epsilon log integral e^((psi(y) - c(x,y))\/epsilon) dif nu(y), quad psi(y) = -epsilon log integral e^((phi(x) - c(x,y))\/epsilon) dif mu(x). $ <eq:schrodinger-system>
The right-hand sides are *soft $c$-transforms*: the softmin $-epsilon log integral e^(-(dot)\/epsilon)$ replaces the infimum in $phi^c$. As $epsilon -> 0$ they converge to the $c$-transform, and @eq:eot-dual to @eq:ot.

*Sinkhorn's algorithm* @sinkhorn1964relationship @cuturi2013sinkhorn alternates the two equations of @eq:schrodinger-system. Equivalently, it is block-coordinate ascent on @eq:eot-dual, or iterative proportional fitting: alternating KL projections of the Gibbs kernel $e^(-c\/epsilon) mu times.o nu$ onto the two marginal constraints. For $m$ points each iteration costs $O(m^2)$ and is embarrassingly parallel; in practice one runs it in the log domain for stability.

== Three properties that matter for TI

#proposition(title: "Smooth potentials and barycentric projection")[
  The Schrödinger potentials defined by @eq:schrodinger-system extend to smooth functions on all of $cal(X)$, and for the quadratic cost
  $ nabla phi(x) = integral nabla_x c(x,y) dif pi^star (y | x) = x - bb(E)_(pi^star) [Y | X = x]. $ <eq:barycentric>
]
Thus $x |-> x - nabla phi(x)$ is the *barycentric projection* of the entropic plan, the entropic analogue of the Brenier map $T = "id" - nabla phi$. It is a statistically efficient estimator of $T$ @pooladian2021entropic.

#proposition(title: "First variation")[
  $(delta OT_epsilon) / (delta mu) (mu, nu) = phi$, the Schrödinger potential (up to an additive constant) @feydy2019interpolating.
]
Combined with @eq:barycentric, this means that the *Wasserstein gradient of $mu |-> OT_epsilon (mu, nu)$ at a particle $x$ points from $x$ to the barycentre of where $x$ is sent*. This is the force that couples the particle clouds in @sec:mfl.

#proposition(title: "Statistics and limits")[
  (i) As $epsilon -> 0$, $OT_epsilon -> OT$ and the optimal plans converge to an optimal transport plan @carlier2017convergence @leonard2012schrodinger. As $epsilon -> infinity$, $pi^star -> mu times.o nu$.
  (ii) Plug-in estimation of $OT_epsilon$ from $n$ samples has error of order $n^(-1\/2)$, independently of the dimension (with constants deteriorating as $epsilon -> 0$) @genevay2019sample @mena2019statistical. Unregularised OT has the $n^(-1\/d)$ rate.
]

#figure(
  image("figures/eot_plans.svg", width: 100%),
  caption: [Entropic plans between two 1-D Gaussian mixtures (log-density). Small $epsilon$: nearly the monotone Brenier map. Large $epsilon$: approaching the independent coupling $mu times.o nu$.],
) <fig:plans>

#connection(title: "Entropy is noise, not just a numerical trick")[
  In machine learning, entropic regularisation is often introduced as a device to make OT smooth, fast and differentiable @cuturi2013sinkhorn. For TI it has a _modelling_ meaning: the entropic plan with $epsilon = tau h$ is exactly the joint law at times $0$ and $h$ of Brownian particles with diffusivity $tau$, conditioned to start from $mu$ and arrive at $nu$. This is Schrödinger's problem, the subject of the next section. The regularisation parameter is therefore not a tuning knob but the _noise level_ of the biological process times the time gap. Moreover, the dimension-free sample complexity (ii) is the statistical reason why entropic methods remain usable with a few hundred cells in high dimension.
]

= The Schrödinger problem <sec:sb>

== Schrödinger's thought experiment

In 1931–32 Schrödinger asked the following question @schrodinger1932theorie @leonard2014survey. Release $N$ independent Brownian particles with diffusivity $epsilon$. At time 0 their empirical distribution is close to $mu$, and at time 1 one observes, to one's surprise, an empirical distribution close to $nu$ instead of the heat-flow prediction. *What is the most likely way the cloud got there?*

Let $R in cal(P)(Omega)$ be the law of one reference Brownian path, and $L^N = 1/N sum_k delta_(X^k)$ the empirical measure of the $N$ paths. By Sanov's theorem, $bb(P)(L^N approx P) approx exp(-N H(P | R))$. Conditioning on the atypical event $\{L^N_0 approx mu, L^N_1 approx nu\}$, the empirical law of paths concentrates on the solution of
$ min_(P in cal(P)(Omega)) H(P | R) quad "s.t." quad P_0 = mu, thick P_1 = nu. $ <eq:dsb>
This is the *dynamic Schrödinger problem*, and its solution the *Schrödinger bridge* @follmer1988random @leonard2014survey.

#connection(title: "The estimator as a conditional large-deviation principle")[
  Read with biological eyes, @eq:dsb is a statement about trajectory inference. Under the null model "cells are independent Brownian particles", the most likely population-level explanation of two observed snapshots is the Schrödinger bridge. The min-entropy estimator of @sec:ti extends exactly this principle to many noisy snapshots.
]

== From dynamic to static: the bridge decomposition

Relative entropy decomposes along the endpoint map: for any $P$,
$ H(P | R) = H(P_(01) | R_(01)) + integral H(P^(x y) | R^(x y)) dif P_(01)(x, y), $ <eq:additivity>
where $P^(x y)$ and $R^(x y)$ are the laws conditioned on $X_0 = x, X_1 = y$. The second term vanishes iff $P$ uses the reference bridges, $P = integral R^(x y) dif P_(01)$. Hence @eq:dsb reduces to the *static Schrödinger problem* $min_(pi in Pi(mu, nu)) H(pi | R_(01))$, and the dynamic solution is obtained by gluing Brownian bridges onto the static plan.

#proposition(title: "Static Schrödinger = entropic OT")[
  Let $R$ be Brownian motion with diffusivity $epsilon$ started from Lebesgue measure, so that $dif R_(01)(x,y) prop e^(-norm(x-y)^2 \/ (2 epsilon)) dif x dif y$. For $pi in Pi(mu, nu)$,
  $ epsilon H(pi | R_(01)) = integral 1/2 norm(x - y)^2 dif pi + epsilon H(pi | mu times.o nu) + epsilon (H(mu) + H(nu)) + "const", $
  with $H(mu) = integral log mu dif mu$. Hence the static Schrödinger problem and $OT_epsilon (mu, nu)$ have the same minimiser.
]
#proof[Write $(dif pi)/(dif R_(01)) = (dif pi)/(dif (mu times.o nu)) dot (dif (mu times.o nu))/(dif x dif y) dot (dif x dif y)/(dif R_(01))$ and integrate the logarithm against $pi$; the middle factor contributes $H(mu) + H(nu)$, the last one $norm(x-y)^2\/(2 epsilon)$ plus a constant.]

On a manifold, or with a non-Gaussian reference, the cost becomes $c_epsilon = -epsilon log p_epsilon$, with $p_epsilon$ the heat kernel. By Varadhan's formula $c_epsilon -> d^2\/2$ as $epsilon -> 0$, which is the cost used in the TI papers.

== Stochastic control and the fluid picture

Let $R^epsilon$ be the law of $dif X_t = sqrt(epsilon) dif B_t$. If $H(P | R^epsilon) < infinity$, then under $P$ the canonical process solves $dif X_t = b_t (X_t) dif t + sqrt(epsilon) dif B_t$ for some drift $b$, and by Girsanov's theorem @follmer1988random
$ H(P | R^epsilon) = H(P_0 | R^epsilon_0) + 1/(2 epsilon) bb(E)_P integral_0^1 norm(b_t (X_t))^2 dif t. $ <eq:girsanov>
The Schrödinger problem is therefore a *stochastic optimal control* problem, which is to steer a diffusion from $mu$ to $nu$ with minimal kinetic energy of the control @daipra1991stochastic @chen2016relation @chen2021stochastic. Replacing the drift by the _current velocity_ $v_t = b_t - epsilon/2 nabla log rho_t$, which satisfies the continuity equation $partial_t rho + nabla dot (rho v) = 0$, and integrating the cross term (it is a time derivative of $integral rho_t log rho_t$) gives the following fluid formulation @chen2016relation @gentil2017analogy.

#theorem(title: "Entropic Benamou–Brenier formula")[
  Up to the additive constant $epsilon/2 (H(nu) - H(mu))$, the value of @eq:dsb times $epsilon$ equals
  $ min_((rho, v)) integral_0^1 integral (1/2 norm(v_t)^2 + epsilon^2/8 norm(nabla log rho_t)^2) dif rho_t dif t quad "s.t." quad partial_t rho_t + nabla dot (rho_t v_t) = 0, thick rho_0 = mu, thick rho_1 = nu. $ <eq:ebb>
]

The second term is $epsilon^2\/8$ times the *Fisher information*, so the Schrödinger bridge is dynamic OT with a penalty on the concentration of the interpolating densities. As $epsilon -> 0$, @eq:ebb $Gamma$-converges to @eq:bb, the Schrödinger bridges converge to the displacement interpolation, and the static plans to the Brenier plan @mikami2004monge @leonard2012schrodinger. Much of displacement convexity theory has an "entropic interpolation" analogue @gentil2017analogy.

#figure(
  image("figures/interpolations.svg", width: 100%),
  caption: [Displacement interpolation vs. Schrödinger bridges between the same 1-D marginals (space-time densities, with a few sample paths). The bridge is a mixture of Brownian bridges glued to the entropic plan, $rho_t = integral cal(N)((1-t)x + t y, thin epsilon t(1-t)) dif pi^star (x,y)$, so it widens in the middle of the interval and its paths cross. As $epsilon -> 0$ it collapses onto the geodesic.],
) <fig:interp>

== Dictionary

#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: left,
    stroke: none,
    inset: (x: 5pt, y: 4pt),
    table.hline(stroke: 0.8pt),
    table.header([], [*Optimal transport* ($epsilon = 0$)], [*Schrödinger / entropic OT* ($epsilon > 0$)]),
    table.hline(stroke: 0.5pt),
    [Optimal plan], [deterministic map $("id", nabla u)_(\#) mu$], [Gibbs plan $e^((phi plus.o psi - c)\/epsilon) mu times.o nu$],
    [Dual operator], [$c$-transform (inf)], [soft $c$-transform (softmin)],
    [Optimality], [Monge–Ampère equation], [Schrödinger system @eq:schrodinger-system],
    [Map], [Brenier map $x - nabla phi(x)$], [barycentric projection $bb(E)[Y | X = x]$],
    [Interpolation], [straight lines, McCann], [Brownian bridges, entropic interpolation],
    [Dynamic form], [Benamou–Brenier @eq:bb], [@eq:ebb (+ Fisher information) / stochastic control],
    [Algorithm], [LP, auction, semi-discrete Newton], [Sinkhorn, $O(m^2)$ per iteration, GPU-friendly],
    [Sample complexity], [$n^(-1\/d)$], [$n^(-1\/2)$ (constants exponential in $1\/epsilon$)],
    [Probabilistic meaning], [—], [most likely evolution of Brownian particles (Sanov)],
    table.hline(stroke: 0.8pt),
  ),
  caption: [Correspondences between OT and the Schrödinger problem.],
) <tab:dictionary>

= Trajectory inference as a multi-marginal Schrödinger problem <sec:ti>

== The population model and why it is a gradient flow

The standard model for a population without births or deaths is the drift-diffusion
$ dif X_t = -nabla Psi(t, X_t) dif t + sqrt(tau) dif B_t, quad X_0 tilde mu_0, $ <eq:sde>
where $Psi$ is the (unknown) *Waddington landscape* @waddington1957strategy @hashimoto2016learning and $tau > 0$ a known diffusivity. Only the gradient part of a drift is identifiable from marginals: a divergence-free component that leaves $rho_t$ invariant changes the trajectories but not the snapshots @weinreb2018fundamental @lavenant2021towards. This is why one assumes a gradient drift.

#connection(title: "Cells follow a Wasserstein gradient flow")[
  By the Fokker–Planck example of @sec:ot, the marginals of @eq:sde are the Wasserstein gradient flow of the free energy $cal(F)_t (rho) = integral Psi(t, dot) dif rho + tau/2 integral rho log rho$ (non-autonomous if $Psi$ depends on $t$). Recovering the dynamics is therefore an *inverse problem for a gradient flow*: given a few points along a gradient-flow curve, recover the energy (the landscape) and the stochastic trajectories. Parametric approaches such as JKOnet @bunne2022proximal take this literally and fit $Psi_theta$ so that JKO steps map each snapshot to the next.
]

== The min-entropy estimator

Let $W^tau$ be the law of (reflected) Brownian motion with diffusivity $tau$. #cite(<lavenant2021towards>, form: "prose") define
$ R^star = argmin_(R in cal(P)(Omega)) 1/lambda sum_(i=1)^T Delta t_i Fit_sigma (R_(t_i), hat(mu)_i) + tau H(R | W^tau), $ <eq:estimator>
with quadrature weights $Delta t_i$, a regularisation $lambda > 0$ and a data-fitting functional, e.g. the negative log-likelihood under Gaussian observation noise used in @chizat2022trajectory:
$ Fit_sigma (rho, hat(mu)) = -integral log (integral e^(-norm(x-y)^2 \/ (2 sigma^2)) dif rho(x)) dif hat(mu)(y). $
This is a *soft-constrained, multi-marginal Schrödinger problem*. With hard constraints $R_(t_i) = mu_i$ it would be the Schrödinger problem with $T$ prescribed marginals; the fitting term relaxes these constraints because the $hat(mu)_i$ are noisy, finite samples.

#theorem(title: [Consistency; #cite(<lavenant2021towards>, form: "prose")])[
  If the ground truth is of the form @eq:sde (on a compact manifold, with suitable regularity), and the measurement times become dense in $[0,1]$, then $R^star -> P$ as $T -> infinity$ followed by $lambda, sigma -> 0$, weakly and almost surely.
]

Why does minimising entropy recover a _drift_? By @eq:girsanov, $tau H(R | W^tau)$ is, up to the initial term, the kinetic energy $1/2 bb(E) integral norm(b_t)^2 dif t$ of the drift. The estimator is thus a *least-action principle*: among all drift-diffusions with noise $tau$ that explain the snapshots, choose the one with the least kinetic energy. For gradient drifts that are densely observed, the least-action explanation is the truth.

== The representer theorem: from paths to $T$ marginals

Because $W^tau$ is Markov, @eq:additivity generalises to a chain. For any $R$ with marginals $mu_i = R_(t_i)$,
$ H(R | W^tau) >= sum_(i=1)^(T-1) H(R_(t_i, t_(i+1)) | p_(tau_i) mu_i times.o mu_(i+1)) + sum_(i=1)^T H(mu_i) + C, quad tau_i = tau (t_(i+1) - t_i), $
with equality iff $R$ is Markov and uses Brownian bridges between the observation times @chizat2022trajectory. Optimising over the couplings with marginals fixed yields entropic OT terms, and we arrive at the *reduced problem* over $bold(mu) = (mu_1, dots, mu_T) in cal(P)(cal(X))^T$:
$ F(bold(mu)) = underbrace(1/lambda sum_i Delta t_i Fit_sigma (mu_i, hat(mu)_i) + sum_(i=1)^(T-1) (OT_(tau_i)(mu_i, mu_(i+1)))/(t_(i+1) - t_i), G(bold(mu))) + tau sum_(i=1)^T H(mu_i). $ <eq:reduced>

#theorem(title: [Representer theorem; #cite(<chizat2022trajectory>, form: "prose")])[
  $F$ is jointly convex with a unique minimiser $bold(mu)^star$, and
  $R^star = integral W^tau (dot | x_1, dots, x_T) dif (gamma_(1,2) gamma_(2,3) dots.c gamma_(T-1,T))(x_1, dots, x_T)$,
  where $gamma_(i,i+1)$ are the optimal entropic plans between $mu_i^star$ and $mu_(i+1)^star$. In words: optimal marginals, chained by Schrödinger couplings, filled in with Brownian bridges.
]

#connection(title: "Three readings of the reduced objective")[
  + *Waddington-OT as a limit.* If the marginals are pinned to the data ($lambda -> 0$ with exact marginals), $F$ reduces to a sum of pairwise entropic OT problems between consecutive snapshots: exactly Waddington-OT @schiebinger2019optimal. WOT is the _interpolating_ special case; the min-entropy estimator is its _smoothing_ (regression) version, which also corrects the marginals.
  + *Multi-marginal OT with a path-structured cost.* Equivalently, $R_(t_1, dots, t_T)$ solves an entropic multi-marginal problem (regularisation $tau$) with cost $sum_i norm(x_(i+1) - x_i)^2 \/ (2 (t_(i+1) - t_i))$. Tree-structured costs keep multi-marginal Sinkhorn tractable @haasler2021multimarginal.
  + *A smoothing spline in Wasserstein space.* $sum_i OT_(tau_i)(mu_i, mu_(i+1)) \/ Delta t$ is a discrete, entropic version of the action $integral norm(v)^2 dif rho dif t$ in @eq:bb. So @eq:reduced is a data term plus a _first-order_ (piecewise geodesic) roughness penalty: a smoothing spline on $(cal(P)_2, W_2)$, blurred by noise. Second-order (cubic) splines on Wasserstein space @benamou2019second @chewi2021fast would penalise acceleration instead. They form a natural, largely unexplored alternative prior for TI.
]

Note the *positive* sign of $tau sum H(mu_i)$ in @eq:reduced: the marginal entropies penalise concentration, which is what makes $F$ jointly convex even though the coupling terms are only separately convex @chizat2022trajectory.

= Two gradient flows: mean-field Langevin dynamics <sec:mfl>

== Two clocks

There are two unrelated "times" in this story, and two gradient flows:

#figure(
  table(
    columns: (auto, 1fr, 1fr),
    align: left, stroke: none, inset: (x: 5pt, y: 4pt),
    table.hline(stroke: 0.8pt),
    table.header([], [*Physical time* $t$], [*Optimisation time* $s$]),
    table.hline(stroke: 0.5pt),
    [Space], [$cal(P)(cal(X))$], [$cal(P)(cal(X))^T$ (one cloud per snapshot)],
    [Energy], [free energy $integral Psi dif rho + tau/2 H(rho)$ (unknown)], [reduced objective $F$ of @eq:reduced (known)],
    [Flow], [Fokker–Planck of the cells], [mean-field Langevin dynamics],
    [Role], [the model to be inferred], [the algorithm that infers it],
    table.hline(stroke: 0.8pt),
  ),
  caption: [The two gradient flows in trajectory inference.],
)

== Mean-field Langevin dynamics

For a functional of the form $F(mu) = G(mu) + tau H(mu)$, the Wasserstein gradient flow @eq:wgf is
$ partial_s mu_s = nabla dot (mu_s nabla (delta G)/(delta mu)[mu_s]) + tau Delta mu_s, $ <eq:mfl>
the law of the McKean–Vlasov diffusion $dif X_s = -nabla (delta G)/(delta mu)[mu_s](X_s) dif s + sqrt(2 tau) dif B_s$, $mu_s = "Law"(X_s)$. This is *mean-field Langevin (MFL) dynamics*. It was introduced to study wide two-layer neural networks @mei2018mean @hu2021mean, and it is simulated with $m$ interacting particles and noisy gradient steps.

For TI, #cite(<chizat2022trajectory>, form: "prose") apply MFL on the product space $cal(P)(cal(X))^T$, one species per snapshot. Using the first variations of @sec:eot, the potential acting on species $i$ is
$ V^((i))[bold(mu)] = (delta Fit)/(delta mu_i) + phi_(i,i+1) / (t_(i+1) - t_i) + psi_(i-1,i) / (t_i - t_(i-1)), $
where $(phi_(i,i+1), psi_(i,i+1))$ are the Schrödinger potentials of $OT_(tau_i)(mu_i, mu_(i+1))$. By @eq:barycentric each particle feels three forces: a soft-min attraction to nearby observations at $t_i$, a pull towards the barycentre of its expected successors at $t_(i+1)$, and a pull towards the barycentre of its expected predecessors at $t_(i-1)$. The $T$ clouds behave like *beads on an elastic chain in Wasserstein space*, kept apart by diffusion and pinned by the data. One step of the particle scheme costs one warm-started Sinkhorn solve per pair of consecutive clouds.

#figure(
  image("figures/mfl_dynamics.png", width: 100%),
  caption: [MFL on the bifurcation example (companion notebook). _Left:_ reduced objective $F$ along the iterations (entropy estimated by $k$-NN). _Right:_ particle clouds at iterations 1, 50, 500, 2500, coloured by physical time $t$.],
) <fig:mfl>

== Why it converges: flat convexity meets Wasserstein dissipation

The convergence theory of MFL is instructive for a gradient-flow audience, because it does *not* rely on displacement convexity. Let $hat(mu) prop e^(-(delta G)/(delta mu)[mu] \/ tau)$ be the *proximal Gibbs measure* of $mu$. Its fixed points are exactly the minimisers of $F$.

#lemma(title: "Entropy sandwich")[
  If $G$ is convex along linear (flat) interpolations $(1-s)mu + s nu$, then
  $F(mu) - min F <= tau H(mu | hat(mu)).$
]
#proof[By flat convexity, $G(nu) >= G(mu) + integral (delta G)/(delta mu)[mu] dif (nu - mu)$. Adding $tau H(nu)$ and minimising over $nu$ with the Gibbs variational principle $min_nu integral V dif nu + tau H(nu) = -tau log integral e^(-V\/tau)$ gives a lower bound on $min F$. Comparing with $F(mu)$ yields exactly $tau H(mu | hat(mu))$.]

Along @eq:mfl, the energy dissipates at the rate of the relative Fisher information,
$ dif/(dif s) F(mu_s) = - integral norm(nabla (delta F)/(delta mu))^2 dif mu_s = -tau^2 I(mu_s | hat(mu)_s). $
If the proximal Gibbs measures satisfy a *log-Sobolev inequality* with a uniform constant $C_"LS"$, i.e. $H(mu | hat(mu)) <= 1/(2 C_"LS") I(mu | hat(mu))$, then combining the two displays gives
$ dif/(dif s) (F(mu_s) - min F) <= -2 C_"LS" tau (F(mu_s) - min F), $
hence $F(mu_s) - min F <= e^(-2 C_"LS" tau s) (F(mu_0) - min F)$.
This is the argument of @nitanda2022convex @chizat2022mean.

#theorem(title: [Global convergence for TI; #cite(<chizat2022trajectory>, form: "prose")])[
  On the torus, with an extra noise $epsilon > 0$ (MFL at temperature $tau + epsilon$), the dynamics converges exponentially to the unique minimiser of $F_epsilon = G + (tau + epsilon) H$, at a rate $C = beta e^(-alpha\/epsilon)$. With a slowly annealed $epsilon_s -> 0$ it converges to $bold(mu)^star$, and thus, by the representer theorem, to the min-entropy estimator.
]

Why the extra noise $epsilon$? In TI the smooth part $G$ is only _separately_ convex in each $mu_i$; it is $G + tau H$ that is jointly convex (the positive marginal entropies make it so). The argument above is therefore applied to the convex functional $G + tau H$ with $epsilon H$ playing the role of the entropy: the temperature in the LSI is $epsilon$, not $tau$. The uniform LSI comes from Holley–Stroock perturbation @holley1987logarithmic. The oscillation of $V^((i))$ is bounded uniformly, by the oscillation bound for Schrödinger potentials and the boundedness of the fitting gradient, so the constant degrades like $e^(-"osc"(V)\/epsilon)$. This is the price of non-convex landscapes, familiar from simulated annealing.

#connection(title: "Convexity in one geometry, dissipation in another")[
  Classical gradient-flow theory @ambrosio2008gradient pairs a metric with convexity _in the same metric_: $lambda$-displacement convexity yields EVI, contraction and rates. Here the objective is convex along _linear_ interpolations (for TI: $G + tau H$, built from soft-min fits and entropic OT terms), and nothing is known about its displacement convexity. MFL pairs flat convexity, which controls the optimality gap (sandwich lemma), with Wasserstein dissipation, which controls its decrease (Fisher information). The two are linked by a functional inequality (LSI), much as in the Otto–Villani theorem linking LSI and transport inequalities @otto2000generalization. This mixed-geometry mechanism is the conceptual core of the guarantee, and it is a natural place to bring other geometries into play (@sec:uot).
]

#figure(
  image("figures/few_samples.png", width: 100%),
  caption: [Grid-free reconstruction from scarce data (companion notebook). _Top:_ one observed cell per intermediate time point. _Bottom:_ 64 cells. MFL particles move freely in $cal(X)$, so information propagates from well-sampled time points through the Schrödinger couplings. A fixed-support discretisation can only place mass on observed points.],
) <fig:few>

= Mass is not conserved: unbalanced transport and growth <sec:uot>

== Growth in the population model

Cells divide and die. With a (state- and time-dependent) net growth rate $g$, the population density satisfies a reaction–drift–diffusion equation
$ partial_t rho = nabla dot (rho nabla Psi) + tau/2 Delta rho + g rho. $ <eq:growth>
Snapshots are normalised, so they show the *relative* abundance of states. A fast-dividing subpopulation looks like a population that cells migrate _into_. Conflating proliferation with migration is the central failure mode that unbalanced methods address.

== Unbalanced optimal transport

Unbalanced OT relaxes mass conservation. In dynamic form, the *Wasserstein–Fisher–Rao* (WFR), or *Hellinger–Kantorovich*, distance @kondratyev2016new @chizat2018interpolating @liero2018optimal allows a source term,
$ "WFR"^2(mu, nu) = min_((rho, v, alpha)) integral_0^1 integral (norm(v_t)^2 + kappa thin alpha_t^2) dif rho_t dif t quad "s.t." quad partial_t rho + nabla dot (rho v) = rho alpha, $
with a length scale $kappa > 0$ (normalisations differ between references). Its static counterpart is the optimal entropy-transport problem
$ min_(gamma in cal(M)_+(cal(X)^2)) integral c dif gamma + rho_"m" H(gamma_1 | mu) + rho_"m" H(gamma_2 | nu), $
in which the hard marginal constraints are replaced by KL penalties of strength $rho_"m"$ @liero2018optimal. With an additional entropic term this is solved by *scaling algorithms*, i.e. Sinkhorn iterations with damped updates @chizat2018scaling. Waddington-OT uses precisely this, with the marginals reweighted by a prior on cellular growth rates @schiebinger2019optimal.

The TI paper @chizat2022trajectory incorporates growth heuristically. Given a prior $g$, the Schrödinger bridge between consecutive clouds uses reweighted marginals, the source $mu_i$ multiplied by $e^(g Delta t\/2)$ and the target $mu_(i+1)$ by $e^(-g Delta t\/2)$, plus KL-relaxed marginal constraints. (The sign convention given here is the one consistent with growth and used in the authors' code.)

#figure(
  image("figures/branching.png", width: 100%),
  caption: [Growth confounds snapshots (companion notebook). The upper branch divides faster, so later snapshots are dominated by it, although 70% of lineages go down. MFL with the growth-reweighted bridges recovers 0.69. Without growth, mass conservation forces a constant branch fraction and the paths over-populate the upper branch.],
) <fig:branching>

== The gradient-flow view of growth

In Otto's picture, the Fisher–Rao (Hellinger) geometry has gradient flows of *birth–death* type,
$ partial_s mu = -mu ((delta F)/(delta mu) - integral (delta F)/(delta mu) dif mu), $
which teleport mass instead of transporting it. The WFR gradient flow sums the transport (Wasserstein) and reaction (Fisher–Rao) parts. For MFL-type objectives, adding birth–death moves to Langevin dynamics is known to accelerate convergence, in particular across energy barriers where pure diffusion is metastable @rotskoff2019global @lu2019accelerating.

#connection(title: "Growth suggests a WFR mean-field Langevin dynamics")[
  Two independent reasons point to the same geometry. _Modelling_: with growth, the marginals of @eq:growth live in $cal(M)_+(cal(X))$, and the natural reference process is branching Brownian motion, whose entropic transport problems are unbalanced. _Optimisation_: the MFL particles of different branches must exchange mass. Under pure Wasserstein dynamics this requires crossing low-density regions, which is exponentially slow in $1\/epsilon$ (the $e^(-alpha\/epsilon)$ in the rate above), whereas Fisher–Rao moves would reweight particles directly. A principled formulation of the min-entropy estimator for branching processes, together with a WFR / Hellinger–Kantorovich gradient flow to compute it and a corresponding "mixed-geometry" convergence analysis, is a natural open problem at the intersection of TI and the theory of unbalanced gradient flows.
]

= The wider landscape <sec:landscape>

The min-entropy estimator is _nonparametric_ and _convex_. Much of the recent literature trades these properties for scalability or structure.

- *Parametric landscapes.* Fit $Psi_theta$ (or a drift network) so that simulated marginals match the snapshots: generative RNNs @hashimoto2016learning, PRESCIENT @yeo2021generative, JKOnet, which learns the energy of a JKO scheme @bunne2022proximal, and TrajectoryNet, a continuous normalising flow regularised by the Benamou–Brenier action @tong2020trajectorynet. These approaches gain a structural prior and the ability to extrapolate, at the cost of non-convex training and no global guarantees. The paper's conclusion suggests learning $Psi$ _jointly_ with $R$ as a hybrid.
- *Neural Schrödinger bridges.* Iterative proportional fitting in path space (diffusion Schrödinger bridges @debortoli2021diffusion) and bridge matching @shi2023diffusion scale the dynamic problem @eq:dsb to high dimension with neural drifts. Flow matching @lipman2023flow with mini-batch OT couplings @tong2024improving is its $epsilon -> 0$, simulation-free relative.
- *Closed forms.* Between Gaussians, both entropic OT and the Schrödinger bridge are explicit @janati2020entropic @bunne2023schrodinger. They are useful as reference processes and as sanity checks.

= Open problems for the OT and gradient-flow community <sec:open>

+ *Statistical rates.* Consistency is known @lavenant2021towards; quantitative rates in terms of $T$, $n_i$, $d$, $lambda$, $sigma$, $tau$ are largely open. The dimension-free sample complexity of entropic OT suggests favourable rates in the number of cells.
+ *Geometry of $F$.* Is the reduced objective (geodesically) convex along entropic or displacement interpolations in $cal(P)(cal(X))^T$? A positive answer would bring EVI-type contraction to TI.
+ *Discrete-time and finite-particle guarantees.* Uniform-in-time propagation of chaos and step-size analysis for the particle scheme with inexact Sinkhorn steps.
+ *Unbalanced min-entropy estimation.* Entropy relative to branching Brownian motion, its reduction to unbalanced Schrödinger bridges, and a WFR mean-field Langevin solver with global guarantees (@sec:uot).
+ *Better priors.* Learning $Psi$ jointly with the path law, non-Brownian reference processes (e.g. Ornstein–Uhlenbeck or learned references), and second-order (spline) priors in Wasserstein space @benamou2019second @chewi2021fast.
+ *Beyond Euclidean state spaces.* Manifold-valued cell states and costs given by learned metrics, where Varadhan's formula links the reference heat kernel to the transport cost.

#v(1em)
#line(length: 100%, stroke: 0.4pt + luma(180))
#text(size: 9pt)[A companion Jupyter notebook implements everything in this handbook from scratch (Sinkhorn, MFL, path sampling, branching) and reproduces the figures; see the repository root.]

#bibliography("refs.bib", title: "References", style: "apa")
