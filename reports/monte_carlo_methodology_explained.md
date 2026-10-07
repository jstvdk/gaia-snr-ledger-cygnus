# The Monte Carlo death ledger, explained from the ground up

*Reference note, written 2026-10-02 so the method can be refreshed later.
All numbers are from the adopted chain **repair_v8** (`scripts/chain.py`), baseline branch = PARSEC,
R_V = 3.1, α = 2.3, coeval (δ = 0), all-explode, unless stated otherwise.
Engine: [scripts/wp7_ledger.py](../scripts/wp7_ledger.py) · constants: [scripts/wp7_ledger_prereg.py](../scripts/wp7_ledger_prereg.py) ·
inputs: [scripts/wp5_fit_imf.py](../scripts/wp5_fit_imf.py), [scripts/wp5_joint_age_fit.py](../scripts/wp5_joint_age_fit.py).*

> Notation. The ledger's columns are called `N_SN`, but every number here is
> $N_{\rm death}$, i.e. $N_{\rm SN}\,|\,\text{all explode}$. It is not an
> observed supernova count.

---

## 0. The idea in one paragraph

The surviving 2–8 M☉ stars fix *how many* stars each subgroup formed. They do not fix *which masses* its massive stars had, and Cyg OB2 is one random draw from the IMF. So the Monte Carlo builds millions of possible versions of each subgroup that are consistent with the data. In each version it:
- decides how many massive stars were born and what masses they had;
- lets the lifetime relation decide which have died;
- dates each death.

Every quantity in the paper is then a count over those versions: the mean number of deaths, its interval, the probability of a recent death, and the death-rate history.

---

## 1. What we know and what we want

**Known, per subgroup:**
- $k$, the IMF normalisation, from the 2–8 M☉ counts;
- $t$, the age;
- $\alpha$, the IMF slope, carried as a branch;
- $\tau(m)$, the lifetime, and the turnoff mass $m_{\rm TO}(t)$, both from the isochrones.

**Wanted:**
1. How many massive stars have died?
2. When did each die?
3. $P(\ge 1\ \text{death})$ and $P(\text{last death} < 100\ \text{kyr})$.
4. The **shape** of the answer, not just its mean. For example, the chance of zero deaths.

## 2. Why the closed-form formula is not enough

The textbook result gives the *expected* number only:

$$
\langle N_{\rm death}\rangle = k\int_{m_{\rm TO}(t)}^{120} m^{-\alpha}\,dm
= k\,\frac{m_{\rm TO}^{\,1-\alpha} - 120^{\,1-\alpha}}{\alpha-1}.
$$

Four things stop it from answering everything in §1:

1. **It is an average, not a distribution.** "4.10 expected deaths in A" means A had 2, or 5, or 7. Questions such as $P(N=0)$ need the whole distribution.
2. **"When was the last death?" is a minimum over a random set of stars.** It depends on whether, by chance, a star was born just above today's turnoff.
3. **The inputs are uncertain, correlated and non-Gaussian.**
   - $k$ and $t$ come from one joint fit: the correlation is −0.13 for A and −0.26 for B.
   - B's age draws pile up at the top of its grid (16–84 %: 3.98–4.26 Myr).
   - Standard error propagation assumes small, Gaussian, independent errors, and none of those holds.
4. **The ingredients combine non-linearly:**
   - a spread of birth times ($\delta$);
   - explodability filters;
   - a turnoff capped at 120 M☉, which switches deaths off completely for C;
   - a sum over three subgroups.

The Monte Carlo is a way to evaluate the resulting multi-dimensional integral numerically by simulating it.

---

## 3. Why deaths are a Poisson process

Given a star's birth mass and the group's age, its death is **not** random: it happens at $\tau(m)$ after birth. The randomness is in how many massive stars were born and at which masses.

- A subgroup formed tens of thousands of stars. Each one had a tiny, independent chance of being born above the turnoff.
- Many independent trials with a small success probability give a Poisson count (the law of rare events). So the number of stars born in any mass interval $[m_1, m_2]$ is
$$
N_{[m_1,m_2]} \sim \mathrm{Poisson}\!\left(k\int_{m_1}^{m_2} m^{-\alpha}\,dm\right).
$$
- Death time is a one-to-one function of mass, so a Poisson scatter of points in mass becomes a Poisson scatter in time, with intensity $R_{\rm SN}(t_{\rm lb})$.
- Counts in disjoint time windows are therefore independent. That gives the closed forms used as sanity checks:

$$
P(\ge 1) = 1 - e^{-N_{\rm death}},\qquad
P(\text{last} < \Delta) = 1-\exp\!\left[-\int_0^{\Delta} R_{\rm SN}\,dt_{\rm lb}\right]\approx 1-e^{-0.8}\approx 0.55,
$$

$$
\text{68 % interval}\approx N\pm\sqrt N = 8.4\pm2.9 \;\Rightarrow\; 5\text{–}11 .
$$

Because $k$ and $t$ are themselves drawn from their posterior, the simulated count is a **mixed Poisson**:

$$
\mathrm{Var}(N) = \mathbb{E}[\mu] + \mathrm{Var}(\mu).
$$

For A this is $4.23 = 4.10 + 0.13$: the uncertainty in $k$ and $t$ adds only about 3 % on top of the IMF-sampling noise.

**When the Poisson assumption could fail:**
- **Fixed mass budget.** If a cloud forms stars until its gas runs out, the counts are slightly less scattered than Poisson.
- **"Optimal sampling".** Weidner & Kroupa argue that a cluster's maximum stellar mass is set by its own mass, which removes most of the randomness. This is debated.
- **Binaries.** Companion masses are correlated, and stars can merge.

---

## 4. The inputs: paired posterior draws of $(k, t)$

The ledger does not start from raw data. For each subgroup × (family, R_V, α) the file `data/processed/wp5_imf_posterior_draws_repair_v8.npz` holds **10,000 paired draws** `k__<key>` and `truth_age_draws__<key>`. They were made in four steps:

1. **Completeness by injection** (itself a Monte Carlo). Synthetic stars of known mass are inserted into the real catalogue with realistic extinction and run through the whole selection pipeline. The recovered fraction gives the completeness $R(\mathrm{obs}\mid m)$ (`wp5_injections*.py`).
2. **The likelihood of the counts.** Members of 2–8 M☉ are binned by mass, each weighted by its membership probability. The expected count in bin $i$ is

$$
\lambda_i = k\int_{\text{bin } i} R(\mathrm{obs}\mid m)\, m^{-\alpha}\,dm ,
$$

   and the observed counts are treated as Poisson with these means.
3. **The posterior for $k$.** A Jeffreys prior $p(k)\propto k^{-1/2}$ makes the conditional posterior of $k$ a Gamma distribution, which is sampled exactly. Uncertainty in the completeness curve enters through a Dirichlet/Beta redraw of the response before each $k$ draw.
4. **Age nodes.** The WP4 colour–magnitude age posterior is discretised into 9 nodes. Each node is reweighted by the marginal likelihood of the counts:

$$
\log \mathcal{L}_j = \sum_i n_i \log r_{ij} - \left(N+\tfrac12\right)\log R_j + \text{const}.
$$

   Each draw first picks a node, then a $k$ consistent with it. **This is why $k$ and $t$ must be drawn as pairs.**

Baseline inputs on repair_v8:

| subgroup | $k$: median [16–84 %] | age: median [16–84 %] (Myr) | $m_{\rm TO}$ | expected deaths $\mu$ |
|---|---|---|---:|---:|
| A | 1692 [1595, 1795] | 4.03 [3.93, 4.10] | 57.3 M☉ | 4.10 |
| B | 1647 [1537, 1761] | 4.26 [3.98, 4.26] | 52.1 M☉ | 4.26 |
| C | 1894 [1767, 2027] | 2.51 [2.47, 2.56] | 120 (capped) | 0.00 |

On the project's PARSEC relation, $\tau(120\,M_\odot) = 2.99$ Myr.

---

## 5. One iteration, step by step (`run_population()`)

One iteration simulates **one possible history of one subgroup**, on one branch.

**Step 1: pick $(k, t)$.** Draw one of the 10,000 pairs at random.

**Step 2: the lowest possible turnoff.** With a birth window of width $\delta$, the oldest stars were born $\delta/2$ before $t$:

$$
m_{\min} = m_{\rm TO}\!\left(t+\tfrac{\delta}{2}\right).
$$

No star below $m_{\min}$ can be dead, so none is simulated. This is exact, not an approximation: a Poisson process restricted to $m > m_{\min}$ is itself a Poisson process, independent of the stars below.

**Step 3: how many were born above $m_{\min}$?**

$$
\mu = k\int_{m_{\min}}^{120} m^{-\alpha}\,dm, \qquad N\sim\mathrm{Poisson}(\mu).
$$

This is where the randomness of star formation enters (`imf_integral`, `rng.poisson`).

**Step 4: their masses, by inverse CDF** (`sample_imf`). With $p = 1-\alpha$, the cumulative fraction of stars below $m$ is

$$
F(m) = \frac{m^{p} - m_{\min}^{p}}{120^{p} - m_{\min}^{p}} .
$$

Set $F(m) = u$, with $u\sim U(0,1)$, and solve:

$$
m = \left[m_{\min}^{p} + u\,\bigl(120^{p} - m_{\min}^{p}\bigr)\right]^{1/p}.
$$

Small $u$ gives masses near $m_{\min}$, where most stars are; $u$ near 1 gives the rare giants near 120 M☉.

**Step 5: birth times.**

$$
t_{\rm birth} = t + \delta\,\bigl(u' - \tfrac12\bigr),\qquad u'\sim U(0,1).
$$

With $\delta = 0$, every star is born at $t$.

**Step 6: dead or alive.** A star is dead if $m > m_{\rm TO}(t_{\rm birth})$. The test uses the same turnoff relation as WP4 and WP6, so the ledger cannot disagree with them because of a different lifetime table (`TurnoffRelation`).

**Step 7: the death date.**

$$
t^{\rm death}_{\rm lb} = t_{\rm birth} - \tau(m),
$$

where $\tau$ is obtained by inverting the same turnoff relation.

**Step 8: explodability as a filter.** "All explode" keeps every death. The threshold scan keeps a death only if $8 \le m < M_{\rm BH}$. These are filters on **the same simulated stars** ("common random numbers"), so any difference between explodability options is physical, not sampling noise.

### Four real iterations (subgroup A, baseline, repair_v8 draws)

```
iter0: k=1777 age=3.882 m_TO=61.4 μ=3.76 → N=4
    m=113.0  died 0.850 Myr ago
    m= 89.9  died 0.629
    m= 85.3  died 0.553
    m= 67.1  died 0.166        ← last death 166 kyr ago: NOT within 100 kyr
iter1: k=1825 age=3.981 m_TO=58.4 μ=4.32 → N=4
    m=105.9  died 0.898
    m= 92.0  died 0.760
    m= 89.2  died 0.717
    m= 69.2  died 0.321        ← no
iter2: k=1672 age=3.882 m_TO=61.4 μ=3.54 → N=1
    m= 72.7  died 0.311        ← no
iter3: k=1741 age=3.882 m_TO=61.4 μ=3.69 → N=3
    m= 85.3  died 0.554
    m= 65.2  died 0.114
    m= 64.0  died 0.080        ← yes: a star just above the turnoff died 80 kyr ago
```

What these show:
- **The same inputs give different histories.** Iterations 0, 2 and 3 share an age, yet produced 4, 1 and 3 deaths. That is IMF-sampling luck, a real physical randomness.
- **Heavier stars die earlier.** This is "mass is a clock", applied star by star.
- **A recent death depends on luck near the turnoff.** $P(\text{last} < 100\ \text{kyr})$ is the fraction of iterations like iteration 3.

---

## 6. Repeat and count

The ledger runs **2,000,000 iterations** per subgroup and branch.

| quantity (baseline, repair_v8) | value | how it is read off |
|---|---|---|
| mean $N_{\rm death}$, association | **8.36** (A 4.10 · B 4.26 · C 0) | average over iterations |
| median, 68 % interval | 8, [5, 11] | 50th, 16th and 84th percentiles |
| $P(\ge 1)$ | 0.9996 | fraction of iterations with $N \ge 1$ |
| $P(\text{last} < 100\ \text{kyr})$ | 0.547 | fraction whose most recent death is under 0.1 Myr ago |
| median time since the last death | 0.088 Myr | median of that minimum |
| $R_{\rm SN}(t_{\rm lb})$ | about 8 per Myr, flat | histogram of all death dates (below) |

The rate history is a normalised histogram of all death dates, in bins of $\Delta t = 0.05$ Myr:

$$
\hat R_{\rm SN}(\text{bin}) = \frac{\#\{\text{deaths in bin}\}}{N_{\rm iter}\,\Delta t},
\qquad \int \hat R_{\rm SN}\,dt_{\rm lb} = \langle N_{\rm death}\rangle .
$$

Every probability is a sample fraction:

$$
\hat P(\text{event}) = \frac{1}{N_{\rm iter}}\sum_{i=1}^{N_{\rm iter}} \mathbf{1}[\text{event in iteration } i].
$$

Distribution of $N$ for subgroup A (200,000 iterations):

| N | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P | 0.018 | 0.070 | 0.142 | 0.190 | **0.191** | 0.157 | 0.108 | 0.064 | 0.034 | 0.016 | 0.007 |

The distribution is right-skewed, with a 1.8 % chance that A has had no deaths. A single "4.10 ± 2" cannot express that.

**The whole association.** The subgroups are simulated separately, then combined iteration by iteration:

$$
N^{(i)}_{\rm assoc} = N^{(i)}_A + N^{(i)}_B + N^{(i)}_C,\qquad
t^{(i)}_{\rm last} = \min\bigl(t^{(i)}_{{\rm last},A},\, t^{(i)}_{{\rm last},B},\, t^{(i)}_{{\rm last},C}\bigr).
$$

Association-level probabilities must be computed on these combined histories, never by combining the separate subgroup numbers.

---

## 7. Which uncertainty lives where

### Inside the Monte Carlo (randomness we can assign probabilities to)
1. **IMF sampling:** the Poisson count and the mass draws. This dominates (§3).
2. **Uncertainty in $k$ and $t$:** the WP5 posterior pairs.
3. **Birth-time spread:** the draws for $\delta > 0$.

### Outside the Monte Carlo (model choices, carried as branches)
- isochrone family: PARSEC, MIST;
- extinction law: $R_V \in \{3.0, 3.1, 3.5\}$;
- IMF slope: $\alpha\in\{2.0, 2.3, 2.6\}$;
- formation duration: $\delta\in\{0,1,2\}$ Myr.

That is 54 branches, and **each branch is its own Monte Carlo.** There is no defensible probability for "PARSEC is right", and drawing branches at random would invent one and hide the disagreement inside a smooth distribution. Since 2026-10-02 the headline slope is $\alpha = 2.3$, with 2.0 and 2.6 reported as sensitivity branches (`provenance/decisions_2026_10_02.json`).

### Not in the Monte Carlo at all
- **Binary evolution:** mass transfer and mergers, which change lifetimes.
- **Realistic explodability:** only all-explode and step thresholds are tested.
- **Distance uncertainty:** the distance is fixed at 1.6245 kpc.
- **Runaways:** they do not enter $N$, because they are living stars below the turnoff and change neither $k$ nor $m_{\rm TO}$.
- **Membership uncertainty:** it enters only through the membership-weighted counts that set $k$.

### The two sensitivity scans reuse the same engine
- **Age scan.** The age is *fixed* at 2.0, 2.25, …, 6.0 Myr for all three subgroups while $k$ is still drawn (`wp7_age_sensitivity_repair_v8.csv`). The result is 0 up to 2.75 Myr, 0.20 at 3.0, 6.47 at 3.5, 12.63 at 4.0 and 37.9 at 6.0. Interpolating, a common age near 3.66 Myr reproduces the resolved 8.36. Note that **each group keeps its own $k$** in this scan: it is not a pooled fit.
- **Black-hole threshold scan.** The same dead stars, filtered with $8 \le m < M_{\rm BH}$. The result is 0 for every threshold at or below 40 M☉ on the baseline, because every star that has died was born above about 52 M☉.

---

## 8. How many iterations: Monte Carlo error is not physical error

The numerical error on a mean shrinks as $1/\sqrt{N_{\rm iter}}$:

$$
\mathrm{SE}\bigl(\langle N\rangle\bigr) \approx \sqrt{\frac{\mathrm{Var}(N)}{N_{\rm iter}}} \approx \sqrt{\frac{8.6}{2\times10^6}} \approx 0.002 .
$$

So 8.36 is numerically good to about ±0.002, while the physical 68 % range is 5–11 and the branch range is wider still.

**The pre-registered convergence test (L6) failed as written.** It required every result to change by less than a fixed *percentage* when the iteration count is doubled. Record from the repair_v7 scan, `provenance/wp7_convergence_scan.json`:

| iterations | worst relative drift, all cells | cells with ≥ 0.5 deaths | association $N$ |
|---:|---:|---:|---:|
| 40,000 | 11.6 % | 1.06 % | 8.453 |
| 400,000 | 3.40 % | 0.28 % | 8.437 |
| 2,000,000 | 2.45 % | 0.155 % | 8.435 |

It plateaus because some cells expect around 0.003 deaths, where a tiny absolute change is a large percentage. The quantities that matter are stable to thousandths of a death. Per project rules, the gate stays recorded as **failed**. The repair_v8 doubling check is in `tables/wp7_convergence_repair_v8.csv`.

---

## 9. Code map

| piece | where |
|---|---|
| $\int m^{-\alpha}dm$, vectorised | `wp7_ledger.imf_integral` |
| inverse-CDF IMF draw | `wp7_ledger.sample_imf` |
| $m_{\rm TO}(t)$ and its inverse $\tau(m)$ | `wp7_ledger.TurnoffRelation` (wraps `wp6_mass_extension_decision.turnoff_mass`) |
| one batch of iterations | `wp7_ledger.run_population` |
| summary statistics | `wp7_ledger.summarize` |
| branch loop, association sums, age and BH scans | `wp7_ledger.main` |
| fixed constants: iterations, δ grid, scan grids, thresholds | `wp7_ledger_prereg.py` |
| the $(k, t)$ posterior | `wp5_fit_imf.py` (Gamma + Dirichlet) and `wp5_joint_age_fit.py` (age nodes) |
| which chain is read | `scripts/chain.py` (`ADOPTED = "repair_v8"`) |

The whole chain reruns with `bash scripts/run_repair_v8_chain.sh`.

---

## 10. What the Monte Carlo inherits, and so cannot fix

The Monte Carlo propagates uncertainty faithfully, but only the uncertainty it is given. Its numbers are **conditional on the ages it receives.**

The main live example is subgroup C (**issue #20**). Its photometric age of 2.51 Myr is what makes C contribute exactly 0 deaths. Phase A (`reports/issue20_phase_a_report.md`) found that C's own spectroscopic stars prefer about 4 Myr instead, and the decision of 2026-10-02 is a Phase A′ with three independent age tests. If C's adopted age changes, every number above changes with it, in a separately pre-registered chain, repair_v9.

## 11. Symbols

| symbol | meaning |
|---|---|
| $\xi(m)=dN/dm=k\,m^{-\alpha}$ | IMF: stars born per unit initial mass |
| $k$ | IMF normalisation (size of the birth event) |
| $\alpha$ | IMF slope; 2.3 is the baseline and headline |
| $m_{\rm TO}(t)$ | turnoff mass: the heaviest star still alive at age $t$, capped at 120 M☉ |
| $\tau(m)$ | lifetime of a star of initial mass $m$ |
| $t$, $t_{\rm birth}$, $t_{\rm lb}$ | subgroup age, a star's birth time, look-back time |
| $\delta$ | width of the star-formation window (0, 1, 2 Myr) |
| $\mu$ | expected number of stars born above $m_{\min}$ in one iteration |
| $R(\mathrm{obs}\mid m)$ | completeness: the probability that a star of mass $m$ is catalogued |
| $R_{\rm SN}(t_{\rm lb})$ | death rate against look-back time |
| $M_{\rm BH}$ | direct-collapse threshold of the explodability scan |
| $N_{\rm iter}$ | number of Monte Carlo iterations (2,000,000) |
