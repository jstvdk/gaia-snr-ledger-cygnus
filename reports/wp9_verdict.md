# WP9 — The verdict

*Pre-registered in [wp9_verdict_prereg.json](../provenance/wp9_verdict_prereg.json)
before the number was computed. Executed on the `repair_v7` chain, 2,000,000
iterations. Machine-readable record:
[wp9_verdict_execution.json](../provenance/wp9_verdict_execution.json).*

**All three predictions passed. The framing rule was applied mechanically and
returns: regular article.**

---

## 1. The verdict sentence

> **Cyg OB2 plausibly supplied a supernova of the age, progenitor type and
> location that Härer et al. 2025 require, with probability 0.32–0.74 depending
> on the assumed IMF slope (median 0.53). The verdict is not resolved by the
> data: it hinges almost entirely on the high-mass IMF slope, being supported
> (P ≥ 0.59) on every α = 2.0 branch and not supported (P ≤ 0.47) on every
> α = 2.3 branch. Under Härer's own permissive age allowance of "a few hundred
> kyr" the scenario is supported on every branch (P = 0.73–0.85).**

The energy condition is **not** included in that probability — see §3.

## 2. How it decomposes

P_verdict = P(C1 age) × P(C3 progenitor type) × P(C4 in-situ), over 36 headline
branches (α ∈ {2.0, 2.3}, all-explode).

| term | requirement | value |
|---|---|---|
| **C1** age | ≥1 SN in 10–100 kyr | **0.379 – 0.861** (median 0.625) |
| **C3** progenitor type | progenitor > 30 M☉ → stripped, type Ib/c | **1.000** |
| **C4** location | in-situ, from the runaway bound | **0.854** |
| **P_verdict** | | **0.323 – 0.736** (median 0.533) |

**C3 = 1.000 exactly.** Every supernova in the ledger came from a star above
52 M☉ — far above the ~30 M☉ envelope-stripping threshold — so *every* event
Cyg OB2 has produced should have been a stripped-envelope type Ib/c. This is
precisely the progenitor channel Härer's model requires, and it is not a
coincidence: the association is young enough that only its most massive stars
have died. Pre-registered as V1, passed.

## 3. What is deliberately *not* in the number

Härer's fourth condition is an explosion energy of **3–5 × 10⁵¹ erg**, well above
the canonical 10⁵¹ erg. **This project measures progenitor masses, not explosion
energies, and refuses to invent a probability for them.**

So P_verdict is properly read as:

> the probability that a supernova of the right **age**, **type** and
> **location** was available, **given** that such an event can reach the
> required energy.

What the ledger *can* say is that the energy requirement and the progenitor
requirement point the same way: a 3–5 × 10⁵¹ erg event needs an energetic
stripped-envelope progenitor, and every supernova in the ledger had one.

## 4. The verdict hinges on one thing

| axis | spread in P_verdict |
|---|---:|
| **IMF slope α** | **0.264** |
| R_V | 0.065 |
| isochrone family | 0.045 |
| formation duration | 0.043 |

And the split is exact:

| | P_verdict | branches ≥ 0.5 |
|---|---|---|
| α = 2.0 | 0.593 – 0.735 | **18 of 18** |
| α = 2.3 | 0.323 – 0.474 | **0 of 18** |

**Every α = 2.0 branch supports the scenario and every α = 2.3 branch does not.**
The boundary at 0.5 falls precisely between the two slopes. That is why the rule
returns INCONCLUSIVE, and it identifies exactly what future work must settle.

Note this is *not* simply inherited from N_SN, where α carries 91.9% of the
variance. The verdict is a probability that saturates — C1 cannot exceed 1 — so
the α = 2.0 branches, which roughly triple the supernova rate, gain far less in
verdict than they do in N_SN. The pre-registration flagged this in advance as
something to look for.

## 5. The framing decision, applied mechanically

The rule was fixed in `paper1_execution_plan.md` at planning time and
operationalized in the pre-registration before the number was known:

| | |
|---|---|
| SUPPORTED | P ≥ 0.5 on **every** headline branch |
| DISFAVOURED | P ≤ 0.1 on **every** headline branch |
| INCONCLUSIVE | otherwise |

18 of 36 branches sit above 0.5 and none below 0.1. **Outcome: INCONCLUSIVE →
regular-article framing.** Lead with the ledger, the constraints and the DR4
forecast; the verdict is one result among several rather than the headline.

**This is not a failure.** A Letter would have claimed a resolution the data do
not support.

## 6. The permissive window changes the answer — and that is worth stating

Härer's conclusion allows an event "within the last few hundred kyr". Under the
pre-registered permissive window (10–500 kyr):

| | P_verdict, permissive |
|---|---|
| α = 2.0 | 0.843 – 0.854 |
| α = 2.3 | 0.727 – 0.841 |
| **all headline branches** | **0.727 – 0.854** |

**Every branch clears 0.5.** So the scenario is *supported* against Härer's own
permissive reading and *unresolved* against their preferred ~50 kyr.

Both windows were pre-registered before either was computed. The strict window
governs the framing decision because it was declared primary; the permissive
result is reported because Härer state it themselves and it is the more
favourable reading of their own claim.

## 7. Excluded branches, reported not deleted

| branch | verdict | why excluded |
|---|---|---|
| **α = 2.6** | P = 0.142 – 0.250 | WP5's calibration-window χ² alone rejects it — best in 1 of 18 cells, worst median. This does **not** spend the census closure |
| **islands explodability** | identically 0 | WP8: PSR J2032+4127 is a neutron star inside the association; the branch predicts zero supernovae |

Both are in the sensitivity table. Neither was deleted.

## 8. Devil's advocate

*The gate requires the strongest referee objection, written as a hostile expert
would raise it, with concessions where they are due.*

**Objection 1 — "You assume a single power-law IMF from 2 to 120 M☉. Your slope
is constrained at 2–8 M☉ but N_SN depends on the slope above 52 M☉. If the IMF
steepens at high mass your budget collapses."**

*The strongest objection, and it has a real answer.* The census closure test is
exactly a test of that extrapolation: `k` is fitted from 2–8 M☉ counts, and the
>8 M☉ census never enters the likelihood. At α = 2.3 the prediction matches the
observed massive-star count to **6.7%**, with 18 of 18 cells' closing slopes
inside the carried grid. A single power law from 2 to 120 M☉ is therefore not
assumed — it is tested, and it holds to 6.7%. **This is the payoff of not
spending the census to select α**; had we done so, this answer would have been
circular.

**Objection 2 — "Your in-situ fraction is an upper bound, so your verdict is
optimistic."**

*Conceded, and we raise it ourselves.* The runaway count is a lower bound twice
over: the traceback is 2D, so line-of-sight ejections are missed, and the ±8° box
caps recovery at 44 km/s over 5 Myr. A higher true runaway fraction lowers C4 and
therefore lowers P_verdict. The effect is bounded — C4 would have to fall from
0.854 to below 0.60 to move any α = 2.0 branch under 0.5 — but the direction is
against us and is stated.

**Objection 3 — "The pulsar does not prove massive stars explode. A
binary-stripped lower-mass progenitor gives the same neutron star."**

*Fully conceded.* This is recorded as tension T2 in WP8 and cannot be resolved
with available data. What the pulsar does establish is narrower and still
sufficient here: **at least one successful explosion occurred**, which excludes
the branch where the budget is identically zero. Every surviving reading gives a
non-zero budget.

**Objection 4 — "CygOB2-C contributes 0 or 7 supernovae depending on which
isochrone family you pick. That is not a measurement."**

*Conceded.* C's fitted age straddles the age at which the turnoff crosses the
120 M☉ IMF ceiling — 2.52 Myr under PARSEC, 3.16 Myr under MIST. It is reported
as a branch range and never marginalized. Its effect on the *verdict* is small
(family spread 0.045) because the verdict saturates, but its effect on N_SN is
large and is disclosed.

**Objection 5 — "P ≈ 0.53 is a coin flip. You have measured nothing."**

*Rejected, and here is the test.* A coin flip returns 0.5 to every question. This
analysis does not: it returns 0.055 for a 7 kyr window, 0.53 for 10–100 kyr, and
0.84 for the permissive window. Against an ignorance baseline — a last supernova
uniformly distributed over the 1.3 Myr active period — P(within 100 kyr) would
be 0.077. We report **0.55, seven times higher**, because the supernova rate is
measured rather than assumed.

**Objection 6 — "α = 2.6 was excluded after you saw it fail."**

*Rejected on the record.* All three slopes were pre-registered as carried
branches before any injection ran, the exclusion criterion is the WP5
calibration-window χ² fixed as a gate statistic long before WP9, and α = 2.6's
verdict is reported here anyway (0.142–0.250). Nothing was deleted.

## 9. What Gaia DR4 would settle

| improvement | what it fixes | expected effect |
|---|---|---|
| **radial velocities** for bright OB stars | the traceback becomes 3D, so runaways stop being a lower bound | tightens C4, currently an upper bound (Objection 2) |
| **better astrometry + photometry** | subgroup ages, especially CygOB2-C's 2.52-vs-3.16 Myr family disagreement | removes the 0.045 family spread and the 0–7 SNe ambiguity in C |
| **deeper, cleaner membership** | the 49 unlabelled massive members and 27 orphan anchors gain subgroup assignments | makes the association-wide census route well-posed, which WP7 showed it currently is not |

**What DR4 will *not* fix is the dominant term.** The verdict hinges on the
high-mass IMF slope (spread 0.264 against 0.045 for everything astrometric).
Settling α = 2.0 versus 2.3 needs a larger or deeper massive-star census, or a
second association analysed the same way — not better astrometry on this one.

## 10. Predictions, scored

| | statement | outcome |
|---|---|---|
| **V1** | C3 > 0.95 on every headline branch | **PASS** — C3 = 1.000 |
| **V2** | the verdict is driven mainly by C1, the age term | **PASS** |
| **V3** | the outcome is INCONCLUSIVE → regular article | **PASS** |

V3 was a prediction about the project's own conclusion, made before the number
existed, and it held.

## 11. Gate

| criterion | status |
|---|---|
| the verdict sentence exists | **PASS** — §1 |
| it is branch-annotated | **PASS** — the α = 2.0 / 2.3 split is stated in the sentence itself |
| it survived a devil's-advocate pass | **PASS** — §8, six objections, three conceded |
