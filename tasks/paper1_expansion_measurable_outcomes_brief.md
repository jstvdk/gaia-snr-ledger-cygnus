# Paper 1 expansion — the measurable outcomes, and what the census can and cannot add

**Date:** 2026-09-16
**Status:** brief. Answers three questions from the 16 Sep discussion, then fixes the
measurable outcomes any expansion must deliver. Nothing here is a product: the two
pilots in §2 were run from a scratch directory, wrote nothing into the repository,
and are quoted only as order-of-magnitude disclosed prior knowledge, exactly as WP13
§2 discloses its own.
**Governing documents:** the 12 Aug novelty audit; `tasks/wp13_pooled_vs_resolved_ablation_brief.md`.

---

## 1. Answers

### 1.1 What a time-resolved wind-versus-supernova budget concludes (question 1)

**The present-day wind power of Cyg OB2 is not a Gaia measurement.** The pilot
(§2.1) sums the Vink et al. (2001) recipe over the 241 labelled members above
8 M☉ that sit on their subgroup's isochrone and gets **3.6 × 10³⁷ erg s⁻¹**
(Björklund et al. 2021 would give 2–3× less). The literature values —
1.5–2.9 × 10³⁸ (Menchiari+24, App. A), 2–3 × 10³⁸ (Vieu+24b, quoted by Härer+25)
— are carried by the three Wolf-Rayet stars, WR 144 alone being ≥ 9 × 10³⁷
(Sander+19). All three are already in `wp1_spectroscopic_anchors` as orphan
anchors, with literature parameters. So a "Gaia census wind luminosity" would
reproduce the literature number by construction and add nothing to it. **Proposal 1
as previously stated is withdrawn in that form.**

What the census *does* add is the time axis, and there the ledger and the wind
history are the same object. A star that has died spent its last ~0.3–0.5 Myr as a
WR star at 10³⁷–10³⁸ erg s⁻¹, i.e. it deposited **10⁵⁰–10⁵¹ erg of wind** just before
depositing **~10⁵¹ erg of supernova**. The dead stars therefore contribute
N_death × (1–2) × 10⁵¹ erg. The living O stars contribute 3.6 × 10³⁷ × 4 Myr ≈
4.5 × 10⁵¹ erg in total; the living WR stars perhaps 5 × 10⁵⁰ so far. On the
baseline (N_death = 8.43):

| channel | energy, erg | share |
|---|---:|---:|
| supernovae of dead stars | 8.4 × 10⁵¹ | ~0.4 |
| WR-phase winds of the same dead stars | 1–8 × 10⁵¹ | ~0.3 |
| living O stars, integrated over 4 Myr | 4.5 × 10⁵¹ | ~0.2 |
| living WR stars so far | ~0.5 × 10⁵¹ | ~0.03 |

**Conclusion A. Roughly two thirds to three quarters of all the mechanical energy
ever injected into the Cygnus superbubble came from stars that no longer exist**
(baseline; about half on the α = 2.6 branches, more on α = 2.0). Every hydrodynamic
and transport model of the region (Menchiari+24, Vieu+24b, Härer+25, Li+26) drives
the bubble with the *present* census at constant power and no supernovae. The ledger
says the source term they use is the minority channel. This is a statement about
the bubble's energy content and thermal state (density, pressure, cavity size),
which is the input Härer's hadronic model is most sensitive to (n = 0.05 cm⁻³).

**Conclusion B. Cyg OB2 is at the wind-to-supernova transition now, and which side
it is on is set by α.** Present wind power 2–3 × 10³⁸ against a time-averaged
supernova power of rate × 10⁵¹ erg: 2.5 × 10³⁸ at the baseline 8.0 Myr⁻¹, ~6 × 10³⁸
on α = 2.0 branches, ~0.8 × 10³⁸ on α = 2.6. Menchiari's "wind-dominated at 3 Myr,
supernova-dominated at 5 Myr" becomes a per-branch probability instead of an
age choice.

**The cosmic-ray connection.** The community's efficiency is η = W_CR / E_available,
and the denominator is where the census enters. Over the window that matters for
the 10 TeV–PeV component — the CR residence time, ≲ 80 kyr in Härer's fit — the
available energy is **bimodal**: winds alone give 2.5 × 10³⁸ × 80 kyr ≈ 6 × 10⁵⁰ erg;
a supernova adds (1–5) × 10⁵¹. Härer's model needs 3 × 10⁵⁰ erg in protons. Without
a supernova in the window that is a 40–50 % wind efficiency, which nobody defends;
with one it is 6–10 %, which is standard. The ledger gives the probability of the
second case: P(death within 80 kyr) ≈ 0.47 on the baseline, and the existing WP12
cocoon score already folds in location and progenitor class (C4 = 0.854 in situ,
0.58–0.72 across scenarios). So the deliverable is not "η_CR = x %" — L_γ is theirs
and L_w is the WR stars' — but **the denominator as a distribution over the
residence time, with the probability that it is supernova-dominated**, per branch.
That is a table any transport paper can consume, and none of them has it.

### 1.2 The recency point, in plain words (question 2)

The death count is the *area* under the death-rate curve. The probability of an
event in the last 100 kyr is the *height* of that curve today. A single-age model
can be tuned to the same area (3.66 Myr reproduces 8.43), but it then has a
different height, for two reasons that add:

1. Its first death is later (turnoff 120 M☉ at ~2.75 Myr, so the history spans
   ~0.9 Myr instead of the resolved 1.30 Myr). The same 8.4 deaths squeezed into a
   shorter interval means a higher rate today: ~12 Myr⁻¹ against 8, hence
   P(< 100 kyr) 0.70 against 0.55.
2. In the pooled model the richest subgroup, C (k_C = 1,894, the largest
   normalization), contributes to today's rate. In the resolved model it
   contributes nothing, because its turnoff is above 120 M☉. The pooled model
   attributes current deaths to a population that has had none.

No single age between 3.25 and 6 Myr gets the height right (the scan sits at
0.70–0.77 throughout). **So a fixed-age estimate can be correct on the count and
still overstate the recent-event probability by ~0.15 absolute, ~30 % relative,
and the recent-event probability is the quantity the gamma-ray argument uses.**
That is WP13 criterion T3, already frozen at |ΔP| ≥ 0.10. Nothing new has to be
built; the point is what T3 is *for*.

### 1.3 What the massive stars alone give (question 4)

The census has, by weighted count, 290 members above 8 M☉ (hard count), 76 above
20, 21 above 30, 3 above 40; 150 carry spectroscopic anchor masses. Four things
come out of them without the 2–8 M☉ machinery:

**(i) The slope, per subgroup — the parameter the whole verdict hinges on.**
Wright+15 measured Γ = 1.39 ± 0.19 (α = 2.39) above 20 M☉ for the whole
association. The closure test is the same measurement per subgroup with Gaia
membership, and it already exists: closing α = 2.34 (A), 2.25 (B), 2.06 (C), with
per-cell 68 % intervals of about ± 0.05. At α = 2.0 the massive stars of A and B are
over-predicted by 30–43 % (closure 0.57–0.74); at α = 2.3 C is under-predicted by
36–46 % (closure 1.36–1.46) while A and B close. **The census itself splits the
headline set: α = 2.0 is excluded for A and B, α = 2.3 is excluded for C.** The
project ruled that closure is a validation and may not choose the headline α; that
rule stands. But a clearly labelled row "N_death | closure-constrained α" is a
different thing from retuning the headline, and it is the number the community
will use. On the baseline family it lands at A ≈ 4, B ≈ 4, C = 0 either way, because
C's zero is a turnoff fact not a slope fact — which is itself worth saying.

**(ii) A direct age bound on C from its most massive living member.** C contains
a 61 M☉ star (photometric, 52–68) and three spectroscopic anchors at 37–40 M☉.
A living 61 M☉ star means C's turnoff is above 61 M☉, i.e. C is younger than the
61 M☉ lifetime (~3.5 Myr, PARSEC), independently of the CMD fit that gives 2.52.
Cheap, model-light, and it is the "C sits below the first-death boundary" claim
tested on one star.

**(iii) A massive-star-only death count, the Fuchs+06 / Wright+15 route.**
N_dead = N_obs(> 20 M☉) × ∫_TO^120 m^−α / ∫_20^TO m^−α. With the weighted labelled
counts and baseline turnoffs: A ≈ 2.0, B ≈ 3.2, C = 0, before the 28 unlabelled
and orphan stars above 20 M☉ are attributed (pro-rata they add 2–3). Consistent
with 8.43 inside the branch spread. This is the check that the count does not
depend on whether the normalization comes from intermediate-mass or massive
stars — a robustness statement reviewers will want, not a result.

**(iv) The top-heavy-C question.** C's closing slope 2.06 is 1.7σ from Wright's
association-wide 2.39 ± 0.19 and 0.3 below A's. Three explanations, each with a
test on existing products: labelling (`wp2_label_stability_per_star.csv` gives
per-star seed fractions — re-run closure with the > 20 M☉ stars of C reassigned
by their stability), mass segregation (radial concentration of the > 20 M☉ stars
of C against its 2–8 M☉ stars), age-conditional completeness (the WP5 response at
C's age). If it survives all three it is a sub-association IMF result; if not, the
systematic is named and the closure excess is explained.

## 2. Pilots (scratch, disclosed, not products)

### 2.1 Wind sum
Vink+01 on `wp4_mass_posteriors_repair_v5` baseline masses, PARSEC R_V = 3.1
isochrone at the subgroup age, membership-weighted, no orphans, no WR, no
unlabelled members: A 7.1 × 10³⁶, B 7.9 × 10³⁶, C 2.1 × 10³⁷, total 3.6 × 10³⁷
erg s⁻¹; 60 % of it from stars above 30 M☉; **58 % of the O-star wind power sits in
C**, the subgroup with no deaths.

### 2.2 Geometry
Membership-weighted centroids (l, b): A (80.24, 0.94), B (79.63, 0.91),
C (80.07, 0.66); rms radii 0.27–0.40°. A–C separation 0.33°, B–C 0.50°, i.e.
9–14 pc at 1.62 kpc. The four most massive living stars (37–61 M☉) are all in C
at b = 0.13–0.35. This is too small against a 50 pc cocoon and a degrees-wide
LHAASO bubble to be a decisive morphological prediction; it is recorded so that it
is not proposed as one.

## 3. The measurable outcomes

Each row is a quantity with a definition, the number it replaces, the consumer,
and the sentence that would falsify it. These are what any expansion must deliver;
anything that does not feed one of them is out of scope.

| id | quantity | replaces | consumer | falsifiable form |
|---|---|---|---|---|
| **O1** | α per subgroup from the massive stars (closing slope with interval), and N_death conditional on it | Wright+15's single Γ = 1.39 ± 0.19; Menchiari's assumed Kroupa | every Cyg OB2 population model; IMF variation literature | "α_A = α_B = α_C" rejected or not at a stated level after the three §1.3(iv) systematics; "α = 2.0 excluded in A and B at > 3σ by closure" |
| **O2** | E_inj(t): cumulative mechanical energy by channel (SN, dead-star WR winds, living winds) and lookback time, per branch; scalars E_tot, f_SN, f_dead | constant L_w × age, no SNe (all four transport/hydro papers) | hydro and CR-transport modellers; superbubble size and pressure arguments | "f_dead > 0.5 on every headline branch"; bubble radius from E_tot vs the HI/X-ray cavity (weak, factor 1.25 for 3× energy) |
| **O3** | available energy over the CR residence time τ ∈ {30, 80, 100, 300 kyr, 1 Myr} as a distribution, with P(SN-dominated \| τ) | Härer's assumed 50 kyr, 3–5 × 10⁵¹ erg event; Menchiari's wind-only denominator | transport papers that quote η_CR | "P(SN within 80 kyr) ∈ [0.3, 0.9] across headline branches, 0.47 baseline"; the wind-only alternative requires η ≥ 0.4 |
| **O4** | ΔP(last < 100 kyr) between M0 and M1 at matched count | fixed-age counts used to argue recency | anyone applying Menchiari-style counts to other associations | WP13 T3, frozen: \|ΔP\| ≥ 0.10 with paired 95 % interval excluding zero |
| **O5** | complete ²⁶Al forecast: SN + dead-star WR winds + living WR stars, decay-weighted; ⁶⁰Fe unchanged | WP11's SN-only ²⁶Al (10–42 % of the complex flux) | INTEGRAL/SPI now (3.9 ± 1.1 × 10⁻⁵ ph cm⁻² s⁻¹ complex-wide), COSI | "the census-predicted ²⁶Al flux of Cyg OB2 is within a factor 2 of the SPI complex flux"; ⁶⁰Fe split at α stays 18/0 |
| **O6** | age bound on C from its most massive living star | the CMD-fit age alone | internal check; reviewers | "C < 3.5 Myr independently of the isochrone fit" |

O1 is the one that changes the verdict. O2–O3 are the ones the high-energy
community would cite. O4 is already briefed. O5 is the one with an existing
measurement to test against and an instrument on the way. O6 is a paragraph.

## 4. Expansion

| package | delivers | new machinery | effort |
|---|---|---|---|
| **WP13** (briefed) | O4 | M0 pooled fit, ledger-level injections | ~3 weeks incl. manuscript |
| **WP14** massive-star mass function | O1, O6, §1.3(iii) | closure re-run under relabelling; radial concentration test; age-conditional response read from WP5 products; PDMF fit > 20 M☉ per subgroup with anchors; N_death \| closure-α row | 1–2 weeks |
| **WP15** energy injection history and isotopes | O2, O3, O5 | per-star wind recipe (Vink and Björklund arms) on isochrone tracks through time; WR-phase energy per dead progenitor mass from the same yield tables WP11 froze (LC06 wind columns); WR anchors from literature; extension of `wp11_isotope_forecast` to the wind channel; the τ-window table from `wp7_rsn_curves` | 2–3 weeks |

**Manuscript shape if all three pass:** "A Gaia DR3 census of Cygnus OB2: mass
function, stellar-death history and energy injection". Sections: census and
closure (O1 leads); death history with the pooled comparison (O4); energy budget
and the CR denominator (O2, O3); isotope forecast (O5); pulsar and cocoon score
stay as applications. The count is a table row. Menchiari is the fixed-age,
wind-only comparison in two places.

**Order.** WP14 first: it is cheapest, and if the C slope dissolves into a
labelling systematic the closure story and O1 change before O2 is built on them.
Then WP13, then WP15.

## 5. Rules carried over

As WP12 and WP13: every input through a hash-verified `frozen()` resolver;
thresholds for O1's equality test, O3's window set and O5's factor-2 band written
into a prereg script and committed before any number is read; the §2 pilots are
part of the disclosed prior knowledge; failed criteria stay failed; the headline
α remains unchosen by closure — the closure-conditional row is labelled as such.
