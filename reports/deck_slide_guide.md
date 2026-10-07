# Talk companion — "Cygnus OB2 Death History" (13 slides)

*For the deck at https://claude.ai/artifact/3EWLoXBDL97KkSKGEs5xHM (private; share it from the page's Share menu).
State as of 6 October 2026: chain repair_v8; issue #19 closed; issues #20 and #21 open.
Option (b) has been decided: redesign the age fit, rebuild the chain as repair_v9, then WP13
(`tasks/stage_repair_v9_wp4_age_redesign_and_wp13.md`).
The talk has three parts: **what I did → why it is not new → what changed and what is next.**
Each slide below gives what is on it, a short script to say, and the questions to expect.
The same script is in each slide's speaker notes. Roughly 13–16 minutes in total.*

---

## Part 1 — What I did

### 1. Cover — "How many massive stars in Cygnus OB2 have already died, and when?"
*Subtitle: what I did, why it is not new yet, and what comes next.*

**Say:** "One question drives the project: how many massive stars in Cygnus OB2 have already died, and when. Three parts: what I did, with the figures; why the main result is not new; and what changed this week and what comes next."

### 2. The problem — "The heaviest stars are gone, but how many, and when?"
**On the slide:** three cards (why it matters / what was missing / what we deliver) and the caveat that a death here is not an observed supernova.

**Say:** "The Cygnus region glows in gamma-rays. Whether that comes from the winds of living stars or from a recent supernova depends on how many massive stars have died, and how recently. Earlier work treated the association as one population with one age; Gaia lets us split it. One caveat throughout: a death in this model is not an observed supernova."

### 3. The idea — "For stars born together, mass is a clock"
**On the slide:** three typeset equations: count the survivors (2–8 M☉ fixes $k$); find who has died (everything above the turnoff mass); date each death ($t - \tau(m)$). The IMF $\xi(m)=k\,m^{-\alpha}$ runs along the bottom.

**Say:** "Stars between 2 and 8 solar masses are all still alive and Gaia sees them, so counting them fixes the size of the birth event. Lifetimes fall steeply with mass, so everything above the turnoff has died. And because the stars were born together, each death happened one lifetime after birth, so mass orders the deaths in time."

**Expect, "why 120 M☉?":** "An assumed upper limit to stellar birth masses, near the most massive stars seen in the Galaxy. It is not fitted. It sets when the first death can happen. Between 100 and 150 it moves the baseline from about 7 to 10 deaths."

### 4. Methodology — "Five steps from Gaia data to a death ledger"
**On the slide:** five step cards (membership, dust and ages, birth count, closure test, death ledger), then the branch product: 2 families × 3 dust laws × 3 IMF slopes × 3 formation durations = 54.

**Say:** "Pick the members and split them into three groups. Correct each star for dust and fit ages. Count the 2 to 8 solar-mass stars to fix each group's size. Check that prediction against massive stars never used in the fit. Then a Monte Carlo, two million simulated histories per branch, dates every death. Model choices run side by side as 54 branches, never averaged, and every check was written down before its result."

**Expect, "why a Monte Carlo?":** "Cyg OB2 is one random draw from the IMF. We want distributions, like the chance of zero deaths or of a recent one, not just averages." The full explanation is in `reports/monte_carlo_methodology_explained.md`.

---

## Part 2 — What came out (figures)

### 5. Membership — "Three groups that share the sky but move differently"
**On the slide:** Figure 2 (sky on the left, proper motions on the right), with the caption: 1,392 members; 1,331 in subgroups (A 476 · B 426 · C 429).

**Say:** "The groups overlap on the sky but move differently, by about 3 to 4 km/s. I found them with a Gaussian mixture on position and motion; three is the only number of groups that comes out the same from every random start. Parallax decides who is a member, but it cannot split the groups: the association sits at one distance, and its depth is smaller than one star's parallax error."

**Expect, "how sure is each label?":** "The split is very reproducible, but individual A/C labels are fuzzy. About half the stars are assigned with more than 90 % confidence. Propagating that is part of the next chain."

### 6. Ages — "Two ways to date a group: they disagree for C"
**On the slide:**
- the colour–magnitude diagrams: A, B, C and the 61 unlabelled spectroscopic stars;
- dots coloured by mass: blue 2–8 M☉ fixes $k$, orange above 8 M☉ is the living census;
- the isochrone drawn as a mass ruler;
- cards with two ages per group, photometric and spectroscopic: A 4.01 and 3.2–3.6; B 4.09 and only 4 stars; C **2.52** and 3.6–4.0 Myr.

**Say:** "Ages come two ways. Photometric: Gaia colour and brightness for about 1,300 stars, after estimating each star's dust. Spectroscopic: about 50 stars whose temperature is measured directly from their spectra. For A they roughly agree. For C they do not: 2.5 against 3.6 to 4. C's bright stars alone agree with the spectra; only its faint stars make it look young."

**Expect, "which is right?":** "The spectra, most likely. Slide 11 explains why."

### 7. Mass function — "The massive-star census closes near a Salpeter-like slope"
**On the slide:** Figure 7, with the caption: the census closes at slope A 2.29 · B 2.25 · C 2.05 (association 2.19); the headline uses α = 2.3.

**Say:** "This checks the extrapolation. From the 2 to 8 solar-mass count I predict the living massive stars, then compare with the census. It closes near the standard slope, so that is the headline; 2.0 and 2.6 are shown as sensitivity cases."

### 8. Death history — "About eight deaths, all in the last 1.3 Myr"
**On the slide:** Figure 8 (the death rate in time, stacked A and B; every branch on the right), with the caption: 8.36 deaths (A 4.10 · B 4.26 · C 0), 68 % range 5–11, P(a death in the last 100 kyr) ≈ 0.55.

**Say:** "On the baseline, about eight deaths, from A and B. The first was about 1.3 million years ago, then about eight per million years since, which is why the chance of one in the last 100 thousand years is about a half. These numbers inherit the photometric ages, so they will move."

### 9. What drives the number — "Assumptions, not noise, set the uncertainty"
**On the slide:** Figure 9 (deaths against an assumed age; explosions against a black-hole threshold), with the caption listing four levers: IMF slope, age, explodability, upper mass limit.

**Say:** "Nothing dies before about 3 Myr, then about twelve deaths per Myr of age, so the ages matter enormously. If stars above 30 solar masses collapse without exploding, there are no explosions at all. The IMF slope moves the count by a factor of ten across 2.0 to 2.6. Model choices dominate; the random scatter is small."

### 10. A forecast for COSI — "The ⁶⁰Fe gamma-ray line could measure the IMF slope"
**On the slide:** the predicted ⁶⁰Fe 1173 keV line flux for every model branch (orange α = 2.0, blue α = 2.3) under three published supernova yield sets, against COSI's 2-year 3σ narrow-line sensitivity of 3 × 10⁻⁶ ph cm⁻² s⁻¹. The caption gives the split on the LC06 NL yields: 18 of 18 α = 2.0 branches above the limit, 0 of 18 α = 2.3 branches.

**Say:** "Supernovae make radioactive iron-60, which emits gamma-ray lines at 1173 and 1332 keV. For each simulated supernova I take the yield for its own progenitor mass, let it decay to today, and convert to a flux. On the first published yield set, every slope-2.0 branch is above COSI's limit and every slope-2.3 branch below. So COSI could measure the IMF slope, the one parameter Gaia cannot settle."

**The caveats, to state up front:**
- A second yield set from the same paper (LC06 Langer), differing only in the Wolf–Rayet mass loss, is about 30 times lower.
- The most recent set (LC18 recommended) gives zero, because every one of these stars collapses directly to a black hole.
- INTEGRAL/SPI's existing upper limit (1.6 × 10⁻⁵, both lines combined) is within about 5 % of the richest α = 2.0 branch (1.51 × 10⁻⁵).
- These fluxes inherit the ages: if C is older, it adds supernovae and raises them.

**Numbers on repair_v8, LC06 NL:**

| | median F(1173 keV) | median ²⁶Al F(1809 keV) |
|---|---:|---:|
| α = 2.0 | 5.6 × 10⁻⁶ | 1.4 × 10⁻⁵ |
| α = 2.3 | 2.1 × 10⁻⁶ | 5.3 × 10⁻⁶ |

**Expect, "what is COSI?":** "NASA's Compton Spectrometer and Imager, a small gamma-ray survey telescope for 0.2–5 MeV, due to launch in 2027." The source is WP11 Part B: `reports/wp11_isotope_forecast.md` and `tables/wp11_isotope_forecast_repair_v8.csv`.

---

## Part 3 — Why it is not new, what changed, what is next

### 11. Context — "Why this is not new"
**On the slide:** five earlier studies (Wright 2015; Berlanas 2019/20; Knödlseder 2002 and Martin 2010; Menchiari 2024; Sukhbold and Ertl), and the conclusion line.

**Say:** "Several deaths, their strong age dependence, and the fact that only the most massive stars have died were all known. Menchiari and collaborators get 7 ± 2.5 at 3 Myr. The one angle that could have been new is that subgroups with different ages change the history, because young C would have no deaths yet. That needs C to be young."

### 12. What changed — "C is probably not young"
**On the slide:** two big numbers: photometric **2.5 Myr** (from faint stars, under a wrong model of the dust errors) against spectra and bright stars **3.2–4 Myr** (as old as A). If C is about 4 Myr, it adds about 3–5 deaths and the subgroup story weakens.

**Say:** "Three independent age tests on C did not settle it formally, but they found why the methods disagree. The photometric fit treats each star's dust uncertainty as two independent errors in brightness and colour. Real dust moves a star along one fixed direction. With that corrected, the faint part of the fit gives impossible ages for every group, and that faint part is exactly where C's young age comes from. The trustworthy measurements say C is about as old as A."

**Expect:**
- "Is the whole analysis wrong?" "No. The method is fine; the age inputs need a redesign, and every number downstream gets recomputed on it."
- "Was anything else fixed?" "Yes, a stale input earlier this week moved the baseline from 8.43 to 8.36."

### 13. Next steps — "Fix the ages, rebuild, then test what is new"
**On the slide:** four steps (redesign the age fit → rebuild the chain as repair_v9 → pooled versus resolved → decide the paper), and three outcomes, with **"equivalent"** marked as likely.

**Say:** "First, redesign the age fit. Propagate the dust uncertainty properly, by integrating over each star's full dust estimate. Date each group from its bright stars, and require agreement with the spectroscopic ages. Then rebuild every downstream number. Then the decisive test: one age for all stars against one per group, with the criteria frozen before looking. If C is as old as A, the likely outcome is 'equivalent': the paper becomes a careful Gaia census and death ledger, with a methods lesson about dating hot massive stars. If C turns out genuinely young, it becomes a study of its own."

**Expect, "how long?":** "The age redesign takes a few days, rebuilding a day, and the ablation through to a manuscript a few weeks."
