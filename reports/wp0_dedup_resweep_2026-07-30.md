# S1 — WP0 dedup sweep, re-run 2026-07-30

*Item **S1** of [pre_wp10_assessment_brief.md](../tasks/pre_wp10_assessment_brief.md).
WP0's own instruction is that the dedup sweep is repeated cheaply mid-project
and again immediately before submission. The previous sweep is dated
2026-07-20 ([wp0_dedup_report.md](../wp0_dedup_report.md)). Machine-readable
record: [wp0_dedup_resweep_execution.json](../provenance/wp0_dedup_resweep_execution.json).*

> **Verdict: unchanged. No direct duplicate.** Four new adjacent papers must be
> cited; none performs a stellar-demographic supernova history of Cyg OB2.
> **This sweep does not satisfy the pre-submission requirement** — it must be
> run again immediately before arXiv posting, since the gap between now and
> submission is exactly the window it cannot cover.

---

## 1. What was searched

**arXiv API**, sorted by submission date, seven queries, everything since
2025-07-01 inspected: `"Cygnus OB2"` (all fields), `Cygnus`+`"supernova
history"`, `"OB association"`+`supernova`+`Gaia`, `"star cluster"`+`"supernova
rate"`+`Gaia`, `IMF`+`missing`+`"massive stars"`+`association`, `"Cygnus
cocoon"`, `"PSR J2032+4127"`. Fourteen papers in window.

**Citation sweep** via the Semantic Scholar graph API for everything citing
Härer et al. (2025), Martin et al. (2010) and Berlanas et al. (2020). Härer+25
has **3 citing papers**, all interpretation-side; Martin+10 has 26 and
Berlanas+20 has 12, of which 6 are from 2025 onward.

The duplicate criterion is WP0's, unchanged: a paper counts as a duplicate only
if it combines Gaia-era membership, stellar mass/extinction inference, an IMF
normalization or equivalent demographic reconstruction, a past-supernova
history or last-supernova posterior, **and** an explicit test of the Cygnus
γ-ray interpretation.

## 2. Verdict

**No paper meets the criterion.** Every 2026 paper touching the Cygnus PeVatron
works the interpretation side — transport, diffusion, alternative accelerators
— and adopts its stellar inputs from the literature. Nobody has done the
stellar census.

The three papers citing Härer+25 make the point: *Microquasar Cygnus X-3 as the
PeVatron powering the Cygnus Bubble*, *Multimessenger Concordance for the
Cygnus Region as the Source of the Cosmic-Ray Knee*, and *The Cosmic-Ray Knee
as a Local Signature of Nearby PeVatrons*. All three argue about **what powers
the bubble**; none measures **how many supernovae the association has
produced**. That is still the gap this paper fills.

## 3. New papers that must be cited

| paper | why |
|---|---|
| **Microquasar Cygnus X-3 as the PeVatron powering the Cygnus Bubble** (arXiv:2607.07100, 2026-07-08) | **The framing must change slightly.** LHAASO's variable UHE source associated with Cyg X-3 offers a *third* explanation for the bubble, alongside collective winds and a supernova. The introduction currently presents a two-way choice; it should present a three-way one. This does not weaken the paper — it makes the supernova budget a discriminant rather than an assumption — but presenting an outdated dichotomy would be a gift to a referee. |
| **Suppressed diffusion and gamma-ray emission from the Cygnus Bubble** (arXiv:2606.03881, 2026-06-02) | Models the same emission under wind-termination-shock versus central-source acceleration; the central-source option is the one our ledger constrains. |
| **The lack of fast rotators in Cyg OB2. I** (arXiv:2510.15540, 2025-10-17) | **Touches our inputs.** A spectral reclassification of Cyg OB2's B0 population. Our anchor table draws spectral types from Wright+15, Berlanas+19/20 and GOSC; a reclassification of B0 stars could move individual anchor temperatures. Impact is expected to be small — anchors enter through extinction and the HRD, and the supernova budget lives above 34 M☉ — but it should be checked and cited, not ignored. |
| **OB runaway stars originating in the Vel OB1 association** (arXiv:2604.11988, 2026-04-13) | Methodologically adjacent: Gaia runaway traceback in an OB association. Worth citing where our runaway recovery and its false-positive treatment are described. |

Also logged, not required: *Massive Clusters and OB Associations as Output of
Massive Star Formation in Gaia Era* (2025, review context) and *Morphology of
Young Massive Stellar Clusters with Next-Generation IACTs* (arXiv:2509.20150).

## 4. Action taken

The four papers above are added to `manuscript/references.bib` and cited; the
introduction's framing is widened from two candidate power sources to three.
No result changes.

## 5. What this sweep does not do

- **Proceedings and theses.** WP0 flags ICRC/TeVPA proceedings as a place
  in-progress duplicates hide. This sweep covers arXiv and the Semantic Scholar
  citation graph only, and did not scan the 2026 proceedings. Neither did it
  reach ADS directly.
- **The pre-submission window.** The sweep must be repeated immediately before
  arXiv posting. Its value is inversely proportional to how stale it is, and it
  costs minutes.
