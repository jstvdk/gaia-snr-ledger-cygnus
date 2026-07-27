#!/usr/bin/env python3
"""Build the commented audit notebook for the versioned WP3--WP5 repair."""
from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parents[1]
output = ROOT / "notebooks" / "wp3_extinction_repair_and_wp5_regate.ipynb"
nb = nbf.v4.new_notebook()
nb["metadata"]["kernelspec"] = {
    "display_name": "Python (cygob2-gaia)",
    "language": "python",
    "name": "cygob2-gaia",
}
nb["cells"] = [
    nbf.v4.new_markdown_cell(
        "# WP3 extinction repair → WP4 posterior masses → WP5 re-gate\n\n"
        "This notebook audits the versioned `repair_v1` artifacts. It does not "
        "overwrite or recompute the frozen WP3/WP4/WP5 products. The old WP5 "
        "blocking verdict is retained as the repair trigger."
    ),
    nbf.v4.new_markdown_cell(
        "## 1. Load the frozen and repaired summaries\n\n"
        "Purpose: establish row counts and the exact versioned files under test. "
        "The notebook intentionally reads generated products; estimator execution "
        "details and hashes live in `provenance/`."
    ),
    nbf.v4.new_code_cell(
        "from pathlib import Path\n"
        "import json\n"
        "import numpy as np\n"
        "import pandas as pd\n"
        "from IPython.display import display, Image, Markdown\n\n"
        "ROOT = Path('..').resolve()\n"
        "gate = json.loads((ROOT/'provenance/wp3_repair_gate.json').read_text())\n"
        "manifest = json.loads((ROOT/'provenance/wp3_wp5_repair_manifest.json').read_text())\n"
        "display(pd.DataFrame([manifest['component_gates']]))\n"
        "print('overall status:', manifest['status'])"
    ),
    nbf.v4.new_markdown_cell(
        "## 2. Permanent WP3 invariant (F4)\n\n"
        "Thresholds: in the fixed 2–12 M☉ diagnostic bins, every median "
        "`A_V − local 8-anchor map` residual must have absolute value below "
        "0.30 mag. The rank test is applied to the bin medians (not thousands "
        "of individual stars, which would turn negligible effects into tiny "
        "p-values) and must have p ≥ 0.05. These thresholds were set in the "
        "repair brief before this run."
    ),
    nbf.v4.new_code_cell(
        "f4 = pd.read_csv(ROOT/'tables/wp3_extinction_mass_invariant_repair_v1.csv')\n"
        "display(f4)\n"
        "assert gate['f4']['pass']\n"
        "print(gate['f4'])"
    ),
    nbf.v4.new_markdown_cell(
        "## 3. WP4 posterior quality and massive-star recensus\n\n"
        "The acceptance checks guard against discrete isochrone snapping and "
        "require the previously empty 2.5–3.2 M☉ interval to be populated. "
        "`N(>8 M☉)` is reported both from posterior medians and as the sum of "
        "per-star posterior probabilities; WP6 must not consume it while the "
        "WP5 gate is blocked."
    ),
    nbf.v4.new_code_cell(
        "mass = pd.read_parquet(ROOT/'data/processed/wp4_mass_posteriors_repair_v1.parquet')\n"
        "print('largest exact fraction:', mass.mass_baseline.value_counts(normalize=True).max())\n"
        "print('N 2.5–3.2:', mass.mass_baseline.between(2.5,3.2).sum())\n"
        "print('median-mass N>8:', (mass.mass_baseline>8).sum())\n"
        "print('posterior E[N>8]:', mass.mass_baseline_p_gt8.sum())\n"
        "assert gate['wp4']['pass']"
    ),
    nbf.v4.new_markdown_cell(
        "## 4. Honest WP5 response and unchanged residual gate\n\n"
        "Injected six-band photometry passed through the actual repaired F1/F2 "
        "estimators. Baseline response resolution is 64 mass draws per recovered "
        "injection. The gate remains χ² p≥0.01, residual-trend p≥0.05, and "
        "maximum |Pearson residual|≤3 for every subgroup. No criterion is relaxed."
    ),
    nbf.v4.new_code_cell(
        "baseline = pd.read_parquet(ROOT/'data/processed/wp5_imf_normalization_repair_v1.parquet')\n"
        "baseline = baseline.query(\"family == 'PARSEC' and R_V == 3.1 and alpha == 2.3\")\n"
        "display(baseline[['subgroup','poisson_chi_square_p','residual_trend_p',\n"
        "                  'max_abs_pearson_residual','residual_gate_pass']])\n"
        "assert not gate['wp5']['all_subgroups_pass']\n"
        "print('Expected blocker:', gate['blocking_reason'])"
    ),
    nbf.v4.new_markdown_cell(
        "## 5. Completeness and mass sanity\n\n"
        "No absolute 95% edge is asserted when the bright plateau itself is "
        "below 95%. Relative-to-plateau edges remain diagnostic only. The "
        "association mass must remain within a factor two of 16,500 M☉."
    ),
    nbf.v4.new_code_cell(
        "display(pd.read_csv(ROOT/'tables/wp5_completeness_baseline_repair_v1.csv'))\n"
        "print(gate['wp5']['association_mass_median_Msun'],\n"
        "      gate['wp5']['association_mass_lo68_Msun'],\n"
        "      gate['wp5']['association_mass_hi68_Msun'])\n"
        "assert gate['wp5']['within_factor_two_literature']"
    ),
    nbf.v4.new_markdown_cell(
        "## 6. Visual gate summary and decision\n\n"
        "The left panel shows that the original mass-correlated extinction "
        "failure is removed. The right panel shows why the repair is still "
        "blocked: the converged CygOB2-C baseline response exceeds the fixed "
        "3σ residual boundary."
    ),
    nbf.v4.new_code_cell(
        "display(Image(filename=str(ROOT/'figures/wp3_repair/wp3_wp5_repair_gates_repair_v1.png')))\n"
        "display(Markdown((ROOT/'reports/WP3_WP5_REPAIR_REPORT_repair_v1.md').read_text()))"
    ),
]
nbf.write(nb, output)
print(output)
