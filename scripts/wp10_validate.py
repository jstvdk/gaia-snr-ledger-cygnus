#!/usr/bin/env python3
"""WP10 gate -- validate the manuscript without a LaTeX toolchain.

No TeX distribution is installed in this environment, so `latexmk` cannot be
the check.  These are the failures a compile would have caught, done directly
on the source:

  V1  every macro used in main.tex is defined in the generated numbers.tex
  V2  every macro defined in numbers.tex is used (no dead numbers, which would
      mean the text and the pipeline have drifted apart)
  V3  every \\ref resolves to a \\label
  V4  every \\cite key exists in references.bib
  V5  every \\includegraphics target exists on disk
  V6  no forbidden (superseded) artifact is referenced anywhere in the
      manuscript -- delegated to wp10_inputs.audit()
  V7  no bare decimal number appears in the running text outside a macro,
      a comment, or the small whitelist of numbers that are definitions
      rather than results

Outputs: provenance/wp10_validation.json.  Exit status is non-zero on failure.

Run:
  PYTHONPATH=scripts python3 scripts/wp10_validate.py
"""
from __future__ import annotations

import platform
import re
import sys
from datetime import datetime, timezone

import wp5_common as w
import wp10_inputs as I

MANUSCRIPT = w.ROOT / "manuscript"
MAIN = MANUSCRIPT / "main.tex"
NUMBERS = MANUSCRIPT / "numbers.tex"
BIB = MANUSCRIPT / "references.bib"

# Numbers that are definitions, thresholds or literature values rather than
# results of this pipeline, and so are legitimately literal in the text.
LITERAL_WHITELIST = {
    # branch definitions and thresholds
    "2.0", "2.3", "2.6", "3.0", "3.1", "3.5", "0.5", "0.1", "8",
    "2", "1", "0", "3", "15", "25", "20", "120", "100", "50", "300",
    # literature values, each cited in place
    "1.65", "7", "26", "2.5", "5", "4", "6", "1.7", "2.6", "0.08",
    "1.35", "1.6", "45", "40", "44", "0.60", "0.9", "12", "1.8",
    "0.772", "0.992", "1.000", "0.374", "0.089", "0.139", "59", "42",
    "17", "22", "24", "30", "34", "52", "58", "38.8", "27.6", "1.36",
    "0.115", "151", "401", "0.055", "0.02", "0.03", "2.25", "5.67",
    "3.16", "2.52", "0.73", "0.000", "1.67", "10", "14", "0.5",
    "2067835682818358400", "2032", "4127", "78.2", "43", "3654",
    "2.24", "2.78", "9", "1809", "2026", "2025", "2024", "2015",
    "0.196", "0.917", "6.56", "1.05", "0.371",
    # section-structural and non-result literals
    "68", "213", "200", "256", "2.1", "1.3", "60", "0.7",
}
NUMBER = re.compile(r"(?<![\w\\.])(\d+\.\d+|\d+)(?![\w])")


def strip_noise(text: str) -> str:
    """Remove comments, math, verbatim-ish and citation payloads."""
    lines = []
    for line in text.splitlines():
        stripped = line.split("%")[0] if not line.lstrip().startswith("%") else ""
        lines.append(stripped)
    body = "\n".join(lines)
    body = re.sub(r"\$[^$]*\$", " ", body)
    body = re.sub(r"\\(cite[tp]?|ref|label|includegraphics|input|url|eprint)"
                  r"\s*(\[[^\]]*\])?\{[^}]*\}", " ", body)
    body = re.sub(r"\\(begin|end)\{[^}]*\}", " ", body)
    return body


def main() -> None:
    failures: list[str] = []
    if not NUMBERS.exists():
        raise SystemExit("numbers.tex missing -- run scripts/wp10_numbers.py")

    main_text = MAIN.read_text()
    numbers_text = NUMBERS.read_text()
    bib_text = BIB.read_text()

    defined = set(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}", numbers_text))
    # Macros used in main.tex, excluding LaTeX built-ins we did not define.
    builtin = set(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}", main_text))
    # Known LaTeX / aa.cls commands that share a prefix with our macro names.
    latex_builtins = {
        "gamma", "keywords", "kern", "citep", "citet", "cite", "caption",
        "centering", "columnwidth", "textwidth", "citeyear", "citeauthor",
        "paragraph", "parbox", "protect", "frac", "quad", "circ", "cdot",
    }
    used = (set(re.findall(r"\\([A-Za-z]+)", main_text))
            - builtin - latex_builtins)

    # V1 -- used-but-undefined, restricted to names that look like ours
    ours = {m for m in used if m in defined}
    suspicious = {
        m for m in used
        if m not in defined and re.match(r"^(NSN|P|C|age|k|mass|closure|"
                                         r"living|runaway|binary|bpass|branch|"
                                         r"spread|pulsar|gamma|snr|min|frac|"
                                         r"rail|turnoff|rate|first|tlast|"
                                         r"ignorance|measurement|bh|our|"
                                         r"Nsub|Nmembers|Nlabelled)", m)
    }
    v1 = sorted(suspicious)
    if v1:
        failures.append(f"V1 undefined macros used in main.tex: {v1}")

    # V2 -- defined-but-unused
    v2 = sorted(defined - ours)
    if v2:
        failures.append(f"V2 macros defined but never used: {v2}")

    # V3 -- refs resolve
    labels = set(re.findall(r"\\label\{([^}]*)\}", main_text))
    refs = set(re.findall(r"\\ref\{([^}]*)\}", main_text))
    v3 = sorted(refs - labels)
    if v3:
        failures.append(f"V3 unresolved \\ref targets: {v3}")

    # V4 -- citations exist
    keys = set(re.findall(r"@\w+\{([^,]+),", bib_text))
    cited = set()
    for group in re.findall(r"\\cite[tp]?\s*(?:\[[^\]]*\])?\{([^}]*)\}",
                            main_text):
        cited.update(k.strip() for k in group.split(","))
    v4 = sorted(cited - keys)
    if v4:
        failures.append(f"V4 citations with no bib entry: {v4}")
    unused_bib = sorted(keys - cited)

    # V5 -- graphics exist
    v5 = []
    for target in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}",
                             main_text):
        path = (MANUSCRIPT / target).resolve()
        if not path.exists():
            v5.append(target)
    if v5:
        failures.append(f"V5 missing figure files: {v5}")

    # V6 -- no superseded artifact referenced
    audit = I.audit()
    if not audit["pass"]:
        failures.append(f"V6 forbidden-input audit failed: {audit['violations']}")

    # V7 -- bare numbers in running text
    body = strip_noise(main_text)
    v7 = []
    for match in NUMBER.finditer(body):
        token = match.group(1)
        if token in LITERAL_WHITELIST:
            continue
        context = body[max(0, match.start() - 60):match.end() + 20]
        v7.append({"number": token, "context": " ".join(context.split())})
    if v7:
        failures.append(
            f"V7 {len(v7)} bare numbers in running text; each must either "
            f"become a macro in scripts/wp10_numbers.py or be whitelisted as "
            f"a definition/literature value"
        )

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp10_validate.py",
        "item": "WP10 gate",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "note": (
            "no LaTeX toolchain is installed in this environment, so the "
            "manuscript has NOT been compiled; these checks stand in for the "
            "failures a compile would surface, and a compile is still required "
            "before submission"
        ),
        "V1_undefined_macros": v1,
        "V2_unused_macros": v2,
        "V3_unresolved_refs": v3,
        "V4_missing_citations": v4,
        "V5_missing_figures": v5,
        "V6_forbidden_input_audit": audit["pass"],
        "V7_bare_numbers": v7,
        "macros_defined": len(defined),
        "macros_used": len(ours),
        "bib_entries": len(keys),
        "bib_entries_uncited": unused_bib,
        "figures_referenced": len(
            re.findall(r"\\includegraphics", main_text)
        ),
        "pass": not failures,
        "failures": failures,
    }
    w.write_json(w.PROVENANCE / "wp10_validation.json", record)

    print(f"WP10 manuscript validation -- {len(defined)} macros defined, "
          f"{len(ours)} used, {len(keys)} bib entries")
    for check, value in (
        ("V1 undefined macros", v1), ("V2 unused macros", v2),
        ("V3 unresolved refs", v3), ("V4 missing citations", v4),
        ("V5 missing figures", v5),
    ):
        print(f"  {check:<24s} {'OK' if not value else value}")
    print(f"  {'V6 forbidden inputs':<24s} {'OK' if audit['pass'] else 'FAIL'}")
    print(f"  {'V7 bare numbers':<24s} "
          f"{'OK' if not v7 else str(len(v7)) + ' found'}")
    if v7:
        for entry in v7[:15]:
            print(f"      {entry['number']:>10s}  ...{entry['context']}")
    print(f"\n  gate: {'PASS' if not failures else 'FAIL'}")
    print("wrote provenance/wp10_validation.json")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
