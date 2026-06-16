# Line-by-line correctness audit of the reversibility cascade + all variants

**Trigger:** after the variant-3.1 sign-inversion bug (a directional call read
the eQuilibrator reversibility-index sign backwards), audit *every* variant /
config knob the same way to ensure no other bug of that class — or any other —
remains.
**Scope:** `scripts/reversibility_lib.py` (the cascade + `ReversibilityConfig`)
and `scripts/variant_catalog.py` (the 14 variant cfg builders), checked against
the upstream parity reference
`ModelSEEDDatabase/Scripts/Thermodynamics/Estimate_Reaction_Reversibility.py`
and the intended behavior in `REACTION_REVERSIBILITY_HEURISTICS_REVIEW.md`.
**Verdict: CLEAN.** The only bug was the already-fixed 3.1 sign inversion.
No further bugs were found by any of the three independent methods below.
**Generated:** 2026-06-16.

---

## Method — three independent layers

1. **Empirical battery (output vs ground truth).** The same test that exposed
   3.1: run the cascade for every variant on all 56,012 MSDB reactions and
   check, among the directional (`>`/`<`) calls, agreement with the ΔG′° sign;
   plus monotonicity of the reversible-count (`=`) deltas vs the documented
   intent of each knob. A sign-inversion bug shows as ~0% sign-agreement (3.1
   buggy = 0%); a correct knob stays at the cascade's natural ~95–97%.
2. **Manual line-by-line read.** Every branch of `estimate_one`,
   `_walk_stoichiometry`, `_stored_bounds`, `_low_energy_points`,
   `_abc_transporter_decision`, `_is_atp_synthase`, and each `_*_cfg` builder,
   read against the upstream module and the sign convention ΔG<0 ⇒ `>`.
3. **Adversarial workflow.** 15 knobs × (1 auditor + 1 independent refuter) =
   30 agents. Auditors hunt for bugs line-by-line; verifiers try to *refute*
   each verdict from first principles. Result: 0 bugs, 0 suspicious, 0
   auditor/verifier disagreements — all 15 knobs `correct` by both.

---

## Layer 1 — empirical battery (all 14 variants, EQ level)

`>agree%` / `<agree%` = fraction of `>` calls with ΔG′°<0 / `<` calls with
ΔG′°>0. `Δ=` = change in reversible-count vs baseline.

| variant | `>` | `<` | `=` | `?` | >agree% | <agree% | Δ= vs base | intent check |
|---|--:|--:|--:|--:|--:|--:|--:|---|
| baseline | 9751 | 1957 | 13202 | 31102 | 95.9 | 95.2 | +0 | reference |
| 3.1 | 10823 | 2201 | 11886 | 31102 | **96.4** | **95.8** | −1316 | resolves `=`→dir; sign now healthy (was 0%) |
| 3.3 | 10627 | 2141 | 12142 | 31102 | 95.9 | 95.4 | −1060 | Bennett conc shifts reference ΔG |
| 3.3_wide | 8379 | 1533 | 14998 | 31102 | 96.3 | 96.7 | +1796 | wider window → fewer directional ✓ |
| 3.5 | 9453 | 1951 | 13506 | 31102 | 96.2 | 95.2 | +304 | k·σ band, net widening ✓ |
| 3.5_wide | 9311 | 1851 | 13748 | 31102 | 96.4 | 95.8 | +546 | wider err → fewer directional ✓ |
| 3.6 | 6869 | 1769 | 16272 | 31102 | 96.7 | 94.9 | +3070 | drop points rule → more `=` ✓ |
| 3.7 | 9751 | 1957 | 13202 | 31102 | 95.9 | 95.2 | +0 | CO₂ override unreachable w/o H3 → no-op ✓ |
| 3.10_tight | 9758 | 1957 | 13195 | 31102 | 95.9 | 95.2 | −7 | tighter band → fewer `=` ✓ |
| 3.10_loose | 9137 | 1875 | 13898 | 31102 | 96.4 | 95.0 | +696 | looser band → more `=` ✓ |
| H1 | 9751 | 1957 | 6680 | 37624 | 95.9 | 95.2 | −6522 | bare-default `=`→`?` (?+6522) ✓ |
| H2 | 9751 | 1956 | 13203 | 31102 | 95.9 | 95.2 | +1 | O₂/H₂ 1µM override (1 rxn) ✓ |
| H3 | 10562 | 2079 | 12269 | 31102 | 95.3 | 94.9 | −933 | ABC transporters forced `>` ✓ |
| H4 | 11220 | 2291 | 11399 | 31102 | 96.6 | 95.8 | −1803 | 3.1+3.5+Bennett stack ✓ |

Every variant's directional sign-agreement is in the same 95–97% band as
baseline (the residual ~4–5% is the legitimate concentration-window effect of
heuristic 1, present in baseline too). No variant is an outlier; a surviving
3.1-class inversion would have stood out at ~0%. All `Δ=` signs and magnitudes
match the knob's documented intent.

---

## Layer 2/3 — per-knob verdicts

All 15 audit units verified **correct** by manual read and by both workflow
agents (auditor + adversarial verifier). Key property checked per knob:

| knob (variant) | property verified | verdict |
|---|---|:--:|
| `ln_ri_by_rxn` (3.1) | NEGATIVE index → `>` (matches ΔG sign); threshold \|ln\|>6.9; fires after heuristic-1/ATPS/ABCT, before mMdeltaG band; inert at default | ✅ fixed & correct |
| `sigma_band_k` (3.5) | band = k·σ; only assigns `=`; precedence over mm_band; σ=0 edge handled | ✅ |
| `sigma_bounds_k` (3.5_wide) | err = σ·k applied ±symmetrically in `_stored_bounds`; widening → fewer directional | ✅ |
| `mm_band` (3.10) | symmetric \|mMdeltaG\|≤band; tight=1.0 fewer `=`, loose=4.0 more `=` | ✅ |
| `per_met_conc_range` + `per_met_conc` (3.3) | keyed on real `rgt["compound"]`; coeff-sign split correct; `BENNETT_2009_MEAN` = geomean of range | ✅ |
| `cell_min`/`cell_max` (3.3_wide) | default window widening → more `=`; matches +1796 | ✅ |
| `low_energy_cpds` (3.6) | empty tuple → cpds loop contributes 0; points·mMdeltaG>2 sign correct (mMdeltaG<0→`>`); only removes directional calls | ✅ |
| `apply_special_conc` + `co2_local_conc` (3.7) | override unreachable in default shadow path → False is a genuine no-op | ✅ |
| `default_direction` (H1) | only the bare fallthrough uses it; no other return path; `=`→`?` only | ✅ |
| `fix_low_local_conc` (H2) | O₂/H₂ 1e-6 override reachable independent of H3 (elif keys on real compound); verified shift = ln(1000) | ✅ |
| `fix_phosphates_shadow` (H3) | fixed path tests `cpd in PHOSPHATE_IDS` + sets real `local_cpd`; default path faithfully reproduces the upstream shadow bug; ABC sign ATP coeff<0→`>` | ✅ |
| `p_forward_threshold` | Gaussian-tail z-scores correct: z_f=(−band−mMdeltaG)/σ→P(fwd)→`>`; z_r=(−band+mMdeltaG)/σ→P(rev)→`<` | ✅ |
| `_H4_cfg` composite | sets ln_ri + sigma_band_k=1.96 + Bennett range/mean; ln_ri fires before sigma_band | ✅ |
| cascade branch order + baseline parity | matches upstream order exactly; new branches (ln_ri, p_forward) inert at default; reproduces MSDB byte-for-byte | ✅ |
| `_stored_bounds` + `_walk` sign math | stored_max pairs pdt_max+rct_min (+err) = max ΔG′; stored_min pairs pdt_min+rct_max (−err) = min ΔG′; no extreme-pairing swap | ✅ |

---

## Deeper checks worth recording

- **Baseline parity is structurally guaranteed.** `ReversibilityConfig()` sets
  `ln_ri_by_rxn=None`, `p_forward_threshold=None`, `sigma_band_k=None`,
  `sigma_bounds_k=None`, `mm_band=2.0`, `default_direction="="`,
  `apply_special_conc=True`, `fix_*_shadow=False`. With those, every new branch
  is skipped and the cascade walks the same path as upstream
  (`stored_max<0→">"` → `stored_min>0→"<"` → ATPS → ABCT → `|mMdeltaG|≤2→"="`
  → points rule → `"="`). The 3.1 fix only changed the comparison *inside* the
  `ln_ri` branch, which the default path never enters — so the fix cannot have
  perturbed the baseline, and the documented byte-for-byte MSDB reproduction
  still holds.

- **The shadow bugs are faithfully preserved, not new bugs.** In the default
  path the `for cpd in PHOSPHATE_IDS` loop leaves `local_cpd = cpd00012` (PPi)
  for every reagent, so (a) the PROTON/WATER skip never fires and (b) the
  CO₂/LOW_LOCAL special-concentration overrides are unreachable — exactly as
  upstream. `fix_phosphates_shadow=True` (H3) and `fix_low_local_conc=True`
  (H2) are the *opt-in repairs*. This was scrutinized specifically because it
  *looks* like a wrong-variable bug; it is a deliberate, documented
  reproduction of upstream behavior and is gated behind the H2/H3 knobs.

- **H2 is independent of H3.** Even with the phosphate shadow unfixed, the
  O₂/H₂ override is reachable because its branch is an `elif` keyed on the real
  `rgt["compound"]` (reached because the shadowed `local_cpd` ≠ CO₂).
  Empirically confirmed: the override shifts an O₂(−1) reaction's `rgt_sum` by
  exactly `ln(1e-3/1e-6) = 6.9078`, and H2 flips the 1 reaction the manifest
  records.

---

## Conclusion

The reversibility cascade and all 14 variant configurations are correct. The
single bug in this subsystem — the variant-3.1 ln(reversibility-index) sign
inversion — was found, fixed (`reversibility_lib.py`: `">" if ln_ri < 0 else
"<"`), and the whole pipeline re-run (see
`VARIANT_3.1_BREAKAGE_INVESTIGATION.md`). Three independent audit methods
(empirical sign/monotonicity battery, manual line-by-line read, and a 30-agent
adversarial workflow) found no further bugs.

### Reproduce
```bash
cd core_models_analysis/
# empirical battery: sign-agreement + monotonicity per variant
#   (uses thermo_variants/*/Estimated_Reaction_Reversibility_Report_EQ.txt
#    + ModelSEEDDatabase MetaNetX_Reaction_Energies.tbl for ΔG′°)
# adversarial audit workflow: 15 knobs × (audit + verify)
```
