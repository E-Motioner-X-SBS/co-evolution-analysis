# Evolution Log — K-map ↔ 3D-Structure Campaign (§17 self-evolution)

## What I tried / what broke / what I learned

| Cycle | Tried | Broke | Learned | Spec change |
|---|---|---|---|---|
| 1 | Draft runner in one pass | Invalid decorator line, flat-dict aggregation silently dropping binary metrics, mapqm OFF-vs-DC truth semantics | Complex multi-stage drivers need a written invariant list BEFORE coding; nested record layouts must be explicit in the aggregation contract | plans/11b gained "truth-table semantics" block |
| 2 | nohup background launches | Harness kills process group on tool timeout — jobs died mid-cache | Use `setsid … < /dev/null & disown` for anything long-running; make every stage resumable by design | cache made per-target resumable |
| 3 | Dual relaunch overlap | Two cache workers raced on same npz (non-atomic save) | Validate-after-write is mandatory for concurrent-written artifacts; single-writer discipline | post-race full-integrity re-validation added (0 corrupt) |
| 4 | Trusting remembered literature numbers ("PSICOV L/5 = 0.44") | Our measurement said 0.727 — 1.6σ conflict with memory | Web-verify EVERY external number against the primary source; ours matched Table-1 exactly (Δ=0.003). Memory of headlines is not a citation | S1 check added to skepticism battery; correction banners policy |
| 5 | Gray enrichment taken at face value (p=1e-77) | Enrichment could be generic chemistry-clustering, not Gray-specific | Adversarial control (code permutation) separates the two: z=2.51 for the specific encoding; honest two-layer reading recorded | permutation control became standard for any encoding claim |
| 6 | Lean proofs written generically first | Instance-synthesis failures (Decidable for nested quantifiers), fragile simp on bit-arith | In this repo's house style, prefer Bool-mirror + native_decide over Prop-quantifier decide; prove arithmetic corollaries with mul_add_div/mod + omega, never simpa-on-xor | ContactCircuits final form uses Bool mirrors throughout |
| 7 | Patch-on-patch edits to the Lean file | File drifted; compiler saw content I thought was replaced | After >2 failed patches, do a full clean rewrite and re-read from disk before compiling | applied immediately |

## Anti-patterns discovered
- **Remembered-number citation**: always re-fetch the paper's actual table.
- **Silent metric-band mixing**: single-number tables hide short-range inflation; always stratify.
- **Decorator/placeholder leftovers in long files**: syntax-check imports before running.
- **Background jobs without setsid**: harness timeout = job death.
- **Non-atomic writes to shared paths**: tmp+rename or single-writer only.

## Strengths confirmed
- Adversarial test-first library development caught 3 real bugs pre-launch.
- Resumable stage design turned job-kills into non-events.
- Literature-calibrated ground truth (identity mapping) eliminated an entire error class.
- Determinism (fixed seeds everywhere) made the ×5 verification trivial to automate.

## Mutation history
| From | To | Trigger | Outcome |
|---|---|---|---|
| Prop-level decide for circuit assertions | Bool-mirror + native_decide | instance-synthesis failures | clean compile |
| Single-number PSICOV comparison | band-stratified + web-reconciled | 0.72-vs-0.44 conflict | corrected docs, calibrated pipeline |
| Spike-only Gray test | 150-protein powered test + permutation control | p=0.345 ambiguity | decisive p=1e-77 + specificity z=2.51 |
