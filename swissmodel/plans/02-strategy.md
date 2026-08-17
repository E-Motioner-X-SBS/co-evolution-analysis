# 02 — Strategy: SwissModel template matching for the Omicron dataset

## Chosen strategy
1. **Primary test: guest web path** (no account needed) — submit seq via `POST /interactive`, poll
   `/interactive/{pid}/templates/`, read `automodel_poll` JSON + embedded `globalTemplateData`.
   Chosen because: no credentials available; returns exactly the template-match list the user asked about.
2. **coreapi path documented + code-ready** (`--mode coreapi`) for when a token is available
   (`SWISSMODEL_API_TOKEN`); endpoint contract fully verified from OpenAPI spec + official wrapper.
3. Deliverable: `swissmodel_template_search.py` in `co-evolution-analysis/` (reusable for any of the 1,299 sequences).

## Alternatives rejected
- Reverse-engineering the minified interactive JS app → unnecessary; `globalTemplateData` is
  server-embedded and parseable.
- Anonymous `template_identification` endpoint → **404 verified**, removed.
- SWISS-MODEL Repository API → returns pre-computed models by UniProt AC, not template matching.

## Validation plan
- [x] crambin (46 aa) end-to-end: submit → poll → 9 templates parsed with full metrics
- [ ] Spike seq1 (WRU87367.1, 1269 aa): project 8XmGU6, search running
- [ ] Cross-check: expected Spike templates ≈ known SARS-CoV-2 Spike PDB structures (6vxx, 7a98, 7bnm, ...)
