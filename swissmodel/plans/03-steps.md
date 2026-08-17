# 03 — Steps: SwissModel CoreAPI test

## STEP 1: Learn API (DONE — research archived in 01-research.md)
Contract verified from live sources:
- `POST /automodel/` — body `{"target_sequences": <str|list>, "project_title": <str>}`, header `Authorization: Token <t>`.
  - 202 = queued; 200 = cached complete; 4xx = error.
- `GET /project/{id}/models/summary/` — status + models[]
- `GET /project/{id}/models/full-details/` — richer per-model payload (templates, scores)
- Auth: token from /account page; rate limit 200/min.

## STEP 2: Check auth availability
- [x] env / config search → no token found
- [ ] anonymous POST test → observe status code

## STEP 3: Submit one Omicron sequence (WRU87367.1, 1269 aa)
- Input: /tmp/opencode/seq1.fasta (gap-stripped, verified 1269 aa)
- If anonymous rejected → BLOCKER: need token (report exact error + how to get token)
- If accepted → record project_id

## STEP 4: Poll until COMPLETED/FAILED
- Poll full-details every 15 s, log status transitions, cap at 30 min

## STEP 5: Parse templates
- Extract: template SMTL ids, seq identity, coverage, GMQE, method, resolution
- Write summary JSON + report

## STEP 6: Verify & report
- Report which templates matched the Omicron Spike sequence; note runtime; archive decisions
