# 00 — Understanding: SwissModel CoreAPI test for Omicron Spike sequences

## Task
1. Learn how to use the SWISS-MODEL Modelling API (`https://swissmodel.expasy.org/coreapi/`, docs at `/api-docs/` + `/docs/help#modelling_api`).
2. Test whether template matches are returned for ONE FASTA sequence from the 1,299-seq Omicron Spike alignment used in the co-evolution analysis.

## Constraints / facts (verified this session)
- API is token-authenticated: `Authorization: Token <token>`; token from https://swissmodel.expasy.org/account (or POST /api-token-auth/ with username/password — we have no account).
- `POST /automodel/` body `{"target_sequences": str-or-list, "project_title": ...}` → 201/202 (202 = job queued, 200 = cached complete).
- Poll `GET /project/{id}/models/summary/` and `.../full-details/` until status COMPLETED|FAILED; sleep 10–17 s between polls.
- Models JSON carries: template PDB/SMTL id, seq identity, coverage, GMQE, QSFE, coordinates_url, modelcif_url.
- Rate limits: 200/min rapid, 10000/6h prolonged → 429 on excess.
- Spike sequence: full-length ~1270 aa; automodel will run BLAST+HHblits template search against SMTL + AlphaFold DB, then model building (minutes).

## Success criteria
- [ ] Correct API call sequence documented (submit → poll → parse templates)
- [ ] One real Omicron sequence (WRU87367.1, 1269 aa gap-stripped) submitted
- [ ] Template list with PDB IDs + metrics retrieved (or a clear, verified blocker with exact error)

## Open questions
- Is an account/token available? (no token found in env/config so far)
- Does /automodel/ accept anonymous requests? (test empirically)
- How long does a 1270-aa automodel job take? (poll)
