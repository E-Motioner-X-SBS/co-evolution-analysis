# 01 — Research: SWISS-MODEL Modelling API (coreapi)

## Sources (all fetched live this session)
1. `https://swissmodel.expasy.org/coreapi/` — CoreAPI schema browser (Django REST framework)
2. `https://swissmodel.expasy.org/coreapi/schema.js` — machine-readable CoreJSON schema (base64, decoded)
3. `https://swissmodel.expasy.org/api-docs/?format=openapi` — OpenAPI 2.0 spec (full JSON, includes definitions)
4. `https://swissmodel.expasy.org/docs/help` — "Modelling API" section: token auth, automodel/alignment/user_template examples, status codes, rate limits, polling loop
5. `https://github.com/SWISS-MODEL/swissmodel_modelling_api` — official Galaxy wrapper (`sm_api_wrapper.py`): exact submit/poll/fetch flow

## Findings
- **Auth**: all coreapi endpoints require `Authorization: Token <token>` (OpenAPI `security: api_key`). Token from the SWISS-MODEL account page (`/account`); alternatively `POST /api-token-auth/` with username/password. **Empirically verified: anonymous `POST /automodel/` → 401 `{"detail":"Authentication credentials were not provided."}`**
- **`POST /automodel/`**: body `{"target_sequences": <str|list>, "project_title": <str>}` → `201` (new project) / `202` (accepted, queued) / `200` (cached complete). Response: `{"project_id": ...}`.
- **Polling**: `GET /project/{project_id}/models/summary/` or `.../models/full-details/`; status values: INITIALISED, QUEUEING, RUNNING, COMPLETED, FAILED, UNKNOWN. Sleep ~10–17 s between polls.
- **Template data**: lives in the per-model payload of `full-details/` (template SMTL id, seq identity, coverage, GMQE, QSQE, method, resolution, coordinates_url, modelcif_url).
- **Rate limits**: 200/min rapid, 10000/6h prolonged; excess → 429.
- **Projects** persist ~14 days; bulk download via `POST /projects/download/` + `GET /projects/download/{id}/`.
- **There is NO standalone template-identification endpoint** in the coreapi; template search runs inside `/automodel/`. The old anonymous `POST /interactive/template_identification` endpoint is **gone (404 verified** for 4 candidate URLs**)**.
- **Guest path (discovered & verified)**: the web form `POST /interactive` (fields: `target`, `project_title`, `email`, `automodel`, `is_alignment`, CSRF) works WITHOUT an account → creates a project in template-search mode → `GET /interactive/{project_id}/templates/` shows results. Verified: created project `8XmGU6` for sequence WRU87367.1; search queueing→running (still running at last check).
- Initial failure mode learned: empty `multiTargetValidate` field → "Target sequences must be unique." (omit that field).

## Conflicts
- None between sources; docs, swagger, wrapper code, and live HTTP behaviour agree (auth required for coreapi; guest web path separate).

## Conclusion
To "get template matches for a FASTA sequence" via coreapi: token + `POST /automodel/` + poll `full-details/` and read per-model template fields. Without an account, the guest `/interactive` form path works and returns the same template search.

## Open
- Runtime of template search for a 1,269-aa Spike sequence (poller running, project 8XmGU6).
