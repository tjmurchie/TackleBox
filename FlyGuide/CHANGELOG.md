# FlyGuide Changelog

## 1.3.1 — 2026-08-31 (guidecheck.sh: fix hard crash + silent column corruption)

Two real bugs found via a live ~10,000-species PalaeoSCOPE Phase B run, both in
`guidecheck.sh`:

1. **Hard crash on a transient bad NCBI response.** NCBI EUtils occasionally returns a
   non-JSON body under sustained load (a transient rate-limit/error page rather than the
   documented JSON response) -- `eutils()` piped this straight into `jq` with no
   validation at all, and `set -euo pipefail` then killed the ENTIRE script on `jq`'s own
   parse error ("Invalid numeric literal at line 1, column 10"), losing every
   already-processed row along with it. Hit twice in a row on the real run (species #3,
   then again at #37 on retry). `eutils()` now validates each response is real JSON
   (`jq empty`) before returning it, retrying up to `GUIDECHECK_EUTILS_MAX_ATTEMPTS`
   times (default 3, `GUIDECHECK_EUTILS_RETRY_DELAY` seconds apart, default 1) and
   degrading to a safe empty JSON object (`{}`) after exhausting retries rather than ever
   propagating bad data downstream -- every caller's own `jq` filter already treats a
   missing/null field as "no data" (NO_TAXID / NA / a zero count), so this degrades one
   query to "nothing found" instead of crashing the whole run.
2. **Silent column corruption for almost every real match.** `matched_name`/`rank` were
   extracted via `get_tax_summary "$taxid" | awk -F'\t' '{print $1, $2}'` -- but `awk`'s
   `print $1, $2` re-joins the two TSV fields with its own OFS (a plain space), which
   destroys the original tab boundary before `read -r matched_name rank` ever sees it.
   Any real binomial `matched_name` (e.g. "Mammuthus primigenius" -- nearly every real
   successful match, not an edge case) has its own internal space, so it was silently
   split into `matched_name="Mammuthus"` / `rank="primigenius species"` instead of the
   correct `matched_name="Mammuthus primigenius"` / `rank="species"`. Fixed by reading the
   TSV line directly with `IFS=$'\t' read -r matched_name rank`, preserving the real
   tab-delimited field boundary regardless of internal spaces.

New test suite: `tests/test_guidecheck.py` (guidecheck.sh's first automated test coverage
of any kind), using a fake `curl` on PATH to deterministically simulate NCBI responses
without touching the real network. 10 passed (was 6), covers both fixes plus the
one-bad-species-doesn't-kill-the-batch guarantee.
