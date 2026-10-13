"""pytest view of a Wright-Fisher validation run.

The simulations take tens of minutes, so they are run with suite.py and this file checks the saved
results:

    python tests/wf_validation/suite.py --outdir tests/wf_validation/runs/<name>
    WFV_RESULTS=tests/wf_validation/runs/<name>/results.json pytest tests/wf_validation/test_wf_validation.py

Without WFV_RESULTS the simulation checks are skipped (the fast theory tests in test_theory.py still run).
"""

import json, os
import pytest

RESULTS = os.environ.get("WFV_RESULTS")
CHECKS = []
if RESULTS:
    with open(RESULTS) as fh:
        CHECKS = [c for c in json.load(fh)["checks"] if c.get("pass") is not None]


@pytest.mark.skipif(not RESULTS, reason="set WFV_RESULTS to a results.json written by suite.py")
@pytest.mark.parametrize("check", CHECKS, ids=[c["check"] for c in CHECKS])
def test_simulation_matches_theory(check):
    detail = {k: v for k, v in check.items() if k not in ("check", "pass")}
    assert check["pass"], f"{check['check']}: {detail}"
