"""Mock-only provider lifecycle. No commercial CLI or live model calls."""
from __future__ import annotations

import json
import subprocess
import sys

SCENARIOS = ("success", "timeout", "error", "empty_transcript", "missing_final", "invalid_output")


def mock_provider(scenario: str, method: str = "cis") -> dict:
    if scenario not in SCENARIOS or method not in ("cis", "krylov"):
        raise ValueError("Unknown mock scenario or method")
    # An actual child process tests timeout, return-code and transcript boundaries.
    script = {
        "success": "import json; print(json.dumps({'type':'candidate','method':%r})); print(json.dumps({'type':'result','status':'completed'}))" % method,
        "timeout": "import time; time.sleep(5)",
        "error": "import sys; print('mock execution error',file=sys.stderr); sys.exit(7)",
        "empty_transcript": "pass",
        "missing_final": "print('{\"type\":\"candidate\",\"method\":\"cis\"}')",
        "invalid_output": "print('invalid JSON')",
    }[scenario]
    result = {"provider": "mock", "model_id": None, "cost_usd": None,
              "scenario": scenario, "live_agent_enabled": False}
    try:
        proc = subprocess.run([sys.executable, "-c", script], capture_output=True,
                              text=True, timeout=0.05 if scenario == "timeout" else 10,
                              env={"PYTHONNOUSERSITE": "1"})
    except subprocess.TimeoutExpired:
        return {**result, "termination": "timeout", "returncode": None, "events": [], "transcript": "", "stderr": ""}
    result.update(returncode=proc.returncode, transcript=proc.stdout, stderr=proc.stderr)
    events = []
    if proc.returncode:
        status = "execution_error"
    elif not proc.stdout.strip():
        status = "empty_transcript"
    else:
        try:
            events = [json.loads(line) for line in proc.stdout.splitlines()]
            status = "completed" if events[-1].get("type") == "result" else "missing_final_event"
        except (ValueError, AttributeError):
            status = "invalid_output"
    return {**result, "termination": status, "events": events}
