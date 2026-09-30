"""Engineering entry points for gray PBTE exploration; no numerical kernels."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "studies/gray_exploration/study.yaml"
PROFILE = ROOT / "studies/gray_exploration/profiles/legacy_t01_smoke.yaml"
REGISTRY = ROOT / "methods/registry.yaml"
BASELINE = ROOT / "docs/refactor/baseline_manifest.json"
