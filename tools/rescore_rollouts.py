#!/usr/bin/env python3
"""Re-grade legacy archives into a NEW evaluation directory.

Example (potentially expensive; no model calls):
    python tools/rescore_rollouts.py experiments/results/t01-square-F3 \
        --task tasks/t01-square --sandbox-root /tmp/pbte-rescore \
        --trials 5 --out runs/f3-rescore-001

--in-place-copy only admits an explicit writable copy under workspaces/ or
outside this repository. Archived-code execution is not a validated OS sandbox.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from exploration.records import (identifier, new_directory, now, safe_extract,
                                 sha256, tree_hashes, write_json)


def grade(task: Path, env: Path, out_json: Path, timeout: float = 1800) -> dict:
    try:
        proc = subprocess.run(
            [sys.executable, str(task / 'verifier/checks.py'), '--env', str(env),
             '--json', str(out_json)], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {'gate': False, 'score': None, 'families': {}, 'termination': 'timeout'}
    out_json.with_suffix('.stdout.log').write_text(proc.stdout)
    out_json.with_suffix('.stderr.log').write_text(proc.stderr)
    if not out_json.exists():
        return {'gate': False, 'score': None, 'families': {}, 'termination': 'verifier_error',
                'returncode': proc.returncode, 'error': (proc.stderr or proc.stdout)[-2000:]}
    try:
        result = json.loads(out_json.read_text())
    except ValueError:
        return {'gate': False, 'score': None, 'families': {}, 'termination': 'invalid_output'}
    result['termination'] = ('passed' if result.get('gate') else 'algorithm_failed') if proc.returncode in (0, 1) else 'verifier_error'
    if proc.returncode not in (0, 1):
        result['gate'] = False
        result['score'] = None
    return result


def row_from(res: dict) -> dict:
    cells = [c for f in res.get('families', {}).values() for c in f['cells']]
    return {'gate': bool(res.get('gate')), 'score': res.get('score'),
            'families': {k: (f['gate'], f['score']) for k, f in res.get('families', {}).items()},
            'cells_passed': sum(bool(c.get('passed')) for c in cells), 'cells': len(cells),
            'termination': res.get('termination', 'unknown')}


def check_writable_copy(results: Path) -> None:
    if results == ROOT or (ROOT in results.parents and results.relative_to(ROOT).parts[0] != 'workspaces'):
        raise ValueError('In-place rescore refuses protected evidence; prepare a new workspaces/ copy')
    for p in results.rglob('*'):
        if p.is_symlink() or (p.is_file() and p.stat().st_nlink != 1):
            raise ValueError('Writable copy must contain no symbolic or hard links')


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('results', type=Path)
    ap.add_argument('--task', type=Path, required=True)
    ap.add_argument('--sandbox-root', type=Path, required=True)
    ap.add_argument('--trials', type=int, nargs='+', default=None)
    ap.add_argument('--oracle', action='store_true')
    ap.add_argument('--timeout', type=float, default=1800)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--out', type=Path)
    mode.add_argument('--in-place-copy', action='store_true')
    args = ap.parse_args(argv)
    try:
        task, results = args.task.resolve(), args.results.resolve()
        if not results.is_dir() or not (task / 'verifier/checks.py').is_file():
            raise ValueError('Missing results directory or verifier')
        scratch = args.sandbox_root.expanduser().resolve()
        if scratch == ROOT or ROOT in scratch.parents:
            raise ValueError('Legacy archived-code execution scratch must be outside the repository')
        if args.timeout <= 0:
            raise ValueError('timeout must be positive')
        trials = args.trials if args.trials is not None else sorted(
            int(p.name.split('_')[1]) for p in results.glob('trial_*_environment.tar.gz'))
        if len(set(trials)) != len(trials) or any(i < 1 for i in trials):
            raise ValueError('Trial IDs must be positive and unique')
        if not trials and not args.oracle:
            raise ValueError('No candidate archives; no evaluation performed')
        for i in trials:
            if not (results / f'trial_{i:02d}_environment.tar.gz').is_file():
                raise ValueError(f'Trial {i}: missing archive')
        before = tree_hashes(results)
        if args.in_place_copy:
            check_writable_copy(results)
            target = results
            evidence = Path(tempfile.mkdtemp(prefix='rescore-eval-', dir=results))
        else:
            target = args.out.resolve()
            if target == results or results in target.parents or target in results.parents:
                raise ValueError('New evaluation output must be disjoint from source evidence')
            identifier(target.name)
            target = new_directory(target)
            evidence = target
        scratch.mkdir(parents=True, exist_ok=True)
        summary_file = results / 'summary.json'
        old = json.loads(summary_file.read_text()) if summary_file.exists() else {'trials': []}
        old_rows = {r['trial']: r for r in old['trials']}
        rows = []
        eval_id = evidence.name
        manifest = {'schema_version': 1, 'eval_id': eval_id, 'started_at': now(),
                    'source': str(results), 'source_files': before,
                    'verifier_files': tree_hashes(task / 'verifier'),
                    'reference_files': tree_hashes(task / 'reference'),
                    'protocol_id': 'legacy_rescore_explicit_verifier_hashes', 'os_isolation': 'not_validated',
                    'counts_as_new_research_evidence': False, 'in_place_copy': args.in_place_copy}
        def evaluate(env, filename):
            score = target / filename
            if args.in_place_copy and score.exists():
                (evidence / (filename + '.previous')).write_bytes(score.read_bytes())
                score.unlink()
            result = grade(task, env, score, args.timeout)
            # Keep checker's raw output; lifecycle classification is a separate record.
            if not score.exists():
                write_json(score, result)
            write_json(evidence / (filename + '.execution.json'), result)
            return result
        for i in trials:
            archive = results / f'trial_{i:02d}_environment.tar.gz'
            with tempfile.TemporaryDirectory(prefix=f'{eval_id}-trial-{i:02d}-', dir=scratch) as td:
                box = Path(td) / 'recovered'
                safe_extract(archive, box, sha256(archive))
                if not (box / 'environment/pybte/driver.py').is_file():
                    raise ValueError(f'Trial {i}: missing candidate entrypoint')
                result = evaluate(box / 'environment', f'trial_{i:02d}_score.json')
                rows.append({'trial': i, **row_from(result), 'source_score': old_rows.get(i, {}).get('score'),
                             'archive_sha256': sha256(archive), 'eval_id': eval_id})
        if args.oracle:
            manifest['oracle'] = row_from(evaluate(task / 'oracle/solution', 'oracle_score.json'))
        summary = {'eval_id': eval_id, 'source_summary': str(summary_file) if summary_file.exists() else None,
                   'trials': rows, 'solved': sum(r['gate'] for r in rows), 'new_protocol_score': None}
        if args.in_place_copy and summary_file.exists():
            (evidence / 'summary.previous.json').write_bytes(summary_file.read_bytes())
            updated = dict(old_rows)
            for row in rows:
                updated[row['trial']] = {**updated.get(row['trial'], {}), **row}
            summary['trials'] = [updated[k] for k in sorted(updated)]
            summary['solved'] = sum(bool(r.get('gate')) for r in summary['trials'] if not r.get('voided'))
            summary_file.write_text(json.dumps(summary, indent=2) + '\n')
        else:
            write_json(target / 'summary.json', summary)
        manifest.update(finished_at=now(), source_unchanged=tree_hashes(results) == before)
        if not args.in_place_copy and not manifest['source_unchanged']:
            raise ValueError('Source evidence changed during rescore')
        manifest['output_files'] = tree_hashes(target)
        write_json(evidence / 'evaluation_manifest.json', manifest)
        (evidence / 'report.md').write_text(
            f"# {eval_id}\n\nLegacy re-evaluation; source unchanged: {manifest['source_unchanged']}.\n"
            'New-protocol score: null. Compare source_score to score in summary.json.\n'
            'OS isolation: not_validated. No model requests made.\n')
        print(json.dumps(summary, indent=2))
        return 0 if all(r['gate'] for r in rows) and manifest.get('oracle', {'gate': True})['gate'] else 1
    except (OSError, ValueError, KeyError) as exc:
        if 'evidence' in locals() and evidence.is_dir():
            failure = evidence / 'failure.json'
            if not failure.exists():
                write_json(failure, {'termination': 'infrastructure_error', 'error': str(exc), 'time': now()})
        print(f'rescore: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
