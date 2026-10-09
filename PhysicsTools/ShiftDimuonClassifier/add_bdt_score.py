#!/usr/bin/env python3
"""Add a provisional Float64 score to every existing dimuon, without a cut.

The input Nano file is copied byte-for-byte before a new Events branch is
appended. Original branch values and non-Events key payloads are verified.
Only PREDICTOR_BRANCHES enter scoring; truth and weights are copied unchanged.
This command currently accepts explicitly declared simulation inputs only.
"""
import argparse
from array import array
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time

import awkward as ak
import numpy as np
import uproot

from features import PREDICTOR_BRANCHES, extract_reco
from portable_inference import file_digest, load_model, predict


DEFAULT_BRANCH = 'ShiftDimuonVertex_bdtScore'


def open_root(path):
    return uproot.open(path, handler=uproot.source.file.MultithreadedFileSource, num_workers=1)


def scores_for_arrays(model, arrays):
    """Keep every candidate in exact event/pair order; mass is not read."""
    reco = extract_reco(arrays, include_nuisance=False)
    scores = predict(model, reco['X'])
    counts = np.asarray(arrays['nShiftDimuonVertex'], dtype=np.int64)
    if (counts < 0).any() or int(counts.sum()) != len(scores):
        raise ValueError('Score count differs from retained dimuons')
    offsets = np.concatenate(([0], np.cumsum(counts)))
    expected_events = np.repeat(np.arange(len(counts)), counts)
    expected_pairs = np.concatenate([np.arange(n) for n in counts]) if len(scores) else np.empty(0, dtype=int)
    np.testing.assert_array_equal(reco['event_index'], expected_events)
    np.testing.assert_array_equal(reco['pair_index'], expected_pairs)
    if not np.isfinite(scores).all() or ((scores < 0) | (scores > 1)).any():
        raise ValueError('Invalid BDT scores')
    return [scores[offsets[i]:offsets[i+1]] for i in range(len(counts))], int((~np.isfinite(reco['X'])).sum())


def array_digest(values):
    """Preserve dtype, jagged offsets, values and NaN bit patterns."""
    form, length, buffers = ak.to_buffers(ak.Array(values))
    result = hashlib.sha256((form.to_json() + ':' + str(length)).encode())
    for key, value in sorted(buffers.items()):
        result.update(key.encode())
        result.update(np.ascontiguousarray(value).tobytes())
    return result.hexdigest()


def key_payload_digest(root, name):
    key = root.key(name)
    data = root.file.source.chunk(key.fSeekKey + key.fKeylen, key.fSeekKey + key.fNbytes).raw_data
    return hashlib.sha256(data).hexdigest()


def streamer_schema(value):
    """ROOT rewrites TObject bookkeeping bits; keep all type/schema fields."""
    if isinstance(value, dict):
        return {key: streamer_schema(item) for key,item in value.items() if key != 'fBits'}
    if isinstance(value, list):
        return [streamer_schema(item) for item in value]
    return value


def verify_original_content(source, output, branch, metadata_key, model, chunk_events):
    """Check the complete original schema/content, then recompute every score."""
    with open_root(source) as original, open_root(output) as scored:
        streamer_count = 0
        for name, versions in original.file.streamers.items():
            for version, info in versions.items():
                replacement = scored.file.streamers.get(name, {}).get(version)
                if replacement is None or streamer_schema(info.tojson()) != streamer_schema(replacement.tojson()):
                    raise ValueError('Scoring changed an original ROOT streamer schema: ' + name)
                streamer_count += 1
        if set(scored.keys(cycle=False)) != set(original.keys(cycle=False)) | {metadata_key}:
            raise ValueError('Scoring changed the original ROOT key inventory')
        for name in original.keys(cycle=False):
            if name != 'Events' and key_payload_digest(original, name) != key_payload_digest(scored, name):
                raise ValueError('Scoring changed a non-Events ROOT key: ' + name)
        before, after = original['Events'], scored['Events']
        if before.num_entries != after.num_entries or set(after.keys()) != set(before.keys()) | {branch}:
            raise ValueError('Scoring changed the event count or original branch inventory')
        if after[branch].typename != 'double[]':
            raise ValueError('BDT score must be a variable-length Float64 branch')
        for name in before.keys():
            if before[name].typename != after[name].typename:
                raise ValueError('Scoring changed an original branch type: ' + name)
        for start in range(0, before.num_entries, chunk_events):
            stop = min(before.num_entries, start + chunk_events)
            a = before.arrays(entry_start=start, entry_stop=stop, library='ak', how=dict)
            b = after.arrays(list(before.keys()), entry_start=start, entry_stop=stop, library='ak', how=dict)
            for name in a:
                if array_digest(a[name]) != array_digest(b[name]):
                    raise ValueError('Scoring changed original values: ' + name)
            expected, _ = scores_for_arrays(model, {name: a[name] for name in PREDICTOR_BRANCHES})
            saved = after[branch].array(entry_start=start, entry_stop=stop, library='ak')
            for index, values in enumerate(expected):
                np.testing.assert_array_equal(np.asarray(saved[index]), values)
        return dict(original_keys=len(original.keys(cycle=False)),
                    original_event_branches=len(before.keys()), events=before.num_entries,
                    original_streamer_schemas=streamer_count, streamer_schemas_equal=True,
                    excluded_streamer_bookkeeping_fields=['fBits'],
                    all_original_content_equal=True, all_scores_recomputed_equal=True,
                    score_type='double[]')


def add_score(source, output, model_path, *, sample_kind, branch=DEFAULT_BRANCH,
              chunk_events=1000, receipt_path=None):
    if sample_kind != 'simulation':
        raise ValueError('This production is explicitly simulation only')
    if type(chunk_events) is not int or chunk_events < 1:
        raise ValueError('Chunk event count must be positive')
    if not re.fullmatch(r'ShiftDimuonVertex_[A-Za-z][A-Za-z0-9_]*', branch):
        raise ValueError('Invalid dimuon score branch name')
    source, output, model_path = map(lambda p: Path(p).resolve(), (source, output, model_path))
    receipt_path = Path(receipt_path).resolve() if receipt_path else output.with_suffix('.bdt.json')
    partial = output.with_name(output.name + '.incomplete')
    if source == output or source == partial or source == receipt_path:
        raise ValueError('Scored output and receipt must be distinct from the input')
    if any(p.exists() for p in (output, partial, receipt_path)):
        raise FileExistsError('Refuse to overwrite existing scored output, incomplete file or receipt')
    model = load_model(model_path)
    metadata_key = branch + 'Metadata'
    started = time.monotonic()
    source_sha256, model_sha256 = file_digest(source), file_digest(model_path)
    with open_root(source) as original:
        if 'Events' not in original or original['Events'].classname != 'TTree':
            raise ValueError('Input requires a Nano Events TTree')
        tree = original['Events']
        if branch in tree or metadata_key in original:
            raise ValueError('Input already contains this BDT score or its metadata')
        missing = set(PREDICTOR_BRANCHES) - set(tree.keys())
        if missing:
            raise ValueError('Missing reconstructed predictors: ' + ', '.join(sorted(missing)))
        counts = np.asarray(tree['nShiftDimuonVertex'].array(library='np'), dtype=np.int64)
        if (counts < 0).any():
            raise ValueError('Negative retained pair count')
        n_events, max_pairs = tree.num_entries, max(1, int(counts.max()) if len(counts) else 1)
    output.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    report = dict(schema='shift-dimuon-bdt-score-receipt-v1', complete=False,
                  source=str(source), source_sha256=source_sha256,
                  output=str(output), model=str(model_path), model_sha256=model_sha256,
                  model_name=model['provenance'].get('model_name'), sample_kind=sample_kind,
                  branch=branch, score_storage='Float64', selection_applied=False,
                  physics_ready=False, normalization_transfer_validated=False,
                  events=n_events, retained_pairs=int(counts.sum()), zero_pair_events=int((counts == 0).sum()),
                  threshold_metadata=model['threshold_metadata'],
                  score_interpretation=model['score_interpretation'], nonfinite_feature_elements=0,
                  scoring_branch_allowlist=list(PREDICTOR_BRANCHES), scoring_requires_gen=False,
                  scorer_source_sha256=file_digest(Path(__file__)))
    with receipt_path.open('x') as stream:
        stream.write(json.dumps(report, indent=2) + '\n')
    root = None
    try:
        with source.open('rb') as old, partial.open('xb') as new:
            shutil.copyfileobj(old, new, length=1024 * 1024)
        if file_digest(partial) != source_sha256:
            raise ValueError('Input changed while copying')
        import ROOT
        root = ROOT.TFile(str(partial), 'UPDATE')
        if not root or root.IsZombie() or root.TestBit(ROOT.TFile.kRecovered):
            raise ValueError('Cannot update the distinct Nano copy')
        tree = root.Get('Events')
        counter = tree.GetBranch('nShiftDimuonVertex')
        buffer = array('d', [0.] * max_pairs)
        new_branch = tree.Branch(branch, buffer, branch + '[nShiftDimuonVertex]/D')
        global_event = 0
        with open_root(source) as original:
            for arrays in original['Events'].iterate(list(PREDICTOR_BRANCHES), step_size=chunk_events,
                                                    library='ak', how=dict):
                event_scores, imputed = scores_for_arrays(model, arrays)
                report['nonfinite_feature_elements'] += imputed
                for values in event_scores:
                    if counter.GetEntry(global_event) < 0 or int(counter.GetLeaf('nShiftDimuonVertex').GetValue()) != len(values):
                        raise ValueError('ROOT counter differs from the scored pair length')
                    for index, value in enumerate(values):
                        buffer[index] = float(value)
                    if new_branch.Fill() < 0:
                        raise ValueError('ROOT failed to fill a score branch entry')
                    global_event += 1
        if global_event != n_events or new_branch.GetEntries() != n_events or tree.GetEntries() != n_events:
            raise ValueError('Score writing changed the event count')
        root.cd()
        tree.Write('', ROOT.TObject.kOverwrite)
        metadata = {key: report[key] for key in
                    ('schema','source','source_sha256','model_name','model_sha256','sample_kind',
                     'branch','score_storage','selection_applied','physics_ready',
                     'normalization_transfer_validated','events','retained_pairs',
                     'threshold_metadata','score_interpretation','scoring_requires_gen',
                     'scorer_source_sha256')}
        ROOT.TObjString(json.dumps(metadata, sort_keys=True, allow_nan=False)).Write(metadata_key)
        root.Close()
        root = None
        report['verification'] = verify_original_content(source, partial, branch, metadata_key, model, chunk_events)
        if file_digest(source) != source_sha256 or file_digest(model_path) != model_sha256:
            raise ValueError('Input or frozen model changed during scoring')
        if output.exists():
            raise FileExistsError('Output appeared during scoring; do not overwrite')
        os.rename(partial, output)
        report.update(complete=True, output_sha256=file_digest(output), elapsed_seconds=time.monotonic()-started)
    except Exception as error:
        report.update(error=repr(error), retained_incomplete_file=str(partial))
        raise
    finally:
        if root:
            root.Close()
        receipt_path.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--sample-kind', choices=('simulation',), required=True)
    parser.add_argument('--score-branch', default=DEFAULT_BRANCH)
    parser.add_argument('--chunk-events', type=int, default=1000)
    parser.add_argument('--receipt', type=Path)
    args = parser.parse_args()
    report = add_score(args.input, args.output, args.model, sample_kind=args.sample_kind,
                       branch=args.score_branch, chunk_events=args.chunk_events, receipt_path=args.receipt)
    print(json.dumps({key: report[key] for key in ('complete','events','retained_pairs','output','model_sha256',
                                                  'selection_applied','elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    main()
