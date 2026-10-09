"""Restore genuinely absent synthetic ROOT descriptors without hiding changes."""

from pathlib import Path
import struct
import sys
import tempfile
import unittest

import uproot
from uproot.writing._cascade import RawStreamerInfo, RawTListOfStrings

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from add_bdt_score import (
    DEFAULT_BRANCH, add_score, array_digest, file_digest, key_payload_digest,
    open_root, restore_missing_streamers, streamer_schema, verify_original_content,
)
from test_add_bdt_score import fixture_arrays, write_fixture, write_model


def rewrite_streamer_list(path, *, remove=(), change_checksum=None):
    """Deliberately corrupt only a test file's actual serialized descriptor list.

    This uses the same pinned Uproot cascading writer as the restoration code.
    Removing a descriptor does not touch keys, event baskets or saved scores.
    """
    with open_root(path) as root:
        infos = [info for versions in root.file.streamers.values() for info in versions.values()]
        rules = list(root.file.streamer_rules)
    with uproot.update(path) as writable:
        node = writable.file._cascading.streamers
        streamers = []
        # Native ROOT shares class tags between descriptors. Rebuild each
        # descriptor independently so removing one does not invalidate the
        # remaining class-tag offsets. The production repair only appends.
        for info in infos:
            identity = (info.name, info.class_version)
            if identity in remove:
                continue
            if identity == change_checksum:
                info._members["fCheckSum"] ^= 1
                # Uproot preserves the original serialized prefix; its final
                # eight bytes contain the checksum and class version.
                info._serialization = info._serialization[:-8] + struct.pack(
                    ">Ii", info.member("fCheckSum"), info.class_version)
            serialized = uproot.serialization.serialize_object_any(info, info.name) + b"\x00"
            streamers.append(RawStreamerInfo(0, serialized, info.name, info.class_version))
        if rules:
            rule_list = uproot.writing.identify.to_TList(
                [uproot.writing.identify.to_TObjString(rule) for rule in rules], name="listOfRules")
            serialized = uproot.serialization.serialize_object_any(rule_list) + b"\x00"
            streamers.append(RawTListOfStrings(0, serialized))
        node._rawstreamers = streamers
        node._lookup = {(info.name, info.class_version) for info in streamers if hasattr(info, "name")}
        node._dependencies = (node.freesegments, node.key, node.header, *streamers)
        node.write(writable.file.sink)
        writable.file.sink.flush()


def schemas(path):
    with open_root(path) as root:
        return {(name, version): streamer_schema(info.tojson())
                for name, versions in root.file.streamers.items() for version, info in versions.items()}


def content(path):
    """Record every event branch and every non-Events compressed key payload."""
    with open_root(path) as root:
        return {
            "branches": {name: array_digest(values) for name, values in
                         root["Events"].arrays(library="ak", how=dict).items()},
            "keys": {name: key_payload_digest(root, name)
                     for name in root.keys(cycle=False) if name != "Events"},
            "types": {name: root["Events"][name].typename for name in root["Events"].keys()},
        }


class OriginalStreamerRestorationTest(unittest.TestCase):
    def setUp(self):
        import ROOT
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.base = Path(directory.name)
        self.source = self.base / "source.root"
        self.scored = self.base / "scored.root"
        self.model_path = self.base / "model.json"
        write_fixture(self.source, fixture_arrays())
        # A histogram is scientifically harmless synthetic metadata. Its
        # descriptor is independent of Events and gives a real missing class
        # version without experiment dictionaries or production payloads.
        root = ROOT.TFile(str(self.source), "UPDATE")
        histogram = ROOT.TH1D("SyntheticStreamerFixture", "", 4, 0, 4)
        histogram.Fill(1.5, -2.75)
        histogram.Write()
        root.Close()
        self.model = write_model(self.model_path)
        self.report = add_score(self.source, self.scored, self.model_path,
                                sample_kind="simulation", chunk_events=1)
        self.assertTrue(self.report["complete"])
        original_schemas = schemas(self.source)
        self.missing = next(identity for identity in original_schemas if identity[0] == "TH1D")
        self.changed = next(identity for identity in original_schemas if identity[0] == "TNamed")

    def verify_scored_file(self):
        return verify_original_content(self.source, self.scored, DEFAULT_BRANCH,
                                       DEFAULT_BRANCH + "Metadata", self.model, 1)

    def test_real_missing_class_version_is_restored_with_identical_schema_and_content(self):
        before = content(self.scored)
        source_hash = file_digest(self.source)
        original_schemas = schemas(self.source)
        rewrite_streamer_list(self.scored, remove={self.missing})
        self.assertNotIn(self.missing, schemas(self.scored))
        self.assertEqual(content(self.scored), before)
        with self.assertRaisesRegex(ValueError, "streamer schema"):
            self.verify_scored_file()

        restored = restore_missing_streamers(self.source, self.scored)
        self.assertEqual(restored, [{"name": self.missing[0], "class_version": self.missing[1]}])
        for identity, schema in original_schemas.items():
            self.assertEqual(schemas(self.scored)[identity], schema)
        self.assertEqual(content(self.scored), before)
        self.assertEqual(file_digest(self.source), source_hash)
        verification = self.verify_scored_file()
        self.assertTrue(verification["streamer_schemas_equal"])
        self.assertTrue(verification["all_original_content_equal"])
        self.assertTrue(verification["all_scores_recomputed_equal"])
        self.assertEqual(verification["original_streamer_schemas"], len(original_schemas))

    def test_changed_existing_checksum_is_rejected_without_modifying_file(self):
        rewrite_streamer_list(self.scored, change_checksum=self.changed)
        self.assertNotEqual(schemas(self.source)[self.changed], schemas(self.scored)[self.changed])
        before = file_digest(self.scored)
        with self.assertRaisesRegex(ValueError, "changed an original ROOT streamer schema"):
            restore_missing_streamers(self.source, self.scored)
        self.assertEqual(file_digest(self.scored), before)
        with self.assertRaisesRegex(ValueError, "streamer schema"):
            self.verify_scored_file()

    def test_schema_change_prevents_restoring_any_other_missing_descriptor(self):
        rewrite_streamer_list(self.scored, remove={self.missing}, change_checksum=self.changed)
        before = file_digest(self.scored)
        with self.assertRaisesRegex(ValueError, "streamer schema"):
            restore_missing_streamers(self.source, self.scored)
        self.assertEqual(file_digest(self.scored), before)
        self.assertNotIn(self.missing, schemas(self.scored))

    def test_no_missing_descriptor_is_a_byte_identical_no_op(self):
        before = file_digest(self.scored)
        self.assertEqual(restore_missing_streamers(self.source, self.scored), [])
        self.assertEqual(file_digest(self.scored), before)
        self.assertTrue(self.verify_scored_file()["all_scores_recomputed_equal"])

    def test_restoration_is_idempotent_and_saved_double_scores_are_unchanged(self):
        with open_root(self.scored) as root:
            saved = root["Events"][DEFAULT_BRANCH].array(library="ak")
            score_digest = array_digest(saved)
            self.assertEqual(root["Events"][DEFAULT_BRANCH].typename, "double[]")
        rewrite_streamer_list(self.scored, remove={self.missing})
        restore_missing_streamers(self.source, self.scored)
        restored_hash = file_digest(self.scored)
        self.assertEqual(restore_missing_streamers(self.source, self.scored), [])
        self.assertEqual(file_digest(self.scored), restored_hash)
        with open_root(self.scored) as root:
            self.assertEqual(array_digest(root["Events"][DEFAULT_BRANCH].array(library="ak")), score_digest)
        self.assertTrue(self.verify_scored_file()["all_scores_recomputed_equal"])

    def test_missing_descriptor_restores_without_root_directory_key_bootstrap(self):
        # Real hadd inputs can have an obsolete/unreadable TKey at fBEGIN,
        # while ROOT's directory header and every actual object remain valid.
        # The high-level update writer wrongly bootstraps that obsolete key.
        rewrite_streamer_list(self.scored, remove={self.missing})
        original_content = content(self.scored)
        for path in (self.source, self.scored):
            with open_root(path) as root:
                begin = root.file.fBEGIN
            with path.open("r+b") as stream:
                stream.seek(begin + 4)
                stream.write(b"\x00\x00")  # Unsupported TKey version, not event data.
        self.assertEqual(content(self.scored), original_content)
        with self.assertRaises(ValueError):
            with uproot.update(self.scored):
                pass
        source_hash = file_digest(self.source)
        with self.scored.open("rb") as stream:
            stream.seek(100)
            directory_prefix = stream.read(64)

        restored = restore_missing_streamers(self.source, self.scored)
        self.assertEqual(restored, [{"name": self.missing[0], "class_version": self.missing[1]}])
        self.assertEqual(file_digest(self.source), source_hash)
        self.assertEqual(content(self.scored), original_content)
        with self.scored.open("rb") as stream:
            stream.seek(100)
            self.assertEqual(stream.read(64), directory_prefix)
        verification = self.verify_scored_file()
        self.assertTrue(verification["streamer_schemas_equal"])
        self.assertTrue(verification["all_original_content_equal"])
        self.assertTrue(verification["all_scores_recomputed_equal"])


if __name__ == "__main__":
    unittest.main()
