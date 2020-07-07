from typing import Dict
import unittest

from Bio.SeqFeature import FeatureLocation

from antismash.common import path
from antismash.common.secmet.test.helpers import DummyCDS
from antismash.common.test.helpers import DummyRecord
from antismash.detection.dmz_cds import dmz_cds


def set_dummy_with_cdses(cds_features: Dict[str, FeatureLocation]) -> DummyRecord:
    cds_dummies = []
    for locus_tag, cds_location in cds_features.items():
        cds_dummy = DummyCDS(locus_tag=locus_tag, start=cds_location.start, end=cds_location.end)
        cds_dummies.append(cds_dummy)
    return DummyRecord(features=cds_dummies)

# TODO: fix test to match

class TestBlastMarker(unittest.TestCase):

    def test_get_cds_lengths(self):
        fake_cds = {"fake1": FeatureLocation(5, 13), "fake2": FeatureLocation(0, 4)}
        fake_record = set_dummy_with_cdses(fake_cds)
        fake_lengths = dmz_cds.get_cds_lengths(fake_record)
        assert "fake1" in fake_lengths and "fake2" in fake_lengths
        assert fake_lengths["fake1"] == 8
        assert fake_lengths["fake2"] == 4

    def test_get_top_splitter(self):
        blast_result = """breakygene	RBAM_RS09105	98.74	398	5	0	1	398	1	398	0.0	815"""
        dummy_splitter = {"RBAM_RS09105": FeatureLocation(1, 398)}
        current_record = set_dummy_with_cdses(dummy_splitter)
        top_splitter_default = dmz_cds.get_top_splitter_name(blast_result, current_record)
        assert top_splitter_default == "RBAM_RS09105"
        # cutoff levels
        top_splitter_higher_id = dmz_cds.get_top_splitter_name(blast_result, current_record, min_perc_identity=99)
        assert not top_splitter_higher_id
        # coverage cutoff
        fake_result_low_coverage = """breakygene	RBAM_RS09105	98.74	20	5	0	1	20	1	20	0.0	815"""
        top_splitter_too_low = dmz_cds.get_top_splitter_name(fake_result_low_coverage, current_record)
        assert not top_splitter_too_low

    def test_get_top_splitter_fail(self):
        broken_result = """breakygene	RBAM_RS09105	98.74	398	5	0	1	398	1	398	0.0"""
        dummy_splitter = {"RBAM_RS09105": FeatureLocation(1, 398)}
        current_record = set_dummy_with_cdses(dummy_splitter)
        with self.assertRaisesRegex(ValueError, "Malformed blast hit"):
            dmz_cds.get_top_splitter_name(broken_result, current_record)
        result_not_in_record = """breakygene	RBAM_RS09106	98.74	398	5	0	1	398	1	398	0.0	815"""
        with self.assertRaisesRegex(KeyError, "Blast hit RBAM_RS09106 not found in genome"):
            dmz_cds.get_top_splitter_name(result_not_in_record, current_record)

    def test_splitter_locations(self):
        feature_1 = DummyCDS(start=5, end=13)
        feature_2 = DummyCDS(start=105, end=113)
        # split on a single CDS
        location_1 = dmz_cds.get_split_locations([feature_1])
        assert location_1 == [5]
        # split on multiple CDSs
        location_1_2 = dmz_cds.get_split_locations([feature_1, feature_2])
        assert location_1_2 == [5, 105]
        location_some_false = dmz_cds.get_split_locations([feature_1, False])
        assert location_some_false == [5]


    def test_splitter_missing_hits(self):
        location_no_hits = dmz_cds.get_split_locations([])
        assert not location_no_hits
