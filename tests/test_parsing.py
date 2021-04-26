""" Tests for v2 parser.

.. codeauthor:: Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL
"""

from pathlib import Path

import awkward as ak
import pytest

from jetscape_analysis.analysis.reader import parse_ascii

@pytest.mark.parametrize(
    "header_version",
    [1, 2],
    ids=["Header v1", "Header v2"]
)
@pytest.mark.parametrize(
    "events_per_chunk",
    [
        5, 16, 50, 5000,
    ], ids=["Multiple, divisible: 5", "Multiple, indivisible: 16", "Equal: 50", "Larger: 5000"]
)
def test_parsing(header_version: int, events_per_chunk: int) -> None:
    here = Path(__file__).parent
    input_filename = here / "parsing" / f"final_state_hadrons_header_v{header_version}.dat"

    # Remap column names that have been renamed.
    rename_columns = {
        "hydro_event_id": "event_ID",
    }

    for i, arrays in enumerate(parse_ascii.read(filename=input_filename, events_per_chunk=events_per_chunk, parser="pandas")):
        # Reduce to the minimum required data.
        arrays = parse_ascii.full_events_to_only_necessary_columns_E_px_py_pz(arrays)

        # We limit the depth of the zip to ensure that we can write the parquet successfully.
        # (parquet can't handle lists of structs at the moment). Later, we'll recreate this
        # structure fully zipped together.
        ak.zip(
            dict(zip(ak.fields(arrays), ak.unzip(arrays))),
            depth_limit = 1
        )

        # Create the reference arrays by checking out the parser v1 (e477e0277fa560f9aba82310c02da8177e61c9e4), setting
        # the chunk size in skim_ascii, and then calling:
        # $ python jetscape_analysis/analysis/reader/skim_ascii.py -i tests/parsing/final_state_hadrons_header_v1.dat -o tests/parsing/events_per_chunk_50/parser_v1_header_v1/test.parquet
        # NOTE: The final state hadron files won't exist when you check out that branch, so
        #       it's best to copy them for your existing branch.
        reference_arrays = ak.from_parquet(
            Path(f"{here}/parsing/events_per_chunk_{events_per_chunk}/parser_v1_header_v1/test_{i:02}.parquet")
        )
        # There are more fields in v2 than in the reference arrays (v1), so only take those
        # that are present in reference for comparison.
        # NOTE: We have to compare the fields one-by-one because the shapes of the fields
        #       are different, and apparently don't broadcast nicely with `__eq__`
        for field in ak.fields(reference_arrays):
            new_field = rename_columns.get(field, field)
            assert ak.all(reference_arrays[field] == arrays[new_field])

        # Check for cross section if header v2
        if header_version == 2:
            assert "cross_section" in ak.fields(arrays)
            assert "cross_section_error" in ak.fields(arrays)

