""" Compare new final state * writer with the more traditional test_out.dat -> FinalState* conversion

"""

from pathlib import Path

import awkward as ak
import numpy as np
import IPython

def compare(filename: Path) -> None:

    new_output = ak.from_parquet(filename)
    reference_output = ak.from_parquet(filename.parent / (str(filename.name).replace("_00.parquet", "") + "_reference_00.parquet"))

    for field in ak.fields(new_output):
        print(f"Testing field: {field}")
        np.testing.assert_allclose(np.asarray(ak.flatten(new_output[field], axis=None)), np.asarray(ak.flatten(reference_output[field], axis=None)), atol=1e-3, rtol=5e-5)

if __name__ == "__main__":
    print("*** hadrons ***")
    compare(filename=Path("test_out_0001_final_state_hadrons_00.parquet"))
    print("*** partons ***")
    compare(filename=Path("test_out_0001_final_state_partons_00.parquet"))
