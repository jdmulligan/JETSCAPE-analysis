#!/usr/bin/env python3

""" Parse JETSCAPE ascii input files in chunks.

.. codeauthor:: Raymond Ehlers
"""

import logging
import os
import re
import typing
from pathlib import Path
from typing import Any, Generator, Iterator, Iterable, List, Optional, Sequence, Union, Tuple, TextIO
from typing_extensions import Literal

import awkward as ak
import attr
import numpy as np


logger = logging.getLogger(__name__)

@attr.s
class CrossSection:
    value: float = attr.ib()
    error: float = attr.ib()


def _retrieve_list_line_of_file(f: TextIO, chunk_size: int = 100) -> str:
    """ Retrieve the last line of the file.

    From: https://stackoverflow.com/a/7167316/12907985

    Args:
        f: File-like object.
        chunk_size: Size of step to read backwards into the file. Default: 100.

    Returns:
        Last line of file, assuming it's found.
    """
    last_line = ""
    while True:
        # We grab chunks from the end of the file towards the beginning until we
        # get a new line
        f.seek(-len(last_line) - chunk_size, os.SEEK_END)
        chunk = f.read(chunk_size)

        if not chunk:
            # The whole file is one big line
            return last_line

        if not last_line and chunk.endswith('\n'):
            # Ignore the trailing newline at the end of the file (but include it
            # in the output).
            last_line = '\n'
            chunk = chunk[:-1]

        nl_pos = chunk.rfind('\n')
        # What's being searched for will have to be modified if you are searching
        # files with non-unix line endings.

        last_line = chunk[nl_pos + 1:] + last_line

        if nl_pos == -1:
            # The whole chunk is part of the last line.
            continue

        return last_line


def _extract_x_sec_and_error(f: TextIO, chunk_size: int = 100) -> Optional[CrossSection]:
    """ Extract cross section and error from EOF.

    Args:
        f: File-like object.
        chunk_size: Size of step to read backwards into the file. Default: 100.

    Returns:
        Cross section and error, if found.
    """
    # Retrieve the last line of the file.
    last_line = _retrieve_list_line_of_file(f=f, chunk_size=chunk_size)
    # Move the file back to the start to reset it for later use.
    f.seek(0)

    if last_line.startswith("# sigmaGen"):
        # The latter two arguments are dummy arguments. This way, we can use one parser for comments.
        return _parse_line_for_header_with_weight(line=last_line, n_events=0, events_per_chunk=100)

    return None


def _parse_line_for_header(line: str, n_events: int, events_per_chunk: int) -> Tuple[bool, Optional[Any]]:
    """ Parse line as appropriate.

    If it's just a standard particle line, we just pass it on. However, if it's a header, we parse
    the header, as well as perform both an event count check. This way, we can split files without
    having to worry about getting only part of an event (which would happen if we just trivially
    split on lines).

    Note:
        We don't even have to yield the line back because we don't ever modify it.

    Args:
        line: Line to be parsed.
        n_events: Number of events processed so far.
        events_per_chunk: Number of events per chunk.
    Returns:
        Whether we've reached the desired number of events and should stop this block, any information parsed from the header.
        The header info is None if it's not a header line.
    """
    time_to_stop = False
    header_info = None
    if line.startswith("#"):
        # We've found a header line.

        # We need to be able to chunk the files into something smaller and more manageable.
        # Therefore, when we hit the target, we provide a signal that it's time to end the block.
        # This is a bit awkward because we don't know that an event has ended until we get to the
        # next event. However, we check for exact agreement because we don't increment the number
        # of events until after we find a header. Basically, it all works out here.
        #if n_events > 0 and n_events % events_per_chunk == 0:
        if n_events == events_per_chunk:
            logger.debug("Hit end of chunk - time to stop!")
            time_to_stop = True

        # Parse the header string.
        # As of 9 October 2020, the formatting isn't really right. This should be fixed in JS.
        # Due to this formatting issue:
        # - We ignore all of the column names.
        # - We only parse the numbers:
        #   1. Event plane angle. float, potentially in scientific notation.
        #   3. (hydro?) event number. int
        #   3. Number of particles. int. (This wasn't clear, originally)
        #
        # It's unclear where the hydro ID is stored.
        #
        # For now, we don't construct any objects to contain the information because
        # it's not worth the computing time - we're not really using this information...
        header_values = line.split()
        # Validation
        if len(header_values) != 9:
            raise ValueError(f"Parsing of header failed: {header_values}")
        # The second value has to be clever because the file is improperly spaced for parsing...
        # It's of the form: "Event5ID", where all we care about is the event number (in this case, 5)
        header_info = [float(header_values[1]), int(header_values[2][5:-2]), int(header_values[3])]

    return time_to_stop, header_info

def _parse_line_for_header_with_weight(line: str, n_events: int, events_per_chunk: int) -> Tuple[bool, Optional[Any]]:
    """ Parse line as appropriate.

    If it's just a standard particle line, we just pass it on. However, if it's a header, we parse
    the header, as well as perform both an event count check. This way, we can split files without
    having to worry about getting only part of an event (which would happen if we just trivially
    split on lines).

    Note:
        We don't even have to yield the line back because we don't ever modify it.

    Args:
        line: Line to be parsed.
        n_events: Number of events processed so far.
        events_per_chunk: Number of events per chunk.
    Returns:
        Whether we've reached the desired number of events and should stop this block, any information parsed from the header.
        The header info is None if it's not a header line.
    """
    time_to_stop = False
    header_info = None
    if line.startswith("#"):
        # We've found a comment line. This is either a header, or a summary of the cross section and error.

        # We need to be able to chunk the files into something smaller and more manageable.
        # Therefore, when we hit the target, we provide a signal that it's time to end the block.
        # This is a bit awkward because we don't know that an event has ended until we get to the
        # next event. However, we check for exact agreement because we don't increment the number
        # of events until after we find a header. Basically, it all works out here.
        #if n_events > 0 and n_events % events_per_chunk == 0:
        if n_events == events_per_chunk:
            logger.debug("Hit end of chunk - time to stop!")
            time_to_stop = True

        # Parse the string.
        # For now, we don't construct any objects to organize the information because I don't think
        # it's worth the computing time to construct the object. This may change as we add more info.
        values = line.split("\t")
        # Compare by length first so we can short circuit immediately if it doesn't match, which should
        # save some string comparisons.
        if len(values) == 5 and values[1] == "sigmaGen":
            ###################################
            # Cross section info specification:
            ###################################
            # The cross section info is formatted as follows, with each entry separated by a `\t` character:
            # # sigmaGen 182.423 sigmaErr 11.234
            # 0 1        2       3        4
            #
            header_info = [
                CrossSection(value=float(values[2]), error=float(values[4]))
            ]

            # This also means that we've hit the end of the file. Let's note it as such.
            logger.debug("Hit end of file - time to stop!")
            time_to_stop = True

        elif len(values) == 19 and values[1] == "Event":
            #######################
            # Header specification:
            #######################
            # As of 20 April 2021, the formatting of the header has been improved.
            # This function was developed to parse it.
            # The header is defined as follows, with each entry separated by a `\t` character:
            # `# Event 1 weight 0.129547 EPangle 0.0116446 N_hadrons 236 | N  pid status E Px Py Pz Eta Phi`
            #  0 1     2 3      4        5       6         7         8   9 10 11  12    13 14 15 16 17  18
            #
            # NOTE: Everything after the "|" is just documentation for the particle entries stored below.
            #
            header_info = [
                int(values[2]),     # Event number
                float(values[4]),   # Event weight
                float(values[6]),   # EP angle
                int(values[8]),     # Number of particles
            ]
        else:
            raise ValueError(f"Parsing of comment line failed: {values}")

    return time_to_stop, header_info


def read_events_in_chunks(filename: Union[Path, str], events_per_chunk: int = int(1e5)) -> Iterator[Tuple[Iterator[str], List[int], List[Any], List[bool]]]:
    """ Read events in chunks from stored JETSCAPE ASCII files.

    Users are encouraged to use `read(...)`.

   This provides the underlying implementation, and in principle could be used directly. However,
   it's missing many useful features that are implemented elsewhere.

    Args:
        filename: Path to the file.
        events_per_chunk: Number of events to store in each chunk. Default: 1e5.
    Returns:
        Chunks generator. When this generator is consumed, it will generate lines from the file
            until it hits the number of events mark.
    """
    # Validation
    filename = Path(filename)

    with open(filename, "r") as f:
        # Setup
        # This is just for convenience.
        return_count = 0
        # This is used to pass a header to the next chunk. This is necessary because we don't know an event
        # is over until we already get the header for the next event. We could keep that line and reparse,
        # but there's no need to parse a file twice.
        keep_header_for_next_chunk = None
        # If we run out of events, we want to notify the user so they can avoid an assertion failure
        # if they check that the expected number of events were returned.
        # NOTE: We need a mutable object so we can modify it in the closure, but the changes will be
        #       passed back through the other return values. So we use a list. This is a dirty, dirty hack...
        reached_end_of_file = []

        # Define an iterator so we can increment it in different locations in the code.
        # Fine to use if it the entire file fits in memory.
        #read_lines = iter(f.readlines())
        # Use this if the file doesn't fit in memory (fairly likely for these type of files)
        read_lines = iter(f)

        # TODO: I think it may be possible to simplify this substantially...?

        for line in read_lines:
            # Setup
            # We keep track of the location of where to split each event. That way, we can come back later
            # and split the 2D numpy array into an awkward array with a jagged structure.
            event_split_index: List[int] = []
            # Store the event header info, to be returned alongside the particles and event split index.
            event_header_info = []
            # If we've kept around a header from a previous chunk, then store that for this iteration.
            if keep_header_for_next_chunk:
                event_header_info.append(keep_header_for_next_chunk)
                # Now that we've stored it, reset it to ensure that it doesn't cause problems for future iterations.
                keep_header_for_next_chunk = None

            # If we're run out of events in the parsing, then we should be done and we should never
            # make it back here.
            if reached_end_of_file:
                raise ValueError("Reached end of file was set, but we're iterating again. This shouldn't occur.")

            def _inner(line: str, kept_header: bool) -> Iterator[str]:
                """ Closure to generate a chunk of events.

                The closure ensures access to the same generator used to access the file.

                This closure is inspired by: https://stackoverflow.com/a/24527424/12907985 .

                Args:
                    kept_header: If True, it means that we kept the header from a previous chunk.
                Returns:
                    Generator yielding the lines of the file.
                """
                # Offset to account for the header lines that we have parsed, and therefore how much
                # we must offset to properly align particles with events. By default, we expect to have
                # one, which accounts for the header that we'll parse over at the end after finishing an event.
                header_offset = 1

                # Anything that is returned from this function will be consumed by np.loadtxt, so we can't
                # directly return any value. Instead, we have to make this nonlocal so we can set it here,
                # and the result will be accessible outside during the next chunk.
                nonlocal reached_end_of_file
                nonlocal keep_header_for_next_chunk
                # If we already have a header, then we already have an event, so we need to increment immediately.
                # NOTE: Together with storing with handling the header in the first line a few lines below, we're
                #       effectively 1-indexing n_events.
                n_events = 0
                n_particles = 0
                if kept_header:
                    n_events += 1
                    # Always reset n particles on a new event (doesn't matter here, but done for constituency).
                    n_particles = 0

                # Handle the first line from the generator.
                _, header_info = _parse_line_for_header(line, n_events, events_per_chunk=events_per_chunk)
                yield line
                # We always increment after yielding.
                # Instead of defining the variable here, we account for it in the enumeration below by
                # starting at 1.

                # If we come across a header immediately (should only happen for the first line of the file),
                # we note the new event, and store the header info.
                # NOTE: Together with incrementing n_events above, we're effectively 1-indexing n_events.
                if header_info:
                    #logger.debug(f"special case header, {line}")

                    # Now account for the new header.
                    n_events += 1
                    # Always reset n particles on a new event.
                    n_particles = 0
                    # We add to the header offset because now we need to account for this initial header in
                    # addition to the later one that we'll parse over.
                    header_offset += 1
                    # Store the header info
                    event_header_info.append(header_info)
                    # NOTE: We don't record any info here for event_split_index because this is line 0, and
                    #       it would try to split on both sides of it, leading to an empty first event.
                else:
                    n_particles += 1

                # Handle additional lines
                # Start at one to account for the first land already being handled.
                for line_count, local_line in enumerate(read_lines, start=1):
                    time_to_stop, header_info = _parse_line_for_header(local_line, n_events, events_per_chunk=events_per_chunk)
                    line_count += 1

                    # A new header signals a new event. It needs some careful handling.
                    if header_info:
                        #logger.debug(f"standard case header: {local_line}")
                        # Cross check that the number of particles extracted from the header matches the number that we think we've extracted.
                        # NOTE: We have to grab the previous header, since the header_info contains the header for the next event.
                        assert event_header_info[-1][2] == n_particles, f"Expected {event_header_info[-1][2]} particles from header, but found {n_particles}"

                        # Now account for the new header.
                        n_events += 1
                        # Always reset n particles on a new event.
                        n_particles = 0
                        # If it's just some event in the middle of the chunk, then we just store the header and event split information.
                        if not time_to_stop:
                            event_header_info.append(header_info)
                            # We need to account for header lines.
                            # Since the header line will be skipped by loadtxt, we need to account for that with:
                            # We need to subtract to account for previous event header lines by subtracting the number of events so far.
                            # We also need to account for other headers that we have seen in this round of parsing.
                            event_split_index.append(line_count - len(event_split_index) - header_offset)
                        else:
                            # If we're about to end this chunk, we need to hold on to the most recent header
                            # (which signaled that we're ready to end the chunk). We'll hold onto it until we
                            # look at the next chunk.
                            keep_header_for_next_chunk = header_info
                            #print(f"header_info: {header_info}, keep_header_for_next_chunk: {keep_header_for_next_chunk}")
                    else:
                        n_particles += 1

                    # Regardless of the header status, we should always yield the line so it can be handled downstream.
                    yield local_line

                    # Finally, if it's time to end the chunk, we need to fully break the loop. We'll pick
                    # up in the next chunk from the first line of the new event (with the header info that's
                    # stored above).
                    if time_to_stop:
                        #print(f"event_split_index len: {len(event_split_index)} - {event_split_index}")
                        break

                # If we've broken out of iteration, but we it's not yet time to stop, that means our
                # iterator ended. That means we've reached the end of the file and should let the user know.
                if not time_to_stop:
                    reached_end_of_file.append(True)

            # Yield the generator for the chunk, along with useful information.
            # NOTE: When we pass these lists, they're empty. This only works because lists are mutable, and thus
            #       our changes are passed on.
            yield _inner(line=line, kept_header=len(event_header_info) > 0), event_split_index, event_header_info, reached_end_of_file

            # Keep track of what's going on. This is basically a debugging tool.
            return_count += 1
            #print(f"return_count: {return_count}")
            #print(f"keep_header_for_next_chunk: {keep_header_for_next_chunk}")

        #print(f"return_count: {return_count}")

    # If we've gotten here, that means we've finally exhausted the file. There's nothing else to do!


class FileLikeGenerator:
    """ Wrapper class to make a generator look like a file.

    Pandas requires passing a filename or a file-like object, but we handle the generator externally
    so we can find each event boundary, parse the headers, chunk, etc. Consequently, we need to make
    this generator appear as if it's a file.

    Based on https://stackoverflow.com/a/18916457/12907985

    Args:
        g: Generator to be wrapped.
    """
    def __init__(self, g: Iterator[str]):
        self.g = g

    def read(self, n: int = 0) -> Any:
        """ Read method is required by pandas. """
        try:
            return next(self.g)
        except StopIteration:
            return ''

    def __iter__(self) -> Iterator[str]:
        """ Iteration is required by pandas. """
        return self.g


def _parse_with_pandas(chunk_generator: Iterator[str]) -> np.ndarray:
    """ Parse the lines with `pandas.read_csv`

    `read_csv` uses a compiled c parser. As of 6 October 2020, it is tested to be the fastest option.

    Args:
        chunk_generator: Generator of chunks of the input file for parsing.
    Returns:
        Array of the particles.
    """
    # Delayed import so we only take the import time if necessary.
    import pandas as pd

    return pd.read_csv(
        FileLikeGenerator(chunk_generator),
        names=["particle_index", "particle_ID", "status", "E", "px", "py", "pz", "eta", "phi"],
        header=None,
        comment="#",
        sep="\s+",
        # Converting to numpy makes the dtype conversion moot.
        #dtype={
        #    "particle_index": np.int32, "particle_ID": np.int32, "status": np.int8,
        #    "E": np.float32, "px": np.float32, "py": np.float32, "pz": np.float32,
        #    "eta": np.float32, "phi": np.float32
        #},
        # We can reduce columns to save a little time reading.
        # However, it makes little difference, and makes it less general. So we disable it for now.
        #usecols=["particle_ID", "status", "E", "px", "py", "eta", "phi"],
    # NOTE: It's important that we convert to numpy before splitting. Otherwise, it will return columns names,
    #       which will break the header indexing and therefore the conversion to awkward.
    ).to_numpy()


def _parse_with_python(chunk_generator: Iterator[str]) -> np.ndarray:
    """ Parse the lines with python.

    We have this as an option because np.loadtxt is surprisingly slow.

    Args:
        chunk_generator: Generator of chunks of the input file for parsing.
    Returns:
        Array of the particles.
    """
    particles = []
    for p in chunk_generator:
        if not p.startswith("#"):
            particles.append(np.array(p.rstrip("\n").split(), dtype=np.float64))
    return np.stack(particles)


def _parse_with_numpy(chunk_generator: Iterator[str]) -> np.ndarray:
    """ Parse the lines with numpy.

    Unfortunately, this option is surprisingly, presumably because it has so many options.
    Pure python appears to be about 2x faster. So we keep this as an option for the future,
    but it is not used by default.

    Args:
        chunk_generator: Generator of chunks of the input file for parsing.
    Returns:
        Array of the particles.
    """
    return np.loadtxt(chunk_generator)


def read(filename: Union[Path, str], events_per_chunk: int, parser: str = "pandas") -> Optional[Iterator[ak.Array]]:
    """ Read a JETSCAPE ascii output file in chunks.

    This is the main user function. We read in chunks to keep the memory usage manageable.

    Note:
        We store the data in the smallest possible types that can still encompass their range.

    Args:
        filename: Filename of the ascii file.
        events_per_chunk: Number of events to provide in each chunk.
        parser: Name of the parser to use. Default: `pandas`, which uses `pandas.read_csv`. It uses
            compiled c, and seems to be the fastest available option. Other options: ["python", "numpy"].
    Returns:
        Generator of an array of events_per_chunk events.
    """
    # Validation
    filename = Path(filename)

    # Setup
    parsing_function_map = {
        "pandas": _parse_with_pandas,
        "python": _parse_with_python,
        "numpy": _parse_with_numpy,
    }
    parsing_function = parsing_function_map[parser]

    # Read the file, creating chunks of events.
    for chunk_generator, event_split_index, event_header_info, reached_end_of_file in read_events_in_chunks(filename=filename, events_per_chunk=events_per_chunk):
        # Give a notification just in case the parsing is slow...
        logger.debug("New chunk")

        # Parse the file and create the awkward event structure.
        array_with_events = ak.Array(
            np.split(
                parsing_function(chunk_generator), event_split_index
            )
        )

        # Cross check that everything is in order and was parsed correctly.
        if events_per_chunk > 0:
            # NOTE: The second check handles the case where n_events evenly divides by events_per_chunk.
            if not reached_end_of_file or len(event_split_index) == events_per_chunk - 1:
                assert len(event_split_index) == events_per_chunk - 1
                assert len(event_header_info) == events_per_chunk
            else:
                logger.warning(f"Requested {events_per_chunk} events, but only {len(event_header_info)} are available because we hit the end of the file.")
                assert len(event_split_index) < events_per_chunk - 1
                assert len(event_header_info) < events_per_chunk

        #print(len(event_split_index))
        #print(f"hadrons: {hadrons}")
        #print(f"array_with_events: {array_with_events}")
        #print(ak.type(array_with_events))
        #print(f"Event header info: {event_header_info}")
        #import IPython; IPython.embed()
        event_plane_angles = [x[0] for x in event_header_info]
        hydro_event_id = [x[1] for x in event_header_info]

        # Convert to the desired structure for our awkward array.
        yield ak.zip(
            {
                # TODO: Does the conversion add any real computation time?
                "event_plane_angle": ak.values_astype(event_plane_angles, np.float32),
                "hydro_event_id": ak.values_astype(hydro_event_id, np.uint16),
                "particle_ID": ak.values_astype(array_with_events[:, :, 1], np.int32),
                # Status is only a couple of numbers, but it's not always 0. It identifies recoils (1?) and holes (-1?)
                "status": ak.values_astype(array_with_events[:, :, 2], np.int8),
                "E": ak.values_astype(array_with_events[:, :, 3], np.float32),
                "px": ak.values_astype(array_with_events[:, :, 4], np.float32),
                "py": ak.values_astype(array_with_events[:, :, 5], np.float32),
                "pz": ak.values_astype(array_with_events[:, :, 6], np.float32),
                # Skip these because we're going to be working with four vectors anyway, so it shouldn't be a
                # big deal to recalculate them, especially compare to the added storage space.
                "eta": ak.values_astype(array_with_events[:, :, 7], np.float32),
                "phi": ak.values_astype(array_with_events[:, :, 8], np.float32),
            },
            depth_limit=1,
        )

    #import IPython; IPython.embed()


def full_events_to_only_necessary_columns_pt_eta_phi(arrays: ak.Array) -> ak.Array:
    return ak.zip(
        {
            "particle_ID": arrays["particle_ID"],
            "status": arrays["status"],
            "pt": np.sqrt(arrays["px"] ** 2 + arrays["py"] ** 2),
            "eta": arrays["eta"],
            "phi": arrays["phi"],
        },
    )


def full_events_to_only_necessary_columns_E_px_py_pz(arrays: ak.Array) -> ak.Array:
    return ak.zip(
        {
            "event_plane_angle": arrays["event_plane_angle"],
            "hydro_event_id": arrays["hydro_event_id"],
            "particle_ID": arrays["particle_ID"],
            "status": arrays["status"],
            "E": arrays["E"],
            "px": arrays["px"],
            "py": arrays["py"],
            "pz": arrays["pz"],
        }, depth_limit=1
    )


def parse_to_parquet(base_output_filename: Union[Path, str], store_only_necessary_columns: bool,
                     input_filename: Union[Path, str], events_per_chunk: int, parser: str = "pandas",
                     max_chunks: int = -1, compression: str = "zstd", compression_level: Optional[int] = None) -> Iterator[ak.Array]:
    """ Parse the JETSCAPE ASCII and convert it to parquet, (potentially) storing only the minimum necessary columns.

    Args:
        base_output_filename: Basic output filename. Should include the entire path.
        store_only_necessary_columns: If True, store only the necessary columns, rather than all of them.
        input_filename: Filename of the input JETSCAPE ASCII file.
        events_per_chunk: Number of events to be read per chunk.
        parser: Name of the parser. Default: "pandas".
        max_chunks: Maximum number of chunks to read. Default: -1.
        compression: Compression algorithm for parquet. Default: "zstd". Options include: ["snappy", "gzip", "ztsd"].
            "gzip" is slightly better for storage, but slower. See the compression tests and parquet docs for more.
        compression_level: Compression level for parquet. Default: `None`, which lets parquet choose the best value.
    Returns:
        None. The parsed events are stored in parquet files.
    """
    # Validation
    base_output_filename = Path(base_output_filename)
    # Setup the base output directory
    base_output_filename.parent.mkdir(parents=True, exist_ok=True)

    for i, arrays in enumerate(read(filename=input_filename, events_per_chunk=events_per_chunk, parser=parser)):
        # Reduce to the minimum required data.
        if store_only_necessary_columns:
            arrays = full_events_to_only_necessary_columns_E_px_py_pz(arrays)

        # We limit the depth of the zip to ensure that we can write the parquet successfully.
        # (parquet can't handle lists of structs at the moment). Later, we'll recreate this
        # structure fully zipped together.
        ak.zip(
            dict(zip(ak.fields(arrays), ak.unzip(arrays))),
            depth_limit = 1
        )

        # Parquet with zlib seems to do about the same as ascii tar.gz when we drop unneeded columns.
        # And it should load much faster!
        if events_per_chunk > 0:
            suffix = base_output_filename.suffix
            output_filename = (base_output_filename.parent / f"{base_output_filename.stem}_{i:02}").with_suffix(suffix)
        else:
            output_filename = base_output_filename
        ak.to_parquet(
            arrays, output_filename,
            compression=compression, compression_level=compression_level,
            explode_records=False,
        )

        # Break now so we don't have to read the next chunk.
        if (i + 1) == max_chunks:
            break


if __name__ == "__main__":
    #read(filename="final_state_hadrons.dat", events_per_chunk=-1, base_output_filename="skim/jetscape.parquet")
    for pt_hat_range in ["7_9", "20_25", "50_55", "100_110", "250_260", "500_550", "900_1000"]:
        print(f"Processing pt hat range: {pt_hat_range}")
        directory_name = "OutputFile_Type5_qhatA10_B0_5020_PbPb_0-10_0.30_2.0_1"
        filename = f"JetscapeHadronListBin{pt_hat_range}"
        parse_to_parquet(
            base_output_filename=f"skim/{filename}.parquet",
            store_only_necessary_columns=True,
            input_filename=f"/alf/data/rehlers/jetscape/osiris/AAPaperData/MATTER_LBT_RunningAlphaS_Q2qhat/{directory_name}/{filename}_test.out",
            events_per_chunk=20,
            #max_chunks=3,
        )
