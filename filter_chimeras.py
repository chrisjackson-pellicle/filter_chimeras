#!/usr/bin/env python


"""
Filter Chimeras with Genome Reference (Captus or HybPiper)

This script processes target-capture assembly data to identify and filter out chimeric sequences
by mapping them to a reference genome.

Use ``--captus_folder`` for Captus output, or ``--hybpiper_folder`` for HybPiper output.

The script performs the following main tasks:
1. Creates a mapping reference from a high-quality genome
2. Parses target gene files to extract gene IDs
3. Maps assembly-derived contig fragments to the reference genome
4. Identifies chimeric sequences based on mapping patterns
5. Filters out chimeric sequences based on configurable thresholds
6. Generates reports and visualizations of the results

########################################################################################################################
Additional information:

This script is designed to work with Captus or HybPiper assembly data and requires:
- A high-quality reference genome in FASTA format
- A target file in FASTA format as used for the assembly runs
- A parent folder containing per-sample folders (``--captus_folder`` or ``--hybpiper_folder``)

The script uses BBmap.sh for genome mapping.
########################################################################################################################
"""

import logging
import sys
import argparse
import os
import socket
import glob
from concurrent.futures import as_completed
from pebble import ProcessPool
from multiprocessing import Manager
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import datetime
import itertools
import subprocess
import traceback
import time
import threading
from tqdm import tqdm
from concurrent.futures import as_completed
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import defaultdict
import datetime
import itertools
import subprocess
import traceback
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_size import Fixed
import seaborn as sns
import itertools
import pandas as pd
from collections import Counter
import numpy as np

import multiprocessing

# Linux-specific multiprocessing configuration to prevent freezing
if sys.platform.startswith('linux'):
    # Use 'spawn' method on Linux to avoid issues with fork
    multiprocessing.set_start_method('spawn', force=True)
    multiprocessing.set_executable(sys.executable)


########################################################################################################################
########################################################################################################################
# Get current working directory and host name:

cwd = os.getcwd()
host = socket.gethostname()

# Initialize logger
logger = None

# Set version
__version__ = '0.0.3'


class LogManager:
    """Global log manager that handles automatic cleanup."""

    def __init__(self):
        self.logger = None
        self.log_queue = None
        self.log_listener = None

    def setup(self, name, log_file, log_directory=None, **kwargs):
        """Set up the global logger."""
        self.logger, self.log_queue, self.log_listener = setup_queue_logger(
            name, log_file, log_directory, **kwargs
        )
        return self.logger, self.log_queue, self.log_listener

    def cleanup(self, timeout=2.0):
        """Clean up the global logger."""
        if self.log_listener:
            self.log_listener.flush(timeout=timeout)
            self.log_listener.stop()
            self.log_listener.join()

    def handle_error(self, error, tb, operation_name, id=None):
        """Handle errors with automatic cleanup."""
        if self.logger:
            if id:
                self.logger.error(f'{"[ERROR]:":15} {operation_name} failed for ID {id}')
            else:
                self.logger.error(f'{"[ERROR]:":15} {operation_name} failed')

            if isinstance(error, subprocess.CalledProcessError):
                self.logger.error(f'{"":15} Command: {" ".join(error.cmd)}')
                self.logger.error(f'{"":15} Return code: {error.returncode}')
                self.logger.error(f'{"":15} STDOUT: {error.stdout}')
                self.logger.error(f'{"":15} STDERR: {error.stderr}')
                self.logger.error(f'{"":15} Full traceback:')
                self.logger.error(tb)

            elif isinstance(error, Exception):
                self.logger.error(f'{"":15} Exception: {error}')
                self.logger.error(f'{"":15} Full traceback:')
                self.logger.error(tb)

        self.cleanup()
        sys.exit(1)


# Global instance
log_manager = LogManager()


class QueueHandler(logging.Handler):
    """A handler that sends log records to a queue for multiprocessing-safe logging.

    This handler is designed to work with multiprocessing environments where
    direct logging to files or console from worker processes can cause issues.
    Instead, it sends log records to a shared queue that is consumed by a
    listener thread in the main process.

    Attributes:
        queue: A multiprocessing-safe queue for sending log records
    """

    def __init__(self, queue):
        """Initialize the QueueHandler.

        Args:
            queue: A multiprocessing-safe queue (e.g., from multiprocessing.Manager)
        """
        super().__init__()
        self.queue = queue

    def emit(self, record):
        """Send the log record to the queue.

        Args:
            record: The log record to send
        """
        try:
            self.queue.put_nowait(record)
        except Exception:
            self.handleError(record)


class QueueListener(threading.Thread):
    """A thread that listens for log records from a queue and forwards them to handlers.

    This listener runs in the main process and consumes log records from a shared
    queue, forwarding them to appropriate handlers (file, console, etc.) based on
    their log levels. This ensures thread-safe logging in multiprocessing environments.

    Attributes:
        queue: A multiprocessing-safe queue for receiving log records
        handlers: List of logging handlers to forward records to
        _stop_event: Threading event to signal when to stop listening
    """

    def __init__(self, queue, *handlers):
        """Initialize the QueueListener.

        Args:
            queue: A multiprocessing-safe queue (e.g., from multiprocessing.Manager)
            *handlers: Variable number of logging handlers to forward records to
        """
        super().__init__()
        self.queue = queue
        self.handlers = handlers
        self._stop_event = threading.Event()

    def run(self):
        """Main loop that listens for log records and forwards them to handlers.

        Continuously polls the queue for log records and forwards them to
        appropriate handlers based on log level. Stops when a sentinel value
        (None) is received or when stop() is called.
        """
        while not self._stop_event.is_set():
            try:
                record = self.queue.get(timeout=1.0)
                if record is None:  # Sentinel to stop
                    break
                for handler in self.handlers:
                    # Only forward to handler if the record level meets the handler's level
                    if record.levelno >= handler.level:
                        handler.handle(record)
            except Exception:
                continue

    def stop(self):
        """Stop the listener thread.

        Sets the stop event and sends a sentinel value to the queue to ensure
        the thread exits cleanly.
        """
        self._stop_event.set()
        self.queue.put(None)  # Sentinel to stop the thread

    def flush(self, timeout=5.0):
        """Wait for all queued log records to be processed.

        Args:
            timeout (float): Maximum time to wait in seconds

        Returns:
            bool: True if queue was flushed successfully, False if timeout occurred
        """
        start_time = time.time()
        while not self.queue.empty() and (time.time() - start_time) < timeout:
            time.sleep(0.01)  # Small sleep to avoid busy waiting

        return self.queue.empty()


def setup_queue_logger(name, log_file, log_directory=None, console_level=logging.INFO,
                       file_level=logging.DEBUG, logger_object_level=logging.DEBUG):
    """Set up a queue-based logger for multiprocessing-safe logging.

    Creates a logger that uses a queue-based system for thread-safe logging
    in multiprocessing environments. The logger sends records to a shared
    queue, which are then consumed by a listener thread that forwards them
    to file and console handlers.

    Args:
        name (str): Name for the logger instance
        log_file (str): Filename for log file
        log_directory (str, optional): Directory for log files. If None, uses current directory
        console_level: Logging level for console output. Defaults to logging.INFO
        file_level: Logging level for file output. Defaults to logging.DEBUG
        logger_object_level: Logging level for logger object. Defaults to logging.DEBUG

    Returns:
        tuple: A tuple containing:
            - logger: Configured logger instance
            - queue: Multiprocessing-safe queue for log records
            - listener: Thread that consumes and forwards log records
    """
    # Get date and time string for log filename
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Create log directory if supplied
    if log_directory:
        if not os.path.exists(log_directory):
            os.makedirs(log_directory)
        log_file_name = f'{log_directory}/{log_file}_{date_and_time}.log'
    else:
        log_file_name = f'{log_file}_{date_and_time}.log'

    # Create handlers
    file_handler = logging.FileHandler(log_file_name, mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter(
        '%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_format)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Create manager and shared queue
    manager = Manager()
    queue = manager.Queue()
    queue_handler = QueueHandler(queue)

    # Set up logger
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)
    logger_object.addHandler(queue_handler)

    # Create and start listener
    listener = QueueListener(queue, file_handler, console_handler)
    listener.start()

    return logger_object, queue, listener


def setup_worker_logger(name, queue):
    """Set up a logger for worker processes that sends records to a shared queue.

    Creates a logger specifically designed for worker processes in multiprocessing
    environments. This logger sends all log records to a shared queue instead of
    directly to handlers, ensuring thread-safe logging across process boundaries.

    Args:
        name (str): Logger name (typically __name__ from the calling module)
        queue: Multiprocessing-safe queue for sending log records

    Returns:
        logging.Logger: Configured logger that sends records to the shared queue
    """
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logging.DEBUG)

    # Remove any existing handlers to avoid duplicates
    for handler in logger_object.handlers[:]:
        logger_object.removeHandler(handler)

    # Add queue handler
    queue_handler = QueueHandler(queue)
    logger_object.addHandler(queue_handler)

    return logger_object


def setup_log_and_report_directories(args):
    """Set up log and report directories and add them to args.

    Creates the log and report directories under the output directory structure
    and sets them as attributes on the args object for later use.

    Args:
        args: argparse.Namespace object containing output_directory attribute

    Returns:
        None: Modifies args object in place by adding log_directory and report_directory attributes
    """
    log_directory = os.path.join(args.output_directory, '00_logs_and_reports/logs')
    report_directory = os.path.join(args.output_directory, '00_logs_and_reports/reports')
    createfolder(log_directory)
    createfolder(report_directory)
    args.log_directory = log_directory
    args.report_directory = report_directory


def log_separator(logger, length=100):
    """
    Log a separator line with consistent formatting.

    Args:
        logger: Logger instance for logging messages
        length (int): Length of the separator line (default: 100)
    """
    logger.info(f'{" " * 11}{"-" * length}')


def format_elapsed_time(elapsed_time):
    """
    Format elapsed time in seconds to hours, minutes, and seconds format.

    Args:
        elapsed_time (float): Elapsed time in seconds

    Returns:
        tuple: (hours, minutes, seconds, formatted_string)
    """
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)
    formatted_string = f'{hours}h {minutes}m {seconds}s ({elapsed_time:.2f} seconds)'

    return hours, minutes, seconds, formatted_string


def print_arguments(args, logger, __version__):
    """Prints the arguments to screen and log.

    Args:
        args (argparse.Namespace): Parsed command line arguments containing input file paths
            and other configuration options.
        logger (logging.Logger): Logger object.
        __version__ (str): Version of the script.
    """

    logger.info(f'{"[INFO]:":10} SCRIPT_NAME version {__version__} was called with these arguments:\n')

    for parameter, value in args.__dict__.items():
        if parameter not in ['func', 'logger', 'log_directory', 'report_directory', 'subcommand_name']:
            logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')


def exit_program():
    """Exit the program with a clean shutdown.

    Args:
        None
    """
    log_manager.cleanup()
    sys.exit(1)


def log_completion_time(start_time, logger=None, label='Completed'):
    """
    Log or print a standardized completion-time message given a start time.

    Args:
        start_time (float): Start timestamp from time.time().
        logger (logging.Logger, optional): Logger to write to; falls back to print if None.
        label (str): Message label preceding the duration (default: 'Completed').
    """
    try:
        elapsed_seconds = time.time() - start_time
        hours, minutes, seconds, _ = format_elapsed_time(elapsed_seconds)
        duration_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
        message = f'{"[INFO]:":10} {label} in {duration_str} ({elapsed_seconds:.2f} seconds)'
        if logger:
            logger.info(message)
        else:
            # If no logger, print a simplified message without the fixed-width tag
            print(f"{label} in {duration_str} ({elapsed_seconds:.2f} seconds)")
    except Exception:
        # Swallow any timing/logging errors to avoid masking real failures
        pass


def _setup_output_directories(args):
    """Set up output directories for the fungi_search pipeline.

    Args:
        args (argparse.Namespace): Parsed command line arguments containing output_directory

    Returns:
        dict: Dictionary mapping directory names to their full paths
    """
    outdirs = {
        'parent_output_directory': args.output_directory,
        'report_directory': args.report_directory,
        'log_directory': args.log_directory,
        'mapped_seqs': os.path.join(args.output_directory, '01_target_capture_seqs_mapped'),
        'unmapped_seqs': os.path.join(args.output_directory, '02_target_capture_seqs_unmapped'),
        'non_chimeric_target_capture_seqs': os.path.join(args.output_directory, f'03_non_chimeric_target_capture_seqs_{args.min_samples_threshold}'),
        'all_target_capture_seqs': os.path.join(args.output_directory, '04_all_target_capture_seqs'),
    }

    return outdirs


class FutureTracker:
    """Class to track Pebble ProcessPool for cleanup on interruption."""

    def __init__(self):
        self.pool = None

    def set_pool(self, pool):
        """Set the Pebble ProcessPool for cleanup."""
        self.pool = pool

    def cancel_all(self):
        """Stop the pool and terminate all running processes.
        
        Pebble's ProcessPool.stop() terminates all running processes immediately,
        unlike concurrent.futures which can only cancel pending (not running) tasks.
        """
        if self.pool:
            print(f"\n{'':10} Stopping pool and terminating all active processes...")
            self.pool.stop()
            self.pool.join()
            self.pool = None


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure

    :param str directory: name of directory to create
    """

    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
        return directory
    except OSError:
        logger.info(f'Error: Creating directory: {directory}')


def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes
    """
    # Check if file exist and is not empty
    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def create_genome_mapping_reference(args):
    """
    Create a mapping reference for a genome using BBmap.sh.

    This function creates a mapping reference for a genome using BBmap.sh.
    If the reference already exists, it skips the creation step.

    Args:
        args (argparse.Namespace): Parsed command line arguments containing reference_genome_fasta and bbmap_Xmx

    Raises:
        ValueError: If there is an issue running BBmap.sh
    """

    genome_basename = os.path.basename(args.reference_genome_fasta)

    log_separator(logger)
    logger.info(f'{"[INFO]:":10} Creating mapping reference for genome "{genome_basename}" using BBmap.sh...')

    expected_file = f'{args.output_directory}/ref/genome/1/summary.txt'

    try:
        assert file_exists_and_not_empty(expected_file)
        logger.info(f'{"[INFO]:":10} BBmap.sh reference exists, skipping!')

    except AssertionError:
        try:
            command = (f'bbmap.sh -Xmx{args.bbmap_Xmx}g '
                       f'k=13 '
                       f'threads={args.threads} '
                       f'path={args.output_directory} '
                       f'ref={args.reference_genome_fasta}')

            result = subprocess.run(command,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=True,
                                    shell=True)

            logger.debug(f'BBmap.sh check_returncode() is: {result.check_returncode()}')
            logger.debug(f'BBmap.sh stdout is: {result.stdout}')
            logger.debug(f'BBmap.sh stderr is: {result.stderr}')

            return f'{args.output_directory}/ref'

        except subprocess.CalledProcessError as exc:
            logger.error(f'BBmap.sh FAILED. Output is: {exc}')
            logger.error(f'BBmap.sh stdout is: {exc.stdout}')
            logger.error(f'BBmap.sh stderr is: {exc.stderr}')

            raise ValueError('There was an issue running BBmap.sh. Check input files!')


def parse_target_file(target_file_fasta):
    """
    Parse a target FASTA file to extract gene IDs.

    This function reads a FASTA file and extracts gene IDs from sequence names.
    The gene ID is assumed to be the last part of the sequence name after the last hyphen.

    Args:
        target_file_fasta (str): Path to the target FASTA file

    Returns:
        list: A sorted list of unique gene IDs extracted from the FASTA file
    """

    target_file_basename = os.path.basename(target_file_fasta)

    log_separator(logger)
    logger.info(f'{"[INFO]:":10} Parsing target file "{target_file_basename}" to recover gene IDs...')

    gene_ids = set()

    for seq in SeqIO.parse(target_file_fasta, 'fasta'):
        gene_id = seq.name.split('-')[-1]
        gene_ids.add(gene_id)

    if len(gene_ids) == 0:
        logger.error(f'{"[ERROR]:":10} No gene IDs found in target file "{target_file_basename}"!')
        log_manager.cleanup()
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Number of unique gene IDs: {len(gene_ids)}')

    return sorted(list(gene_ids))


def _parse_ref_protein_coord_interval(ref_coord_str):
    """Parse ``start-end`` into ``(start, end)`` integers, or None if not parseable."""
    ref_coord_str = ref_coord_str.strip()
    if '-' not in ref_coord_str:
        return None
    left, right = ref_coord_str.split('-', 1)
    try:
        return int(left.strip()), int(right.strip())
    except ValueError:
        return None


def _intervals_merge_under_wiggle(s1, e1, s2, e2, wiggle):
    """True if closed intervals [s1,e1] and [s2,e2] overlap after expansion by ``wiggle`` on both sides."""
    a1, b1 = s1 - wiggle, e1 + wiggle
    a2, b2 = s2 - wiggle, e2 + wiggle
    return not (b1 < a2 or b2 < a1)


def _cluster_ref_coord_strings(coord_strings, wiggle):
    labels = list(coord_strings)
    n = len(labels)
    parent = list(range(n))

    def find(i):
        if parent[i] != i:
            parent[i] = find(parent[i])
        return parent[i]

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    parsed = [_parse_ref_protein_coord_interval(lab) for lab in labels]
    for i in range(n):
        for j in range(i + 1, n):
            if parsed[i] is None or parsed[j] is None:
                continue
            s1, e1 = parsed[i]
            s2, e2 = parsed[j]
            if _intervals_merge_under_wiggle(s1, e1, s2, e2, wiggle):
                union(i, j)

    clusters = defaultdict(list)
    for idx, lab in enumerate(labels):
        clusters[find(idx)].append(lab)
    return list(clusters.values())


def _sort_key_ref_coord_label(lab):
    p = _parse_ref_protein_coord_interval(lab)
    if p is None:
        return (1 << 60, 0, lab)
    return (p[0], p[1], lab)


def _filename_coord_suffix_for_cluster(member_labels):
    ordered = sorted(set(member_labels), key=_sort_key_ref_coord_label)
    return '_'.join(ordered)


def _group_introns_by_ref_protein_coords(records, intron_ref_coords_wiggle):
    str_to_recs = defaultdict(list)
    for rec in records:
        if 'ref_protein_coords:' not in rec.description:
            continue
        ref_coords = rec.description.split('ref_protein_coords:')[-1].strip()
        str_to_recs[ref_coords].append(rec)

    if not str_to_recs:
        return []

    if intron_ref_coords_wiggle < 0:
        clusters = [[lab] for lab in str_to_recs.keys()]
    else:
        clusters = _cluster_ref_coord_strings(str_to_recs.keys(), intron_ref_coords_wiggle)

    out = []
    for members in clusters:
        rec_list = []
        for lab in members:
            rec_list.extend(str_to_recs[lab])
        suffix = _filename_coord_suffix_for_cluster(members)
        out.append((suffix, rec_list))
    return out


def sample_name_from_folder(assembly_source, sample_folder):
    if assembly_source == 'captus':
        return os.path.basename(sample_folder).split('__')[0]
    if assembly_source == 'hybpiper':
        return hybpiper_sample_name_from_folder(sample_folder)
    raise ValueError(f'Unknown assembly_source: {assembly_source!r}')


def load_sample_nuc_seqrecords_by_gene_paralog(assembly_source, sample_path, sample_name):
    gene_paralog = {}
    if assembly_source == 'captus':
        combined_gene_fasta_file = f'{sample_path}/01_coding_NUC/NUC_coding_NT.fna'
        if not os.path.isfile(combined_gene_fasta_file):
            logger.warning(f'{"[WARNING]:":10} Missing Captus NUC file: {combined_gene_fasta_file}')
            return gene_paralog
        for seq in SeqIO.parse(combined_gene_fasta_file, 'fasta'):
            name_split = seq.name.split('__')
            gene = name_split[1]
            if len(name_split) == 3:
                paralog_number = int(name_split[2])
                gene_paralog.setdefault(gene, {})[f'paralog_{int(paralog_number) + 1}'] = seq
            else:
                gene_paralog.setdefault(gene, {})['paralog_1'] = seq
        return gene_paralog
    if assembly_source == 'hybpiper':
        logger.debug(f'HybPiper NUC load not implemented for sample {sample_name!r}; returning empty dict.')
        return gene_paralog
    raise ValueError(f'Unknown assembly_source: {assembly_source!r}')


def _read_gene_sample_pivot(tsv_path, value_name):
    """Load a gene × sample TSV and return a pivoted matrix suitable for seaborn heatmaps."""
    df = pd.read_csv(tsv_path, sep='\t')
    df = df.melt(id_vars=['gene'], var_name='sample_id', value_name=value_name)
    return df.pivot(index='gene', columns='sample_id', values=value_name)


def _non_chimeric_paralog_pivot_capped(tsv_path, cap):
    """Read non-chimeric paralog count TSV, cap counts at ``cap``, return pivot (gene × sample).

    Rows are alphabetical by gene name (pandas' default ``pivot`` behaviour) so this heatmap aligns
    row-for-row with the matrix produced by ``_read_gene_sample_pivot`` (used by heatmaps 04 and 07).
    """
    df = pd.read_csv(tsv_path, sep='\t')
    sample_cols = df.columns[1:]
    capped = df.loc[:, sample_cols]
    df.loc[:, sample_cols] = capped.where((capped < cap) | capped.isna(), cap)
    melted = df.melt(id_vars=['gene'], var_name='sample_id', value_name='non_chimeric_paralog_count')
    return melted.pivot(index='gene', columns='sample_id', values='non_chimeric_paralog_count')


def _sorted_sample_names_from_chimera_dict(chimera_count_dict):
    names = {s for _g, sd in chimera_count_dict.items() for s in sd}
    return sorted(names)


def _make_cell_sized_figure(df, cell_size_inches=0.30):
    """Create a (fig, ax) whose dimensions scale with ``df.shape`` so each heatmap cell is ``cell_size_inches`` inches.

    Adds ~4 in of width budget for the colorbar + y-tick labels and ~2 in of height for the x-tick
    labels, with a 6-in floor so very small heatmaps remain legible. Returns the figure, axes, and
    the resolved ``(fig_w, fig_h)`` for downstream figure-fraction calculations.
    """
    n_rows, n_cols = df.shape
    fig_w = max(6.0, n_cols * cell_size_inches + 4.0)
    fig_h = max(6.0, n_rows * cell_size_inches + 2.0)
    sns.set_style('ticks')
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    return fig, ax, fig_w, fig_h


def _attach_cell_sized_cbar(fig, ax, label, cell_size_inches, *, ticks=None, ticklabels=None,
                            integer_ticks=False):
    """Attach a manual colorbar whose width is exactly one heatmap cell.

    The colorbar title is placed on the LEFT side, snug against the bar (``labelpad=2``), so swatch
    space to the right of the bar stays unambiguous. ``ticks`` is forwarded to ``fig.colorbar``;
    ``ticklabels`` (if given) replaces the auto-generated tick text via ``set_yticklabels``;
    ``integer_ticks=True`` installs a ``MaxNLocator(integer=True)`` for continuous integer data.
    Returns the cbar's Axes so callers can position swatches relative to it.
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=Fixed(cell_size_inches), pad=Fixed(cell_size_inches))
    mappable = ax.collections[0]
    cbar = fig.colorbar(mappable, cax=cax, extend='neither', ticks=ticks)
    if ticklabels is not None:
        cbar.ax.set_yticklabels(ticklabels)
    if integer_ticks:
        cbar.ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    cbar.ax.yaxis.set_label_position('left')
    cbar.set_label(label, rotation=90, labelpad=2)
    return cax


def _build_sample_to_gene_seqrecord_dict(assembly_source, sample_paths):
    """Return ``sample_name -> gene -> paralog -> SeqRecord`` for every assembly output folder."""
    sample_to_gene_seqrecord_dict = defaultdict(lambda: defaultdict(dict))
    for sample_path in sample_paths:
        sample_name = sample_name_from_folder(assembly_source, sample_path)
        gene_paralog = load_sample_nuc_seqrecords_by_gene_paralog(assembly_source, sample_path, sample_name)
        for gene, paralog_dict in gene_paralog.items():
            sample_to_gene_seqrecord_dict[sample_name][gene].update(paralog_dict)
    return sample_to_gene_seqrecord_dict


def _classify_paralog_sequences(combined_samples_dict, sample_to_gene_seqrecord_dict,
                                combined_intron_seqrecord_dict, non_chimera_statuses, not_tested_statuses):
    """Bucket every (sample, gene, paralog) into non-chimeric / chimeric / not-tested, with running totals.

    Returns ``(non_chimeric_seqs_dict, non_chimeric_introns_dict, unfiltered_seqs_dict,
    chimera_count_dict, counts)`` where ``counts`` is a small dict of five totals
    (``all``, ``non_chimeric``, ``chimeric``, ``not_tested``, ``filtered_out``).
    """
    non_chimeric_seqs_dict = defaultdict(lambda: defaultdict(list))
    non_chimeric_introns_dict = defaultdict(lambda: defaultdict(list))
    unfiltered_seqs_dict = defaultdict(lambda: defaultdict(list))
    chimera_count_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    counts = {'all': 0, 'non_chimeric': 0, 'chimeric': 0, 'not_tested': 0, 'filtered_out': 0}

    for sample, genes_dict in sorted(combined_samples_dict.items()):
        for gene_name, paralog_dict in genes_dict.items():
            # zero all counters
            chimera_count_dict[gene_name][sample]['non_chimera_count'] = 0
            chimera_count_dict[gene_name][sample]['chimera_count'] = 0
            chimera_count_dict[gene_name][sample]['not_tested'] = 0

            for paralog_name, chimera_text in paralog_dict.items():
                counts['all'] += 1
                seqrecord = sample_to_gene_seqrecord_dict[sample][gene_name][paralog_name]
                unfiltered_seqs_dict[gene_name][sample].append(seqrecord)

                if chimera_text in non_chimera_statuses:
                    non_chimeric_seqs_dict[gene_name][sample].append(seqrecord)
                    chimera_count_dict[gene_name][sample]['non_chimera_count'] += 1
                    counts['non_chimeric'] += 1
                    intron_seqrecord_list = (
                        combined_intron_seqrecord_dict
                        .get(sample, {})
                        .get(gene_name, {})
                        .get(paralog_name, []))
                    for intron_seqrecord in intron_seqrecord_list:
                        non_chimeric_introns_dict[gene_name][sample].append(intron_seqrecord)
                elif chimera_text in not_tested_statuses:
                    chimera_count_dict[gene_name][sample]['not_tested'] += 1
                    counts['filtered_out'] += 1
                    counts['not_tested'] += 1
                else:
                    chimera_count_dict[gene_name][sample]['chimera_count'] += 1
                    counts['chimeric'] += 1
                    counts['filtered_out'] += 1

    return non_chimeric_seqs_dict, non_chimeric_introns_dict, unfiltered_seqs_dict, chimera_count_dict, counts


def _write_chimera_count_reports(outdir_reports, chimera_count_dict, sample_names_sorted,
                                 sample_names_joined):
    """Write ``01_chimera_report_...tsv`` (raw N/C counts) and ``02_chimera_report_heatmap_data.tsv``.

    Returns the heatmap-data TSV path so the 03 heatmap helper can read it back.
    """
    outfile_chimera_report = os.path.join(
        outdir_reports, '01_chimera_report_non_chimera_count_vs_chimera_count.tsv')
    outfile_chimera_report_heatmap = os.path.join(
        outdir_reports, '02_chimera_report_heatmap_data.tsv')

    log_separator(logger)
    logger.info(f'{"[INFO]:":10} Writing report: {outfile_chimera_report}')
    logger.info(f'{"[INFO]:":10} Writing report: {outfile_chimera_report_heatmap}')

    with (open(outfile_chimera_report, 'w') as chimera_report_handle,
          open(outfile_chimera_report_heatmap, 'w') as heatmap_data_handle):
        chimera_report_handle.write(f'gene\t{sample_names_joined}\n')
        heatmap_data_handle.write(f'gene\t{sample_names_joined}\n')

        for gene, sample_dict in sorted(chimera_count_dict.items()):
            gene_row = [gene]
            gene_row_heatmap_data = [gene]
            for sample in sample_names_sorted:
                sample_count_dict = sample_dict[sample]
                if len(sample_count_dict) == 0:
                    logger.debug(f'No assembly sequences for gene {gene} sample {sample}!')
                    gene_row.append('no_seqs')
                    gene_row_heatmap_data.append('NaN')
                else:
                    non_chimera_count = sample_count_dict['non_chimera_count']
                    chimera_count = sample_count_dict['chimera_count']
                    gene_row.append(f'{non_chimera_count}/{chimera_count}')
                    total_count = non_chimera_count + chimera_count
                    if total_count == 0:  # i.e. no _tested_ sequences found for this gene/sample combination
                        logger.debug(f'For sample {sample} and gene {gene}, no _tested_ sequences found!')
                        gene_row_heatmap_data.append('2')
                    else:
                        gene_row_heatmap_data.append(str(chimera_count / total_count))
            chimera_report_handle.write('\t'.join(gene_row) + '\n')
            heatmap_data_handle.write('\t'.join(gene_row_heatmap_data) + '\n')

    return outfile_chimera_report_heatmap


def _filter_min_paralogs_per_sample(non_chimeric_seqs_dict, min_num_paralogs_per_sample):
    """Return ``gene -> sample -> seqrecord_list`` keeping only samples with >= ``min_num_paralogs_per_sample``."""
    filtered = defaultdict(dict)
    for gene_name, sample_dict in sorted(non_chimeric_seqs_dict.items()):
        for sample, seqrecord_list in sample_dict.items():
            if len(seqrecord_list) >= min_num_paralogs_per_sample:
                filtered[gene_name][sample] = seqrecord_list
    return filtered


def _write_single_seq_genes_report(outdir_reports, non_chimeric_seqs_dict, single_seq_ratio=0.90):
    """Write ``04_genes_with_one_seq_per_sample_in_90_percent_samples.tsv``.

    Lists genes whose share of single-sequence samples (i.e. samples that recovered exactly one
    non-chimeric paralog) is at least ``single_seq_ratio`` (default 0.90).
    """
    single_seq_genes_dict = {}
    for gene_name, sample_dict in sorted(non_chimeric_seqs_dict.items()):
        gene_number_of_seqs = [len(seqs) for seqs in sample_dict.values()]
        if not gene_number_of_seqs:
            continue
        number_of_samples = len(gene_number_of_seqs)
        samples_with_only_one_seq = sum(1 for v in gene_number_of_seqs if v == 1)
        ratio = samples_with_only_one_seq / number_of_samples
        if ratio >= single_seq_ratio:
            logger.debug(
                f'{"[DEBUG]:":10} For gene {gene_name}, the number of samples with a single sequence, which is also '
                f'non-chimeric: {samples_with_only_one_seq}. The total number of samples is: {number_of_samples}. '
                f'The ratio of single sequences is: {ratio}')
            single_seq_genes_dict[gene_name] = (samples_with_only_one_seq, number_of_samples, ratio)

    outfile = os.path.join(outdir_reports, '04_genes_with_one_seq_per_sample_in_90_percent_samples.tsv')
    logger.info(f'{"[INFO]:":10} Writing report: {outfile}')
    with open(outfile, 'w') as fh:
        fh.write('gene_name\tnumber_samples_with_single_non_chimeric_seq\ttotal_samples\tratio_of_single_seqs\n')
        for gene_name, (single, total, ratio) in single_seq_genes_dict.items():
            fh.write(f'{gene_name}\t{single}\t{total}\t{ratio}\n')


def _apply_sample_threshold_and_write_fastas(outdir_non_chimeras, outdir_reports, filtered_dict,
                                             non_chimeric_introns_dict, min_threshold,
                                             min_samples_threshold, min_num_paralogs_per_sample,
                                             total_number_of_samples, intron_ref_coords_wiggle):
    """For every gene clearing ``min_threshold``, write ``<gene>.fasta`` (+ grouped intron FASTAs).

    Genes that don't clear the threshold are recorded in
    ``05_genes_not_recovered_after_sample_threshold_<min_samples_threshold>_filtering.tsv`` and their
    samples are stamped with a count of 0 in the returned ``recovered_dict`` so the downstream
    paralog-count TSV/heatmaps render those rows as NaN cells (lightgrey ``no_seqs``).

    Returns ``recovered_dict`` (``gene -> sample -> non_chimeric_paralog_count``).
    """
    recovered_dict = defaultdict(lambda: defaultdict(int))
    genes_not_recovered_dict = {}

    for gene, sample_dict in filtered_dict.items():
        gene_seqs_to_write = []
        gene_introns_to_write = []
        if len(sample_dict) >= min_threshold:
            for sample, seqrecord_list in sample_dict.items():
                for seq in seqrecord_list:
                    gene_seqs_to_write.append(seq)
                    recovered_dict[gene][sample] += 1
                intron_seqrecord_list = non_chimeric_introns_dict.get(gene, {}).get(sample, [])
                gene_introns_to_write.extend(intron_seqrecord_list)

            gene_outfile = os.path.join(outdir_non_chimeras, f'{gene}.fasta')
            with open(gene_outfile, 'w') as fh:
                SeqIO.write(gene_seqs_to_write, fh, 'fasta')

            if gene_introns_to_write:
                outdir_non_chimeras_introns = createfolder(os.path.join(outdir_non_chimeras, 'introns'))
                for coord_suffix, rec_list in _group_introns_by_ref_protein_coords(
                        gene_introns_to_write, intron_ref_coords_wiggle):
                    intron_path = os.path.join(
                        outdir_non_chimeras_introns, f'{gene}_ref_protein_coords_{coord_suffix}.fasta')
                    with open(intron_path, 'w') as fh:
                        SeqIO.write(rec_list, fh, 'fasta')
        else:
            logger.debug(
                f'Not recovering sequences from any samples for gene {gene}, as there are fewer than '
                f'{min_threshold} of {total_number_of_samples} samples containing '
                f'{min_num_paralogs_per_sample} or more sequences!')
            genes_not_recovered_dict[gene] = len(sample_dict)
            for sample in sample_dict:
                recovered_dict[gene][sample] = 0

    if genes_not_recovered_dict:
        genes_path = os.path.join(
            outdir_reports,
            f'05_genes_not_recovered_after_sample_threshold_{min_samples_threshold}_filtering.tsv')
        with open(genes_path, 'w') as fh:
            fh.write(
                'min_samples_threshold_fraction\ttotal_number_of_samples\tmin_threshold_number\n'
                f'{min_samples_threshold}\t{total_number_of_samples}\t{min_threshold}\n\n'
                'gene_name\tnumber_samples_with_paralogs\t\n')
            for gene_name, n in genes_not_recovered_dict.items():
                fh.write(f'{gene_name}\t{n}\t\n')

    return recovered_dict


def _write_paralog_count_report(outdir_reports, min_samples_threshold, recovered_dict,
                                sample_names_sorted, sample_names_joined):
    """Write ``06_non_chimeric_paralogs_after_sample_threshold_<min_samples_threshold>_filtering.tsv``.

    Returns the path so the downstream heatmaps (07/08/09) can pivot it for plotting.
    """
    outfile = os.path.join(
        outdir_reports,
        f'06_non_chimeric_paralogs_after_sample_threshold_{min_samples_threshold}_filtering.tsv')
    logger.info(f'{"[INFO]:":10} Writing report: {outfile}')

    with open(outfile, 'w') as fh:
        fh.write(f'gene\t{sample_names_joined}\n')
        for gene, sample_dict in sorted(recovered_dict.items()):
            gene_row = [gene]
            for sample_name in sample_names_sorted:
                if sample_name not in sample_dict:
                    logger.debug(
                        f'No post-threshold non-chimeric paralog count for gene {gene!r} sample {sample_name!r} '
                        f'(writing NaN); sample may be below min_num_paralogs_per_sample for this gene, '
                        f'or absent from the per-gene recovery set while still present in chimera stats for other genes.')
                    gene_row.append('NaN')
                else:
                    count = sample_dict[sample_name]
                    gene_row.append('NaN' if count == 0 else str(count))
            fh.write('\t'.join(gene_row) + '\n')

    return outfile


def _write_all_seqs_fastas(outdir_all_seqs, unfiltered_seqs_dict):
    """Write ``<gene>.all.fasta`` for every gene (every recovered paralog, chimeric or not)."""
    for gene, sample_dict in unfiltered_seqs_dict.items():
        gene_seqs_to_write = []
        for seqrecord_list in sample_dict.values():
            gene_seqs_to_write.extend(seqrecord_list)
        gene_outfile = os.path.join(outdir_all_seqs, f'{gene}.all.fasta')
        with open(gene_outfile, 'w') as fh:
            SeqIO.write(gene_seqs_to_write, fh, 'fasta')


def _add_cbar_swatches(fig, cax, swatches, cell_size_inches, fig_w, fig_h,
                       *, gap_inches=0.9, inter_swatch_gap_cells=0.5):
    """Place a vertical stack of single-cell legend swatches to the right of a colorbar.

    ``swatches`` is an iterable of ``(facecolor, label_text)`` pairs. Each swatch is a Rectangle
    of exactly ``cell_size_inches × cell_size_inches`` with the label half a cell to its right.
    The top swatch is top-aligned with the colorbar, each subsequent swatch sits one cell-height
    plus ``inter_swatch_gap_cells`` of a cell below the previous one. ``gap_inches`` is the
    horizontal clearance between the cbar tick labels and the swatch column (default 0.9 in).
    """
    fig.canvas.draw()  # commits the divider layout so cax.get_position() is meaningful
    cax_pos = cax.get_position()
    cell_w_frac = cell_size_inches / fig_w
    cell_h_frac = cell_size_inches / fig_h
    gap_frac = gap_inches / fig_w
    swatch_left = cax_pos.x1 + gap_frac
    inter_gap_frac = (cell_size_inches * inter_swatch_gap_cells) / fig_h
    label_x = swatch_left + cell_w_frac + (cell_size_inches * 0.5) / fig_w

    swatch_bottom = cax_pos.y1 - cell_h_frac  # top-align first swatch with cbar top
    for facecolor, label in swatches:
        sw_ax = fig.add_axes([swatch_left, swatch_bottom, cell_w_frac, cell_h_frac])
        sw_ax.set_xlim(0, 1)
        sw_ax.set_ylim(0, 1)
        sw_ax.set_xticks([])
        sw_ax.set_yticks([])
        sw_ax.add_patch(Rectangle((0, 0), 1, 1, facecolor=facecolor, edgecolor='black', linewidth=0.5))
        fig.text(label_x, swatch_bottom + cell_h_frac / 2, label, ha='left', va='center', fontsize=8)
        swatch_bottom -= cell_h_frac + inter_gap_frac


def write_sequences_and_reports(outdirs,
                                assembly_source,
                                sample_parent_folder,
                                combined_samples_dict,
                                combined_intron_seqrecord_dict,
                                min_samples_threshold,
                                min_num_paralogs_per_sample,
                                intron_ref_coords_wiggle,
                                non_chimeric_paralog_max_count):
    """
    Build chimera/non-chimera reports and write filtered FASTA outputs after per-sample mapping.
    """
    non_chimera_statuses = frozenset({'no_chimera_from_multi_hit', 'no_chimera_as_single_contig_hit'})
    not_tested_statuses = frozenset({
        'no_detectable_chimera_as_single_contig_hit_after_length_filtering',
        'no_detectable_chimera_as_no_seqs_left_after_length_filtering',
        'no_chimera_test_performed_as_too_little_seq_length_after_length_filtering',
        'no_chimera_test_performed_as_too_few_contigs_after_length_filtering',
        'unknown_repeated_subject',
    })

    outdir_reports = createfolder(outdirs['report_directory'])
    outdir_non_chimeras = createfolder(outdirs['non_chimeric_target_capture_seqs'])
    outdir_all_seqs = createfolder(outdirs['all_target_capture_seqs'])
    sample_paths = sorted(glob.glob(f'{sample_parent_folder}/*'))

    #------------------------------------------------------------------------------------------------------------------
    # Calculate thresholds for writing filtered fasta files:
    #------------------------------------------------------------------------------------------------------------------
    total_number_of_samples = len(combined_samples_dict)
    min_threshold = min_samples_threshold * total_number_of_samples
    log_separator(logger)
    logger.info(f'{"[INFO]:":10} Total number of samples: {total_number_of_samples}')
    logger.info(
        f'{"[INFO]:":10} Per-gene sample-count threshold for writing non-chimeric FASTAs: '
        f'{min_threshold:.2f} (--min_samples_threshold {min_samples_threshold} x '
        f'{total_number_of_samples} samples).')
    logger.info(
        f'{"":10} A gene is only written when at least this many '
        f'samples each recover >= {min_num_paralogs_per_sample} non-chimeric paralog(s) '
        f'(--min_num_paralogs_per_sample).')

    #------------------------------------------------------------------------------------------------------------------
    # Classify every (sample, gene, paralog) and accumulate per-gene/per-sample counters:
    #------------------------------------------------------------------------------------------------------------------
    sample_to_gene_seqrecord_dict = _build_sample_to_gene_seqrecord_dict(assembly_source, sample_paths)

    (non_chimeric_seqs_dict, non_chimeric_introns_dict, unfiltered_seqs_dict,
     chimera_count_dict, counts) = _classify_paralog_sequences(
        combined_samples_dict, sample_to_gene_seqrecord_dict, combined_intron_seqrecord_dict,
        non_chimera_statuses, not_tested_statuses)

    logger.info(f'{"[INFO]:":10} Number of genes: {len(chimera_count_dict)}')
    logger.info(f'{"[INFO]:":10} Count of all sequences: {counts["all"]}')
    logger.info(f'{"[INFO]:":10} Count of non-chimeric sequences: {counts["non_chimeric"]}')
    logger.info(f'{"[INFO]:":10} Count of chimeric sequences: {counts["chimeric"]}')
    logger.info(f'{"[INFO]:":10} Count of sequences not tested (did not pass length or contig number filtering): {counts["not_tested"]}')
    logger.info(f'{"[INFO]:":10} Count of filtered out sequences (chimeras and not tested): {counts["filtered_out"]}')

    #------------------------------------------------------------------------------------------------------------------
    # Write the chimera count reports (01 raw counts, 02 chimera-proportion for heatmap 03):
    #------------------------------------------------------------------------------------------------------------------
    sample_names_sorted = _sorted_sample_names_from_chimera_dict(chimera_count_dict)
    sample_names_joined = '\t'.join(sample_names_sorted)
    outfile_chimera_report_heatmap = _write_chimera_count_reports(
        outdir_reports, chimera_count_dict, sample_names_sorted, sample_names_joined)

    #------------------------------------------------------------------------------------------------------------------
    # Heatmap 03: chimera proportion per (gene, sample), with chartreuse over-cells for "not tested":
    #------------------------------------------------------------------------------------------------------------------
    heatmap_name = os.path.join(
        outdir_reports, f'03_chimera_proportion_heatmap_{assembly_source}.png')
    logger.info(f'{"[INFO]:":10} Writing heatmap: {heatmap_name}')
    df = _read_gene_sample_pivot(outfile_chimera_report_heatmap, 'chimera_proportion')

    # Special-case the "no seqs tested" sentinel (value 2) by routing it to set_over / set_under
    # so the actual data range stays 0..1 (real chimera proportions).
    over_colour = 'chartreuse'
    cmap = matplotlib.colormaps['OrRd'].copy()
    cmap.set_over(over_colour)
    cmap.set_under(over_colour)
    df_for_scaling = df.replace(2, np.nan)
    vmin = df_for_scaling.min().min()
    vmax = df_for_scaling.max().max()
    if pd.isna(vmin) or pd.isna(vmax):
        vmin, vmax = 0, 1

    cell_size_inches = 0.30
    fig, ax, fig_w, fig_h = _make_cell_sized_figure(df, cell_size_inches)
    sns.heatmap(df, vmin=vmin, vmax=vmax, cmap=cmap, xticklabels=1, yticklabels=1, square=True,
                linewidth=0.5, cbar=False, ax=ax)
    ax.set_facecolor('lightgrey')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    cax = _attach_cell_sized_cbar(fig, ax, 'Fraction of sequences that are chimeric', cell_size_inches)
    _add_cbar_swatches(
        fig, cax,
        [(over_colour, 'No seqs tested'), ('lightgrey', 'No seqs recovered')],
        cell_size_inches, fig_w, fig_h,
    )

    plt.savefig(heatmap_name, dpi=300, bbox_inches='tight')
    plt.close(fig)

    #------------------------------------------------------------------------------------------------------------------
    # Filter to samples with >= min_num_paralogs_per_sample non-chimeric paralogs, then report "single-seq" genes:
    #------------------------------------------------------------------------------------------------------------------
    filtered_dict = _filter_min_paralogs_per_sample(non_chimeric_seqs_dict, min_num_paralogs_per_sample)

    #------------------------------------------------------------------------------------------------------------------
    # Write a report of genes where >=90% of samples had a single non-chimeric sequence after filtering:
    #------------------------------------------------------------------------------------------------------------------
    _write_single_seq_genes_report(outdir_reports, non_chimeric_seqs_dict)

    #------------------------------------------------------------------------------------------------------------------
    # Apply per-gene sample-count threshold: write recovered FASTAs (+ grouped introns) and the
    # "genes not recovered" TSV; return per-(gene, sample) recovered paralog counts for heatmaps 07/08/09.
    #------------------------------------------------------------------------------------------------------------------
    recovered_dict = _apply_sample_threshold_and_write_fastas(
        outdir_non_chimeras, outdir_reports, filtered_dict, non_chimeric_introns_dict,
        min_threshold, min_samples_threshold, min_num_paralogs_per_sample, total_number_of_samples,
        intron_ref_coords_wiggle)

    #------------------------------------------------------------------------------------------------------------------
    # Write the 06 paralog-count TSV (used by heatmaps 07, 08, 09):
    #------------------------------------------------------------------------------------------------------------------
    outfile_non_chimeric_paralogs_report = _write_paralog_count_report(
        outdir_reports, min_samples_threshold, recovered_dict, sample_names_sorted, sample_names_joined)

    #------------------------------------------------------------------------------------------------------------------
    # Heatmap 07: raw paralog count per (gene, sample) after threshold filtering:
    #------------------------------------------------------------------------------------------------------------------
    outfile_heatmap_name = os.path.join(
        outdir_reports,
        f'07_non_chimeric_paralog_count_heatmap_{assembly_source}_after_sample_threshold_'
        f'{min_samples_threshold}_filtering.png')
    logger.info(f'{"[INFO]:":10} Writing heatmap: {outfile_heatmap_name}')
    df = _read_gene_sample_pivot(outfile_non_chimeric_paralogs_report, 'non_chimeric_paralog_count')

    cell_size_inches = 0.30
    fig, ax, fig_w, fig_h = _make_cell_sized_figure(df, cell_size_inches)
    sns.heatmap(df, vmin=0, cmap=matplotlib.colormaps['OrRd'].copy(), xticklabels=1, yticklabels=1,
                square=True, linewidth=0.5, cbar=False, ax=ax)
    ax.set_facecolor('lightgrey')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    cax = _attach_cell_sized_cbar(
        fig, ax, 'number of non-chimeric paralogs', cell_size_inches, integer_ticks=True)
    _add_cbar_swatches(fig, cax, [('lightgrey', 'no_seqs')], cell_size_inches, fig_w, fig_h)

    plt.savefig(outfile_heatmap_name, dpi=300, bbox_inches='tight')
    plt.close(fig)

    #------------------------------------------------------------------------------------------------------------------
    # Write a heatmap for sample vs gene showing number of non-chimeric paralogs after sample threshold filtering,
    # with any paralog count above non_chimeric_paralog_max_count standardised to the cap [STEPPED COLOURS]:
    #------------------------------------------------------------------------------------------------------------------
    cap = int(non_chimeric_paralog_max_count)
    outfile_heatmap_name = os.path.join(
        outdir_reports,
        f'08_non_chimeric_paralog_count_heatmap_{assembly_source}_after_sample_threshold_'
        f'{min_samples_threshold}_filtering_max_count_{cap}_stepped.png')
    logger.info(f'{"[INFO]:":10} Writing heatmap: {outfile_heatmap_name}')
    df = _non_chimeric_paralog_pivot_capped(outfile_non_chimeric_paralogs_report, cap)

    # Discrete cmap: tab10 (cap <= 10) or viridis sampled at `cap` points (cap > 10, because tab10
    # only has 10 categorical colours). Shift vmin/vmax by 0.5 so integer ticks land at band centres.
    if cap <= 10:
        bucket_colours = matplotlib.colormaps['tab10'].colors[:cap]
    else:
        logger.info(
            f'{"":10} cap > 10 ({cap}); sampling viridis at {cap} points instead of categorical tab10.')
        bucket_colours = matplotlib.colormaps['viridis'](np.linspace(0, 1, cap))
    cmap = matplotlib.colors.ListedColormap(bucket_colours)

    cell_size_inches = 0.30
    fig, ax, fig_w, fig_h = _make_cell_sized_figure(df, cell_size_inches)
    sns.heatmap(df, vmin=0.5, vmax=cap + 0.5, cmap=cmap, xticklabels=1, yticklabels=1, square=True,
                linewidth=0.5, cbar=False, ax=ax)
    ax.set_facecolor('lightgrey')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    cax = _attach_cell_sized_cbar(
        fig, ax, f'non-chimeric paralog count (capped at {cap})', cell_size_inches,
        ticks=range(1, cap + 1),
        ticklabels=[str(i) if i < cap else f'{cap}+' for i in range(1, cap + 1)],
    )
    _add_cbar_swatches(fig, cax, [('lightgrey', 'no_seqs')], cell_size_inches, fig_w, fig_h)

    plt.savefig(outfile_heatmap_name, dpi=300, bbox_inches='tight')
    plt.close(fig)

    #------------------------------------------------------------------------------------------------------------------
    # Write a heatmap for sample vs gene showing number of non-chimeric paralogs after sample threshold filtering,
    # with any paralog count above non_chimeric_paralog_max_count standardised to the cap [CONTINUOUS COLOURS]:
    #------------------------------------------------------------------------------------------------------------------
    outfile_heatmap_name = os.path.join(
        outdir_reports,
        f'09_non_chimeric_paralog_count_heatmap_{assembly_source}_after_sample_threshold_'
        f'{min_samples_threshold}_filtering_max_count_{cap}_continuous.png')
    logger.info(f'{"[INFO]:":10} Writing heatmap: {outfile_heatmap_name}')
    df = _non_chimeric_paralog_pivot_capped(outfile_non_chimeric_paralogs_report, cap)

    cell_size_inches = 0.30
    fig, ax, fig_w, fig_h = _make_cell_sized_figure(df, cell_size_inches)
    # viridis: perceptually uniform, dark-purple low end unambiguous vs the lightgrey no_seqs swatch.
    sns.heatmap(df, vmin=1, vmax=cap, cmap=matplotlib.colormaps['viridis'].copy(), xticklabels=1,
                yticklabels=1, square=True, linewidth=0.5, cbar=False, ax=ax)
    ax.set_facecolor('lightgrey')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    cax = _attach_cell_sized_cbar(
        fig, ax, f'non-chimeric paralog count (capped at {cap})', cell_size_inches, integer_ticks=True)
    _add_cbar_swatches(fig, cax, [('lightgrey', 'no_seqs')], cell_size_inches, fig_w, fig_h)

    plt.savefig(outfile_heatmap_name, dpi=300, bbox_inches='tight')
    plt.close(fig)

    #------------------------------------------------------------------------------------------------------------------
    # Write per-gene FASTAs containing every recovered paralog (chimeric and non-chimeric) for comparison:
    #------------------------------------------------------------------------------------------------------------------
    _write_all_seqs_fastas(outdir_all_seqs, unfiltered_seqs_dict)


def map_target_capture_sequences_to_genome(outdirs,
                                           log_queue,
                                           gene_list,
                                           reference_genome_fasta,
                                           assembly_source,
                                           sample_parent_folder,
                                           future_tracker,
                                           min_seq_length=75,
                                           min_length_percentage=0.80,
                                           min_contig_number_percentage=0.80,
                                           mappacbio_opts=None,
                                           pool=1,
                                           threads=1,
                                           ):
    """
    Map stitched assembly contigs to a reference genome using multiprocessing (Captus or HybPiper).

    Args:
        outdirs (dict): Dictionary of output directories: 'report_directory', 'log_directory', 'mapped_seqs', 'unmapped_seqs', 'non_chimeric_target_capture_seqs', 'all_target_capture_seqs'
        log_queue (Queue): Queue for logging messages
        gene_list (list): List of gene IDs to process
        reference_genome_fasta (str): Path to the reference genome FASTA file
        assembly_source (str): 'captus' or 'hybpiper'
        sample_parent_folder (str): Parent folder containing per-sample assembly output folders
        future_tracker (FutureTracker, optional): Future tracker for cleanup on interruption. Defaults to None.
        min_seq_length (int, optional): Minimum sequence length to map. Defaults to 75.
        min_length_percentage (float): Minimum percentage of the paralog sequence length remaining after
                                      filtering contig hits via <min_seq_length> for chimera detection to be performed.
                                      Otherwise, no sequence is retained for this sample/gene. Defaults to 0.80.
        min_contig_number_percentage (float): Minimum percentage of the total number of contig hits remaining for a
                                             given paralog after filtering contig hits via <min_seq_length> for
                                             chimera detection to be performed. Otherwise, no sequence is retained
                                             for this sample/gene. Defaults to 0.80.
        mappacbio_opts (dict): Tuning options for mapPacBio.sh. Contains keys  ``xmx``, ``k``, ``fastareadlen``, ``maxindel``, ``minid``.
        pool (int, optional): Number of processes to run concurrently. Defaults to 1.
        threads (int, optional): Number of threads to use for each concurrent process. Defaults to 1.
    """

    # Track wall-clock runtime for completion message
    start_time = time.time()

    current_time = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    logger.info(f'{"[INFO]:":10} Mapping {assembly_source} assembly seqs to genome {os.path.basename(reference_genome_fasta)} using MapPacBio.sh at {current_time}...')

    # Setup output directories:
    outdir_mapped_seqs = createfolder(outdirs['mapped_seqs'])
    outdir_unmapped_seqs = createfolder(outdirs['unmapped_seqs'])

    # Get input data
    input_data = sorted(glob.glob(f'{sample_parent_folder}/*'))
    total_records = len(input_data)

    # Initialize dictionaries to store results
    combined_samples_dict = dict()
    combined_intron_seqrecord_dict = dict()

    worker = map_captus_stitched_contigs if assembly_source == 'captus' else map_hybpiper_stitched_contigs

    with tqdm(total=total_records, desc=f"{'':10} Mapping samples sequences to genome with {pool} workers",
              unit="samples", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]') as progress_bar:

        with ProcessPool(max_workers=pool) as pool:

            # Register pool with future tracker for cleanup on interruption
            future_tracker.set_pool(pool)

            # Submit all batch jobs
            futures = {}
            for count, input_item in enumerate(input_data, 1):
                future = pool.schedule(worker,
                                       args=(input_item,
                                             log_queue,
                                             gene_list,
                                             min_seq_length,
                                             min_length_percentage,
                                             min_contig_number_percentage,
                                             outdirs['parent_output_directory'],
                                             outdir_mapped_seqs,
                                             outdir_unmapped_seqs,
                                             mappacbio_opts,
                                             threads),
                                            )

                futures[future] = (count, len(input_data))

            for future in as_completed(futures):

        

                count, total_records = futures[future]
                progress_bar.update(1)

                success, result = future.result()

                if success:
                    sample_name, sample_dict, sample_intron_seqrecord_dict = result
                    combined_samples_dict[sample_name] = sample_dict
                    combined_intron_seqrecord_dict[sample_name] = sample_intron_seqrecord_dict

                else:
                    logger.error(f'{"[ERROR]:":10} Failed to with item {count}: {result[0]}. See log file for full traceback.')
                    logger.debug(f'{"":10} Full traceback: {result[1]}')

            # Clear pool reference after normal completion
            future_tracker.set_pool(None)
    
    # Log completion
    log_completion_time(start_time, logger, label="map_target_capture_sequences_to_genome completed")

    return combined_samples_dict, combined_intron_seqrecord_dict

def map_captus_stitched_contigs(sample_folder,
                                log_queue,      
                                gene_list,
                                min_seq_length,
                                min_length_percentage,
                                min_contig_number_percentage,
                                parent_output_directory,
                                outdir_mapped_seqs,
                                outdir_unmapped_seqs,
                                mappacbio_opts,
                                threads):
    """
    Process one Captus sample: build exon/intron records, then run shared chimera mapping.

    Returns:
        tuple: (sample_name, sample_dict, gene_intron_seq_records_dict picklable plain dict)
    """

    try:
        worker_logger = setup_worker_logger(__name__, log_queue)

        worker_logger.debug(f'{"[DEBUG]:":10} Mapping Captus sample: {sample_folder}...')

        # Build gene seq and intron records
        sample_name, gene_seq_records_dict, gene_intron_seq_records_dict = (
            captus_build_gene_seq_and_intron_records(sample_folder, gene_list))

        sample_dict, plain_introns = map_sample_stitched_contigs_core(
            worker_logger,
            sample_name,
            gene_seq_records_dict,
            gene_intron_seq_records_dict,
            min_seq_length,
            min_length_percentage,
            min_contig_number_percentage,
            parent_output_directory,
            outdir_mapped_seqs,
            outdir_unmapped_seqs,
            mappacbio_opts,
            threads=threads)

        return True, (sample_name, sample_dict, plain_introns)

    except Exception as error:
        return False, (error, traceback.format_exc())


def hybpiper_sample_name_from_folder(sample_folder):
    """TODO(HybPiper): Derive sample id from HybPiper directory layout."""
    return os.path.basename(sample_folder.rstrip(os.sep))


def hybpiper_build_gene_seq_and_intron_dicts_stub(sample_folder, gene_list):
    """
    TODO(HybPiper): Populate ``gene_seq_records_dict`` and intron records like ``captus_build_gene_seq_and_intron_records``.

    Returns:
        tuple: (sample_name, gene_seq_records_dict, gene_intron_seq_records_dict)

    Each maps gene -> ``paralog_N`` -> list of exon/coding SeqRecords (one entry per separately mapped contig
    fragment), and gene -> ``paralog_N`` -> list of intron SeqRecords (optional).
    """
    sample_name = hybpiper_sample_name_from_folder(sample_folder)
    logger.warning(f'{"[WARNING]:":10} HybPiper sequence extraction is not implemented; skipping sample {sample_name!r}.')
    return sample_name, {}, {}


def captus_build_gene_seq_and_intron_records(sample_folder, gene_list):
    """
    From Captus ``recovery_stats.tsv`` and MEGAHIT hit contigs, build per-paralog exon fragments and introns.

    Returns:
        tuple: (sample_name, gene_seq_records_dict, gene_intron_seq_records_dict)
    """
    sample_basename = os.path.basename(sample_folder)
    sample_name = sample_basename.split('__')[0]

    #------------------------------------------------------------------------------------------------------------------
    # Read in the recovery stats file to get MEGAHIT contig coordinates and seqs, and reference protein coordinates:
    #------------------------------------------------------------------------------------------------------------------
    gene_stats_dict = defaultdict(list)
    with open(f'{sample_folder}/06_assembly_annotated/{sample_name}_recovery_stats.tsv') as recovery_stats_handle:
        lines = iter(recovery_stats_handle.readlines())
        for line in lines:
            if line.split('\t')[1] != 'NUC':
                continue  # bit of a hack to skip the variable number of header lines for Captus v1.6.1 vs earlier versions
            gene = line.split('\t')[2]

            if gene not in gene_list:
                continue

            megahit_contigs = [item.strip() for item in line.split('\t')[22].split(';')]
            megahit_contigs_strands = [item.strip() for item in line.split('\t')[23].split(';')]
            megahit_contigs_coords = [item.strip() for item in line.split('\t')[24].split(';')]
            ref_protein_coords = [item.strip() for item in line.split('\t')[4].split(';')]

            zipped_stats = list(zip(megahit_contigs,
                                    megahit_contigs_strands,
                                    megahit_contigs_coords,
                                    ref_protein_coords))
            gene_stats_dict[gene].append(zipped_stats)

    all_raw_megahit_contigs_dict = SeqIO.to_dict(SeqIO.parse(f'{sample_folder}/06_assembly_annotated/'
                                                             f'{sample_name}_hit_contigs.fasta', 'fasta'))

    #------------------------------------------------------------------------------------------------------------------
    # Build a dictionary of constituent coding seqrecords for each sample, for each gene/paralog seq (i.e. one 
    # seqrecord for each MEGAHIT contig used to build a given paralog sequence):
    #------------------------------------------------------------------------------------------------------------------
    gene_seq_records_dict = defaultdict(lambda: defaultdict(list))

    for gene, paralog_stats_list in gene_stats_dict.items():
        for count, paralog_stats in enumerate(paralog_stats_list, 1):

            for node_stats_tuple in paralog_stats:
                node_name, strand, coords, ref_protein_coords = node_stats_tuple
                coords_list = coords.split(',')
                ref_protein_coords_list = ref_protein_coords.split(',')
                raw_megahit_contig = all_raw_megahit_contigs_dict[node_name]

                if strand == '+':
                    exon_regions = ''

                    for coord, ref_protein_coord in zip(coords_list, ref_protein_coords_list):
                        start, end = coord.split('-')
                        exon_seq = raw_megahit_contig.seq[int(start) - 1: int(end)]
                        exon_regions = exon_regions + exon_seq

                    seqrecord = SeqRecord(seq=exon_regions,
                                          name=raw_megahit_contig.name,
                                          id=raw_megahit_contig.id,
                                          description=raw_megahit_contig.description)

                    gene_seq_records_dict[gene][f'paralog_{count}'].append(seqrecord)

                else:
                    exon_regions = ''

                    for coord in coords_list:
                        start, end = coord.split('-')
                        exon_seq = raw_megahit_contig.seq[int(start) - 1: int(end)]
                        exon_seq = exon_seq.reverse_complement()
                        exon_regions = exon_regions + exon_seq

                    seqrecord = SeqRecord(seq=exon_regions,
                                          name=raw_megahit_contig.name,
                                          id=raw_megahit_contig.id,
                                          description=raw_megahit_contig.description)

                    gene_seq_records_dict[gene][f'paralog_{count}'].append(seqrecord)

    #------------------------------------------------------------------------------------------------------------------
    # Build a dictionary of intron seqrecords for each sample, for each gene/paralog seq (i.e. a list of seqrecords, 
    # one for each intron found in the MEGAHIT contigs used to build a given paralog sequence):
    #------------------------------------------------------------------------------------------------------------------
    gene_intron_seq_records_dict = defaultdict(lambda: defaultdict(list))

    for gene, paralog_stats_list in gene_stats_dict.items():
        for count, paralog_stats in enumerate(paralog_stats_list, 1):

            for node_stats_tuple in paralog_stats:
                node_name, strand, coords, ref_protein_coords = node_stats_tuple
                coords_list = coords.split(',')
                ref_protein_coords_list = ref_protein_coords.split(',')
                raw_megahit_contig = all_raw_megahit_contigs_dict[node_name]

                coord_pairs = []
                for coord, ref_protein_coord in zip(coords_list, ref_protein_coords_list):
                    if not coord:
                        continue
                    start_str, end_str = coord.split('-')
                    ref_protein_start_str, ref_protein_end_str = ref_protein_coord.split('-')
                    coord_pairs.append((int(start_str), int(end_str), int(ref_protein_start_str), int(ref_protein_end_str)))

                if len(coord_pairs) > 1:
                    prefix = f'{sample_name}__{gene}__{count - 1:02d}'
                    coord_pairs_sorted = sorted(coord_pairs, key=lambda x: x[0])
                    for intron_idx, ((start1, end1, ref_start1, ref_end1), (start2, end2, ref_start2, ref_end2)) in enumerate(
                            zip(coord_pairs_sorted, coord_pairs_sorted[1:])):
                        intron_start = end1
                        intron_end = start2 - 1
                        if intron_end <= intron_start:
                            continue
                        intron_seq = raw_megahit_contig.seq[intron_start: intron_end]
                        if strand == '-':
                            intron_seq = intron_seq.reverse_complement()

                        description = raw_megahit_contig.description
                        intron_ref_start = min(ref_end1, ref_end2)
                        intron_ref_end = max(ref_start1, ref_start2)
                        ref_coord_str = f'{intron_ref_start}-{intron_ref_end}'
                        description = f'{description} ref_protein_coords:{ref_coord_str}'

                        intron_id = f'{prefix}_{raw_megahit_contig.id}_i{intron_idx}'
                        intron_seqrecord = SeqRecord(seq=intron_seq,
                                                     name=intron_id,
                                                     id=intron_id,
                                                     description=description)
                        gene_intron_seq_records_dict[gene][f'paralog_{count}'].append(intron_seqrecord)

    return sample_name, gene_seq_records_dict, gene_intron_seq_records_dict


def map_hybpiper_stitched_contigs(sample_folder,
                                  log_queue,
                                  gene_list,
                                  reference_genome_fasta,
                                  assembly_source,
                                  sample_parent_folder,
                                  min_seq_length,
                                  min_length_percentage,
                                  min_contig_number_percentage,
                                  parent_output_directory,
                                  outdir_mapped_seqs,
                                  outdir_unmapped_seqs,
                                  bbmap_Xmx,
                                  mappacbio_opts,
                                  threads):
    """Process one HybPiper sample; same worker contract as ``map_captus_stitched_contigs``."""
    try:
        worker_logger = setup_worker_logger(__name__, log_queue)
        worker_logger.debug(f'{"[DEBUG]:":10} Mapping HybPiper sample: {sample_folder}...')

        sample_name, gseq, gint = hybpiper_build_gene_seq_and_intron_dicts_stub(sample_folder, gene_list)

        gene_seq_records_dict = defaultdict(lambda: defaultdict(list))
        for gene, paralog_map in (gseq or {}).items():
            for paralog, recs in paralog_map.items():
                if isinstance(recs, list):
                    gene_seq_records_dict[gene][paralog].extend(recs)
                else:
                    gene_seq_records_dict[gene][paralog].append(recs)

        gene_intron_seq_records_dict = defaultdict(lambda: defaultdict(list))
        for gene, paralog_map in (gint or {}).items():
            for paralog, recs in paralog_map.items():
                if isinstance(recs, list):
                    gene_intron_seq_records_dict[gene][paralog].extend(recs)
                else:
                    gene_intron_seq_records_dict[gene][paralog].append(recs)

        sample_dict, plain_introns = map_sample_stitched_contigs_core(
            worker_logger,
            sample_name,
            gene_seq_records_dict,
            gene_intron_seq_records_dict,
            min_seq_length,
            min_length_percentage,
            min_contig_number_percentage,
            parent_output_directory,
            outdir_mapped_seqs,
            outdir_unmapped_seqs,
            bbmap_Xmx,
            mappacbio_opts,
            threads=threads)
        return True, (sample_name, sample_dict, plain_introns)

    except Exception as error:
        return False, (error, traceback.format_exc())


def _paralog_fragment_stats(node_list, min_seq_length):
    """Return (n_fragments, total_bp, n_fragments_ge_min, sum_bp_ge_min).

    Call after ``_strip_ns_inplace_on_fragments`` so lengths exclude N padding.
    """
    lengths = [len(n) for n in node_list]
    n_frag = len(lengths)
    total_bp = sum(lengths)
    ok = [L for L in lengths if L >= min_seq_length]
    n_ok = len(ok)
    sum_ok = sum(ok)
    return n_frag, total_bp, n_ok, sum_ok


def _strip_ns_inplace_on_fragments(node_list):
    """Remove N from each exon fragment in place (before length gates and FASTA output)."""
    for seq in node_list:
        seq.seq = seq.seq.replace('N', '')


def _write_unmapped_paralog_fasta(sample_name, gene, paralog, exon_records, intron_records, status, outdir_unmapped):
    """Write per-paralog ``*_scipio_hits.fasta`` (+ introns) for paralogs that bypass genome mapping.

    Files land under ``{outdir_unmapped}/{status}/{sample_name}/{gene}/{paralog}_scipio_hits.fasta`` so
    each chimera-status bucket can be inspected independently. ``exon_records`` may be a single
    ``SeqRecord`` or a list/iterable of ``SeqRecord`` (``SeqIO.write`` accepts both).
    """
    status_dir = createfolder(f'{outdir_unmapped}/{status}/{sample_name}/{gene}')
    exon_path = f'{status_dir}/{paralog}_scipio_hits.fasta'
    with open(exon_path, 'w') as fh:
        SeqIO.write(exon_records, fh, 'fasta')
    if intron_records:
        intron_path = f'{status_dir}/{paralog}_scipio_hits_introns.fasta'
        with open(intron_path, 'w') as fh:
            SeqIO.write(intron_records, fh, 'fasta')


def _filter_fragments_by_min_length(node_list, min_seq_length, worker_logger, sample_name, gene, paralog):
    """Return fragments with length >= min_seq_length (call after N-stripping)."""
    filtered = []
    for seq in node_list:
        if len(seq) < min_seq_length:
            worker_logger.debug(
                f'Sample {sample_name} gene {gene} paralog {paralog} sequence {seq.name} is less than '
                f'min_seq_length = {min_seq_length} bases long, sequence will not be written or mapped.')
        else:
            filtered.append(seq)
    return filtered


def _run_map_pac_bio(outfile_gene, expected_alignment, threads, parent_output_directory, worker_logger,
                     mappacbio_opts):
    """Run mapPacBio.sh unless a non-empty SAM already exists.

    Args:
        outfile_gene: Per-paralog query FASTA passed as ``in=``.
        expected_alignment: SAM path passed as ``out=`` and checked for non-empty completion.
        threads: Number of mapper threads (``threads=``).
        parent_output_directory: Reference index location (``path=``).
        worker_logger: Logger for this worker process.
        bbmap_Xmx: BBMap/Java heap in GB (``Xmx=...g``).
        mappacbio_opts (dict): Tuning options for ``mapPacBio.sh``. Must contain
            ``k``, ``fastareadlen``, ``maxindel``, ``minid``, ``local``, ``ambiguous`` (populated in
            ``main`` from the ``--mappacbio_*`` CLI flags). ``noheader=t`` is fixed because the
            downstream SAM parser splits on tabs and assumes headerless rows.
    """
    if file_exists_and_not_empty(expected_alignment):
        worker_logger.debug(f'Expected file {expected_alignment} exists, skipping...')
        return

    command = (
        f'mapPacBio.sh '
        f'Xmx={mappacbio_opts["xmx"]}g '
        f'in={outfile_gene} '
        f'k={mappacbio_opts["k"]} '
        f'fastareadlen={mappacbio_opts["fastareadlen"]} '
        f'maxindel={mappacbio_opts["maxindel"]} '
        f'minid={mappacbio_opts["minid"]} '
        f'local=true '
        f'ambiguous=all '
        f'noheader=t '
        f'threads={threads} '
        f'path={parent_output_directory} '
        f'out={expected_alignment}')
    try:
        result = subprocess.run(command,
                                universal_newlines=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                check=True,
                                shell=True)
        worker_logger.debug(f'mapPacBio.sh check_returncode() is: {result.check_returncode()}')
        worker_logger.debug(f'mapPacBio.sh stdout is: {result.stdout}')
        worker_logger.debug(f'mapPacBio.sh stderr is: {result.stderr}')

    except subprocess.CalledProcessError as exc:
        worker_logger.error(f'mapPacBio.sh FAILED. Output is: {exc}')
        worker_logger.error(f'mapPacBio.sh stdout is: {exc.stdout}')
        worker_logger.error(f'mapPacBio.sh stderr is: {exc.stderr}')
        raise ValueError('There was an issue running mapPacBio.sh. Check input files!') from exc

    if not file_exists_and_not_empty(expected_alignment):
        raise ValueError(f'mapPacBio.sh did not produce expected alignment: {expected_alignment}')


def _multihit_status_from_sam(worker_logger, sample_name, gene, paralog,
                             expected_alignment, outfile_gene):
    """
    Parse SAM (query -> reference names) and classify stitched fragments.

    Returns:
        One of: chimera_from_multi_hit, no_chimera_from_multi_hit, unknown_repeated_subject
    """
    query_to_subjects = defaultdict(list)
    with open(expected_alignment) as sam_handle:
        for line in sam_handle:
            stats = line.split('\t')
            query_name = stats[0].split(',')[0]
            subject_name = stats[2]
            query_to_subjects[query_name].append(subject_name)

    query_names_in_fasta = [seq.name for seq in SeqIO.parse(outfile_gene, 'fasta')]
    query_counts = Counter(query_names_in_fasta)

    for name in query_names_in_fasta:
        if name not in query_to_subjects:
            worker_logger.warning(
                f'{"[WARNING]:":10} No mapping result in file {expected_alignment} for query {name}.')
            raise AssertionError(f'No SAM rows for query {name!r}')

    sorted_by_n_subjects = dict(sorted(query_to_subjects.items(), key=lambda x: len(x[1])))
    first = dict(itertools.islice(sorted_by_n_subjects.items(), 1))
    rest = dict(itertools.islice(sorted_by_n_subjects.items(), 1, len(sorted_by_n_subjects) + 1))

    chimera_found = False
    unknown_repeated_subject = False

    for query, subject_list in first.items():
        query_count = query_counts[query]
        if query_count == 0:
            raise AssertionError(f'query_count for {query!r}')

        if query_count > 1:
            worker_logger.debug(f'query_count: {query_count}')

        for subject_contig_name, cnt in Counter(subject_list).items():
            if cnt > query_count:
                worker_logger.debug(
                    f'{"[WARNING]:":10} Sample {sample_name}, gene {gene}: subject {subject_contig_name} '
                    f'occurs more ({cnt}) than the number of times the query {query} occurs in the '
                    f'"{paralog}_scipio_hits.fasta" file ({query_count}). Sequence will be flagged as '
                    f'"unknown_repeated_subject".')
                unknown_repeated_subject = True

        for subject in subject_list:
            for remaining_query, remaining_subject_list in rest.items():
                remaining_query_count = query_counts[remaining_query]
                if remaining_query_count == 0:
                    raise AssertionError(f'remaining_query_count for {remaining_query!r}')
                if remaining_query_count > 1:
                    worker_logger.debug(f'remaining_query_count: {remaining_query_count}')
                rem_counter = Counter(remaining_subject_list)

                if subject not in remaining_subject_list:
                    chimera_found = True
                else:
                    rem_subj_count = rem_counter[subject]
                    if rem_subj_count == 0:
                        raise AssertionError(f'subject count for {subject!r}')
                    if rem_subj_count > remaining_query_count:
                        worker_logger.debug(
                            f'{"[WARNING]:":10} Sample {sample_name}, gene {gene}: subject {subject} '
                            f'occurs more ({rem_subj_count}) than the number of times the query {remaining_query} '
                            f'occurs in the "{paralog}_scipio_hits.fasta" file ({remaining_query_count}). Sequence '
                            f'will be flagged as "unknown_repeated_subject".')
                        unknown_repeated_subject = True

    if unknown_repeated_subject:
        return 'unknown_repeated_subject'
    if chimera_found:
        return 'chimera_from_multi_hit'
    return 'no_chimera_from_multi_hit'


def map_sample_stitched_contigs_core(worker_logger,
                                     sample_name,
                                     gene_seq_records_dict,
                                     gene_intron_seq_records_dict,
                                     min_seq_length,
                                     min_length_percentage,
                                     min_contig_number_percentage,
                                     parent_output_directory,
                                     outdir,
                                     outdir_unmapped,
                                     mappacbio_opts,
                                     threads=1):
    """
    Shared chimera logic: strip Ns from fragments, apply length gates, mapPacBio when needed, classify multi-hit.

    Nucleotide ``N`` is removed from each exon fragment **before** length-based gates
    (``_paralog_fragment_stats`` / ``min_length_percentage`` / ``min_contig_number_percentage``) so those
    thresholds apply to effective sequence length.

    Args:
        worker_logger: Logger for this worker process (queue handler).
        sample_name: Short sample label for output paths.
        gene_seq_records_dict: gene -> paralog -> list of SeqRecord (one per separately mapped fragment).
        gene_intron_seq_records_dict: gene -> paralog -> list of intron SeqRecords.
        parent_output_directory: Passed to mapPacBio ``path=`` (reference index location).
        outdir: Base directory for per-sample mapping outputs (mapped paralogs).
        outdir_unmapped: Base directory for paralogs that bypass mapping; ``*_scipio_hits.fasta`` (+ introns)
            are written under ``{outdir_unmapped}/{status}/{sample_name}/{gene}/`` per chimera-status bucket.
        mappacbio_opts (dict): Tuning options forwarded to ``_run_map_pac_bio``. 
        threads: Number of threads to use for mapping.
    """
    sample_dict = {}

    for gene, paralogs_dict in gene_seq_records_dict.items():
        for paralog, node_list in paralogs_dict.items():
            introns = gene_intron_seq_records_dict.get(gene, {}).get(paralog, [])

            _strip_ns_inplace_on_fragments(node_list)
            n_frag, total_bp, n_ok, sum_ok = _paralog_fragment_stats(node_list, min_seq_length)
            passes_length = (sum_ok / total_bp) >= min_length_percentage
            passes_contig = (n_ok / n_frag) >= min_contig_number_percentage

            pd = sample_dict.setdefault(gene, {})

            if n_frag == 1:
                worker_logger.debug(
                    f'Only one seq for gene {gene} paralog {paralog} sample {sample_name}, so no chimera from '
                    f'stitching.')
                status = 'no_chimera_as_single_contig_hit'
                pd[paralog] = status
                _write_unmapped_paralog_fasta(
                    sample_name, gene, paralog, [node_list[0]], introns, status, outdir_unmapped)
                continue

            if n_frag > 1 and n_ok == 1:
                worker_logger.debug(
                    f'Only one seq for gene {gene} paralog {paralog} sample {sample_name} after filtering out '
                    f'seqs < {min_seq_length}, so no _detectable_ chimera from stitching.')
                status = 'no_detectable_chimera_as_single_contig_hit_after_length_filtering'
                pd[paralog] = status
                _write_unmapped_paralog_fasta(
                    sample_name, gene, paralog, node_list, introns, status, outdir_unmapped)
                continue

            if n_frag > 1 and n_ok == 0:
                worker_logger.debug(
                    f'No seqs left for gene {gene} paralog {paralog} sample {sample_name} after filtering out '
                    f'seqs < {min_seq_length}, so no _detectable_ chimera from stitching.')
                status = 'no_detectable_chimera_as_no_seqs_left_after_length_filtering'
                pd[paralog] = status
                _write_unmapped_paralog_fasta(
                    sample_name, gene, paralog, node_list, introns, status, outdir_unmapped)
                continue

            if not passes_length:
                worker_logger.debug(
                    f'Too little of the paralog sequence remaining for gene {gene} paralog {paralog} sample '
                    f'{sample_name} after filtering out seqs < {min_seq_length}, so no chimera detection performed.')
                status = 'no_chimera_test_performed_as_too_little_seq_length_after_length_filtering'
                pd[paralog] = status
                _write_unmapped_paralog_fasta(
                    sample_name, gene, paralog, node_list, introns, status, outdir_unmapped)
                continue

            if not passes_contig:
                worker_logger.debug(
                    f'Too few contigs remaining for gene {gene} paralog {paralog} sample {sample_name} after '
                    f'filtering out seqs < {min_seq_length}, so no chimera detection performed.')
                status = 'no_chimera_test_performed_as_too_few_contigs_after_length_filtering'
                pd[paralog] = status
                _write_unmapped_paralog_fasta(
                    sample_name, gene, paralog, node_list, introns, status, outdir_unmapped)
                continue

            # Multi-fragment: write FASTA, map, classify from SAM
            gene_dir = createfolder(f'{outdir}/{sample_name}/{gene}')
            outfile_gene = f'{gene_dir}/{paralog}_scipio_hits.fasta'
            seqs_filtered = _filter_fragments_by_min_length(
                node_list, min_seq_length, worker_logger, sample_name, gene, paralog)

            with open(outfile_gene, 'w') as out_fh:
                SeqIO.write(seqs_filtered, out_fh, 'fasta')

            if introns:
                intron_path = f'{gene_dir}/{paralog}_scipio_hits_introns.fasta'
                with open(intron_path, 'w') as out_fh:
                    SeqIO.write(introns, out_fh, 'fasta')

            expected_sam = f'{outdir}/{sample_name}/{gene}/{paralog}_mapping.sam'
            _run_map_pac_bio(outfile_gene, expected_sam, threads, parent_output_directory, worker_logger, mappacbio_opts)

            pd[paralog] = _multihit_status_from_sam(
                worker_logger, sample_name, gene, paralog, expected_sam, outfile_gene)

    gene_intron_plain = {
        g: {p: list(recs) for p, recs in par_dict.items()}
        for g, par_dict in gene_intron_seq_records_dict.items()
    }

    return sample_dict, gene_intron_plain


########################################################################################################################
# argparse options:
########################################################################################################################
def parse_arguments():
    """
    Parse command line arguments for the script.

    This function sets up the command line argument parser and defines all the
    required and optional arguments for the script.

    Returns:
        argparse.Namespace: An object containing the parsed command line arguments
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', '-v',
                        action='version',
                        version=f'filter_chimeras.py version {__version__}',
                        help='Print the script version number.')
    asm_group = parser.add_mutually_exclusive_group(required=True)
    asm_group.add_argument('--captus_folder',
                           type=str,
                           default=None,
                           help='Parent folder containing Captus sample output folders.')
    asm_group.add_argument('--hybpiper_folder',
                           type=str,
                           default=None,
                           help='Parent folder containing HybPiper sample output folders.')
    parser.add_argument('reference_genome_fasta',
                        type=str,
                        help='High quality reference genome for mapping target-capture contigs.')
    parser.add_argument('target_file_fasta',
                        type=str,
                        help='Target FASTA used for the assembly runs (Captus or HybPiper).')
    parser.add_argument('--output_directory',
                        type=str,
                        default='output_directory',
                        help='Root folder for logs, reports, mapped intermediates, and filtered FASTAs. Default is: '
                             '%(default)s')
    parser.add_argument('--min_seq_length',
                        type=int,
                        default=75,
                        help='Minimum length for a contig-derived sequence to be mapped. Default is: %(default)s')
    parser.add_argument('--min_length_percentage',
                        type=float,
                        default=0.80,
                        help='The minimum percentage of the paralog sequence length retained after '
                             'filtering contig hits via <min_seq_length> for mapping and chimera detection to be '
                             'performed. Default is: %(default)s')
    parser.add_argument('--min_contig_number_percentage',
                        type=float,
                        default=0.80,
                        help='The minimum percentage of the total number of contig hits remaining for a given paralog '
                             'after filtering contig hits via <min_seq_length> for mapping and chimera detection to be '
                             'performed. Default is: %(default)s')
    parser.add_argument('--min_samples_threshold',
                        type=float,
                        default=0.75,
                        help='For a given gene, the minimum percentage of total samples to have '
                             '>= <min_num_paralogs_per_sample> non-chimeric sequences for sequences to be written to '
                             'file. Useful to increase stringency of alignment sample occupancy, at the cost of fewer '
                             'gene alignments. Default is: %(default)s')
    parser.add_argument('--min_num_paralogs_per_sample',
                        type=int,
                        default=0,
                        help='For a given gene for a given sample, the minimum number of non-chimeric paralog '
                             'sequences recovered for the sequences to be written to file. This can be useful when you '
                             'expect paralogs to be present (e.g. due to polyploidy), and want to skip genes that '
                             'might have hidden paralogy due to sequence or assembly issues. Default is: %(default)s')
    parser.add_argument('--non_chimeric_paralog_max_count',
                        type=int,
                        default=10,
                        help='Cap value used by the discrete non-chimeric paralog count heatmaps. Per-(gene, sample) '
                             'non-chimeric paralog counts above this value are collapsed into a single "<cap>+" bucket '
                             'so the categorical legend stays readable. Values <= 10 use the categorical tab10 palette; '
                             'values > 10 sample the viridis sequential colormap at that many points. Default is: '
                             '%(default)s')
    parser.add_argument('--pool',
                        type=int,
                        default=1,
                        help='Number of samples to run concurrently. Default is: %(default)s')
    parser.add_argument('--threads',
                        type=int,
                        default=1,
                        help='Number of threads to use for mapping for a given sample. Default is: %(default)s')
    parser.add_argument('--intron_ref_coords_wiggle',
                        type=int,
                        default=-1,
                        help='When writing non-chimeric intron FASTAs under introns/, merge ref_protein_coords groups '
                             'whose query intervals overlap or touch after expanding each end by this many aa. '
                             'Use 0 to merge adjacent junctions (e.g. 82-83 with 83-84). Use -1 (default) for legacy '
                             'behaviour: one file per exact coord string only. Merged groups include every coord label '
                             'in the filename, sorted (e.g. gene_ref_protein_coords_82-83_83-84.fasta).')
    parser.add_argument('--bbmap_Xmx',
                        type=str,
                        default=30,
                        help='Amount of memory to allocate to BBmap.sh in GB. Default is: %(default)s')

    mappacbio_group = parser.add_argument_group('mapPacBio.sh tuning options')
    mappacbio_group.add_argument('--mappacbio_k',
                                 type=int,
                                 default=13,
                                 help='mapPacBio.sh k-mer size (k=). Default is: %(default)s')
    mappacbio_group.add_argument('--mappacbio_fastareadlen',
                                 type=int,
                                 default=6000,
                                 help='mapPacBio.sh fastareadlen= (max FASTA read length). Default is: %(default)s')
    mappacbio_group.add_argument('--mappacbio_maxindel',
                                 type=int,
                                 default=100000,
                                 help='mapPacBio.sh maxindel= (max indel length tolerated). Default is: %(default)s')
    mappacbio_group.add_argument('--mappacbio_minid',
                                 type=float,
                                 default=0.76,
                                 help='mapPacBio.sh minid= (minimum identity to keep an alignment). '
                                      'Default is: %(default)s')

    results = parser.parse_args()

    return results


########################################################################################################################
def main():
    args = parse_arguments()

    # Track wall-clock runtime for completion message
    start_time = time.time()

    # Create future tracker for cleanup
    future_tracker = FutureTracker()

    try:
        global logger, log_queue, log_listener

        # Set up log and report directories
        setup_log_and_report_directories(args)

        # Set up global logger
        logger, log_queue, log_listener = log_manager.setup(
            __name__, 'filter_chimeras', log_directory=args.log_directory
        )

        # Print arguments to screen and log:
        print_arguments(args, logger, __version__)

        # Set up output directories:
        outdirs = _setup_output_directories(args)

        # Create genome mapping reference
        create_genome_mapping_reference(args)

        # Parse target file
        gene_list = parse_target_file(args.target_file_fasta)

        # Get assembly source and sample parent folder
        if args.captus_folder:
            assembly_source = 'captus'
            sample_parent_folder = args.captus_folder
        else:
            assembly_source = 'hybpiper'
            sample_parent_folder = args.hybpiper_folder

        # Map target capture sequences to genome
        combined_samples_dict, combined_intron_seqrecord_dict = map_target_capture_sequences_to_genome(
                                               outdirs,
                                               log_queue,
                                               gene_list, 
                                               args.reference_genome_fasta,
                                               assembly_source,
                                               sample_parent_folder,
                                               future_tracker,
                                               args.min_seq_length,
                                               args.min_length_percentage,
                                               args.min_contig_number_percentage,
                                               {'xmx': args.bbmap_Xmx,
                                                'k': args.mappacbio_k,
                                                'fastareadlen': args.mappacbio_fastareadlen,
                                                'maxindel': args.mappacbio_maxindel,
                                                'minid': args.mappacbio_minid,
                                                },
                                               args.pool,
                                               args.threads,
                                               )

        # Write non-chimeric sequences and reports to file:
        write_sequences_and_reports(
            outdirs,
            assembly_source,
            sample_parent_folder,
            combined_samples_dict,
            combined_intron_seqrecord_dict,
            args.min_samples_threshold,
            args.min_num_paralogs_per_sample,
            args.intron_ref_coords_wiggle,
            args.non_chimeric_paralog_max_count,
        )

    except KeyboardInterrupt:
        print("\n" + "=" * 80)
        print("INTERRUPTED: User requested termination (Ctrl+C)")  
        print("=" * 80)

        if 'logger' in globals() and logger:
            logger.info(f'{"[INFO]:":10} Pipeline interrupted by user (KeyboardInterrupt)')
            log_separator(logger)
            log_completion_time(start_time, logger, label=__file__ + " INTERRUPTED")

        # Cancel any active futures
        future_tracker.cancel_all()

        print(f'\n{" " * 10}Cleaning up processes and exiting...')
        log_manager.cleanup()
        sys.exit(0)

    except Exception as e:
        log_manager.handle_error(e, traceback.format_exc(), "main()")

    finally:
        # Log total completion time before cleaning up the logger
        log_separator(logger)
        log_completion_time(start_time, logger if ('logger' in globals() and logger) else None,
                                  label=__file__ + " completed")

        log_manager.cleanup()

########################################################################################################################
if __name__ == '__main__':
    main()
