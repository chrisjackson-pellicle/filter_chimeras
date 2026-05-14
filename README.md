# filter_chimeras — reference-guided chimera filtering for Captus or HybPiper

`filter_chimeras.py` screens target-capture sequences assembled by [`Captus`](https://github.com/edgardomortiz/captus) or [`HybPiper`](https://github.com/mossmatters/HybPiper) for chimeric stitched genes by mapping each constituent stitched sequence to a reference genome and testing whether they map to the same reference sequence(s) (chromosome-level contig, scaffold, or haplotigs).

---

## Why this exists

Phylogenetic reconstruction using sequences assembled from target-capture data often assumes that each sequence represents a single genomic locus. However, when assembling sequences from target-capture data for multi-exonic genes, `Captus` and `HybPiper` often stitch together sequences extracted from multiple contigs; these contigs usually represent individual exons, and are assembled by MEGAHIT (`Captus`) or SPAdes (`HybPiper`). This stitching is appropriate when all contigs truly come from one locus, but this is not always the case. For genes with multiple paralogs, exons from different paralogs can be artificially stitched together. The result is a chimera: a sequence that looks like one gene but is not a single orthologous locus. This can occur due to e.g. differential recovery of exons from different paralogs during target-capture, variable assembly of exons from different paralogs during assembly of target-capture data, or variable sequence identify between exons from different paralogs and sequences in the target-file used for sequence extraction. 

Such chimeras often cannot be inferred from the primary FASTA sequence alone. This script uses an independent reference assembly: if stitched-together sequences derived from different MEGAHIT/SPAdes contigs for a given target-capture-assembled sequence map to different reference targets, the sequence is flagged as a likely chimera.

---


## Processing pipeline

Running `filter_chimeras.py` does the following in order.

**1. Logging and reference index**  
Runs `bbmap.sh` to create a reference genome indexed under `<output_directory>/ref/genome/1/`.

**2. Gene list from the target FASTA**  
Parses the target file FASTA used in the target-capture assembly to recover a list of gene IDs to process.

**3. Per-sample chimera screen**  
For each folder under `--captus_folder` or `--hybpiper_folder`:

- **Reconstruct** stitched coding fragments (and introns where applicable) from that sample’s assembly outputs — for `Captus`, this data is extracted from `06_assembly_annotated/<sample>_recovery_stats.tsv` plus `<sample>_hit_contigs.fasta`; `HybPiper` uses XXX.


- **Gate** each stitched gene/paralog with `--min_seq_length`, `--min_length_percentage`, and `--min_contig_number_percentage`. Failing the gate skips full multi-hit testing and the sequence is flagged as 'not tested' rather than chimera vs non-chimera from mapping. See [Figure_1](https://github.com/chrisjackson-pellicle/filter_chimeras/wiki/Pre%E2%80%90map-filtering-of-stitched-and-non%E2%80%90stiched-sequences).


- **Map** eligible multi-fragment paralogs with `mapPacBio.sh` against the indexed reference; read the SAM and compare reference sequence names across hits. Agreement → non-chimera; disagreement → chimera; ambiguous repeated subjects → `unknown_repeated_subject` (and related paths). See [Figure_2](https://github.com/chrisjackson-pellicle/filter_chimeras/wiki/Mapping-and-chimera-detection).


- **Write** per-sample intermediates for **mapped** paralogs (multi-fragment query FASTAs, intron FASTAs, SAMs) to `01_target_capture_seqs_mapped/<sample>/<gene>/`.


- **Write** per-paralog `*_scipio_hits.fasta` (and `*_scipio_hits_introns.fasta` when introns are present) for **unmapped** paralogs to `02_target_capture_seqs_unmapped/<status>/<sample>/<gene>/`. One subfolder per gate/status: `no_chimera_as_single_contig_hit`, `no_detectable_chimera_as_single_contig_hit_after_length_filtering`, `no_detectable_chimera_as_no_seqs_left_after_length_filtering`, `no_chimera_test_performed_as_too_little_seq_length_after_length_filtering`, `no_chimera_test_performed_as_too_few_contigs_after_length_filtering`.

**4. Reports and filtered FASTAs**  

- Writes a numbered series of TSVs and heatmaps under `00_logs_and_reports/reports/`:
  - `01_chimera_report_non_chimera_count_vs_chimera_count.tsv` — raw `non_chimera/chimera` count per (gene, sample).
  - `02_chimera_report_heatmap_data.tsv` — the chimera-proportion matrix used to generate heatmap 03.
  - `03_chimera_proportion_heatmap_<assembly_source>.png` — chimera proportion heatmap.
  - `04_genes_with_one_seq_per_sample_in_90_percent_samples.tsv` — genes where >=90% of samples recovered a single non-chimeric paralog.
  - `05_genes_not_recovered_after_sample_threshold_<min_samples_threshold>_filtering.tsv` — genes dropped because fewer than `min_samples_threshold` × N samples cleared `--min_num_paralogs_per_sample`.
  - `06_non_chimeric_paralogs_after_sample_threshold_<min_samples_threshold>_filtering.tsv` — per-(gene, sample) non-chimeric paralog count for the surviving genes.
  - `07_non_chimeric_paralog_count_heatmap_<assembly_source>_after_sample_threshold_<min_samples_threshold>_filtering.png` — continuous heatmap of the counts in `06*`.
  - `08_non_chimeric_paralog_count_heatmap_<assembly_source>_after_sample_threshold_<min_samples_threshold>_filtering_max_count_<cap>_stepped.png` — categorical heatmap of the counts in `06*`. Count >= `--non_chimeric_paralog_max_count` (default: 10) collapse to the top bucket.
  - `09_non_chimeric_paralog_count_heatmap_<assembly_source>_after_sample_threshold_<min_samples_threshold>_filtering_max_count_<cap>_continuous.png` — continuous version of the same capped matrix.
- Writes non-chimeric sequences to FASTA files. Applies `--min_samples_threshold` and `--min_num_paralogs_per_sample` so only genes with enough qualifying samples get per-gene exon FASTAs written to `03_non_chimeric_target_capture_seqs_<min_samples_threshold>/`; matching introns go under `introns/` when present, with `--intron_ref_coords_wiggle` controlling merging of intron records across samples based on their coordinates with respect to the protein reference used by the target-capture software.
- Writes per-gene `.all.fasta` files to `04_all_target_capture_seqs/`; these contain all recovered coding sequences for that gene (chimeric and non-chimeric) for downstream comparison.

---

## Assumptions and limitations

- The reference should be close enough that exon-sized pieces from your assembly usually map to one primary reference sequence name when they are truly co-locus, and to different names when pieces come from different genomic regions. Fragmented references, mis-joined reference contigs, or naming that splits one biological locus across several FASTA records will weaken interpretation.

---

## Requirements

- **Python** 
- Python modules: **biopython**, **pandas**, **numpy**, **matplotlib**, **seaborn**, **pebble**, **tqdm**  
- **BBTools** on `PATH`:  
  - `bbmap.sh` — index build 
  - `mapPacBio.sh` — align per-contig query FASTAs to the indexed reference

---

## Install

The recommended way to get everything (Python, the listed Python modules, and `bbmap.sh` / `mapPacBio.sh` from BBTools) is a dedicated `conda` environment using the `conda-forge` and `bioconda` channels. If you don't already have a conda distribution installed, the lightweight [Miniforge](https://github.com/conda-forge/miniforge) installer is recommended.

**1. Create the environment**

```bash
conda create -n filter_chimeras -c conda-forge -c bioconda \
  python=3.11 \
  biopython pandas numpy matplotlib seaborn pebble tqdm \
  bbmap
```

The `bbmap` package from bioconda installs the full BBTools suite, which includes both `bbmap.sh` and `mapPacBio.sh`.

**2. Activate it**

```bash
conda activate filter_chimeras
```

**3. Verify**

```bash
bbmap.sh --version
mapPacBio.sh --version
python -c "import Bio, pandas, numpy, matplotlib, seaborn, pebble, tqdm; print('python deps OK')"
```

**4. Get the script**

```bash
git clone https://github.com/<your-fork-or-upstream>/filter_chimeras.git
cd filter_chimeras
python filter_chimeras.py --help
```

---

## Usage

Provide either `--captus_folder` or `--hybpiper_folder` (mutually exclusive), plus the reference genome and target FASTA used for the assemblies.

**Captus example**

```bash
python filter_chimeras.py \
  --captus_folder path/to/captus_sample_parent_folder \
  reference_genome.fasta \
  targets_used_in_assembly.fasta \
  [--output_directory my_run_output] \
  [options]
```

**HybPiper example**

```bash
python filter_chimeras.py \
  --hybpiper_folder path/to/hybpiper_sample_parent_folder \
  reference_genome.fasta \
  targets_used_in_assembly.fasta \
  [--output_directory my_run_output] \
  [options]
```

**Full options**

```
usage: filter_chimeras.py [-h] [--version] (--captus_folder CAPTUS_FOLDER |
                          --hybpiper_folder HYBPIPER_FOLDER)
                          [--output_directory OUTPUT_DIRECTORY]
                          [--min_seq_length MIN_SEQ_LENGTH]
                          [--min_length_percentage MIN_LENGTH_PERCENTAGE]
                          [--min_contig_number_percentage MIN_CONTIG_NUMBER_PERCENTAGE]
                          [--min_samples_threshold MIN_SAMPLES_THRESHOLD]
                          [--min_num_paralogs_per_sample MIN_NUM_PARALOGS_PER_SAMPLE]
                          [--non_chimeric_paralog_max_count NON_CHIMERIC_PARALOG_MAX_COUNT]
                          [--pool POOL] [--threads THREADS]
                          [--intron_ref_coords_wiggle INTRON_REF_COORDS_WIGGLE]
                          [--bbmap_Xmx BBMAP_XMX] [--mappacbio_k MAPPACBIO_K]
                          [--mappacbio_fastareadlen MAPPACBIO_FASTAREADLEN]
                          [--mappacbio_maxindel MAPPACBIO_MAXINDEL]
                          [--mappacbio_minid MAPPACBIO_MINID]
                          reference_genome_fasta target_file_fasta

positional arguments:
  reference_genome_fasta
                        High quality reference genome for mapping target-
                        capture contigs.
  target_file_fasta     Target FASTA used for the assembly runs (Captus or
                        HybPiper).

positional arguments:
  reference_genome_fasta
                        High quality reference genome for mapping target-
                        capture contigs.
  target_file_fasta     Target FASTA used for the assembly runs (Captus or
                        HybPiper).

options:
  -h, --help            show this help message and exit
  --version, -v         Print the script version number.
  --captus_folder CAPTUS_FOLDER
                        Parent folder containing Captus sample output folders.
  --hybpiper_folder HYBPIPER_FOLDER
                        Parent folder containing HybPiper sample output
                        folders.
  --output_directory OUTPUT_DIRECTORY
                        Root folder for logs, reports, mapped intermediates,
                        and filtered FASTAs. Default is: output_directory
  --min_seq_length MIN_SEQ_LENGTH
                        Minimum length for a contig-derived sequence to be
                        mapped. Default is: 75
  --min_length_percentage MIN_LENGTH_PERCENTAGE
                        The minimum percentage of the paralog sequence length
                        retained after filtering contig hits via
                        <min_seq_length> for chimera detection to be
                        performed. Default is: 0.8
  --min_contig_number_percentage MIN_CONTIG_NUMBER_PERCENTAGE
                        The minimum percentage of the total number of contig
                        hits remaining for a given paralog after filtering
                        contig hits via <min_seq_length> for chimera detection
                        to be performed.Default is: 0.8
  --min_samples_threshold MIN_SAMPLES_THRESHOLD
                        For a given gene, the minimum percentage of total
                        samples to have >= <min_num_paralogs_per_sample> non-
                        chimeric sequences for sequences to be written to
                        file. Useful to increase stringency of alignment
                        sample occupancy, at the cost of fewer gene
                        alignments. Default is: 0.75
  --min_num_paralogs_per_sample MIN_NUM_PARALOGS_PER_SAMPLE
                        For a given gene for a given sample, the minimum
                        number of non-chimeric paralog sequences recovered for
                        the sequences to be written to file. This can be
                        useful when you expect paralogs to be present (e.g.
                        due to polyploidy), and want to skip genes that might
                        have hidden paralogy due to sequence or assembly
                        issues. Default is: 0
  --non_chimeric_paralog_max_count NON_CHIMERIC_PARALOG_MAX_COUNT
                        Cap value used by the discrete non-chimeric paralog
                        count heatmaps. Per-(gene, sample) non-chimeric
                        paralog counts above this value are collapsed into a
                        single "<cap>+" bucket so the categorical legend stays
                        readable. Values <= 10 use the categorical tab10
                        palette; values > 10 sample the viridis sequential
                        colormap at that many points. Default is: 10
  --pool POOL           Number of samples to run concurrently. Default is: 1
  --threads THREADS     Number of threads to use for mapping for a given
                        sample. Default is: 1
  --intron_ref_coords_wiggle INTRON_REF_COORDS_WIGGLE
                        When writing non-chimeric intron FASTAs under
                        introns/, merge ref_protein_coords groups whose query
                        intervals overlap or touch after expanding each end by
                        this many aa. Use 0 to merge adjacent junctions (e.g.
                        82-83 with 83-84). Use -1 (default) for legacy
                        behaviour: one file per exact coord string only.
                        Merged groups include every coord label in the
                        filename, sorted (e.g.
                        gene_ref_protein_coords_82-83_83-84.fasta).
  --bbmap_Xmx BBMAP_XMX
                        Amount of memory to allocate to BBmap.sh in GB.
                        Default is: 30

mapPacBio.sh tuning options:
  --mappacbio_k MAPPACBIO_K
                        mapPacBio.sh k-mer size (k=). Default is: 13
  --mappacbio_fastareadlen MAPPACBIO_FASTAREADLEN
                        mapPacBio.sh fastareadlen= (max FASTA read length).
                        Default is: 6000
  --mappacbio_maxindel MAPPACBIO_MAXINDEL
                        mapPacBio.sh maxindel= (max indel length tolerated).
                        Default is: 100000
  --mappacbio_minid MAPPACBIO_MINID
                        mapPacBio.sh minid= (minimum identity to keep an
                        alignment). Default is: 0.76
```


---

## Output

| Path | Contents |
|------|----------|
| `00_logs_and_reports/reports/` | Numbered TSV and PNG reports (`01_..tsv` through `09_..png`) — see the "Reports and filtered FASTAs" section above for what each file contains. |
| `00_logs_and_reports/logs/` | Per-run log files. |
| `01_target_capture_seqs_mapped/<sample>/<gene>/` | Per-paralog query FASTAs (`<paralog>_scipio_hits.fasta`), per-paralog intron FASTAs (`<paralog>_scipio_hits_introns.fasta`), and `mapPacBio.sh` SAMs for multi-fragment paralogs that went through mapping. |
| `02_target_capture_seqs_unmapped/<status>/<sample>/<gene>/` | Per-paralog `<paralog>_scipio_hits.fasta` (+ `_introns.fasta` when present) for paralogs that bypassed mapping; one subfolder per status (e.g. `no_chimera_as_single_contig_hit/`, `no_detectable_chimera_as_single_contig_hit_after_length_filtering/`, `no_detectable_chimera_as_no_seqs_left_after_length_filtering/`, `no_chimera_test_performed_as_too_little_seq_length_after_length_filtering/`, `no_chimera_test_performed_as_too_few_contigs_after_length_filtering/`). |
| `03_non_chimeric_target_capture_seqs_<min_samples_threshold>/` | One `<gene>.fasta` per gene (non-chimeric exonic sequences passing thresholds), plus an `introns/` subfolder with per-gene intron FASTAs grouped by `ref_protein_coords` (see `--intron_ref_coords_wiggle`). |
| `04_all_target_capture_seqs/` | `<gene>.all.fasta` — all recovered coding sequences per gene regardless of chimera call. |

The BBMap index is written to `<output_dir>/ref/genome/1/`.

---
