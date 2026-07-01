# Rotavirus VP4/VP7 NCBI Retrieval Pipeline

A Python pipeline for retrieving, filtering, and exporting Rotavirus A nucleotide sequences from NCBI GenBank via Biopython Entrez, with Mozambique prioritization and optional MAFFT alignment.

---

## Overview

The pipeline is designed to support phylogenetic analysis of rotavirus VP4 (P-type) and VP7 (G-type) sequences, with a focus on Mozambican strains.

**Core workflow:**

1. Search NCBI Nucleotide (`nuccore`) via Entrez using a broad query
2. Fetch full GenBank records in batches
3. Deduplicate records by accession
4. Extract metadata locally: country, year, G-type, P-type, sequence length
5. Filter by year range, minimum sequence length, and target genotype
6. Select the final set: either a simple priority-first, recency-sorted slice (default), or optional random stratified sampling across countries with a configurable per-country cap and an optional seed for reproducibility
7. Export outputs: FASTA (informative + standard headers), metadata CSVs, accession lists, summary, and optionally Excel and a MAFFT alignment

---

## Requirements

```
biopython
pandas       # optional — for Excel export
openpyxl     # optional — for Excel export
mafft        # optional — for alignment (must be on PATH)
```

Install Python dependencies:

```bash
pip install biopython pandas openpyxl
```

---

## Usage

### Running the script

```bash
python seq_quest.py
```

The script prompts for all parameters at runtime:

| Prompt | Default |
|---|---|
| NCBI email | `youremail@email.com` |
| Organism | `Rotavirus A` |
| Gene/Protein | `VP4` |
| Target genotype | `P6` |
| Priority country | `Mozambique` |
| Start year | `2000` |
| End year | `2025` |
| Minimum sequence length | `1500` |
| Maximum sequences | `500` |
| Max IDs to retrieve from NCBI (retmax) | `10000` |
| Use per-country random sampling? | `n` |
| ↳ Max sequences per country *(if sampling = y)* | `auto` (blank = split remaining slots evenly across countries) |
| ↳ Random seed for reproducibility *(if sampling = y)* | `random` (blank = a different sample each run) |
| Export Excel files too? | `y` |
| Run MAFFT alignment on selected sequences? | `n` |

### Notebook / IPython

There's no `argparse` interface — every parameter is gathered through `input()` prompts in `get_user_inputs()`. This means the script runs the same way whether invoked from a terminal or pasted into a Jupyter cell; the prompts just appear as inline input boxes.

---

## Outputs

Each run creates a timestamped directory under `outputs/`, e.g. `outputs/20250601_143022_vp4_p6/`:

| File | Contents |
|---|---|
| `selected.fasta` | Selected sequences with RVA strain name headers (e.g. `RVA/Human-wt/MOZ/xxx/2018/G9P[6]`); falls back to `accession\|country\|year\|genotype` if no strain name is present; records with unreliable years are flagged with `\|YEAR_UNVERIFIED` |
| `selected_verified_year.fasta` | Same as `selected.fasta` but with fallback-year records excluded entirely (Option A — use this when year accuracy is critical) |
| `selected_standard_header.fasta` | Same sequences with original NCBI-style headers (for tools that require standard format) |
| `selected_metadata.csv` | Metadata for selected sequences |
| `all_metadata.csv` | Metadata for all fetched records (pre-filter) |
| `selected_accessions.txt` | Accession list for selected sequences |
| `all_accessions.txt` | Accession list for all fetched records |
| `summary.txt` | Run parameters + country/genotype/year distribution report |
| `all_metadata.xlsx` | *(optional)* Excel version of all metadata |
| `selected_metadata.xlsx` | *(optional)* Excel version of selected metadata |
| `selected_aligned.fasta` | *(optional)* MAFFT alignment of selected sequences |

A persistent query log is appended to `logs/query_history.csv` after each run.

### Metadata CSV columns

`accession`, `version`, `strain_name`, `country`, `country_source`, `year`, `year_source`, `g_type`, `g_source`, `p_type`, `p_source`, `length`, `is_priority`, `description`

The `*_source` columns indicate which GenBank field each value was extracted from, enabling audit of uncertain extractions. `year_source` values of `annotations_date_fallback` or `annotations_createdate_fallback` indicate the submission date was used as a last resort — these should be treated with caution.

---

## Key design decisions

**Broad search, local filtering.** The Entrez query is kept intentionally broad (organism + gene name only). Genotype filtering happens locally after metadata extraction, which is more reliable than relying on NCBI search term matching across inconsistently annotated records.

**Year extraction priority.** The pipeline uses a strict priority order to avoid conflating the GenBank submission date (often recent) with the true biological collection year:
1. `collection_date` qualifier on the source feature
2. Strain name (rotavirus strains conventionally embed the collection year, e.g. `RVA/Human-wt/MOZ/.../2003/G9P[6]`)
3. Record description
4. `note` qualifier
5. Last resort: `annotations['date']` / `annotations['createdate']` (submission date — flagged in output)

**Country extraction priority.** Tries `country` qualifier → `geo_loc_name` qualifier → text inference from `note`/`strain`/`isolation_source` fields → regex fallback on full record text. Also detects Mozambique via the `/MOZ/` ISO code pattern in strain names.

**Selection mode (standard vs. random sampling).** Two selection modes are available, chosen at runtime:

- **Standard (default, `n`):** Priority-country sequences are taken first (sorted by year descending), followed by all other sequences (also sorted by year descending), trimmed to `max_seqs`. Deterministic and unambiguous — use this when you want the most recent sequences without any randomness.
- **Per-country random sampling (`y`):** Priority-country sequences are shuffled and taken first. Remaining slots are distributed across the other countries present — either an even auto-split of the remaining budget or a fixed `seqs_per_country` cap — and filled by random sampling within each country's pool. Provide a `random_seed` for a fully reproducible draw. Use this when you want geographic breadth rather than just the most recent records.

**FASTA headers and year reliability.** `selected.fasta` uses the RVA strain name from the GenBank `strain` qualifier as the header (e.g. `RVA/Human-wt/MOZ/BRT-001/2018/G9P[6]`), falling back to `accession|country|year|genotype` for records without one. Any record whose year was inferred from the GenBank submission date rather than a true collection date is flagged with `|YEAR_UNVERIFIED` in the header. A companion file `selected_verified_year.fasta` excludes those records entirely — use it when year accuracy is critical for phylogenetic dating. `selected_standard_header.fasta` retains unmodified NCBI headers for tools that validate FASTA header format.

---

## Configuration

Technical defaults are set in the `CONFIG` dict at the top of the script:

```python
CONFIG = {
    "email": "youremail@gmail.com",       # REQUIRED by NCBI
    "api_key": None,                      # Optional NCBI API key (increases rate limit)
    "retmax": 10000,                      # Default max IDs to retrieve — overridable at the prompt
    "batch_size": 200,                    # Records per efetch batch
    "sleep": 0.34,                        # Polite delay between batches (seconds)
    "default_seqs_per_country": None,     # Per-country sampling cap; None = auto-derive
    ...
}
```

An NCBI API key is not required but recommended for large queries — it raises the rate limit from 3 to 10 requests/second.

---

## Downstream use

The selected FASTA files are ready for:

- **MAFFT** multiple sequence alignment (can be triggered within the pipeline)
- **IQ-TREE** / **RAxML** phylogenetic tree inference
- **FigTree** / **iTOL** tree visualization

For phylogenetic work, use `selected_standard_header.fasta` if your alignment tool is strict about header format, and `selected.fasta` for annotated visualization.

---

## Notes

- `normalize_country()` trims GenBank's `"Country: Region"` formatting down to just the country (e.g. `Mozambique: Maputo` → `Mozambique`), and `normalize_text()` lowercases and strips common Portuguese diacritics (ç, ã, á, é, etc.) before comparing country names. `country_matches()` applies an alias table on top of this to catch variants that diacritic-stripping alone cannot resolve — for example, `Moçambique` normalises to `mocambique` (ç→c), which is mapped back to `mozambique` before comparison.
- To switch from VP4 (P-type) to VP7 (G-type), change the gene prompt input to `VP7` — the pipeline auto-detects which genotype letter to filter on.
