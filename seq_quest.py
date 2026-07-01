from __future__ import annotations

import csv
import random
import re
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import Entrez, SeqIO

try:
    import pandas as pd
except Exception:
    pd = None


# ============================================================
# TECHNICAL DEFAULTS
# ============================================================
CONFIG = {
    "email": "youremail@email.com" ,
    "api_key": None,
    "retmax": 10000,
    "batch_size": 200,
    "sleep": 0.34,
    "default_min_year": 2000,
    "default_max_year": 2025,
    "default_min_length": 1500,
    "default_max_seqs": 500,
    # How many sequences to allow per non-priority country at most.
    # Set to None to let the quota be derived automatically from
    # max_seqs and the number of countries found.
    "default_seqs_per_country": None,
}

QUERY_HISTORY_FILE = Path("logs/query_history.csv")


# ============================================================
# INPUTS
# ============================================================
def prompt_text(label: str, default: str) -> str:
    value = input(f"{label} [{default}]: ").strip()
    return value if value else default


def prompt_int(label: str, default: int) -> int:
    value = input(f"{label} [{default}]: ").strip()
    if not value:
        return default
    return int(value)


def prompt_int_or_none(label: str, default_label: str = "auto") -> Optional[int]:
    """Prompt for an integer; return None if left blank (means 'auto')."""
    value = input(f"{label} [{default_label}]: ").strip()
    if not value:
        return None
    return int(value)


def prompt_yes_no(label: str, default: str = "n") -> bool:
    default = default.lower().strip()
    value = input(f"{label} [{default}]: ").strip().lower()
    if not value:
        value = default
    return value in {"y", "yes", "true", "1"}


def get_user_inputs() -> Dict:
    print("\n====================================")
    print("ROTAVIRUS NCBI RETRIEVAL PIPELINE")
    print("====================================\n")

    email = prompt_text("NCBI email", CONFIG["email"])
    organism = prompt_text("Organism", "Rotavirus A")
    gene = prompt_text("Gene/Protein", "VP4")
    target_genotype = prompt_text("Target genotype", "P6")
    priority_country = prompt_text("Priority country", "Mozambique")
    min_year = prompt_int("Start year", CONFIG["default_min_year"])
    max_year = prompt_int("End year", CONFIG["default_max_year"])
    min_length = prompt_int("Minimum sequence length", CONFIG["default_min_length"])
    max_seqs = prompt_int("Maximum sequences", CONFIG["default_max_seqs"])
    retmax = prompt_int("Max IDs to retrieve from NCBI (retmax)", CONFIG["retmax"])

    # --- Random sampling controls (optional) ---
    use_sampling = prompt_yes_no(
        "Use per-country random sampling? (n = take top N filtered sequences by priority then year)",
        "n",
    )
    seqs_per_country = None
    random_seed = None
    if use_sampling:
        print("\n--- Random sampling options (leave blank for defaults) ---")
        seqs_per_country = prompt_int_or_none(
            "Max sequences per country (blank = auto-distribute remaining slots)", "auto"
        )
        random_seed = prompt_int_or_none(
            "Random seed for reproducibility (blank = random each run)", "random"
        )

    export_xlsx = prompt_yes_no("Export Excel files too?", "y")
    run_alignment = prompt_yes_no("Run MAFFT alignment on selected sequences?", "n")

    # --- Input validation ---
    errors = []
    if min_year > max_year:
        errors.append(f"  Start year ({min_year}) is greater than end year ({max_year}).")
    if min_length <= 0:
        errors.append(f"  Minimum sequence length must be positive (got {min_length}).")
    if max_seqs <= 0:
        errors.append(f"  Maximum sequences must be positive (got {max_seqs}).")
    if retmax <= 0:
        errors.append(f"  NCBI retmax must be positive (got {retmax}).")
    if errors:
        print("\n[ERROR] Invalid inputs:")
        for e in errors:
            print(e)
        raise SystemExit(1)

    return {
        "email": email,
        "api_key": CONFIG["api_key"],
        "organism": organism,
        "gene": gene,
        "target_genotype": target_genotype,
        "priority_country": priority_country,
        "min_year": min_year,
        "max_year": max_year,
        "min_length": min_length,
        "max_seqs": max_seqs,
        "retmax": retmax,
        "use_sampling": use_sampling,
        "seqs_per_country": seqs_per_country,
        "random_seed": random_seed,
        "export_xlsx": export_xlsx,
        "run_alignment": run_alignment,
    }


# ============================================================
# HELPERS
# ============================================================
def normalize_text(text: Optional[str]) -> str:
    if not text:
        return ""
    text = text.strip().lower()
    replacements = {
        "ç": "c",
        "ã": "a",
        "á": "a",
        "à": "a",
        "â": "a",
        "é": "e",
        "ê": "e",
        "í": "i",
        "ó": "o",
        "ô": "o",
        "ú": "u",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"\s+", " ", text)
    return text


def safe_join(values) -> str:
    return " | ".join(str(v) for v in values if v is not None)


def normalize_country(country: Optional[str]) -> str:
    """
    Trims GenBank-style country values down to just the country.
    GenBank's 'country'/'geo_loc_name' qualifiers are frequently formatted
    as "Country: Region, Locality" (e.g. "Mozambique: Maputo"), so we split
    on ':' as well as ';,/' to avoid those failing an exact-match comparison
    against the priority country.
    """
    if not country:
        return "Unknown"
    country = country.strip()
    country = re.split(r"[;,/:]", country)[0].strip()
    return country if country else "Unknown"


def country_matches(country: str, priority_country: str) -> bool:
    # Maps normalized variants that diacritic-stripping can't resolve correctly
    # e.g. Moçambique → normalize_text gives "mocambique" ≠ "mozambique"
    COUNTRY_ALIASES: Dict[str, str] = {
        "mocambique": "mozambique",
    }
    def resolve(name: str) -> str:
        norm = normalize_text(name)
        return COUNTRY_ALIASES.get(norm, norm)
    return resolve(country) == resolve(priority_country)


def infer_country_from_text(text: str, priority_country: str) -> Optional[str]:
    text_u = text.upper()
    if "/MOZ/" in text_u:
        return "Mozambique"
    if normalize_text(priority_country) in normalize_text(text):
        return priority_country
    return None


def make_safe_token(text: str) -> str:
    text = normalize_text(text)
    text = re.sub(r"[^a-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unknown"


def extract_year_from_text(text: str) -> Optional[int]:
    if not text:
        return None
    m = re.search(r"(19\d{2}|20\d{2})", text)
    return int(m.group(1)) if m else None


def extract_strain_name(record) -> Optional[str]:
    """
    Extract the rotavirus strain name from the source feature.
    Rotavirus sequences typically carry a strain name in RVA nomenclature,
    e.g. RVA/Human-wt/MOZ/BRT-001/2018/G9P[6].
    Returns None if no strain qualifier is present.
    """
    for feat in getattr(record, "features", []):
        if feat.type != "source":
            continue
        quals = feat.qualifiers
        if "strain" in quals and quals["strain"]:
            return quals["strain"][0].strip()
    return None


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


# ============================================================
# RECORD TEXT
# ============================================================
def get_record_text(record) -> str:
    parts = [
        record.id or "",
        record.name or "",
        record.description or "",
        str(getattr(record, "annotations", {})),
    ]

    for feat in getattr(record, "features", []):
        if feat.type != "source":
            continue
        quals = feat.qualifiers
        for key in [
            "country",
            "geo_loc_name",
            "note",
            "isolation_source",
            "collection_date",
            "host",
            "strain",
            "segment",
        ]:
            if key in quals:
                parts.append(safe_join(quals[key]))

    return " | ".join(parts)


# ============================================================
# DEDUPLICATION
# ============================================================
def deduplicate_records(records: List) -> List:
    seen = set()
    unique = []
    for rec in records:
        acc = rec.id
        if acc in seen:
            continue
        seen.add(acc)
        unique.append(rec)
    return unique


# ============================================================
# COUNTRY EXTRACTION
# ============================================================
def extract_country_with_source(record, priority_country: str) -> Tuple[str, str]:
    # 1) direct qualifiers
    for feat in getattr(record, "features", []):
        if feat.type != "source":
            continue
        quals = feat.qualifiers

        if "country" in quals and quals["country"]:
            return normalize_country(quals["country"][0]), "country"
        if "geo_loc_name" in quals and quals["geo_loc_name"]:
            return normalize_country(quals["geo_loc_name"][0]), "geo_loc_name"

    # 2) source-feature text
    for feat in getattr(record, "features", []):
        if feat.type != "source":
            continue
        quals = feat.qualifiers
        for key in ["note", "strain", "isolation_source"]:
            if key in quals:
                joined = safe_join(quals[key])
                inferred = infer_country_from_text(joined, priority_country)
                if inferred:
                    return inferred, f"{key}_text"

    # 3) description / annotations fallback
    text = get_record_text(record)
    patterns = [
        r'country="([^"]+)"',
        r'geo_loc_name="([^"]+)"',
        r"country=([^;|]+)",
        r"geo_loc_name=([^;|]+)",
        r"country:\s*([^;|]+)",
        r"geo_loc_name:\s*([^;|]+)",
    ]
    for p in patterns:
        m = re.search(p, text, flags=re.IGNORECASE)
        if m:
            return normalize_country(m.group(1).strip()), "text_pattern"

    inferred = infer_country_from_text(text, priority_country)
    if inferred:
        return inferred, "text_inferred"

    return "Unknown", "unknown"


# ============================================================
# YEAR EXTRACTION
# ============================================================
def get_record_year(record) -> Tuple[Optional[int], str]:
    """
    Extracts the biological collection year, in priority order.

    Priority order (most reliable first):
      1. collection_date qualifier on the source feature
      2. strain name (rotavirus strains conventionally embed the year)
      3. record description
      4. note qualifier on the source feature
      5. LAST RESORT: annotations 'date' / 'createdate' (submission date)
    """
    # 1) Explicit collection date
    for feat in getattr(record, "features", []):
        if feat.type != "source":
            continue
        quals = feat.qualifiers
        if "collection_date" in quals:
            year = extract_year_from_text(safe_join(quals["collection_date"]))
            if year:
                return year, "collection_date"

    # 2) Strain name
    for feat in getattr(record, "features", []):
        if feat.type != "source":
            continue
        quals = feat.qualifiers
        if "strain" in quals:
            year = extract_year_from_text(safe_join(quals["strain"]))
            if year:
                return year, "strain"

    # 3) Description
    if record.description:
        year = extract_year_from_text(record.description)
        if year:
            return year, "description"

    # 4) Note qualifier
    for feat in getattr(record, "features", []):
        if feat.type != "source":
            continue
        quals = feat.qualifiers
        if "note" in quals:
            year = extract_year_from_text(safe_join(quals["note"]))
            if year:
                return year, "note"

    # 5) Fallback only: GenBank submission/modification date.
    for key in ["date", "createdate"]:
        val = record.annotations.get(key)
        if val:
            year = extract_year_from_text(str(val))
            if year:
                return year, f"annotations_{key}_fallback"

    return None, "unknown"


# ============================================================
# GENOTYPE EXTRACTION
# ============================================================
def extract_genotype_component(record, letter: str) -> Tuple[Optional[str], str]:
    letter = letter.upper().strip()
    text_sources = []

    for feat in getattr(record, "features", []):
        if feat.type != "source":
            continue
        quals = feat.qualifiers
        for key in ["note", "strain", "host", "isolation_source"]:
            if key in quals:
                text_sources.append((f"source:{key}", safe_join(quals[key])))

    text_sources.append(("record:description", record.description or ""))
    text_sources.append(("record:annotations", str(getattr(record, "annotations", {}))))
    text_sources.append(("record:text", get_record_text(record)))

    patterns = [
        rf"\b{letter}\[(\d+)\]",
        rf"\b{letter}\s*[\[\(]?\s*(\d+)\s*[\]\)]?\b",
        rf"\b{letter}(\d+)\b",
        rf"{letter}\s*-\s*(\d+)",
    ]

    for source_name, text in text_sources:
        for pattern in patterns:
            m = re.search(pattern, text, flags=re.IGNORECASE)
            if m:
                return f"{letter}{m.group(1)}", source_name

    return None, "unknown"


def extract_both_genotypes(record) -> Tuple[Optional[str], str, Optional[str], str]:
    g_type, g_src = extract_genotype_component(record, "G")
    p_type, p_src = extract_genotype_component(record, "P")
    return g_type, g_src, p_type, p_src


def infer_target_letter(gene: str) -> str:
    gene_u = gene.upper().strip()
    if "VP4" in gene_u:
        return "P"
    if "VP7" in gene_u:
        return "G"
    return "P"


def normalize_target_genotype(target: str, target_letter: str) -> str:
    target = target.strip().upper().replace(" ", "")
    if re.fullmatch(r"[GP]\d+", target):
        return target
    if re.fullmatch(r"\d+", target):
        return f"{target_letter.upper()}{target}"
    return target


# ============================================================
# ENTREZ
# ============================================================
def build_query(organism: str, gene: str) -> str:
    org_block = f'"{organism}"[Organism]'
    gene_block = f"({gene}[All Fields])"
    return f"{org_block} AND {gene_block}"


def fetch_ids(query: str, retmax: int) -> List[str]:
    handle = Entrez.esearch(
        db="nuccore",
        term=query,
        retmax=retmax,
        usehistory="n",
    )
    result = Entrez.read(handle)
    handle.close()
    return result.get("IdList", [])


def fetch_records(id_list: List[str], batch_size: int, sleep_s: float):
    total = len(id_list)
    n_batches = (total + batch_size - 1) // batch_size
    for batch_num, start in enumerate(range(0, total, batch_size), 1):
        chunk = id_list[start : start + batch_size]
        end = min(start + batch_size, total)
        print(f"[INFO] Fetching batch {batch_num}/{n_batches} (records {start + 1}–{end} of {total})...")
        handle = Entrez.efetch(
            db="nuccore",
            id=",".join(chunk),
            rettype="gb",
            retmode="text",
        )
        records = list(SeqIO.parse(handle, "genbank"))
        handle.close()

        for rec in records:
            yield rec

        time.sleep(sleep_s)


# ============================================================
# METADATA
# ============================================================
def build_metadata(records: List, priority_country: str) -> List[Dict]:
    rows = []
    for rec in records:
        country, country_source = extract_country_with_source(rec, priority_country)
        year, year_source = get_record_year(rec)
        g_type, g_source, p_type, p_source = extract_both_genotypes(rec)
        strain_name = extract_strain_name(rec)

        rows.append(
            {
                "accession": rec.id,
                "version": rec.id,
                "strain_name": strain_name,
                "country": country,
                "country_source": country_source,
                "year": year,
                "year_source": year_source,
                "g_type": g_type,
                "g_source": g_source,
                "p_type": p_type,
                "p_source": p_source,
                "length": len(rec.seq),
                "description": rec.description or "",
                "record": rec,
            }
        )
    return rows


def report_metadata(rows: List[Dict]) -> str:
    lines = []
    lines.append("===== COUNTRY SUMMARY =====")
    country_counts = Counter(row["country"] for row in rows)
    for country, count in country_counts.most_common(20):
        lines.append(f"{country:35} {count}")

    lines.append("")
    lines.append("===== G TYPE SUMMARY =====")
    g_counts = Counter(row["g_type"] for row in rows)
    for g, count in g_counts.most_common():
        lines.append(f"{str(g):10} {count}")

    lines.append("")
    lines.append("===== P TYPE SUMMARY =====")
    p_counts = Counter(row["p_type"] for row in rows)
    for p, count in p_counts.most_common():
        lines.append(f"{str(p):10} {count}")

    lines.append("")
    lines.append("===== YEAR SUMMARY =====")
    year_counts = Counter(row["year"] for row in rows if row["year"] is not None)
    for year, count in year_counts.most_common(20):
        lines.append(f"{year}: {count}")

    lines.append("")
    lines.append("===== YEAR SOURCE SUMMARY =====")
    lines.append("(how the year was determined - 'annotations_*_fallback' means")
    lines.append(" it came from the GenBank submission date, not a true collection")
    lines.append(" date, and is the least reliable category)")
    year_source_counts = Counter(row["year_source"] for row in rows)
    for source, count in year_source_counts.most_common():
        lines.append(f"{source:30} {count}")

    report = "\n".join(lines)
    print("\n" + report + "\n")
    return report


# ============================================================
# FILTERING
# ============================================================
def filter_rows(
    rows: List[Dict],
    min_year: int,
    max_year: int,
    min_length: int,
    target_genotype: str,
    target_letter: str,
    priority_country: str,
) -> List[Dict]:
    target_genotype = normalize_target_genotype(target_genotype, target_letter)
    filtered = []

    for row in rows:
        year = row["year"]
        if year is None:
            continue
        if not (min_year <= year <= max_year):
            continue
        if row["length"] < min_length:
            continue

        if target_letter == "P":
            candidate = row["p_type"]
        else:
            candidate = row["g_type"]

        if candidate is None:
            continue
        if candidate.upper().strip() != target_genotype.upper().strip():
            continue

        row["is_priority"] = country_matches(row["country"], priority_country)
        filtered.append(row)

    return filtered


# ============================================================
# STANDARD SELECTION (no sampling)
# ============================================================
def sort_and_select(
    rows: List[Dict],
    max_seqs: int,
    priority_country: str,
) -> List[Dict]:
    """
    Simple selection without random sampling:
    priority-country sequences first (sorted by year descending),
    then all other sequences (sorted by year descending),
    trimmed to max_seqs.
    """
    priority = sorted(
        [r for r in rows if r.get("is_priority")],
        key=lambda r: r["year"] or 0,
        reverse=True,
    )
    others = sorted(
        [r for r in rows if not r.get("is_priority")],
        key=lambda r: r["year"] or 0,
        reverse=True,
    )
    return (priority + others)[:max_seqs]


# ============================================================
# RANDOM SAMPLING (OPTIONAL)
# ============================================================
def sample_per_country(
    rows: List[Dict],
    max_seqs: int,
    priority_country: str,
    seqs_per_country: Optional[int] = None,
    random_seed: Optional[int] = None,
) -> List[Dict]:
    """
    Randomly samples sequences so the final selection is drawn from across
    countries rather than always picking the first (most-recent) sequences.

    Strategy:
      1. Priority-country sequences are always taken first (up to max_seqs).
      2. Remaining slots are filled by randomly sampling from each other
         country's pool, respecting seqs_per_country quota.
      3. If seqs_per_country is None (auto), the remaining budget is split
         evenly across the non-priority countries that have sequences.
      4. A random_seed ensures reproducibility when supplied.

    Parameters
    ----------
    rows             : filtered metadata rows (each has 'is_priority' flag)
    max_seqs         : hard ceiling on total sequences selected
    priority_country : country that gets precedence (e.g. "Mozambique")
    seqs_per_country : max sequences to take from each non-priority country;
                       None = auto-derive from remaining budget
    random_seed      : seed for random.sample(); None = non-deterministic
    """
    rng = random.Random(random_seed)  # isolated RNG so global state isn't affected

    # --- Split into priority and others ---
    priority_rows = [r for r in rows if r.get("is_priority")]
    other_rows    = [r for r in rows if not r.get("is_priority")]

    # --- Sample priority country ---
    # Shuffle so we don't always get the same subset when there are more
    # priority sequences than we can fit.
    rng.shuffle(priority_rows)
    selected_priority = priority_rows[:max_seqs]

    remaining_slots = max_seqs - len(selected_priority)
    if remaining_slots <= 0:
        print(
            f"[INFO] Priority country filled all {max_seqs} slots; "
            "no room for other countries."
        )
        return selected_priority

    # --- Group non-priority rows by country ---
    by_country: Dict[str, List[Dict]] = defaultdict(list)
    for row in other_rows:
        by_country[row["country"]].append(row)

    n_countries = len(by_country)
    if n_countries == 0:
        return selected_priority

    # --- Determine per-country quota ---
    if seqs_per_country is None:
        # Distribute remaining slots as evenly as possible across countries
        quota = max(1, remaining_slots // n_countries)
    else:
        quota = seqs_per_country

    print(
        f"[INFO] Sampling up to {quota} sequence(s) per non-priority country "
        f"across {n_countries} countries ({remaining_slots} remaining slots)."
    )

    # --- Random sample from each country ---
    selected_others: List[Dict] = []
    for country, country_rows in sorted(by_country.items()):  # sorted = deterministic order
        k = min(quota, len(country_rows))
        sampled = rng.sample(country_rows, k)
        selected_others.extend(sampled)
        print(
            f"  {country:35}  pool={len(country_rows):4d}  sampled={k}"
        )

    # Trim to remaining budget
    rng.shuffle(selected_others)          # shuffle before trimming to avoid
    selected_others = selected_others[:remaining_slots]  # always favouring A-Z first

    selected = selected_priority + selected_others
    return selected


# ============================================================
# OUTPUTS
# ============================================================
def write_csv(rows: List[Dict], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "accession",
                "version",
                "strain_name",
                "country",
                "country_source",
                "year",
                "year_source",
                "g_type",
                "g_source",
                "p_type",
                "p_source",
                "length",
                "is_priority",
                "description",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "accession": row["accession"],
                    "version": row["version"],
                    "strain_name": row.get("strain_name") or "",
                    "country": row["country"],
                    "country_source": row["country_source"],
                    "year": row["year"],
                    "year_source": row["year_source"],
                    "g_type": row["g_type"],
                    "g_source": row["g_source"],
                    "p_type": row["p_type"],
                    "p_source": row["p_source"],
                    "length": row["length"],
                    "is_priority": row.get("is_priority", False),
                    "description": row["description"],
                }
            )


def write_accessions(rows: List[Dict], path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        for row in rows:
            f.write(row["accession"] + "\n")


FALLBACK_YEAR_SOURCES = {"annotations_date_fallback", "annotations_createdate_fallback"}


def write_fasta(rows: List[Dict], path: Path) -> None:
    """
    Writes selected.fasta with RVA strain names as headers when available.
    Falls back to accession|country|year|genotype for records without a strain name.
    Records whose year came from a GenBank submission-date fallback are flagged
    with |YEAR_UNVERIFIED at the end of the header.
    """
    with path.open("w", encoding="utf-8") as f:
        for row in rows:
            rec = row["record"]
            strain = row.get("strain_name")
            if strain:
                header = strain
            else:
                country = make_safe_token(row["country"])
                year = row["year"] if row["year"] is not None else "NA"
                genotype = row["p_type"] if row["p_type"] else row["g_type"] or "NA"
                header = f"{row['accession']}|{country}|{year}|{genotype}"
            if row.get("year_source") in FALLBACK_YEAR_SOURCES:
                header = f"{header}|YEAR_UNVERIFIED"
            f.write(f">{header}\n")
            seq = str(rec.seq).upper()
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def write_fasta_standard_header(rows: List[Dict], path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        for row in rows:
            rec = row["record"]
            description = (rec.description or "").strip()
            header = f"{rec.id} {description}".strip()
            f.write(f">{header}\n")
            seq = str(rec.seq).upper()
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")



def write_fasta_verified_year(rows: List[Dict], path: Path) -> None:
    """
    Option A: writes a FASTA excluding any record whose year was determined
    from a GenBank submission-date fallback (year_source contains 'fallback').
    These records have unreliable years and should be treated with caution
    in phylogenetic analyses. The main selected.fasta flags them with
    |YEAR_UNVERIFIED; this file removes them entirely.
    """
    verified = [r for r in rows if r.get("year_source") not in FALLBACK_YEAR_SOURCES]
    excluded = len(rows) - len(verified)
    if excluded:
        print(
            f"[INFO] Verified-year FASTA: {excluded} record(s) with fallback year "
            "source excluded (flagged in selected.fasta as |YEAR_UNVERIFIED)."
        )
    else:
        print("[INFO] Verified-year FASTA: no fallback-year records in selection — identical to selected.fasta.")
    write_fasta(verified, path)


def write_summary_file(
    rows_all: List[Dict],
    rows_selected: List[Dict],
    summary_text: str,
    path: Path,
    inputs: Dict,
    query: str,
    total_ids: int,
) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write("ROTAVIRUS RETRIEVAL SUMMARY\n")
        f.write("===========================\n\n")
        f.write(f"Date: {datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"Query: {query}\n")
        f.write(f"Organism: {inputs['organism']}\n")
        f.write(f"Gene/Protein: {inputs['gene']}\n")
        f.write(f"Target genotype: {inputs['target_genotype']}\n")
        f.write(f"Priority country: {inputs['priority_country']}\n")
        f.write(f"Year range: {inputs['min_year']} - {inputs['max_year']}\n")
        f.write(f"Minimum length: {inputs['min_length']}\n")
        f.write(f"Max sequences: {inputs['max_seqs']}\n")
        f.write(f"NCBI retmax: {inputs.get('retmax', CONFIG['retmax'])}\n")
        # NEW: sampling parameters in summary
        f.write(f"Seqs per country: {inputs.get('seqs_per_country', 'auto')}\n")
        seed = inputs.get("random_seed")
        f.write(f"Random seed: {seed if seed is not None else 'random (not fixed)'}\n")
        f.write(f"IDs retrieved: {total_ids}\n")
        f.write(f"Records fetched: {len(rows_all)}\n")
        f.write(f"Selected records: {len(rows_selected)}\n\n")
        f.write(summary_text)
        f.write("\n")


def write_excel(rows: List[Dict], path: Path) -> None:
    if pd is None:
        print("[WARN] pandas is not available; skipping Excel export.")
        return

    df = pd.DataFrame(rows).copy()
    if "record" in df.columns:
        df = df.drop(columns=["record"])
    df.to_excel(path, index=False)


def append_query_history(
    inputs: Dict,
    total_ids: int,
    total_records: int,
    selected_records: int,
    output_dir: Path,
) -> None:
    ensure_dir(QUERY_HISTORY_FILE.parent)

    row = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "organism": inputs["organism"],
        "gene": inputs["gene"],
        "target_genotype": inputs["target_genotype"],
        "priority_country": inputs["priority_country"],
        "min_year": inputs["min_year"],
        "max_year": inputs["max_year"],
        "min_length": inputs["min_length"],
        "max_seqs": inputs["max_seqs"],
        "retmax": inputs.get("retmax", CONFIG["retmax"]),
        "seqs_per_country": inputs.get("seqs_per_country", "auto"),
        "random_seed": inputs.get("random_seed", ""),                 # NEW
        "total_ids": total_ids,
        "total_records": total_records,
        "selected_records": selected_records,
        "output_dir": str(output_dir),
    }

    file_exists = QUERY_HISTORY_FILE.exists()
    with QUERY_HISTORY_FILE.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


# ============================================================
# MAFFT
# ============================================================
def run_mafft(input_fasta: Path, output_fasta: Path) -> bool:
    mafft = shutil.which("mafft")
    if not mafft:
        print("[WARN] MAFFT not found on PATH; alignment skipped.")
        return False

    try:
        with output_fasta.open("w", encoding="utf-8") as out:
            subprocess.run(
                [mafft, "--auto", str(input_fasta)],
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )
        return True
    except subprocess.CalledProcessError as e:
        print("[WARN] MAFFT failed; alignment skipped.")
        if e.stderr:
            print(e.stderr[:1000])
        return False


# ============================================================
# PIPELINE
# ============================================================
def main() -> int:
    inputs = get_user_inputs()

    Entrez.email = inputs["email"]
    Entrez.tool = "rotavirus_metadata_first_pipeline"
    if inputs["api_key"]:
        Entrez.api_key = inputs["api_key"]

    target_letter = infer_target_letter(inputs["gene"])
    target_genotype = normalize_target_genotype(inputs["target_genotype"], target_letter)

    query = build_query(inputs["organism"], inputs["gene"])
    print(f"\n[INFO] Query: {query}")

    total_ids = fetch_ids(query, inputs["retmax"])
    print(f"[INFO] IDs retrieved: {len(total_ids)} (retmax={inputs['retmax']})")
    if not total_ids:
        print("[WARN] No IDs found. Broaden the query or check the gene name.")
        return 1

    records = list(fetch_records(total_ids, CONFIG["batch_size"], CONFIG["sleep"]))
    print(f"[INFO] Records fetched: {len(records)}")

    records = deduplicate_records(records)
    print(f"[INFO] Unique records: {len(records)}")

    metadata = build_metadata(records, inputs["priority_country"])
    report_text = report_metadata(metadata)

    filtered = filter_rows(
        rows=metadata,
        min_year=inputs["min_year"],
        max_year=inputs["max_year"],
        min_length=inputs["min_length"],
        target_genotype=target_genotype,
        target_letter=target_letter,
        priority_country=inputs["priority_country"],
    )
    print(f"[INFO] Records passing filters: {len(filtered)}")

    if inputs.get("use_sampling"):
        seed = inputs.get("random_seed")
        if seed is not None:
            print(f"[INFO] Random seed: {seed} (run is reproducible)")
        else:
            print("[INFO] No random seed set — sample will differ each run")

        selected = sample_per_country(
            rows=filtered,
            max_seqs=inputs["max_seqs"],
            priority_country=inputs["priority_country"],
            seqs_per_country=inputs.get("seqs_per_country"),
            random_seed=seed,
        )
    else:
        print("[INFO] Per-country sampling disabled — selecting top sequences by priority then year.")
        selected = sort_and_select(
            rows=filtered,
            max_seqs=inputs["max_seqs"],
            priority_country=inputs["priority_country"],
        )
    print(f"[INFO] Final selected records: {len(selected)}")

    if not selected:
        print("[WARN] No records matched the filters.")
        return 1

    # Output directory
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    gene_token = make_safe_token(inputs["gene"])
    geno_token = make_safe_token(target_genotype)
    out_dir = ensure_dir(Path("outputs") / f"{ts}_{gene_token}_{geno_token}")

    # Files
    all_csv              = out_dir / "all_metadata.csv"
    selected_csv         = out_dir / "selected_metadata.csv"
    all_acc              = out_dir / "all_accessions.txt"
    selected_acc         = out_dir / "selected_accessions.txt"
    selected_fasta            = out_dir / "selected.fasta"
    selected_fasta_std        = out_dir / "selected_standard_header.fasta"
    selected_fasta_verified   = out_dir / "selected_verified_year.fasta"
    summary_txt          = out_dir / "summary.txt"
    aligned_fasta        = out_dir / "selected_aligned.fasta"
    all_xlsx             = out_dir / "all_metadata.xlsx"
    selected_xlsx        = out_dir / "selected_metadata.xlsx"

    # Write outputs
    write_csv(metadata, all_csv)
    write_csv(selected, selected_csv)
    write_accessions(metadata, all_acc)
    write_accessions(selected, selected_acc)
    write_fasta(selected, selected_fasta)
    write_fasta_standard_header(selected, selected_fasta_std)
    write_fasta_verified_year(selected, selected_fasta_verified)
    write_summary_file(
        rows_all=metadata,
        rows_selected=selected,
        summary_text=report_text,
        path=summary_txt,
        inputs=inputs,
        query=query,
        total_ids=len(total_ids),
    )

    if inputs["export_xlsx"]:
        write_excel(metadata, all_xlsx)
        write_excel(selected, selected_xlsx)

    if inputs["run_alignment"]:
        aligned_ok = run_mafft(selected_fasta, aligned_fasta)
        if aligned_ok:
            print(f"[OK] Alignment written to: {aligned_fasta}")

    append_query_history(
        inputs=inputs,
        total_ids=len(total_ids),
        total_records=len(metadata),
        selected_records=len(selected),
        output_dir=out_dir,
    )

    print(f"\n[OK] Output directory: {out_dir}")
    print(f"[OK] All metadata CSV: {all_csv}")
    print(f"[OK] Selected metadata CSV: {selected_csv}")
    print(f"[OK] All accessions: {all_acc}")
    print(f"[OK] Selected accessions: {selected_acc}")
    print(f"[OK] Selected FASTA (strain name headers): {selected_fasta}")
    print(f"[OK] Selected FASTA (standard headers): {selected_fasta_std}")
    print(f"[OK] Selected FASTA (verified year only): {selected_fasta_verified}")
    print(f"[OK] Summary: {summary_txt}")

    if inputs["export_xlsx"]:
        print(f"[OK] All metadata Excel: {all_xlsx}")
        print(f"[OK] Selected metadata Excel: {selected_xlsx}")

    # Country breakdown of final selection
    print("\n[INFO] Country breakdown of selected sequences:")
    country_tally = Counter(r["country"] for r in selected)
    for country, cnt in country_tally.most_common():
        marker = " ← priority" if country_matches(country, inputs["priority_country"]) else ""
        print(f"  {country:35} {cnt}{marker}")

    moz = [r for r in selected if r.get("is_priority")]
    print(f"\n[INFO] Priority-country sequences among selected: {len(moz)}")
    for r in moz[:20]:
        g = r["g_type"] or "NA"
        p = r["p_type"] or "NA"
        print(f"  {r['accession']} | {r['country']} | {r['year']} | G={g} | P={p}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
