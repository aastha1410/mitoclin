#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════╗
║  MITOCLIN  –  Mitochondrial Genome Analysis Platform         ║
║  Single-file app: Flask backend + embedded frontend +        ║
║  HTML/PDF report generation                                  ║
║                                                              ║
║  Usage:                                                      ║
║    python3 mitoclin_app.py                                   ║
║    then open  http://localhost:5001                          ║
║                                                              ║
║  CLI report only:                                            ║
║    python3 mitoclin_app.py --report-only                     ║
║      --sample S105                                           ║
║      --csv  output_S105/10_final_report_data/S105_final_report.csv ║
║      --cov  output_S105/06_coverage/S105_coverage.txt        ║
║      --outdir output_S105/reports/                           ║
╚══════════════════════════════════════════════════════════════╝
"""

# ══════════════════════════════════════════════════════════════
# STDLIB IMPORTS
# ══════════════════════════════════════════════════════════════
import argparse, base64, csv, io, json, math, os, re, sys
import subprocess, threading, uuid, textwrap
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# ══════════════════════════════════════════════════════════════
# MATPLOTLIB (non-interactive)
# ══════════════════════════════════════════════════════════════
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np

# Suppress MatplotlibDeprecationWarnings from internal savefig kwarg handling
# (orientation, facecolor, edgecolor, bbox_inches_restore) — deprecated since 3.3
import warnings
warnings.filterwarnings(
    "ignore",
    category=matplotlib.MatplotlibDeprecationWarning
)
warnings.filterwarnings(
    "ignore",
    message=".*savefig.*unexpected keyword argument.*",
    category=UserWarning,
)

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
# ══════════════════════════════════════════════════════════════
# FLASK
# ══════════════════════════════════════════════════════════════
from flask import Flask, request, jsonify, send_file, Response

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 20 * 1024 * 1024 * 1024  # 20 GB

# ── CORS helper (no flask-cors needed)
def cors(response):
    response.headers["Access-Control-Allow-Origin"]  = "*"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type"
    response.headers["Access-Control-Allow-Methods"] = "GET,POST,DELETE,OPTIONS"
    return response

@app.after_request
def after(r): return cors(r)

@app.before_request
def preflight():
    if request.method == "OPTIONS":
        return cors(Response(status=200))

# ══════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════
BASE_DIR      = os.path.abspath(os.path.dirname(__file__))
UPLOAD_DIR    = os.path.join(BASE_DIR, "uploads")

# ── Load BRIC-CDFD logo (bric_cdfd_logo.jpeg) as base64 for embedding
def _load_logo_b64(filename="bric_cdfd_logo.jpeg") -> str:
    """Load logo file from BASE_DIR and return base64 data-URI string."""
    logo_path = os.path.join(BASE_DIR, filename)
    if os.path.isfile(logo_path):
        import base64 as _b64
        ext = filename.rsplit(".",1)[-1].lower()
        mime = "image/jpeg" if ext in ("jpg","jpeg") else f"image/{ext}"
        with open(logo_path,"rb") as _f:
            return f"data:{mime};base64,{_b64.b64encode(_f.read()).decode()}"
    # Fallback: inline SVG placeholder (shown if file not found)
    return "data:image/svg+xml;base64," + __import__("base64").b64encode(
        b'''<svg viewBox="0 0 80 80" xmlns="http://www.w3.org/2000/svg">
        <circle cx="40" cy="40" r="38" fill="#1a3a5c"/>
        <text x="40" y="36" text-anchor="middle" fill="white"
              font-size="13" font-weight="bold" font-family="Arial">BRIC</text>
        <text x="40" y="52" text-anchor="middle" fill="white"
              font-size="11" font-weight="bold" font-family="Arial">CDFD</text>
        </svg>'''
    ).decode()

BRIC_CDFD_LOGO_B64 = _load_logo_b64()
OUTPUT_DIR    = BASE_DIR
PIPELINE_SH   = os.path.join(BASE_DIR, "mtdna_pipeline_new.sh")
REFERENCE_DIR = os.path.join(BASE_DIR, "reference")
REPORTS_DIR   = os.path.join(BASE_DIR, "reports")

for d in (UPLOAD_DIR, REFERENCE_DIR, REPORTS_DIR):
    os.makedirs(d, exist_ok=True)

# In-memory job store
JOBS: dict = {}
JOBS_LOCK = threading.Lock()

# ══════════════════════════════════════════════════════════════
# MTDNA CONSTANTS
# ══════════════════════════════════════════════════════════════
MTDNA_LEN = 16569

MT_GENES = [
    (1,576,"D-loop","dloop"),(577,647,"MT-TF","trna"),
    (648,1601,"MT-RNR1","rrna"),(1602,1670,"MT-TV","trna"),
    (1671,3229,"MT-RNR2","rrna"),(3230,3304,"MT-TL1","trna"),
    (3307,4262,"MT-ND1","complex1"),(4263,4331,"MT-TI","trna"),
    (4329,4400,"MT-TQ","trna"),(4402,4469,"MT-TM","trna"),
    (4470,5511,"MT-ND2","complex1"),(5512,5579,"MT-TW","trna"),
    (5587,5655,"MT-TA","trna"),(5657,5729,"MT-TN","trna"),
    (5761,5826,"MT-TC","trna"),(5826,5891,"MT-TY","trna"),
    (5904,7445,"MT-CO1","complex4"),(7446,7514,"MT-TS1","trna"),
    (7518,7585,"MT-TD","trna"),(7586,8269,"MT-CO2","complex4"),
    (8295,8364,"MT-TK","trna"),(8366,8572,"MT-ATP8","complex5"),
    (8527,9207,"MT-ATP6","complex5"),(9207,9990,"MT-CO3","complex4"),
    (9991,10058,"MT-TG","trna"),(10059,10404,"MT-ND3","complex1"),
    (10405,10469,"MT-TR","trna"),(10470,10766,"MT-ND4L","complex1"),
    (10760,12137,"MT-ND4","complex1"),(12138,12206,"MT-TH","trna"),
    (12207,12265,"MT-TS2","trna"),(12266,12336,"MT-TL2","trna"),
    (12337,14148,"MT-ND5","complex1"),(14149,14673,"MT-ND6","complex1"),
    (14674,14742,"MT-TE","trna"),(14747,15887,"MT-CYB","complex3"),
    (15888,15953,"MT-TT","trna"),(15956,16023,"MT-TP","trna"),
    (16024,16569,"D-loop","dloop"),
]

GENE_COLORS = {
    "dloop":"#4a90d9","rrna":"#27ae60","trna":"#e67e22",
    "complex1":"#c0392b","complex3":"#8e44ad",
    "complex4":"#2980b9","complex5":"#16a085",
}
# Pathogenicity from MITOMAP (confirmed/reported/polymorphism) + MITOMASTER
# Scientific classification per MITOMAP & ClinVar mtDNA standards
PATH_COLORS = {
    "Pathogenic":        "#c0392b",   # MITOMAP confirmed pathogenic
    "Likely_Pathogenic": "#e67e22",   # MITOMAP reported pathogenic
    "Uncertain":         "#8e7cc3",   # Unknown significance
    "Benign":            "#27ae60",   # MITOMAP polymorphism / haplogroup-defining
    "NA":                "#7f8c8d",   # Not in MITOMAP
}
PATH_LABELS = {
    "Pathogenic":        "Pathogenic",
    "Likely_Pathogenic": "Likely Pathogenic",
    "Uncertain":         "Uncertain Significance",
    "Benign":            "Benign / Polymorphism",
    "NA":                "Not Annotated",
}
# MITOMAP status → internal key
MITOMAP_STATUS_MAP = {
    "pathogenic":    "Pathogenic",
    "confirmed":     "Pathogenic",
    "reported":      "Likely_Pathogenic",
    "likely":        "Likely_Pathogenic",
    "polymorphism":  "Benign",
    "benign":        "Benign",
    "uncertain":     "Uncertain",
    "unknown":       "Uncertain",
    "na":            "NA",
    "":              "NA",
}

# ══════════════════════════════════════════════════════════════
# DATA CLEANING
# ══════════════════════════════════════════════════════════════

def clean_gene(raw):
    if not raw or raw.strip() in ("—","NA",""): return "—"
    raw = raw.strip().strip('"')
    parts = [p.strip() for p in raw.split("|") if p.strip()]
    for p in parts:
        if re.match(r'^(ND\d\w*|CO[I123]+|CYB|ATP[68]|RNR[12]|TL[12]|'
                    r'TK|TH|TT|TP|TV|TW|TA|TN|TF|TC|TY|TM|TG|TS[12]|'
                    r'TD|TR|TI|TQ|12S|16S|SHLP\d|Humanin|SHMOOSE|'
                    r'ALTND\d|ATPase\d)', p, re.I):
            return p
    non_cr = [p for p in parts if not p.startswith("CR:") and p != "ATT"]
    return (non_cr[0] if non_cr else parts[0])[:14] if parts else "—"

def clean_mutation(raw, pos, ref, alt):
    if not raw or raw.strip().lower() in ("","na","—","transition","transversion","indel"):
        return f"m.{pos}{ref}>{alt}"
    return raw.strip()

def clean_protein(raw):
    if not raw or raw.strip() in ("","NA","—","non-coding"): return "—"
    return raw.strip()

def clean_disease(raw, maxlen=50):
    if not raw or raw.strip() in ("","NA","—"): return "—"
    s = raw.strip()
    return (s[:maxlen-1]+"…") if len(s)>maxlen else s

def parse_haplogroup(raw):
    if not raw: return "Unknown"
    s = raw.strip().strip('"').strip("'")
    return "Unknown" if s.upper() in ("YES","NO","NA","") else s

# ══════════════════════════════════════════════════════════════
# DATA LOADING
# ══════════════════════════════════════════════════════════════

def normalise_pathogenicity(raw):
    """Map MITOMAP/pipeline status to canonical key."""
    if not raw or raw.strip() in ("","—"): return "NA"
    s = raw.strip().lower()
    for k, v in MITOMAP_STATUS_MAP.items():
        if k and k in s: return v
    return "Uncertain"

def load_report(csv_path):
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            for fld in ("Position","Ref_Count","Alt_Count","Coverage_at_Position"):
                try: row[fld] = int(float(row.get(fld,0) or 0))
                except: row[fld] = 0
            for fld in ("Heteroplasmy_pct","Mean_Coverage"):
                try: row[fld] = float(row.get(fld,0) or 0)
                except: row[fld] = 0.0
            # Normalise pathogenicity — handle both old (VUS/Novel) and new (MITOMAP) values
            raw_path = row.get("Pathogenicity","NA")
            # Keep existing valid keys directly
            if raw_path in PATH_COLORS:
                row["_path"] = raw_path
            else:
                row["_path"] = normalise_pathogenicity(raw_path)
            row["_gene"]    = clean_gene(row.get("Gene",""))
            row["_mut"]     = clean_mutation(row.get("Mutation_Type",""), row["Position"],
                                              row.get("Ref","?"), row.get("Alt","?"))
            row["_protein"] = clean_protein(row.get("Protein_Change",""))
            row["_disease"] = clean_disease(row.get("Disease",""), maxlen=55)
            # Haplogroup from new pipeline column (column 10 in new script)
            row["_haplo"]   = parse_haplogroup(row.get("Haplogroup",""))
            # Variant type for hover tooltip
            ref = row.get("Ref",""); alt = row.get("Alt","")
            if len(ref)==1 and len(alt)==1:
                row["_var_type"] = "Transition" if (
                    {ref.upper(),alt.upper()} in [{"A","G"},{"C","T"}]) else "Transversion"
            else:
                row["_var_type"] = "Indel"
            rows.append(row)
    return rows
# ══════════════════════════════════════════════════════════════
# VCF PARSER + FULL ANNOTATION
# Mirrors the pipeline annotation steps:
#   1. Parse VCF → extract position, DP, AF, REF, ALT
#   2. Assign gene/region from MT_GENES coordinate table
#   3. Classify variant type (Transition / Transversion / Indel)
#   4. MITOMAP lookup (if local DB present)
#   5. MITOMASTER API submission for protein/haplogroup annotation
#   6. Normalise pathogenicity using MITOMAP_STATUS_MAP
# ══════════════════════════════════════════════════════════════

def _assign_gene(pos):
    """Return (gene_name, gene_type) for a given rCRS position."""
    for start, end, name, gtype in MT_GENES:
        if start <= pos <= end:
            return name, gtype
    return "—", "dloop"

def _variant_type(ref, alt):
    """Classify SNV as Transition / Transversion, or Indel."""
    if len(ref) != 1 or len(alt) != 1:
        return "Indel"
    transitions = {frozenset({"A","G"}), frozenset({"C","T"})}
    return "Transition" if frozenset({ref.upper(), alt.upper()}) in transitions else "Transversion"

def _load_mitomap_db():
    """Load mitomap_foswiki.csv from the database directory.
    Returns a dict keyed by (pos, ref, alt) → {disease, pathogenicity, protein}.
    Returns empty dict if the file is not present.
    """
    db = {}
    db_path = os.path.join(BASE_DIR, "database", "mitomap_foswiki.csv")
    if not os.path.isfile(db_path):
        return db
    try:
        with open(db_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    pos = int(float(row.get("Position", row.get("Pos", 0)) or 0))
                except Exception:
                    continue
                ref = str(row.get("Ref", row.get("ref", "")) or "").strip().upper()
                alt = str(row.get("Alt", row.get("alt", "")) or "").strip().upper()
                key = (pos, ref, alt)
                db[key] = {
                    "disease":       str(row.get("Disease", row.get("disease", "")) or "").strip(),
                    "pathogenicity": str(row.get("Pathogenicity", row.get("Status", "")) or "").strip(),
                    "protein":       str(row.get("Protein_Change", row.get("AA_Change", "")) or "").strip(),
                }
    except Exception as e:
        app.logger.warning(f"MITOMAP DB load error: {e}")
    return db

def _submit_mitomaster(variants):
    """Submit variants to MITOMASTER API using multipart/form-data.

    Matches the bash pipeline curl command:
      curl -F "file=@input.txt" -F "fileType=snvlist" -F "output=detail"
           https://www.mitomap.org/mitomaster/websrvc.cgi

    MITOMASTER detail output columns (after <br> and CSV->TSV conversion):
      col2=pos, col3=ref, col4=alt, col7=gene, col8=protein,
      col10=haplogroup, col12=disease
    Returns empty dict on any error (pipeline continues without annotation).
    """
    import urllib.request
    result = {}
    if not variants:
        return result
    try:
        # Build snvlist input (tab-separated: sample pos ref var)
        sample_name = "SAMPLE"
        rows = ["sample\tpos\tref\tvar"]
        for v in variants:
            ref = str(v.get("Ref","")).strip()
            alt = str(v.get("Alt","")).strip()
            if not ref or not alt or ref in ("-","?","") or alt in ("-","?",""):
                continue
            rows.append(f"{sample_name}\t{v['Position']}\t{ref}\t{alt}")

        if len(rows) <= 1:
            return result

        file_content = "\n".join(rows).encode("utf-8")

        # Build multipart/form-data body manually (no third-party deps)
        boundary = "MITOCLINbnd42"
        nl = b"\r\n"

        def make_field(name, value):
            return (
                b"--" + boundary.encode() + nl +
                f"Content-Disposition: form-data; name=\"{name}\"".encode() + nl + nl +
                value.encode() + nl
            )

        body = (
            make_field("fileType", "snvlist") +
            make_field("output",   "detail")  +
            b"--" + boundary.encode() + nl +
            b"Content-Disposition: form-data; name=\"file\"; filename=\"variants.txt\"" + nl +
            b"Content-Type: text/plain" + nl + nl +
            file_content + nl +
            b"--" + boundary.encode() + b"--" + nl
        )

        req = urllib.request.Request(
            "https://www.mitomap.org/mitomaster/websrvc.cgi",
            data=body,
            headers={
                "User-Agent":   "Mozilla/5.0 (MITOCLIN/2.1; Research)",
                "Content-Type": f"multipart/form-data; boundary={boundary}",
            },
        )

        with urllib.request.urlopen(req, timeout=120) as resp:
            raw = resp.read().decode("utf-8", errors="replace")

        # Reject HTML error pages
        if not raw.strip() or "<html" in raw.lower()[:300]:
            app.logger.warning("MITOMASTER returned HTML error page or empty response")
            return result

        # Mirror pipeline sed transforms:
        #   s/<br>/|/g    replace <br> field separators
        #   s/,/\t/g      convert CSV -> TSV
        raw_tsv = raw.replace("<br>", "|").replace(",", "\t")
        tsv_lines = [ln for ln in raw_tsv.splitlines() if ln.strip()]
        if not tsv_lines:
            return result

        # Parse: skip header row, then extract per-column values
        # col2=pos, col3=ref, col4=alt, col5=mutation,
        # col7=gene, col8=protein, col10=haplogroup, col12=disease
        for ln in tsv_lines[1:]:
            cols = ln.split("\t")
            if len(cols) < 5:
                continue
            pos_s = cols[1].strip() if len(cols) > 1 else ""
            if not pos_s or pos_s.lower() in ("pos", "position", ""):
                continue
            try:
                pos_int = int(pos_s)
            except ValueError:
                continue
            ref_c  = cols[2].strip().upper()  if len(cols) > 2  else ""
            alt_c  = cols[3].strip().upper()  if len(cols) > 3  else ""
            gene_c = cols[6].strip()          if len(cols) > 6  else "-"
            prot_c = cols[7].strip()          if len(cols) > 7  else "-"
            hap_c  = cols[9].strip()          if len(cols) > 9  else "Unknown"
            dis_c  = cols[11].strip()         if len(cols) > 11 else "-"
            if ref_c and alt_c:
                result[(pos_int, ref_c, alt_c)] = {
                    "gene":    gene_c or "-",
                    "protein": prot_c or "-",
                    "haplo":   hap_c  or "Unknown",
                    "disease": dis_c  or "-",
                }

        app.logger.info(f"MITOMASTER returned {len(result)} annotations")

    except Exception as e:
        app.logger.warning(f"MITOMASTER API error: {type(e).__name__}: {e}")
    return result
def load_vcf(vcf_path):
    """Parse VCF and apply full annotation:
    gene assignment, variant type, MITOMAP lookup, MITOMASTER API.
    """
    raw_variants = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_s, vid, ref, alt, qual, flt, info = parts[:8]

            # Only process PASS variants (or no filter field)
            if flt not in ("PASS", ".", ""):
                continue

            try:
                pos = int(pos_s)
            except ValueError:
                continue

            dp, af = 0, 0.0
            ad_ref, ad_alt = 0, 0

            # Parse INFO field
            for x in info.split(";"):
                if "=" not in x:
                    continue
                k, _, v = x.partition("=")
                k = k.strip(); v = v.strip()
                if k == "DP":
                    try: dp = int(v)
                    except: pass
                elif k == "AF":
                    try: af = float(v.split(",")[0])
                    except: pass

            # Parse FORMAT/SAMPLE fields for AD (allele depth) if present
            if len(parts) >= 10:
                fmt = parts[8].split(":")
                smp = parts[9].split(":")
                fmt_map = {f: i for i, f in enumerate(fmt)}
                if "AD" in fmt_map:
                    try:
                        ads = smp[fmt_map["AD"]].split(",")
                        ad_ref = int(ads[0]); ad_alt = int(ads[1]) if len(ads) > 1 else 0
                        if dp == 0:
                            dp = ad_ref + ad_alt
                        if af == 0.0 and dp > 0:
                            af = ad_alt / dp
                    except Exception:
                        pass
                if "AF" in fmt_map and af == 0.0:
                    try: af = float(smp[fmt_map["AF"]])
                    except: pass
                if "DP" in fmt_map and dp == 0:
                    try: dp = int(smp[fmt_map["DP"]])
                    except: pass

            # Skip if below 5% AF threshold (mirrors pipeline filter)
            if af < 0.05:
                continue

            ref = ref.upper(); alt = alt.upper()
            if ad_ref == 0 and ad_alt == 0:
                ad_ref = int(dp * (1 - af))
                ad_alt = int(dp * af)

            raw_variants.append({
                "Position":              pos,
                "Ref":                   ref,
                "Alt":                   alt,
                "Ref_Count":             ad_ref,
                "Alt_Count":             ad_alt,
                "Coverage_at_Position":  dp,
                "Heteroplasmy_pct":      af,
                "Mean_Coverage":         dp,
            })

    if not raw_variants:
        return []

    # ── Step 1: Gene assignment from MT_GENES coordinates ──
    for v in raw_variants:
        gene_name, gene_type = _assign_gene(v["Position"])
        v["_gene_name"] = gene_name
        v["_gene_type"] = gene_type

    # ── Step 2: MITOMAP local DB lookup ──────────────────
    mitomap_db = _load_mitomap_db()

    # ── Step 3: MITOMASTER API ────────────────────────────
    app.logger.info(f"Submitting {len(raw_variants)} variants to MITOMASTER…")
    mitomaster = _submit_mitomaster(raw_variants)
    app.logger.info(f"MITOMASTER returned {len(mitomaster)} annotations")

    # ── Step 4: Assemble final annotated variant dicts ────
    # Use first variant's haplo as sample haplogroup if MITOMASTER gives it
    sample_haplo = "Unknown"

    variants = []
    for v in raw_variants:
        pos = v["Position"]; ref = v["Ref"]; alt = v["Alt"]
        key = (pos, ref.upper(), alt.upper())

        # Gene from MITOMASTER if available, else from coordinate table
        mm  = mitomaster.get(key, {})
        mit = mitomap_db.get(key, {})

        gene_raw = mm.get("gene") or v["_gene_name"]
        gene     = clean_gene(gene_raw)

        protein_raw = mm.get("protein") or mit.get("protein") or "—"
        protein  = clean_protein(protein_raw)

        # Use MITOMASTER disease if available, fall back to MITOMAP local DB
        disease_raw = mm.get("disease") or mit.get("disease") or "—"
        disease  = clean_disease(disease_raw, maxlen=55)

        path_raw    = mit.get("pathogenicity") or "NA"
        path_key    = normalise_pathogenicity(path_raw)

        haplo = mm.get("haplo") or "Unknown"
        if haplo != "Unknown":
            sample_haplo = haplo

        mut_notation = f"m.{pos}{ref}>{alt}"
        vtype = _variant_type(ref, alt)

        variants.append({
            # Raw numeric fields
            "Position":             pos,
            "Ref":                  ref,
            "Alt":                  alt,
            "Ref_Count":            v["Ref_Count"],
            "Alt_Count":            v["Alt_Count"],
            "Coverage_at_Position": v["Coverage_at_Position"],
            "Heteroplasmy_pct":     v["Heteroplasmy_pct"],
            "Mean_Coverage":        v["Mean_Coverage"],
            # Annotation fields (raw, for tooltips)
            "Gene":                 gene_raw,
            "Mutation_Type":        mut_notation,
            "Protein_Change":       protein_raw,
            "Disease":              disease_raw,
            "Pathogenicity":        path_raw,
            "Haplogroup":           haplo,
            # Cleaned display fields
            "_path":                path_key,
            "_gene":                gene,
            "_mut":                 mut_notation,
            "_protein":             protein,
            "_disease":             disease,
            "_haplo":               haplo,
            "_var_type":            vtype,
        })

    # Backfill sample haplogroup to all variants if resolved
    if sample_haplo != "Unknown":
        for v in variants:
            if v["_haplo"] == "Unknown":
                v["_haplo"] = sample_haplo
                v["Haplogroup"] = sample_haplo

    return variants
    
def load_coverage(cov_path):
    pos, dep = [], []
    with open(cov_path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 3:
                try: pos.append(int(p[1])); dep.append(int(p[2]))
                except: pass
    return pos, dep

def list_references():
    refs = []
    for ext in (".fasta",".fa",".fna"):
        refs += [f.name for f in Path(REFERENCE_DIR).glob(f"*{ext}")]
    return sorted(refs)

# ══════════════════════════════════════════════════════════════
# PLOTS
# ══════════════════════════════════════════════════════════════

def fig_to_b64(fig):
    """Convert matplotlib figure to base64 PNG.
    Sets facecolor/edgecolor on the figure before saving — avoids
    MatplotlibDeprecationWarning on matplotlib >= 3.3 where these
    are no longer valid savefig() keyword arguments.
    """
    buf = io.BytesIO()
    bg = fig.patch.get_facecolor()
    fig.patch.set_edgecolor("none")
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                facecolor=bg, edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

def plot_histogram(variants):
    het = [v["Heteroplasmy_pct"]*100 for v in variants]
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    bins = list(range(0,101,10))
    counts, edges = np.histogram(het, bins=bins)
    clrs = ["#5dade2" if e+5<85 else "#1a5276" for e in edges[:-1]]
    bars = ax.bar(edges[:-1], counts, width=9.5, align="edge",
                  color=clrs, edgecolor="white", linewidth=0.8)
    for b,c in zip(bars,counts):
        if c>0:
            ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.1, str(c),
                    ha="center", va="bottom", fontsize=8,
                    fontweight="bold", color="#2c3e50")
    ax.set_xlabel("Heteroplasmy (%)", fontsize=9)
    ax.set_ylabel("Variant Count", fontsize=9)
    ax.set_title("Heteroplasmy Distribution", fontsize=10, fontweight="bold")
    ax.set_xlim(0,100); ax.set_xticks([0,20,40,60,80,100])
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.spines[["top","right"]].set_visible(False)
    ax.set_facecolor("#f8f9fa"); fig.patch.set_facecolor("white")
    fig.tight_layout(pad=0.8)
    return fig_to_b64(fig)

def plot_circle(variants):
    fig, ax = plt.subplots(figsize=(5.0, 5.2), subplot_kw={"projection":"polar"})
    ax.set_facecolor("white"); fig.patch.set_facecolor("white")
    outer_r, inner_r = 1.0, 0.60
    for start,end,name,gtype in MT_GENES:
        t1 = 2*math.pi*start/MTDNA_LEN
        t2 = 2*math.pi*end/MTDNA_LEN
        thetas = np.linspace(t1,t2,max(3,int((end-start)/10)))
        color  = GENE_COLORS.get(gtype,"#bdc3c7")
        ax.fill_between(thetas,[inner_r]*len(thetas),[outer_r]*len(thetas),
                        color=color,alpha=0.85,linewidth=0)
    tf = np.linspace(0,2*math.pi,300)
    ax.fill_between(tf,0,inner_r-0.02,color="white",zorder=2)
    ax.plot(tf,[outer_r]*300,color="#2c3e50",linewidth=0.8,zorder=3)
    ax.plot(tf,[inner_r]*300,color="#2c3e50",linewidth=0.5,zorder=3)
    for v in variants:
        theta = 2*math.pi*v["Position"]/MTDNA_LEN
        p     = v.get("_path","NA")
        color = PATH_COLORS.get(p,"#95a5a6")
        size  = 80 if p in ("Pathogenic","Likely_Pathogenic") else 22
        ax.scatter(theta, outer_r+0.12, s=size, color=color,
                   zorder=10, edgecolors="white", linewidths=0.5)
    for pos in [1,2000,4000,6000,8000,10000,12000,14000,16000]:
        theta = 2*math.pi*pos/MTDNA_LEN
        ax.text(theta,inner_r-0.10,str(pos),ha="center",va="center",fontsize=5,color="#666")
    large = {"MT-ND1","MT-ND2","MT-ND4","MT-ND5","MT-CO1","MT-CO2",
             "MT-CO3","MT-CYB","MT-ATP6","MT-RNR1","MT-RNR2"}
    for start,end,name,gtype in MT_GENES:
        if name in large:
            theta = 2*math.pi*((start+end)/2)/MTDNA_LEN
            ax.text(theta,(inner_r+outer_r)/2,name.replace("MT-",""),
                    ha="center",va="center",fontsize=5,
                    fontweight="bold",color="white",zorder=5)
    ax.text(0,0,"mtDNA\n16,569 bp",ha="center",va="center",
            fontsize=8,fontweight="bold",color="#2c3e50",zorder=6)
    ax.set_ylim(0,outer_r+0.30); ax.axis("off")
    handles = [mpatches.Patch(color=GENE_COLORS[k], label=l) for k,l in [
        ("rrna","rRNA"),("trna","tRNA"),("complex1","Complex I"),
        ("complex3","Complex III"),("complex4","Complex IV"),
        ("complex5","Complex V"),("dloop","D-loop")]]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5,-0.05),
              ncol=4, fontsize=5.5, frameon=False,
              handlelength=1.0, columnspacing=0.5)
    fig.tight_layout(pad=0.2)
    return fig_to_b64(fig)

def plot_coverage(positions, depths, variants):
    fig, ax = plt.subplots(figsize=(11, 3.2))
    if len(depths)>200:
        w = max(1,len(depths)//300)
        ds = np.convolve(depths,np.ones(w)/w,mode="same")
    else:
        ds = np.array(depths)
    pa = np.array(positions)
    ax.fill_between(pa,ds,alpha=0.15,color="#2980b9",zorder=1)
    ax.plot(pa,ds,color="#1a5276",linewidth=0.7,zorder=2)
    for v in variants:
        pos = v["Position"]; p = v.get("_path","NA")
        color = PATH_COLORS.get(p,"#95a5a6")
        lw = 1.5 if p in ("Pathogenic","Likely_Pathogenic") else 0.6
        ax.axvline(x=pos,color=color,linewidth=lw,alpha=0.7,zorder=3)
    ax.set_xlabel("mtDNA Position",fontsize=9)
    ax.set_ylabel("Read Depth (×)",fontsize=9)
    ax.set_title("Coverage Across Mitochondrial Genome",fontsize=10,fontweight="bold",loc="left")
    ax.set_xlim(1,MTDNA_LEN)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,_:f"{int(x):,}"))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,_:f"{int(x):,}"))
    ax.spines[["top","right"]].set_visible(False)
    ax.set_facecolor("#f8f9fa"); fig.patch.set_facecolor("white")
    fig.tight_layout(pad=0.8)
    return fig_to_b64(fig)

def plot_pie(variants):
    counts = defaultdict(int)
    for v in variants: counts[v.get("_path","NA")] += 1
    order = ["Pathogenic","Likely_Pathogenic","Uncertain","Benign","NA"]
    labels,sizes,clrs = [],[],[]
    for k in order:
        if counts.get(k,0)>0:
            labels.append(f"{PATH_LABELS[k]} ({counts[k]})")
            sizes.append(counts[k]); clrs.append(PATH_COLORS[k])
    fig, ax = plt.subplots(figsize=(4.5, 4.0))
    wedges,texts,autotexts = ax.pie(
        sizes, labels=None, colors=clrs,
        autopct=lambda p: f"{p:.0f}%" if p>=5 else "",
        startangle=90, wedgeprops={"edgecolor":"white","linewidth":1.5},
        pctdistance=0.72)
    for at in autotexts:
        at.set_fontsize(8); at.set_color("white"); at.set_fontweight("bold")
    ax.legend(wedges,labels,loc="lower center",bbox_to_anchor=(0.5,-0.18),
              fontsize=7,frameon=False,ncol=2)
    ax.set_title("Pathogenicity Classification",fontsize=9,fontweight="bold",pad=4)
    fig.patch.set_facecolor("white")
    fig.tight_layout(pad=0.5)
    return fig_to_b64(fig)

# ══════════════════════════════════════════════════════════════
# CLINICAL INTERPRETATION
# ══════════════════════════════════════════════════════════════

def clinical_interp(variants, haplogroup):
    path_v = [v for v in variants if v.get("_path") in ("Pathogenic","Likely_Pathogenic")]
    vus_v  = [v for v in variants if v.get("_path") == "Uncertain"]
    if not path_v and not vus_v:
        return ("No pathogenic or likely pathogenic variants identified. "
                "All detected variants are classified as benign or novel. "
                "Haplogroup-defining variants are consistent with the assigned haplogroup.")
    parts = []
    if path_v:
        names    = " and ".join(f"<em>{v['_mut']}</em>" for v in path_v)
        diseases = ", ".join(v.get("Disease","").strip() for v in path_v
                             if v.get("Disease","").strip() not in ("","—","NA"))
        parts.append(f"Detected pathogenic variant(s) {names}"
                     + (f" associated with {diseases}." if diseases else "."))
        if any(v["Heteroplasmy_pct"]<0.99 for v in path_v):
            parts.append("Heteroplasmy levels suggest potential clinical relevance.")
    if vus_v:
        parts.append(f"Variant(s) of uncertain significance: "
                     + ", ".join(f"<em>{v['_mut']}</em>" for v in vus_v)
                     + ". Clinical correlation recommended.")
    parts.append("Correlation with clinical phenotype is recommended.")
    return " ".join(parts)

# ══════════════════════════════════════════════════════════════
# HTML REPORT TEMPLATE
# ══════════════════════════════════════════════════════════════

BADGE_CSS = {
    "Pathogenic":        "bp",
    "Likely_Pathogenic": "blp",
    "Uncertain":         "bu",
    "Benign":            "bb",
    "NA":                "bna",
}

REPORT_HTML = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>MITOCLIN Report — {sample}</title>
<style>
@import url('https://fonts.googleapis.com/css2?family=Source+Sans+3:wght@300;400;600;700&family=Source+Code+Pro:wght@400;600&display=swap');
*{{box-sizing:border-box;margin:0;padding:0}}
:root{{
  --navy:#0b2135;--blue:#1163a8;--blue2:#1a8fc1;--ice:#d4eaf7;
  --text:#1c2e3e;--muted:#5a7896;--border:#c8dde8;
  --success:#1a6e3a;--warn:#8a5a00;--danger:#8b1a1a;
  --radius:8px;
}}
body{{font-family:'Source Sans 3',sans-serif;font-size:10pt;color:var(--text);
      background:#fff;max-width:1010px;margin:0 auto;padding:22px 28px;line-height:1.5}}
/* HEADER */
.rpt-header{{display:flex;justify-content:space-between;align-items:center;
             padding-bottom:12px;border-bottom:3px solid var(--navy);margin-bottom:12px}}
.cdfd-logo{{display:flex;align-items:center;gap:10px}}
.cdfd-badge{{width:52px;height:52px;border-radius:50%;overflow:hidden;
             display:flex;align-items:center;justify-content:center;
             flex-shrink:0;background:#f0f4f8}}
.cdfd-org{{font-size:9pt;font-weight:700;color:var(--navy)}}
.cdfd-sub{{font-size:7pt;color:var(--muted)}}
.brand{{font-size:20pt;font-weight:700;color:var(--navy);letter-spacing:1px;text-align:right}}
.brand span{{color:var(--blue2)}}
.brand-tag{{font-size:8pt;color:var(--muted);text-align:right}}
.brand-gh{{font-size:7.5pt;color:var(--blue2);text-decoration:none;display:block;text-align:right}}
/* META */
.meta-bar{{display:flex;flex-wrap:wrap;font-size:8pt;color:var(--muted);
           background:#f5f8fc;border:1px solid var(--border);border-radius:var(--radius);
           padding:7px 12px;margin-bottom:12px;gap:0}}
.mi{{padding:0 12px;border-right:1px solid var(--border)}}
.mi:last-child{{border-right:none}}
.mi b{{color:var(--text)}}
/* STAT CARDS */
.stat-row{{display:flex;gap:7px;margin-bottom:12px;flex-wrap:wrap}}
.sc{{flex:1;min-width:90px;border-radius:var(--radius);padding:9px 6px;text-align:center;color:#fff}}
.sv{{font-size:18pt;font-weight:700;line-height:1.1;font-family:'Source Code Pro',monospace}}
.sl{{font-size:7pt;opacity:.85;margin-top:2px}}
.c1{{background:linear-gradient(135deg,#0d5e8a,#1a8fc1)}}
.c2{{background:linear-gradient(135deg,#8b1a1a,#b92d27)}}
.c3{{background:linear-gradient(135deg,#6a4c00,#c47e00)}}
.c4{{background:linear-gradient(135deg,#4a235a,#8e44ad)}}
.c5{{background:linear-gradient(135deg,#1a6e3a,#27ae60)}}
.c6{{background:linear-gradient(135deg,#3d3d3d,#717171)}}
/* SECTIONS */
.sec{{margin:13px 0}}
.sec-title{{font-size:11pt;font-weight:700;color:var(--navy);
            border-bottom:2px solid var(--ice);padding-bottom:3px;margin-bottom:9px;
            display:flex;align-items:center;gap:7px}}
.sec-title::before{{content:'';width:4px;height:16px;background:var(--blue2);
                    border-radius:2px;flex-shrink:0}}
/* QC TABLE */
.qc-table{{width:100%;border-collapse:collapse;font-size:8.5pt;margin-bottom:4px}}
.qc-table th,.qc-table td{{padding:5px 9px;border:1px solid var(--border)}}
.qc-table thead tr{{background:#f0f5fa}}
.qc-table thead th{{font-weight:600;color:var(--navy)}}
.qc-table tbody tr:nth-child(even){{background:#f9fbfd}}
.chip{{display:inline-flex;align-items:center;gap:3px;padding:2px 8px;
       border-radius:10px;font-size:7.5pt;font-weight:600;border:1px solid currentColor}}
.chip-ok{{color:#1a6e3a;background:rgba(26,110,58,0.08)}}
.chip-warn{{color:#8a5a00;background:rgba(138,90,0,0.08)}}
.chip-err{{color:#8b1a1a;background:rgba(139,26,26,0.08)}}
/* INTERP */
.ibox{{background:#edf5fb;border-left:4px solid var(--blue2);
       padding:11px 15px;border-radius:0 var(--radius) var(--radius) 0;
       font-size:9.5pt;line-height:1.65}}
.hap-row{{display:flex;align-items:center;gap:10px;margin-bottom:7px;flex-wrap:wrap}}
.hap-badge{{background:var(--navy);color:#fff;padding:3px 13px;border-radius:20px;
            font-size:9pt;font-weight:700}}
.qc-chips{{display:flex;gap:5px;flex-wrap:wrap;margin-top:8px}}
/* PLOTS */
.plots-row{{display:flex;gap:11px;align-items:flex-start;margin-bottom:9px}}
.pb{{flex:1;text-align:center}}
.pb img,.cov-img img{{width:100%;border-radius:var(--radius);border:1px solid var(--border);
                      transition:box-shadow 0.2s,transform 0.2s}}
.pb img:hover{{box-shadow:0 6px 20px rgba(17,99,168,0.20);transform:scale(1.01)}}
.cov-img img{{margin-bottom:9px}}
.cov-img img:hover{{box-shadow:0 4px 16px rgba(17,99,168,0.18)}}
/* VARIANT TABLE */
.tbl-wrap{{overflow-x:auto;border-radius:var(--radius);border:1px solid var(--border);margin-bottom:4px}}
table{{width:100%;border-collapse:collapse;font-size:7.8pt}}
thead tr{{background:var(--navy);color:#fff}}
thead th{{padding:7px 6px;text-align:left;font-weight:600;font-size:7.5pt;white-space:nowrap}}
tbody tr{{transition:background 0.12s}}
tbody tr:nth-child(even){{background:#f3f8fc}}
tbody tr:hover{{background:#cde4f2!important;cursor:default}}
tbody td{{padding:5px 6px;border-bottom:1px solid #e0eaf2;vertical-align:middle;word-break:break-word}}
/* tooltip */
.has-tip{{position:relative;cursor:help;border-bottom:1px dashed var(--muted)}}
.has-tip .tip{{display:none;position:absolute;bottom:calc(100% + 4px);left:0;z-index:200;
               background:#1e2d3e;color:#fff;padding:6px 10px;border-radius:6px;
               font-size:7.5pt;min-width:160px;max-width:300px;white-space:normal;
               box-shadow:0 4px 14px rgba(0,0,0,0.28);line-height:1.4}}
.has-tip:hover .tip{{display:block}}
/* badges */
.badge{{display:inline-block;padding:2px 8px;border-radius:10px;font-size:7.5pt;
        font-weight:700;color:#fff;white-space:nowrap}}
.bp{{background:#b92d27}}.blp{{background:#d4811a}}.bu{{background:#7b5ea7}}
.bb{{background:#27ae60}}.bna{{background:#7f8c8d}}
/* het bar */
.hw{{background:#dde8f0;border-radius:4px;height:8px;width:58px;
     display:inline-block;vertical-align:middle;overflow:hidden}}
.hf{{height:100%;border-radius:4px;background:linear-gradient(90deg,#1a8fc1,#0d5e8a)}}
.hf.hi{{background:linear-gradient(90deg,#c0392b,#8b1a1a)}}
/* var type */
.vt{{font-size:7pt;padding:1px 5px;border-radius:4px;
     background:#e8f0f8;color:#1163a8;font-family:'Source Code Pro',monospace}}
.vt.trv{{background:#f5eaf8;color:#7b2fa8}}
/* methodology */
.mbox{{background:#f8fafc;border:1px solid var(--border);border-radius:var(--radius);
       padding:10px 14px;font-size:8.5pt;line-height:1.7;color:#444}}
.ref-list{{font-size:7.5pt;color:var(--muted);margin-top:6px;padding-left:16px}}
.ref-list li{{margin-bottom:2px}}
/* FOOTER */
.rpt-footer{{margin-top:16px;border-top:2px solid var(--navy);padding-top:10px}}
.footer-top{{display:flex;justify-content:space-between;flex-wrap:wrap;gap:8px;font-size:7.5pt;color:var(--muted)}}
.team-title{{font-size:8pt;font-weight:700;color:var(--navy);margin-bottom:3px}}
.footer-disc{{max-width:400px;text-align:right;font-size:7.5pt;color:var(--muted)}}
.footer-bottom{{text-align:center;margin-top:8px;font-size:6.5pt;color:#aaa;font-style:italic}}
a{{color:var(--blue2)}}
</style>
</head>
<body>

<!-- HEADER -->
<div class="rpt-header">
  <div class="cdfd-logo">
    <div class="cdfd-badge">
      <img src="{bric_logo}" alt="BRIC-CDFD Logo"
           style="width:48px;height:48px;object-fit:contain;border-radius:50%"/>
    </div>
    <div>
      <div class="cdfd-org">Centre for DNA Fingerprinting and Diagnostics</div>
      <div class="cdfd-sub">Hyderabad, India &nbsp;&middot;&nbsp; Ministry of Science &amp; Technology</div>
    </div>
  </div>
  <div>
    <div class="brand">MITO<span>CLIN</span></div>
    <div class="brand-tag">Mitochondrial Genome Analysis Platform</div>
    <a href="https://github.com/lgi/mitoclin" target="_blank" class="brand-gh">
      &#128279; github.com/lgi/mitoclin
    </a>
  </div>
</div>

<!-- META BAR -->
<div class="meta-bar">
  <span class="mi">Sample: <b>{sample}</b></span>
  <span class="mi">Case: <b>{case_id}</b></span>
  <span class="mi">Date: <b>{date}</b></span>
  <span class="mi">Platform: <b>{platform} (Illumina)</b></span>
  <span class="mi">Reference: <b>rCRS &mdash; NC_012920.1</b></span>
  <span class="mi">Pipeline: <b>MITOCLIN v1.0</b></span>
</div>

<!-- STAT CARDS -->
<div class="stat-row">
  <div class="sc c1"><div class="sv">{mean_cov}</div><div class="sl">Mean Coverage</div></div>
  <div class="sc {c_path}"><div class="sv">{n_path}</div><div class="sl">Pathogenic</div></div>
  <div class="sc c3"><div class="sv">{n_lpath}</div><div class="sl">Likely Pathogenic</div></div>
  <div class="sc c4"><div class="sv">{n_unc}</div><div class="sl">Uncertain Significance</div></div>
  <div class="sc c5"><div class="sv">{n_ben}</div><div class="sl">Benign / Polymorphism</div></div>
  <div class="sc c6"><div class="sv">{n_total}</div><div class="sl">Total Variants</div></div>
</div>

<!-- QC SUMMARY -->
<div class="sec">
  <div class="sec-title">Quality Control Summary</div>
  <table class="qc-table">
    <thead><tr><th>Metric</th><th>Value</th><th>Threshold</th><th>Status</th></tr></thead>
    <tbody>
      <tr><td>Mean Depth of Coverage</td><td><b>{mean_cov_f}&times;</b></td>
          <td>&ge;1,000&times; recommended for heteroplasmy detection</td><td>{cov_status}</td></tr>
      <tr><td>Contamination Score (HaploCheck)</td><td><b>{contam}</b></td>
          <td>&lt;0.02 &mdash; clean sample</td><td>{contam_status}</td></tr>
      <tr><td>Heteroplasmy Detection Threshold</td><td><b>5% AF</b></td>
          <td>Allele fraction &ge;0.05 (GATK recommendation)</td>
          <td><span class="chip chip-ok">&#10003; Applied</span></td></tr>
      <tr><td>Variant Caller</td><td><b>GATK Mutect2 (mtDNA mode)</b></td>
          <td>PASS filter + mitochondria-mode</td>
          <td><span class="chip chip-ok">&#10003; Standard</span></td></tr>
      <tr><td>Total SNVs Detected</td><td><b>{n_total}</b></td>
          <td>Filtered: AF&ge;5%, PASS, SNVs only</td>
          <td><span class="chip chip-ok">&#10003; Filtered</span></td></tr>
    </tbody>
  </table>
</div>

<!-- CLINICAL INTERPRETATION -->
<div class="sec">
  <div class="sec-title">Clinical Interpretation</div>
  <div class="ibox">
    <div class="hap-row">
      <span class="hap-badge">Haplogroup: {haplogroup}</span>
      <span style="font-size:8pt;color:var(--muted)">Assigned by HaploCheck v1.3</span>
      <span style="font-size:8pt;color:var(--muted)">&nbsp;|&nbsp; Contamination: {contam}</span>
    </div>
    {interp}
    <div class="qc-chips">
      {path_chip}{lpath_chip}{unc_chip}{ben_chip}
    </div>
  </div>
</div>

<!-- PLOTS -->
<div class="sec">
  <div class="sec-title">Variant Analysis &nbsp;<span style="font-weight:400;font-size:8.5pt;color:var(--muted)">(hover charts to enlarge)</span></div>
  <div class="plots-row">
    <div class="pb">
      <img src="data:image/png;base64,{img_hist}" alt="Heteroplasmy Distribution"
           title="Heteroplasmy distribution: binned allele fractions (0–100%). Majority of haplogroup-defining variants appear at 100%."/>
    </div>
    <div class="pb">
      <img src="data:image/png;base64,{img_circle}" alt="mtDNA Genome Map"
           title="Circular mitochondrial genome (16,569 bp). Variant positions marked; colour indicates pathogenicity. Gene regions colour-coded by complex."/>
    </div>
    <div class="pb">
      <img src="data:image/png;base64,{img_pie}" alt="Pathogenicity Classification"
           title="MITOMAP-based pathogenicity classification of all detected variants."/>
    </div>
  </div>
</div>

<!-- COVERAGE PLOT -->
<div class="sec cov-img">
  <img src="data:image/png;base64,{img_cov}" alt="Coverage across mitochondrial genome"
       title="Read depth across the 16,569 bp mitochondrial genome. Vertical lines mark variant positions coloured by pathogenicity. Dip regions may indicate difficult-to-map loci."/>
</div>

<!-- VARIANT TABLE -->
<div class="sec">
  <div class="sec-title">Variant Summary &nbsp;<span style="font-weight:400;font-size:8.5pt;color:var(--muted)">{n_total} variants &middot; hover cells for details</span></div>
  <div class="tbl-wrap">
  <table>
    <thead><tr>
      <th>Position (rCRS)</th>
      <th>Gene / Region</th>
      <th>Variant Notation</th>
      <th>Var. Type</th>
      <th>Protein Change</th>
      <th>Heteroplasmy %</th>
      <th>Read Depth</th>
      <th>Ref → Alt</th>
      <th>Pathogenicity</th>
      <th>Disease Association</th>
      <th>Haplogroup</th>
    </tr></thead>
    <tbody>{rows}</tbody>
  </table>
  </div>
  <div style="font-size:7.5pt;color:var(--muted);margin-top:4px">
    <b>Pathogenicity:</b>
    <span style="color:#b92d27">&bull; Pathogenic</span> = MITOMAP confirmed;
    <span style="color:#d4811a">&bull; Likely Pathogenic</span> = MITOMAP reported;
    <span style="color:#7b5ea7">&bull; Uncertain Significance</span> = not classified;
    <span style="color:#27ae60">&bull; Benign / Polymorphism</span> = haplogroup-defining or MITOMAP polymorphism;
    <span style="color:#7f8c8d">&bull; Not Annotated</span> = absent from MITOMAP.
    Heteroplasmy &ge;95% indicates likely homoplasmic status.
  </div>
</div>

<!-- METHODOLOGY -->
<div class="sec">
  <div class="sec-title">Methodology</div>
  <div class="mbox">
    <b>Read QC &amp; Trimming:</b> FastQC for quality assessment; Trim Galore
    (quality &ge;30, minimum length 30 bp, adapter auto-detection).<br/>
    <b>Alignment:</b> BWA-MEM to rCRS (NC_012920.1); GATK MarkDuplicates for PCR duplicate removal.<br/>
    <b>Variant Calling:</b> GATK Mutect2 in mitochondria-mode with allele fraction &ge;5% filter
    and FilterMutectCalls (PASS only, SNVs only, split multiallelic sites with BCFtools norm).<br/>
    <b>Heteroplasmy:</b> Allele fraction = Alt reads / (Ref + Alt reads) from AD FORMAT field.<br/>
    <b>Haplogroup &amp; Contamination:</b> HaploCheck v1.3 using Phylotree Build 17.<br/>
    <b>Annotation:</b> MITOMAP database (Foswiki CSV) for disease association and
    pathogenicity; MITOMASTER API for functional and amino acid annotation.<br/>
    <ul class="ref-list">
      <li>Van Oven M &amp; Kayser M. (2009). Updated comprehensive phylogenetic tree of global human mitochondrial DNA variation. <em>Hum Mutat.</em> 30(2):E386&ndash;394.</li>
      <li>MITOMAP: A Human Mitochondrial Genome Database. <a href="https://www.mitomap.org">mitomap.org</a></li>
      <li>McKenna A et al. (2010). The Genome Analysis Toolkit: A MapReduce framework. <em>Genome Res.</em> 20:1297&ndash;1303.</li>
      <li>Weissensteiner H et al. (2021). Haplocheck: contamination detection in mtDNA studies. <em>Genome Res.</em> 31:1723&ndash;1732.</li>
      <li>Li H &amp; Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Aligner. <em>Bioinformatics.</em> 25(14):1754&ndash;1760.</li>
    </ul>
  </div>
</div>

<!-- FOOTER -->
<div class="rpt-footer">
  <div class="footer-top">
    <div>
      <div class="team-title">&#128101; Analysis Team</div>
      <div>Dr. Ajay Kumar Mahato &nbsp;&middot;&nbsp; Ms. Aastha Yadav</div>
      <div style="color:#aaa;font-size:7pt;margin-top:2px">
        Centre for DNA Fingerprinting &amp; Diagnostics (CDFD), Hyderabad
      </div>
      <div style="margin-top:4px">
        <a href="https://github.com/lgi/mitoclin" target="_blank" style="font-size:7.5pt">
          &#128279; github.com/lgi/mitoclin
        </a>
      </div>
    </div>
    <div class="footer-disc">
      <div><b>&#9888; Research Use Only</b></div>
      <div style="margin-top:3px;line-height:1.5">
        This report is generated for research purposes only and is not
        intended for direct clinical diagnosis without expert clinical review.
        Results should be interpreted in conjunction with clinical phenotype.
      </div>
      <div style="margin-top:5px">
        &copy; 2026 CDFD &amp; MITOCLIN Team. All rights reserved.
      </div>
    </div>
  </div>
  <div class="footer-bottom">
    Generated by MITOCLIN v1.0 &nbsp;|&nbsp; {date} &nbsp;|&nbsp;
    BWA-MEM + GATK Mutect2 + HaploCheck + MITOMAP &nbsp;|&nbsp;
    rCRS (NC_012920.1) &nbsp;|&nbsp; For queries: contact CDFD Genome Informatics
  </div>
</div>

</body></html>
"""


def build_html_rows(variants):
    """Build HTML table rows with correct Ref/Alt alleles (production-safe)."""
    rows = []

    for v in variants:
        p      = v.get("_path", "NA")
        bc     = BADGE_CSS.get(p, "bna")
        lbl    = PATH_LABELS.get(p, p)

        het    = float(v.get("Heteroplasmy_pct", 0))
        hi     = "hi" if het >= 0.5 else ""

        vtype  = v.get("_var_type", "SNV")
        vt_cls = "trv" if "Transversion" in vtype else ""

        ref_c  = int(v.get("Ref_Count", 0))
        alt_c  = int(v.get("Alt_Count", 0))
        cov    = int(v.get("Coverage_at_Position", 0))

        haplo  = v.get("_haplo") or v.get("Haplogroup") or "—"

        gene   = v.get("_gene", "—")
        mut    = v.get("_mut", "")
        prot   = v.get("_protein", "—")
        dis    = v.get("_disease", "—")

        raw_dis  = v.get("Disease", "") or ""
        raw_gene = v.get("Gene", "") or ""

        # Tooltips
        dis_html = (
            f"<span class='has-tip'>{dis}<span class='tip'>{raw_dis}</span></span>"
            if len(raw_dis) > 35 else dis
        )

        gene_html = (
            f"<span class='has-tip'>{gene}<span class='tip'>{raw_gene}</span></span>"
            if "|" in raw_gene else gene
        )

        # --- FIX: Extract REF / ALT safely ---
        ref = v.get("Ref")
        alt = v.get("Alt")

        if not ref or not alt or ref == "?" or alt == "?":
            # fallback from mutation string like m.73A>G
            match = re.search(r'([ACGT])>([ACGT])', mut.upper())
            if match:
                ref, alt = match.group(1), match.group(2)
            else:
                ref, alt = "—", "—"

        # Display — clean alleles only, no read counts
        ref_alt_display = f"{ref}&nbsp;→&nbsp;{alt}"

        rows.append(
            f"<tr>"
            f"<td><b>{v.get('Position','—')}</b></td>"
            f"<td>{gene_html}</td>"
            f"<td><em>{mut}</em></td>"
            f"<td><span class='vt {vt_cls}'>{vtype[:3]}</span></td>"
            f"<td>{prot}</td>"
            f"<td><span class='hw'><span class='hf {hi}' style='width:{het*100:.0f}%'></span></span>&nbsp;<b>{het*100:.1f}%</b></td>"
            f"<td>{cov:,}</td>"
            f"<td style='font-family:monospace;font-size:7pt'>{ref_alt_display}</td>"
            f"<td><span class='badge {bc}'>{lbl}</span></td>"
            f"<td>{dis_html}</td>"
            f"<td>{haplo}</td>"
            f"</tr>"
        )

    return "\n".join(rows)

def make_chips(n_path, n_lpath, n_unc, n_ben):
    chips = []
    if n_path>0:  chips.append(f'<span class="chip" style="color:#b92d27;background:rgba(185,45,39,0.08);border-color:#b92d27">&#9679; {n_path} Pathogenic</span>')
    if n_lpath>0: chips.append(f'<span class="chip" style="color:#d4811a;background:rgba(212,129,26,0.08);border-color:#d4811a">&#9679; {n_lpath} Likely Pathogenic</span>')
    if n_unc>0:   chips.append(f'<span class="chip" style="color:#7b5ea7;background:rgba(123,94,167,0.08);border-color:#7b5ea7">&#9679; {n_unc} Uncertain</span>')
    if n_ben>0:   chips.append(f'<span class="chip" style="color:#1a6e3a;background:rgba(26,110,58,0.08);border-color:#1a6e3a">&#9679; {n_ben} Benign</span>')
    return " ".join(chips)

def make_cov_status(mean_cov):
    if mean_cov>=1000: return '<span class="chip chip-ok">&#10003; Excellent</span>'
    if mean_cov>=100:  return '<span class="chip chip-warn">&#9651; Acceptable</span>'
    return '<span class="chip chip-err">&#9888; Low</span>'

def make_contam_status(contam_str):
    try:
        val = float(str(contam_str).strip('"\''))
        if val<0.02: return '<span class="chip chip-ok">&#10003; Clean</span>'
        if val<0.05: return '<span class="chip chip-warn">&#9651; Borderline</span>'
        return '<span class="chip chip-err">&#9888; Contaminated</span>'
    except: return '<span class="chip">N/A</span>'

def render_html_report(variants, positions, depths,
                       sample, case_id, date_str, platform,
                       haplogroup, contamination, mean_cov,
                       img_hist, img_circle, img_cov, img_pie):
    n_path  = sum(1 for v in variants if v.get("_path")=="Pathogenic")
    n_lpath = sum(1 for v in variants if v.get("_path")=="Likely_Pathogenic")
    n_unc   = sum(1 for v in variants if v.get("_path")=="Uncertain")
    n_ben   = sum(1 for v in variants if v.get("_path")=="Benign")
    chips   = make_chips(n_path, n_lpath, n_unc, n_ben)
    return REPORT_HTML.format(
        sample=sample, case_id=case_id, date=date_str, platform=platform,
        haplogroup=haplogroup, contamination=contamination,
        mean_cov=f"{mean_cov:.0f}\u00d7",
        mean_cov_f=f"{mean_cov:.2f}",
        n_path=n_path, n_lpath=n_lpath, n_unc=n_unc, n_ben=n_ben,
        n_total=len(variants), contam=contamination,
        c_path="c2" if n_path>0 else ("c3" if n_lpath>0 else "c5"),
        cov_status=make_cov_status(mean_cov),
        contam_status=make_contam_status(contamination),
        path_chip=chips, lpath_chip="", unc_chip="", ben_chip="",
        interp=clinical_interp(variants, haplogroup),
        rows=build_html_rows(variants),
        img_hist=img_hist, img_circle=img_circle,
        img_cov=img_cov, img_pie=img_pie,
        bric_logo=BRIC_CDFD_LOGO_B64,
    )

# ══════════════════════════════════════════════════════════════
# PDF REPORT  (ReportLab — no overflow, proper image sizing)
# ══════════════════════════════════════════════════════════════

def render_pdf_report(variants, positions, depths,
                      sample, case_id, date_str, platform,
                      haplogroup, contamination, mean_cov,
                      img_hist, img_circle, img_cov, img_pie,
                      out_path):

    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm, cm
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table,
        TableStyle, HRFlowable, Image as RLImg, KeepTogether
    )
    from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT

    W, H   = A4
    LM = RM = 15*mm
    TM = BM = 14*mm
    AVAIL  = W - LM - RM          # ≈ 165 mm

    C = colors.HexColor

    doc = SimpleDocTemplate(out_path, pagesize=A4,
        leftMargin=LM, rightMargin=RM,
        topMargin=TM, bottomMargin=BM,
        title=f"MITOCLIN Report – {sample}")

    def P(txt, fs=8.5, bold=False, italic=False, color="#1e2733",
          align=TA_LEFT, leading=None, space_before=0, space_after=2,
          bg=None, left_i=0, right_i=0):
        fn = ("Helvetica-BoldOblique" if bold and italic else
              "Helvetica-Bold"   if bold   else
              "Helvetica-Oblique" if italic else "Helvetica")
        kw = dict(fontName=fn, fontSize=fs,
                  textColor=C(color) if isinstance(color,str) else color,
                  alignment=align,
                  leading=leading or (fs*1.35),
                  spaceBefore=space_before, spaceAfter=space_after,
                  leftIndent=left_i, rightIndent=right_i)
        if bg: kw["backColor"] = C(bg) if isinstance(bg,str) else bg
        style = ParagraphStyle("_", parent=getSampleStyleSheet()["Normal"], **kw)
        return Paragraph(txt, style)

    def img_from_b64(b64, w_mm, h_mm):
        """Decode base64 PNG into ReportLab Image.
        Safely handles empty/None b64 (returns a blank Spacer).
        Also strips 'data:image/...;base64,' URI prefix if present.
        """
        if not b64:
            return Spacer(w_mm * mm, h_mm * mm)
        try:
            # Strip data-URI prefix if present (e.g. from BRIC logo)
            raw_b64 = b64.split(",", 1)[1] if b64.startswith("data:") else b64
            buf = io.BytesIO(base64.b64decode(raw_b64))
            buf.seek(0)
            im = RLImg(buf)
            im.drawWidth  = w_mm * mm
            im.drawHeight = h_mm * mm
            return im
        except Exception as _img_err:
            app.logger.warning(f"img_from_b64 failed ({_img_err}) — using blank space")
            return Spacer(w_mm * mm, h_mm * mm)

    n_path  = sum(1 for v in variants if v.get("_path")=="Pathogenic")
    n_lpath = sum(1 for v in variants if v.get("_path")=="Likely_Pathogenic")
    n_unc   = sum(1 for v in variants if v.get("_path")=="Uncertain")
    n_ben   = sum(1 for v in variants if v.get("_path")=="Benign")
    n_vus   = n_unc  # legacy compat

    story = []

    # ── HEADER ──────────────────────────────────────────
    # Logo + title side by side
    # ── HEADER: Logo + Org (left) | MITOCLIN brand (right) — mirrors HTML ──
    logo_path = os.path.join(BASE_DIR, "bric_cdfd_logo.jpeg")
    if os.path.isfile(logo_path):
        logo_img = RLImg(logo_path)
        logo_img.drawWidth  = 14 * mm
        logo_img.drawHeight = 14 * mm
        logo_cell = logo_img
    else:
        logo_cell = P("BRIC<br/>CDFD", fs=8, bold=True, color="#ffffff",
                      bg="#1a3a5c", align=TA_CENTER)

    # Left side: logo + organisation name
    left_col = Table([[
        logo_cell,
        [
            P("Centre for DNA Fingerprinting and Diagnostics",
              fs=8, bold=True, color="#0b2135", leading=10),
            P("Hyderabad, India  ·  Ministry of Science &amp; Technology",
              fs=6.5, color="#5a7896", leading=9),
        ]
    ]], colWidths=[16*mm, AVAIL - 16*mm - 70*mm])
    left_col.setStyle(TableStyle([
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("LEFTPADDING",  (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (-1,-1), 4),
        ("GRID", (0,0), (-1,-1), 0, colors.white),
    ]))

    # Right side: MITOCLIN brand block
    right_col = Table([[
        P("MITO<font color=\"#1a8fc1\">CLIN</font>",
          fs=18, bold=True, color="#0b2135", align=TA_RIGHT, leading=20),
    ], [
        P("Mitochondrial Genome Analysis Platform",
          fs=7, color="#5a7896", align=TA_RIGHT, leading=9),
    ], [
        P("github.com/lgi/mitoclin",
          fs=7, color="#1a8fc1", align=TA_RIGHT, leading=9),
    ]], colWidths=[70*mm])
    right_col.setStyle(TableStyle([
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("LEFTPADDING",  (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (-1,-1), 0),
        ("GRID", (0,0), (-1,-1), 0, colors.white),
    ]))

    hdr_tbl = Table([[left_col, right_col]],
                    colWidths=[AVAIL - 70*mm, 70*mm])
    hdr_tbl.setStyle(TableStyle([
        ("VALIGN",       (0,0), (-1,-1), "MIDDLE"),
        ("LEFTPADDING",  (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (-1,-1), 0),
        ("GRID",         (0,0), (-1,-1), 0, colors.white),
    ]))
    story.append(hdr_tbl)
    story.append(HRFlowable(width="100%", thickness=2.5,
                             color=C("#0b2135"), spaceAfter=4))

    # Meta bar: Sample | Case | Date | Platform | Reference | Pipeline
    meta_txt = (
        f"Sample: <b>{sample}</b>  |  Case: <b>{case_id}</b>"
        f"  |  Date: <b>{date_str}</b>"
        f"  |  Platform: <b>{platform} (Illumina)</b>"
        f"  |  Reference: <b>rCRS — NC_012920.1</b>"
        f"  |  Pipeline: <b>MITOCLIN v1.0</b>"
    )
    meta_tbl = Table([[P(meta_txt, fs=7, color="#5a7896", align=TA_LEFT, leading=10)]],
                     colWidths=[AVAIL])
    meta_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), C("#f5f8fc")),
        ("BOX",        (0,0), (-1,-1), 0.4, C("#c8dde8")),
        ("TOPPADDING",    (0,0), (-1,-1), 5),
        ("BOTTOMPADDING", (0,0), (-1,-1), 5),
        ("LEFTPADDING",   (0,0), (-1,-1), 8),
        ("RIGHTPADDING",  (0,0), (-1,-1), 8),
    ]))
    story.append(meta_tbl)
    story.append(Spacer(1, 6))

    # ── STATS CARDS ─────────────────────────────────────
    def stat_pair(val, lbl, bg):
        return [
            P(f"<b>{val}</b>", fs=15, color="#ffffff", align=TA_CENTER, leading=17),
            P(lbl, fs=6.5, color="#ffffff", align=TA_CENTER, leading=9),
        ]

    stat_data = [[
        stat_pair(f"{mean_cov:.0f}\u00d7", "Mean Coverage",          "#0d5e8a"),
        stat_pair(str(n_path),              "Pathogenic",              "#8b1a1a"),
        stat_pair(str(n_lpath),             "Likely Pathogenic",       "#7a4200"),
        stat_pair(str(n_unc),               "Uncertain Significance",  "#4a235a"),
        stat_pair(str(n_ben),               "Benign / Polymorphism",   "#1a6e3a"),
        stat_pair(str(len(variants)),        "Total Variants",          "#3d3d3d"),
    ]]
    cw6 = [AVAIL/6]*6
    stat_tbl = Table(stat_data, colWidths=cw6)
    stat_tbl.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(0,0),C("#0d5e8a")),
        ("BACKGROUND",(1,0),(1,0),C("#8b1a1a")),
        ("BACKGROUND",(2,0),(2,0),C("#7a4200")),
        ("BACKGROUND",(3,0),(3,0),C("#4a235a")),
        ("BACKGROUND",(4,0),(4,0),C("#1a6e3a")),
        ("BACKGROUND",(5,0),(5,0),C("#3d3d3d")),
        ("TOPPADDING",(0,0),(-1,-1),8),("BOTTOMPADDING",(0,0),(-1,-1),8),
        ("LEFTPADDING",(0,0),(-1,-1),3),("RIGHTPADDING",(0,0),(-1,-1),3),
        ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
        ("GRID",(0,0),(-1,-1),0,colors.white),
    ]))
    story.append(stat_tbl)
    story.append(Spacer(1, 7))

    # ── CLINICAL INTERPRETATION ─────────────────────────
    story.append(P("Clinical Interpretation", fs=11, bold=True,
                   color="#0d5e8a", space_before=4, space_after=4))

    interp_raw   = clinical_interp(variants, haplogroup)
    interp_plain = re.sub(r'</?em>','',interp_raw)

    # Haplogroup badge on its own line, then interpretation below — no overlap
    interp_tbl = Table(
        [[P(f"Haplogroup: {haplogroup}", fs=9, bold=True,
            color="#ffffff", bg="#0d5e8a",
            left_i=4, right_i=4, space_after=0)]],
        colWidths=[60*mm]
    )
    interp_tbl.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(0,0),C("#0d5e8a")),
        ("TOPPADDING",(0,0),(-1,-1),5),("BOTTOMPADDING",(0,0),(-1,-1),5),
        ("LEFTPADDING",(0,0),(-1,-1),8),("RIGHTPADDING",(0,0),(-1,-1),8),
        ("ROUNDEDCORNERS",[10]),
    ]))

    # Wrap both inside a background box using a 1-col table
    interp_block = Table(
        [[interp_tbl],
         [P(interp_plain, fs=8.5, leading=12.5, left_i=4, right_i=4,
            color="#1e2733", space_after=0)]],
        colWidths=[AVAIL]
    )
    interp_block.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,-1),C("#edf5fb")),
        ("LINEBEFORE",(0,0),(0,-1),4,C("#1a8fc1")),
        ("TOPPADDING",(0,0),(-1,-1),6),("BOTTOMPADDING",(0,0),(-1,-1),6),
        ("LEFTPADDING",(0,0),(-1,-1),8),("RIGHTPADDING",(0,0),(-1,-1),8),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
    ]))
    story.append(interp_block)
    story.append(Spacer(1, 7))

    # ── VARIANT ANALYSIS PLOTS (3 equal columns) ─────────
    story.append(P("Variant Analysis", fs=11, bold=True,
                   color="#0d5e8a", space_before=2, space_after=4))

    # Each plot: width = AVAIL/3, height proportional
    pw_mm = AVAIL/mm / 3          # width in mm per plot
    ph_mm = pw_mm * 0.78           # height

    plots_row = [[
        img_from_b64(img_hist,   pw_mm, ph_mm),
        img_from_b64(img_circle, pw_mm, ph_mm),
        img_from_b64(img_pie,    pw_mm, ph_mm),
    ]]
    plots_tbl = Table(plots_row, colWidths=[AVAIL/3]*3)
    plots_tbl.setStyle(TableStyle([
        ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
        ("ALIGN", (0,0),(-1,-1),"CENTER"),
        ("GRID",  (0,0),(-1,-1),0,colors.white),
        ("TOPPADDING",(0,0),(-1,-1),0),
        ("BOTTOMPADDING",(0,0),(-1,-1),0),
        ("LEFTPADDING",(0,0),(-1,-1),2),
        ("RIGHTPADDING",(0,0),(-1,-1),2),
    ]))
    story.append(plots_tbl)
    story.append(Spacer(1, 6))

    # ── COVERAGE PLOT ───────────────────────────────────
    story.append(P("Coverage Plot", fs=11, bold=True,
                   color="#0d5e8a", space_before=2, space_after=4))
    cov_h_mm = AVAIL/mm * 0.285
    story.append(img_from_b64(img_cov, AVAIL/mm, cov_h_mm))
    story.append(Spacer(1, 8))

    # ── VARIANT TABLE ───────────────────────────────────
    story.append(P(f"Variant Summary ({len(variants)} variants)",
                   fs=11, bold=True, color="#0d5e8a",
                   space_before=2, space_after=4))

    # ── Cell paragraph helpers (must be defined before QC table and variant table) ──
    def cp(txt, fs=6.5, bold=False, italic=False, clr="#1e2733"):
        fn = ("Helvetica-Bold" if bold else
              "Helvetica-Oblique" if italic else "Helvetica")
        st = ParagraphStyle("_", parent=getSampleStyleSheet()["Normal"],
                            fontName=fn, fontSize=fs, leading=8.5,
                            textColor=C(clr) if isinstance(clr, str) else clr,
                            wordWrap="CJK")
        return Paragraph(txt, st)

    def hp(txt):  # white bold header cell
        return cp(txt, fs=7, bold=True, clr="#ffffff")

    # ── QC SUMMARY TABLE (mirrors HTML QC section) ───────
    story.append(P("Quality Control Summary", fs=10, bold=True,
                   color="#0d5e8a", space_before=4, space_after=4))

    qc_hdr = [hp("Metric"), hp("Value"), hp("Threshold"), hp("Status")]
    cov_status_str = (
        "Excellent (≥1000×)" if mean_cov >= 1000 else
        "Acceptable (≥100×)" if mean_cov >= 100 else
        "Low (<100×)"
    )
    contam_val_str = str(contamination)
    try:
        cv = float(str(contamination).strip('"\''))
        contam_status_str = "Clean (<0.02)" if cv < 0.02 else ("Borderline (<0.05)" if cv < 0.05 else "Contaminated (≥0.05)")
    except Exception:
        contam_status_str = "N/A"

    qc_rows = [
        [cp("Mean Depth of Coverage"), cp(f"{mean_cov:.0f}×", bold=True),
         cp("≥1,000× recommended"), cp(cov_status_str)],
        [cp("Contamination (HaploCheck)"), cp(contam_val_str, bold=True),
         cp("<0.02 — clean sample"), cp(contam_status_str)],
        [cp("Heteroplasmy Threshold"), cp("5% AF", bold=True),
         cp("AF ≥0.05 (GATK recommendation)"), cp("Applied")],
        [cp("Variant Caller"), cp("GATK Mutect2 (mtDNA mode)", bold=True),
         cp("PASS filter + mitochondria-mode"), cp("Standard")],
        [cp("Total SNVs Detected"), cp(str(len(variants)), bold=True),
         cp("Filtered: AF≥5%, PASS, SNVs"), cp("Filtered")],
    ]
    qc_tbl_data = [qc_hdr] + qc_rows
    qc_cw = [AVAIL*w for w in [0.30, 0.18, 0.32, 0.20]]
    qc_tbl = Table(qc_tbl_data, colWidths=qc_cw)
    qc_tcmds = [
        ("BACKGROUND",(0,0),(-1,0),C("#0d5e8a")),
        ("TEXTCOLOR",(0,0),(-1,0),colors.white),
        ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),
        ("FONTSIZE",(0,0),(-1,-1),6.5),
        ("LEADING",(0,0),(-1,-1),8.5),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[C("#f0f5fa"),colors.white]),
        ("GRID",(0,0),(-1,-1),0.3,C("#c8dde8")),
        ("TOPPADDING",(0,0),(-1,-1),4),("BOTTOMPADDING",(0,0),(-1,-1),4),
        ("LEFTPADDING",(0,0),(-1,-1),4),("RIGHTPADDING",(0,0),(-1,-1),4),
        ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
    ]
    qc_tbl.setStyle(TableStyle(qc_tcmds))
    story.append(qc_tbl)
    story.append(Spacer(1, 8))

    story.append(P(f"Variant Table ({len(variants)} variants)",
                   fs=11, bold=True, color="#0d5e8a",
                   space_before=2, space_after=4))

    # ── 11 columns matching HTML exactly ─────────────────
    # Position | Gene | Variant | Var.Type | Protein | Het% | Depth | Ref→Alt | Pathogenicity | Disease | Haplogroup
    cw = [AVAIL*w for w in [0.06, 0.07, 0.10, 0.05, 0.09, 0.06, 0.06, 0.06, 0.12, 0.19, 0.08]]

    tdata = [[
        hp("Position"), hp("Gene"), hp("Variant"), hp("Var.Type"),
        hp("Protein"), hp("Het%"), hp("Depth"), hp("Ref→Alt"),
        hp("Pathogenicity"), hp("Disease"), hp("Haplogroup")
    ]]

    # Pathogenicity → background colour map (keys = PATH_LABELS values)
    PATH_BG_MAP = {
        "Pathogenic":             "#c0392b",
        "Likely Pathogenic":      "#e67e22",
        "Uncertain Significance": "#8e44ad",
        "Benign / Polymorphism":  "#27ae60",
        "Not Annotated":          "#7f8c8d",
    }
    tcmds = [
        ("BACKGROUND",(0,0),(-1,0),C("#0d5e8a")),
        ("TEXTCOLOR",(0,0),(-1,0),colors.white),
        ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),
        ("FONTSIZE",(0,0),(-1,-1),6.5),
        ("LEADING",(0,0),(-1,-1),8.5),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[C("#edf5fb"),colors.white]),
        ("GRID",(0,0),(-1,-1),0.3,C("#c8dff0")),
        ("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3),
        ("LEFTPADDING",(0,0),(-1,-1),3),("RIGHTPADDING",(0,0),(-1,-1),3),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
    ]

    for ri, v in enumerate(variants, start=1):
        path  = v.get("_path","NA")
        label = PATH_LABELS.get(path, path)
        het   = float(v.get("Heteroplasmy_pct", 0))
        cov_p = int(v.get("Coverage_at_Position", 0))
        vtype = v.get("_var_type", "SNV")

        # Ref / Alt — same fallback logic as HTML
        ref = v.get("Ref") or "—"
        alt = v.get("Alt") or "—"
        mut = v.get("_mut","")
        if ref in ("?","") or alt in ("?",""):
            m = re.search(r'([ACGT])>([ACGT])', mut.upper())
            ref, alt = (m.group(1), m.group(2)) if m else ("—","—")

        haplo = v.get("_haplo") or v.get("Haplogroup") or "—"

        row = [
            cp(str(v.get("Position","—")), bold=True),
            cp(v.get("_gene","—")),
            cp(mut, italic=True),
            cp(vtype[:3]),
            cp(v.get("_protein","—")),
            cp(f"{het*100:.1f}%"),
            cp(f"{cov_p:,}"),
            cp(f"{ref}→{alt}"),
            cp(f"<b>{label}</b>", clr="#ffffff"),
            cp(v.get("_disease","—")),
            cp(haplo),
        ]
        tdata.append(row)
        bg = PATH_BG_MAP.get(label)
        if bg:
            tcmds += [
                ("BACKGROUND",(8,ri),(8,ri),C(bg)),
                ("TEXTCOLOR",(8,ri),(8,ri),colors.white),
            ]

    vtbl = Table(tdata, colWidths=cw, repeatRows=1)
    vtbl.setStyle(TableStyle(tcmds))
    story.append(vtbl)
    story.append(Spacer(1, 8))

    # ── METHODOLOGY ─────────────────────────────────────
    story.append(P("Methodology", fs=11, bold=True,
                   color="#0d5e8a", space_before=4, space_after=4))
    story.append(P(
        "<b>Read QC &amp; Trimming:</b> FastQC for quality assessment; Trim Galore "
        "(quality ≥30, minimum length 30 bp, adapter auto-detection). "
        "<b>Alignment:</b> BWA-MEM to rCRS (NC_012920.1); GATK MarkDuplicates for PCR duplicate removal. "
        "<b>Variant Calling:</b> GATK Mutect2 in mitochondria-mode with allele fraction ≥5% filter "
        "and FilterMutectCalls (PASS only, SNVs only, split multiallelic sites with BCFtools norm). "
        "<b>Heteroplasmy:</b> Allele fraction = Alt reads / (Ref + Alt reads) from AD FORMAT field. "
        "<b>Haplogroup &amp; Contamination:</b> HaploCheck v1.3 using Phylotree Build 17. "
        "<b>Annotation:</b> MITOMAP database (Foswiki CSV) for disease association and "
        "pathogenicity; MITOMASTER API for functional and amino acid annotation.",
        fs=8, leading=11.5, bg="#f5f8fc", left_i=5, right_i=5,
        color="#444444", space_after=4
    ))
    # References (matching HTML)
    refs_data = [
        "Van Oven M & Kayser M. (2009). Updated comprehensive phylogenetic tree of global human mitochondrial DNA variation. Hum Mutat. 30(2):E386–394.",
        "MITOMAP: A Human Mitochondrial Genome Database. mitomap.org",
        "McKenna A et al. (2010). The Genome Analysis Toolkit: A MapReduce framework. Genome Res. 20:1297–1303.",
        "Weissensteiner H et al. (2021). Haplocheck: contamination detection in mtDNA studies. Genome Res. 31:1723–1732.",
        "Li H & Durbin R. (2009). Fast and accurate short read alignment with BWA. Bioinformatics. 25(14):1754–1760.",
    ]
    for i, ref in enumerate(refs_data, 1):
        story.append(P(f"{i}. {ref}", fs=6.5, color="#5a7896", leading=9,
                       left_i=8, space_after=1))
    story.append(Spacer(1, 8))

    # ── FOOTER (mirrors HTML footer) ────────────────────
    story.append(HRFlowable(width="100%", thickness=2,
                             color=C("#0b2135"), spaceAfter=6))
    footer_tbl = Table([[
        [
            P("<b>Analysis Team</b>", fs=8, color="#0b2135", leading=11),
            P("Dr. Ajay Kumar Mahato  ·  Ms. Aastha Yadav", fs=7.5, color="#1e2733", leading=10),
            P("Centre for DNA Fingerprinting &amp; Diagnostics (CDFD), Hyderabad",
              fs=6.5, color="#5a7896", leading=9),
            P("github.com/lgi/mitoclin", fs=7, color="#1a8fc1", leading=9),
        ],
        [
            P("<b>⚠ Research Use Only</b>", fs=8, color="#0b2135",
              align=TA_RIGHT, leading=11),
            P("This report is generated for research purposes only and is not "
              "intended for direct clinical diagnosis without expert clinical review. "
              "Results should be interpreted in conjunction with clinical phenotype.",
              fs=7, color="#5a7896", align=TA_RIGHT, leading=10),
            P(f"© 2026 CDFD &amp; MITOCLIN Team. All rights reserved.",
              fs=6.5, color="#5a7896", align=TA_RIGHT, leading=9),
        ],
    ]], colWidths=[AVAIL * 0.5, AVAIL * 0.5])
    footer_tbl.setStyle(TableStyle([
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("LEFTPADDING",  (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (-1,-1), 0),
        ("TOPPADDING",   (0,0), (-1,-1), 0),
        ("BOTTOMPADDING",(0,0), (-1,-1), 4),
        ("GRID", (0,0), (-1,-1), 0, colors.white),
    ]))
    story.append(footer_tbl)
    story.append(P(
        f"Generated by MITOCLIN v1.0  |  {date_str}  |  "
        "BWA-MEM + GATK Mutect2 + HaploCheck + MITOMAP  |  rCRS (NC_012920.1)",
        fs=6.5, italic=True, color="#aaaaaa", align=TA_CENTER, space_before=4
    ))

    doc.build(story)

# ══════════════════════════════════════════════════════════════
# REPORT ORCHESTRATOR
# ══════════════════════════════════════════════════════════════

def generate_reports(sample, csv_path, cov_path, case_id,
                     platform, outdir, vcf_path=None):

    os.makedirs(outdir, exist_ok=True)

    # ===============================
    # VCF OR CSV HANDLING
    # ===============================
    if vcf_path:
        app.logger.info(f"Loading VCF: {vcf_path}")
        variants = load_vcf(vcf_path)
        app.logger.info(f"VCF loaded: {len(variants)} variants after AF≥5% + SNV filter")
        if not variants:
            app.logger.warning(
                "No variants passed filters from VCF — check that the file contains "
                "PASS-filtered variants with AF≥5%"
            )

        if variants:
            # Mean coverage from per-variant read depths
            dp_vals  = [v["Coverage_at_Position"] for v in variants if v["Coverage_at_Position"] > 0]
            mean_cov = sum(dp_vals) / len(dp_vals) if dp_vals else 0

            # Build a sparse coverage profile across the full mtDNA genome.
            # We only have depth at variant positions; interpolate a flat
            # baseline so the coverage plot renders the full 16,569 bp axis.
            var_depth_map = {v["Position"]: v["Coverage_at_Position"] for v in variants}
            if mean_cov > 0:
                # Fill all positions with mean; overwrite with actual depths at variants
                positions = list(range(1, MTDNA_LEN + 1))
                depths    = [var_depth_map.get(p, int(mean_cov)) for p in positions]
            else:
                positions = [v["Position"] for v in variants]
                depths    = [v["Coverage_at_Position"] for v in variants]

            # Haplogroup resolved by MITOMASTER inside load_vcf
            haplo_candidates = [v["_haplo"] for v in variants
                                 if v["_haplo"] not in ("Unknown","—","")]
            haplogroup = haplo_candidates[0] if haplo_candidates else "Unknown"
        else:
            positions, depths = [], []
            haplogroup = "Unknown"
            mean_cov   = 0

        contamination = "NA"

    else:
        variants = load_report(csv_path)

        if cov_path:
            positions, depths = load_coverage(cov_path)
        else:
            positions, depths = [], []

        row0 = variants[0] if variants else {}

        haplogroup = row0.get("_haplo","Unknown") or parse_haplogroup(row0.get("Haplogroup","Unknown"))
        contamination = str(row0.get("Contamination","NA")).strip('"').strip("'")
        mean_cov = float(row0.get("Mean_Coverage",0) or 0)

    date_str = datetime.now().strftime("%d %b %Y")

    # ===============================
    # PLOTS (safe)
    # ===============================
    img_hist   = plot_histogram(variants) if variants else ""
    img_circle = plot_circle(variants) if variants else ""
    img_cov    = plot_coverage(positions, depths, variants) if positions else ""
    img_pie    = plot_pie(variants) if variants else ""

    # ===============================
    # OUTPUT PATHS
    # ===============================
    html_path = os.path.join(outdir, f"{sample}_mitoclin_report.html")
    pdf_path  = os.path.join(outdir, f"{sample}_mitoclin_report.pdf")

    # In VCF mode: copy the source VCF into outdir for download
    vcf_out_path = None
    if vcf_path and os.path.isfile(vcf_path):
        import shutil as _shutil
        vcf_out_path = os.path.join(outdir, f"{sample}_annotated.vcf")
        _shutil.copy2(vcf_path, vcf_out_path)

    # ===============================
    # HTML
    # ===============================
    html = render_html_report(
        variants, positions, depths, sample, case_id,
        date_str, platform, haplogroup, contamination,
        mean_cov, img_hist, img_circle, img_cov, img_pie)

    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)

    # ===============================
    # PDF
    # ===============================
    render_pdf_report(
        variants, positions, depths, sample, case_id,
        date_str, platform, haplogroup, contamination,
        mean_cov, img_hist, img_circle, img_cov, img_pie,
        pdf_path)

    return html_path, pdf_path, len(variants), vcf_out_path
    
# ══════════════════════════════════════════════════════════════
# JOB RUNNER
# ══════════════════════════════════════════════════════════════

def _update(jid, **kw):
    with JOBS_LOCK:
        if jid in JOBS:
            JOBS[jid].update(kw)
            JOBS[jid]["updated_at"] = datetime.now().isoformat()

def _run_pipeline(jid, sample, r1, r2, reference,
                  het_thresh, base_qual, map_qual, case_id,
                  platform, email="", seq_mode="PE"):
    _update(jid, status="running", progress=5, message="Pipeline started")
    log_dir  = os.path.join(OUTPUT_DIR, f"output_{sample}", "logs")
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, "pipeline_runtime.log")
    _update(jid, log_file=log_file)

    if not os.path.isfile(PIPELINE_SH):
        _update(jid, status="failed", progress=0,
                message=f"Pipeline script not found: {PIPELINE_SH}")
        return

    cmd = ["bash", PIPELINE_SH, sample, r1, r2, reference]
    prog_map = {
        "FastQC":10,"Trim Galore":20,"BWA MEM":35,"Sort BAM":45,
        "Mark Duplicates":55,"Mutect2":65,"FilterMutectCalls":72,
        "HaploCheck":80,"MITOMASTER":88,"final report":95,
    }
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            line = line.strip()
            if line:
                for kw, pct in prog_map.items():
                    if kw.lower() in line.lower():
                        _update(jid, progress=pct, message=line)
                        break
        proc.wait()

        if proc.returncode == 0:
            # Pipeline writes output_<SAMPLE> in the dir where the script lives
            pipeline_cwd = os.path.dirname(os.path.abspath(PIPELINE_SH))
            search_dirs  = [pipeline_cwd, OUTPUT_DIR, os.getcwd()]
            csv_p = cov_p = rep_dir = None
            for base in search_dirs:
                _csv = os.path.join(base, f"output_{sample}",
                                    "10_final_report_data", f"{sample}_final_report.csv")
                _cov = os.path.join(base, f"output_{sample}",
                                    "06_coverage", f"{sample}_coverage.txt")
                if os.path.isfile(_csv) and os.path.isfile(_cov):
                    csv_p   = _csv
                    cov_p   = _cov
                    rep_dir = os.path.join(base, f"output_{sample}", "reports")
                    break
            if csv_p and cov_p:
                _update(jid, progress=97, message="Generating reports…")
                rep_dir = os.path.join(os.path.dirname(csv_p), "..", "reports", jid[:8])

                # ── Find the filtered VCF produced by the pipeline ──────────
                vcf_p = None
                for base in search_dirs:
                    # Standard pipeline output paths for filtered VCF
                    for vcf_candidate in [
                        os.path.join(base, f"output_{sample}", "07_filtered_vcf",
                                     f"{sample}_filtered_snps.vcf"),
                        os.path.join(base, f"output_{sample}", "07_filtered_vcf",
                                     f"{sample}_filtered.vcf"),
                        os.path.join(base, f"output_{sample}", "07_filtered_vcf",
                                     f"{sample}.vcf"),
                        os.path.join(base, f"output_{sample}", "05_mutect2",
                                     f"{sample}_filtered.vcf"),
                        os.path.join(base, f"output_{sample}", "05_mutect2",
                                     f"{sample}_somatic_filtered.vcf"),
                        os.path.join(base, f"output_{sample}", "05_mutect2",
                                     f"{sample}_somatic.vcf"),
                        os.path.join(base, f"output_{sample}", "07_filtered_vcf",
                                     f"{sample}_norm_filtered.vcf"),
                    ]:
                        if os.path.isfile(vcf_candidate):
                            vcf_p = vcf_candidate
                            break
                    if vcf_p:
                        break
                    # Glob fallback: find any *filtered*.vcf in output dir
                    if not vcf_p:
                        import glob as _glob
                        out_dir = os.path.join(base, f"output_{sample}")
                        patterns = [
                            f"{out_dir}/**/*filtered*.vcf",
                            f"{out_dir}/**/{sample}*.vcf",
                        ]
                        for pat in patterns:
                            found = _glob.glob(pat, recursive=True)
                            if found:
                                vcf_p = found[0]
                                break
                    if vcf_p:
                        break

                html_p, pdf_p, n, _vcf_out = generate_reports(
                    sample, csv_p, cov_p, case_id, platform, rep_dir)
                from datetime import timedelta
                expiry_dt = (datetime.now()+timedelta(days=REPORT_EXPIRY_DAYS)).isoformat()
                html_url = f"{SERVER_BASE}/api/download/{jid}/html"
                pdf_url  = f"{SERVER_BASE}/api/download/{jid}/pdf"
                upd = dict(status="completed", progress=100,
                           message=f"Done — {n} variants",
                           html_report=html_p, pdf_report=pdf_p,
                           expires_at=expiry_dt,
                           completed_at=datetime.now().isoformat())
                if vcf_p:
                    upd["vcf"] = vcf_p
                _update(jid, **upd)
                if email:
                    send_report_email(email, sample, jid, html_url, pdf_url)
            else:
                # Still mark completed so user can manually generate report
                _update(jid, status="completed_no_report", progress=100,
                        message=f"Pipeline done. CSV not auto-found — use Generate Report tab")
        else:
            _update(jid, status="failed", progress=0,
                    message=f"Pipeline exited with code {proc.returncode}")
    except Exception as e:
        _update(jid, status="failed", progress=0, message=str(e))

# ══════════════════════════════════════════════════════════════
# FLASK ROUTES — API
# ══════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════
# GMAIL CONFIGURATION — genomeinformatics.cdfd@gmail.com
# ══════════════════════════════════════════════════════════════
# Steps to enable (one-time setup):
#  1. Go to https://myaccount.google.com/security
#  2. Enable 2-Step Verification (required for App Passwords)
#  3. Go to https://myaccount.google.com/apppasswords
#  4. App name: "MITOCLIN"  → Google generates a 16-char password
#  5. Paste that password as GMAIL_APP_PASSWORD below (no spaces)
#  6. Do NOT use your regular Gmail login password here
# ══════════════════════════════════════════════════════════════

GMAIL_USER = "genomeinformatics.cdfd@gmail.com"

# FIXED: correct env variable usage
GMAIL_APP_PASSWORD = os.environ.get("GMAIL_APP_PASSWORD", "")

# (Optional fallback for testing — uncomment if needed)
# GMAIL_APP_PASSWORD = "abcdefghijklmnop"

SMTP_HOST = "smtp.gmail.com"
SMTP_PORT = 587
SMTP_USER = GMAIL_USER
SMTP_PASS = GMAIL_APP_PASSWORD
SMTP_FROM = f"MITOCLIN Analysis <{GMAIL_USER}>"
SERVER_BASE = os.environ.get("SERVER_BASE", "http://localhost:5001")
REPORT_EXPIRY_DAYS = 7


def email_ready() -> bool:
    return bool(GMAIL_APP_PASSWORD and GMAIL_APP_PASSWORD.strip())


def send_report_email(to_addr, sample, job_id, html_url, pdf_url):
    """Send analysis completion email via Gmail SMTP (7-day expiry)."""
    if not email_ready() or not to_addr:
        print("Email not sent: App Password not configured or no recipient.")
        return False

    try:
        from datetime import timedelta
        exp_date = (datetime.now() + timedelta(days=REPORT_EXPIRY_DAYS)).strftime("%d %b %Y")

        msg = MIMEMultipart("alternative")
        msg["Subject"] = f"[MITOCLIN] Analysis Complete — {sample}"
        msg["From"] = SMTP_FROM
        msg["To"] = to_addr

        html_body = f"""
        <html><body style="font-family:Arial,sans-serif;color:#1e2733;max-width:600px;margin:0 auto">
        <div style="background:#0d2b45;padding:20px;border-radius:8px 8px 0 0;text-align:center">
          <span style="font-size:22px;font-weight:700;color:#fff">MITO<span style="color:#4db3d8">CLIN</span></span>
          <div style="font-size:11px;color:rgba(255,255,255,0.6);margin-top:4px">
            Mitochondrial Genome Analysis Platform
          </div>
        </div>
        <div style="background:#f5f8fc;padding:24px;border:1px solid #d0dde8">
          <h2 style="color:#0d5e8a;margin:0 0 12px">Analysis Complete</h2>
          <p>Your mitochondrial genome analysis for sample <strong>{sample}</strong> has completed.</p>
          <table style="width:100%;border-collapse:collapse;margin:16px 0;font-size:13px">
            <tr style="background:#e8f0f8">
              <td style="padding:8px 12px;font-weight:600">Sample</td>
              <td style="padding:8px 12px">{sample}</td>
            </tr>
            <tr>
              <td style="padding:8px 12px;font-weight:600">Job ID</td>
              <td style="padding:8px 12px;font-family:monospace;font-size:11px">{job_id}</td>
            </tr>
            <tr style="background:#e8f0f8">
              <td style="padding:8px 12px;font-weight:600">Report Expiry</td>
              <td style="padding:8px 12px;color:#c0392b"><strong>{exp_date}</strong> (7 days)</td>
            </tr>
          </table>
          <p>Download your reports below. Links expire on <strong>{exp_date}</strong>.</p>
          <div style="text-align:center;margin:20px 0">
            <a href="{html_url}" style="background:#1a5276;color:#fff;padding:12px 24px;
               border-radius:6px;text-decoration:none;margin:0 8px;font-weight:600">
              Download HTML Report
            </a>
            <a href="{pdf_url}" style="background:#0d5e8a;color:#fff;padding:12px 24px;
               border-radius:6px;text-decoration:none;margin:0 8px;font-weight:600">
              Download PDF Report
            </a>
          </div>
          <p style="font-size:11px;color:#888;border-top:1px solid #dde;padding-top:12px;margin-top:12px">
            This is an automated message from MITOCLIN. Reports are stored for 7 days.<br/>
            For research use only. Not for clinical diagnosis.<br/>
            &copy; 2026 CDFD &amp; MITOCLIN Team. All rights reserved.
          </p>
        </div>
        </html>
        """
        msg.attach(MIMEText(html_body, "html"))

        with smtplib.SMTP(SMTP_HOST, SMTP_PORT, timeout=15) as server:
            server.ehlo()
            server.starttls()
            server.ehlo()
            server.login(SMTP_USER, SMTP_PASS)
            server.sendmail(GMAIL_USER, to_addr, msg.as_string())
        app.logger.info(f"Email sent to {to_addr} for sample {sample}")
        return True
    except smtplib.SMTPAuthenticationError:
        app.logger.error(
            "Gmail auth failed. Ensure 2FA is ON and you used an App Password "
            "(not your regular password). Get one at: "
            "https://myaccount.google.com/apppasswords"
        )
        return False
    except Exception as e:
        app.logger.warning(f"Email send failed: {type(e).__name__}: {e}")
        return False

@app.route("/api/health")
def api_health():
    return jsonify({"status":"ok","time":datetime.now().isoformat(),
                    "pipeline":os.path.isfile(PIPELINE_SH),
                    "email_configured": email_ready(),
                    "email_from": GMAIL_USER})

@app.route("/api/references")
def api_refs():
    return jsonify({"references": list_references()})

@app.route("/api/submit", methods=["POST"])
def api_submit():
    sample    = request.form.get("sample_name","").strip()
    ref_name  = request.form.get("reference","").strip()
    case_id   = request.form.get("case_id","N/A").strip()
    platform  = request.form.get("platform","Illumina").strip()
    het_thr   = request.form.get("het_threshold","5")
    base_qual = request.form.get("base_quality","30")
    map_qual  = request.form.get("map_quality","20")
    seq_mode  = request.form.get("seq_mode","PE")   # PE or SE
    email     = request.form.get("email","").strip()

    if not sample:
        return jsonify({"error":"sample_name required"}), 400
    if not re.match(r'^[A-Za-z0-9_\-]+$', sample):
        return jsonify({"error":"Sample name: letters/numbers/underscore/hyphen only"}), 400

    r1f = request.files.get("r1")
    r2f = request.files.get("r2") if seq_mode == "PE" else None
    if not r1f:
        return jsonify({"error":"R1 FASTQ file required"}), 400
    if seq_mode == "PE" and not r2f:
        return jsonify({"error":"R2 FASTQ file required for paired-end mode"}), 400

    ref_path = os.path.join(REFERENCE_DIR, ref_name) if ref_name else ""
    if ref_name and not os.path.isfile(ref_path):
        return jsonify({"error":f"Reference not found: {ref_name}"}), 404

    jid      = str(uuid.uuid4())
    jup_dir  = os.path.join(UPLOAD_DIR, jid)
    os.makedirs(jup_dir, exist_ok=True)

    from werkzeug.utils import secure_filename
    r1_path = os.path.join(jup_dir, secure_filename(r1f.filename) or "R1.fastq.gz")
    r1f.save(r1_path)
    r2_path = ""
    if r2f:
        r2_path = os.path.join(jup_dir, secure_filename(r2f.filename) or "R2.fastq.gz")
        r2f.save(r2_path)

    from datetime import timedelta
    expiry_dt = (datetime.now() + timedelta(days=REPORT_EXPIRY_DAYS)).isoformat()

    with JOBS_LOCK:
        JOBS[jid] = {
            "job_id":jid,"sample":sample,"status":"queued",
            "progress":0,"message":"Queued",
            "seq_mode":seq_mode,"email":email,
            "created_at":datetime.now().isoformat(),
            "updated_at":datetime.now().isoformat(),
            "expires_at":expiry_dt,
            "html_report":None,"pdf_report":None,
        }

    t = threading.Thread(
        target=_run_pipeline, daemon=True,
        args=(jid, sample, r1_path, r2_path,
              ref_path or "", het_thr, base_qual,
              map_qual, case_id, platform, email, seq_mode))
    t.start()
    return jsonify({"job_id":jid,"sample":sample,"status":"queued"}), 202

@app.route("/api/report-only", methods=["POST"])
def api_report_only():
    """Generate reports from uploaded CSV+COV or VCF file (multipart).
    Accepts:
      - csv_file + cov_file  → CSV pipeline mode
      - vcf_file             → VCF direct analysis mode
    """
    from werkzeug.utils import secure_filename

    sample   = request.form.get("sample","").strip()
    case_id  = request.form.get("case_id","N/A").strip()
    platform = request.form.get("platform","NGS").strip()

    if not sample:
        return jsonify({"error":"sample name required"}), 400

    vcf_file = request.files.get("vcf_file")
    csv_file = request.files.get("csv_file")
    cov_file = request.files.get("cov_file")

    if not vcf_file and not csv_file:
        return jsonify({"error":"Either vcf_file or (csv_file + cov_file) must be provided"}), 400

    tmp_dir = os.path.join(REPORTS_DIR, "tmp_uploads", sample + "_" + str(uuid.uuid4())[:8])
    os.makedirs(tmp_dir, exist_ok=True)

    jid    = str(uuid.uuid4())
    outdir = os.path.join(REPORTS_DIR, f"{sample}_{jid[:8]}")

    try:
        if vcf_file:
            # ── VCF MODE ────────────────────────────────────
            vcf_path = os.path.join(tmp_dir, secure_filename(vcf_file.filename) or f"{sample}.vcf")
            vcf_file.save(vcf_path)
            html_p, pdf_p, n, vcf_out = generate_reports(
                sample, None, None, case_id, platform, outdir,
                vcf_path=vcf_path
            )
            # vcf_out is the copy in outdir; fall back to original upload if copy failed
            final_vcf = vcf_out if (vcf_out and os.path.isfile(vcf_out)) else vcf_path
            with JOBS_LOCK:
                JOBS[jid] = {
                    "job_id": jid, "sample": sample,
                    "status": "completed", "progress": 100,
                    "message": f"Report generated (VCF) — {n} variants",
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat(),
                    "completed_at": datetime.now().isoformat(),
                    "html_report": html_p, "pdf_report": pdf_p,
                    "vcf": final_vcf,
                }
        else:
            # ── CSV + COVERAGE MODE ──────────────────────────
            if not cov_file:
                return jsonify({"error":"cov_file is required when uploading CSV"}), 400
            csv_path = os.path.join(tmp_dir, secure_filename(csv_file.filename) or f"{sample}_final_report.csv")
            cov_path = os.path.join(tmp_dir, secure_filename(cov_file.filename) or f"{sample}_coverage.txt")
            csv_file.save(csv_path)
            cov_file.save(cov_path)
            html_p, pdf_p, n, _vcf_out = generate_reports(
                sample, csv_path, cov_path, case_id, platform, outdir
            )
            with JOBS_LOCK:
                JOBS[jid] = {
                    "job_id": jid, "sample": sample,
                    "status": "completed", "progress": 100,
                    "message": f"Report generated — {n} variants",
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat(),
                    "completed_at": datetime.now().isoformat(),
                    "html_report": html_p, "pdf_report": pdf_p,
                }

        return jsonify({"job_id": jid, "html": html_p, "pdf": pdf_p, "variants": n,
                        "has_vcf": bool(JOBS[jid].get("vcf"))})
    except Exception as e:
        import traceback
        return jsonify({"error": str(e), "detail": traceback.format_exc()}), 500

@app.route("/api/download-file")
def api_download_file():
    """Serve a report file directly by its absolute path (for browser preview)."""
    path = request.args.get("path","")
    if not path or not os.path.isfile(path):
        return jsonify({"error":"file not found"}), 404
    # Security: only serve files inside BASE_DIR
    if not os.path.abspath(path).startswith(BASE_DIR):
        return jsonify({"error":"access denied"}), 403
    fmt  = "pdf" if path.endswith(".pdf") else "html"
    mime = "application/pdf" if fmt=="pdf" else "text/html"
    name = os.path.basename(path)
    return send_file(path, as_attachment=True,
                     download_name=name, mimetype=mime)

@app.route("/api/status/<jid>")
def api_status(jid):
    with JOBS_LOCK:
        j = JOBS.get(jid,{}).copy()
    if not j: return jsonify({"error":"not found"}),404
    # expose vcf presence as boolean (not the full path) for security
    out = {k:v for k,v in j.items() if k not in ("html_report","pdf_report")}
    out["has_vcf"] = bool(j.get("vcf") and os.path.isfile(j.get("vcf","")))
    # keep vcf key for download — frontend checks j.vcf
    if j.get("vcf"):
        out["vcf"] = j["vcf"]
    return jsonify(out)

@app.route("/api/jobs")
def api_jobs():
    with JOBS_LOCK:
        jobs = []
        for j in JOBS.values():
            row = {k:v for k,v in j.items() if k not in ("html_report","pdf_report")}
            # pass vcf flag so frontend can show VCF download button
            if j.get("vcf") and os.path.isfile(j.get("vcf","")):
                row["vcf"] = j["vcf"]
            jobs.append(row)
    return jsonify({"jobs":sorted(jobs, key=lambda x:x.get("created_at",""), reverse=True)})

@app.route("/api/download/<jid>/<fmt>")
def api_download(jid, fmt):
    with JOBS_LOCK:
        j = JOBS.get(jid,{}).copy()
    if not j: return jsonify({"error":"not found"}),404
    if fmt == "vcf":
        path = j.get("vcf")
        mime = "text/plain"
        ext  = "vcf"
    elif fmt == "pdf":
        path = j.get("pdf_report")
        mime = "application/pdf"
        ext  = "pdf"
    else:
        path = j.get("html_report")
        mime = "text/html"
        ext  = "html"
    if not path or not os.path.isfile(path):
        return jsonify({"error":"report not ready"}),404
    return send_file(path, as_attachment=True,
                     download_name=f"{j['sample']}_mitoclin.{ext}",
                     mimetype=mime)

@app.route("/api/logs/<jid>")
def api_logs(jid):
    with JOBS_LOCK:
        j = JOBS.get(jid,{}).copy()
    if not j: return jsonify({"error":"not found"}),404
    lf   = j.get("log_file","")
    tail = request.args.get("tail",80,type=int)
    if not lf or not os.path.isfile(lf):
        return jsonify({"logs":[]})
    with open(lf) as f:
        lines = f.readlines()
    return jsonify({"logs":[l.rstrip() for l in lines[-tail:]]})

@app.route("/download_vcf/<job_id>")
def download_vcf(job_id):
    job = JOBS.get(job_id)
    if not job:
        return jsonify({"error": "Invalid job ID"}), 404
    vcf_path = job.get("vcf")
    if not vcf_path or not os.path.exists(vcf_path):
        return jsonify({"error": "VCF not available"}), 404
    return send_file(vcf_path, as_attachment=True)
    
# ══════════════════════════════════════════════════════════════
# FLASK ROUTE — SERVE FRONTEND (embedded HTML)
# ══════════════════════════════════════════════════════════════

FRONTEND_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>MITOCLIN – Mitochondrial Genome Analysis Platform</title>
<style>
/* ── RESET & FONTS ── */
*{box-sizing:border-box;margin:0;padding:0}
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:ital,wght@0,300;0,400;0,500;0,600;0,700;1,400&family=DM+Mono:wght@400;500&display=swap');

:root{
  --navy:#061d33;
  --navy2:#0a2d4d;
  --blue:#1163a8;
  --blue2:#1a8fc1;
  --blue3:#4db3d8;
  --ice:#d4eaf7;
  --ice2:#eaf4fb;
  --white:#ffffff;
  --text:#1c2e3e;
  --muted:#5a7a96;
  --border:#c2d9ed;
  --success:#1da462;
  --warn:#d4811a;
  --danger:#c0302a;
  --card-shadow:0 4px 24px rgba(6,29,51,0.13);
  --radius:10px;
  --radius-sm:6px;
}

html,body{
  height:100%;font-family:'DM Sans',sans-serif;
  color:var(--text);background:var(--navy);
  overflow:hidden;
}

/* ── ANIMATED BG ── */
.bg-canvas{
  position:fixed;inset:0;z-index:0;
  background:linear-gradient(135deg,#061d33 0%,#0a2d4d 50%,#0d3a60 100%);
  overflow:hidden;
}
.bg-canvas::before{
  content:'';position:absolute;inset:0;
  background-image:
    radial-gradient(ellipse 60% 40% at 20% 30%,rgba(26,143,193,0.12) 0%,transparent 70%),
    radial-gradient(ellipse 40% 60% at 80% 70%,rgba(77,179,216,0.08) 0%,transparent 70%);
}
.dna-strand{
  position:absolute;right:-80px;top:50%;transform:translateY(-50%);
  width:220px;height:90%;opacity:0.04;
}

/* ── LAYOUT ── */
.shell{
  position:relative;z-index:1;
  display:grid;
  grid-template-columns:220px 1fr 260px;
  grid-template-rows:60px 1fr;
  height:100vh;
  gap:0;
}

/* ── TOPBAR ── */
.topbar{
  grid-column:1/-1;
  display:flex;align-items:center;
  padding:0 24px;
  background:rgba(6,18,32,0.85);
  backdrop-filter:blur(12px);
  border-bottom:1px solid rgba(77,179,216,0.15);
  gap:14px;
}
.logo-icon{
  display:flex;align-items:center;gap:10px;
}
.logo-helix{
  width:36px;height:36px;
}
.logo-text{
  font-size:18px;font-weight:700;
  color:var(--white);letter-spacing:0.5px;
}
.logo-text span{color:var(--blue3)}
.logo-sub{
  font-size:11px;color:var(--blue3);
  opacity:0.8;margin-left:4px;
  font-weight:400;
}
.topbar-right{
  margin-left:auto;display:flex;align-items:center;gap:16px;
}
.status-dot{
  width:8px;height:8px;border-radius:50%;
  background:var(--success);
  box-shadow:0 0 6px var(--success);
  animation:pulse-dot 2s infinite;
}
@keyframes pulse-dot{0%,100%{opacity:1}50%{opacity:0.4}}
.status-label{font-size:11px;color:rgba(255,255,255,0.5)}

/* ── SIDEBAR ── */
.sidebar{
  background:rgba(6,18,32,0.75);
  backdrop-filter:blur(10px);
  border-right:1px solid rgba(77,179,216,0.10);
  padding:20px 0;
  overflow-y:auto;
}
.nav-section{
  padding:0 14px;margin-bottom:6px;
  font-size:9px;font-weight:600;letter-spacing:1.5px;
  color:rgba(77,179,216,0.5);text-transform:uppercase;
}
.nav-item{
  display:flex;align-items:center;gap:10px;
  padding:10px 18px;margin:2px 10px;
  border-radius:var(--radius-sm);
  cursor:pointer;
  color:rgba(255,255,255,0.55);
  font-size:13.5px;font-weight:500;
  transition:all 0.18s ease;
  user-select:none;
}
.nav-item:hover{
  background:rgba(77,179,216,0.10);
  color:rgba(255,255,255,0.85);
}
.nav-item.active{
  background:linear-gradient(90deg,rgba(17,99,168,0.5),rgba(26,143,193,0.2));
  color:var(--white);
  border-left:3px solid var(--blue3);
  margin-left:7px;padding-left:15px;
}
.nav-item svg{width:16px;height:16px;flex-shrink:0}
.nav-divider{
  height:1px;background:rgba(77,179,216,0.08);
  margin:10px 18px;
}

/* ── MAIN CONTENT ── */
.main{
  overflow-y:auto;
  padding:20px 22px;
  background:transparent;
}

/* ── PANELS ── */
.panel{display:none}
.panel.active{display:block}

/* ── CARDS ── */
.card{
  background:rgba(255,255,255,0.04);
  border:1px solid rgba(77,179,216,0.15);
  border-radius:var(--radius);
  backdrop-filter:blur(8px);
  margin-bottom:16px;
  overflow:hidden;
}
.card-header{
  padding:14px 18px;
  border-bottom:1px solid rgba(77,179,216,0.10);
  display:flex;align-items:center;gap:10px;
}
.card-header h2{
  font-size:14.5px;font-weight:600;color:var(--white);
}
.card-header p{
  font-size:11px;color:var(--muted);margin-top:2px;
}
.card-body{padding:18px}

/* ── FORM ELEMENTS ── */
.field-label{
  font-size:11.5px;font-weight:600;
  color:rgba(255,255,255,0.7);
  margin-bottom:6px;display:block;
}
.form-input{
  width:100%;padding:10px 13px;
  background:rgba(6,29,51,0.6);
  border:1px solid rgba(77,179,216,0.25);
  border-radius:var(--radius-sm);
  color:var(--white);font-size:13px;
  font-family:'DM Sans',sans-serif;
  outline:none;transition:border-color 0.2s;
}
.form-input:focus{
  border-color:var(--blue3);
  box-shadow:0 0 0 3px rgba(77,179,216,0.12);
}
.form-input::placeholder{color:rgba(255,255,255,0.25)}
.form-select{
  appearance:none;
  background-image:url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='12' height='8'%3E%3Cpath d='M1 1l5 5 5-5' stroke='%234db3d8' stroke-width='1.5' fill='none'/%3E%3C/svg%3E");
  background-repeat:no-repeat;
  background-position:right 12px center;
  padding-right:34px;
  cursor:pointer;
}
.form-grid{
  display:grid;gap:14px;
}
.form-grid-2{grid-template-columns:1fr 1fr}
.form-grid-3{grid-template-columns:1fr 1fr 1fr}

/* ── FILE UPLOAD ── */
.file-zone{
  border:2px dashed rgba(77,179,216,0.3);
  border-radius:var(--radius-sm);
  padding:14px 16px;
  text-align:center;cursor:pointer;
  transition:all 0.2s;
  background:rgba(6,29,51,0.4);
  position:relative;
}
.file-zone:hover,.file-zone.drag{
  border-color:var(--blue3);
  background:rgba(77,179,216,0.06);
}
.file-zone input[type=file]{
  position:absolute;inset:0;opacity:0;cursor:pointer;width:100%;height:100%;
}
.file-icon{
  width:28px;height:28px;margin:0 auto 6px;
  color:var(--blue3);opacity:0.7;
}
.file-zone p{font-size:11.5px;color:var(--muted)}
.file-zone .fname{
  font-size:11px;color:var(--blue3);
  font-family:'DM Mono',monospace;
  margin-top:4px;word-break:break-all;
}

/* ── SECTION DIVIDER ── */
.sec-divider{
  display:flex;align-items:center;gap:10px;
  margin:6px 0 14px;
}
.sec-divider span{
  font-size:11px;font-weight:600;letter-spacing:0.8px;
  color:var(--blue3);text-transform:uppercase;white-space:nowrap;
}
.sec-divider::before,.sec-divider::after{
  content:'';flex:1;height:1px;
  background:rgba(77,179,216,0.2);
}

/* ── BUTTONS ── */
.btn-primary{
  display:flex;align-items:center;justify-content:center;gap:8px;
  background:linear-gradient(135deg,var(--blue),var(--blue2));
  color:#fff;border:none;border-radius:var(--radius-sm);
  padding:12px 28px;font-size:14px;font-weight:600;
  cursor:pointer;transition:all 0.2s;width:100%;
  font-family:'DM Sans',sans-serif;
  box-shadow:0 4px 14px rgba(17,99,168,0.4);
}
.btn-primary:hover:not(:disabled){
  background:linear-gradient(135deg,#1472bd,#22a3d8);
  transform:translateY(-1px);
  box-shadow:0 6px 18px rgba(17,99,168,0.5);
}
.btn-primary:disabled{opacity:0.55;cursor:not-allowed;transform:none}
.btn-secondary{
  display:inline-flex;align-items:center;gap:6px;
  background:rgba(77,179,216,0.12);
  color:var(--blue3);border:1px solid rgba(77,179,216,0.25);
  border-radius:var(--radius-sm);
  padding:8px 16px;font-size:12px;font-weight:500;
  cursor:pointer;transition:all 0.18s;
  font-family:'DM Sans',sans-serif;
}
.btn-secondary:hover{
  background:rgba(77,179,216,0.2);
  border-color:var(--blue3);
}
.btn-sm{padding:5px 12px;font-size:11px}

/* ── STATUS BADGE ── */
.status-row{
  display:flex;align-items:center;gap:8px;
  margin-top:12px;justify-content:center;
}
.status-badge{
  display:inline-flex;align-items:center;gap:5px;
  padding:5px 12px;border-radius:20px;font-size:11.5px;font-weight:600;
}
.sb-ready{background:rgba(29,164,98,0.15);color:#2de08a;border:1px solid rgba(29,164,98,0.3)}
.sb-running{background:rgba(77,179,216,0.15);color:var(--blue3);border:1px solid rgba(77,179,216,0.3)}
.sb-failed{background:rgba(192,48,42,0.15);color:#ef6c66;border:1px solid rgba(192,48,42,0.3)}
.sb-done{background:rgba(29,164,98,0.15);color:#2de08a;border:1px solid rgba(29,164,98,0.3)}

/* ── PROGRESS ── */
.progress-wrap{
  background:rgba(6,29,51,0.6);border-radius:20px;
  height:8px;overflow:hidden;margin:8px 0;
}
.progress-bar{
  height:100%;border-radius:20px;
  background:linear-gradient(90deg,var(--blue),var(--blue3));
  transition:width 0.5s ease;
  box-shadow:0 0 8px rgba(77,179,216,0.5);
}

/* ── JOB TABLE ── */
.job-table{width:100%;border-collapse:collapse}
.job-table th{
  font-size:10px;font-weight:600;letter-spacing:0.8px;
  color:var(--muted);text-transform:uppercase;
  padding:8px 10px;border-bottom:1px solid rgba(77,179,216,0.12);
  text-align:left;
}
.job-table td{
  padding:9px 10px;border-bottom:1px solid rgba(77,179,216,0.06);
  font-size:12.5px;color:rgba(255,255,255,0.8);vertical-align:middle;
}
.job-table tr:hover td{background:rgba(77,179,216,0.04)}

/* ── STAT CARDS (dashboard) ── */
.stat-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-bottom:16px}
.stat-card{
  background:rgba(255,255,255,0.04);
  border:1px solid rgba(77,179,216,0.15);
  border-radius:var(--radius);padding:16px;
}
.stat-card .stat-val{
  font-size:26px;font-weight:700;color:var(--white);
  font-family:'DM Mono',monospace;
}
.stat-card .stat-lbl{
  font-size:11px;color:var(--muted);margin-top:3px;
}
.stat-card .stat-icon{
  float:right;width:38px;height:38px;
  background:rgba(77,179,216,0.10);
  border-radius:8px;display:flex;align-items:center;justify-content:center;
}

/* ── RIGHT PANEL ── */
.right-panel{
  background:rgba(6,18,32,0.65);
  backdrop-filter:blur(10px);
  border-left:1px solid rgba(77,179,216,0.10);
  padding:18px 16px;
  overflow-y:auto;
}
.rp-section{margin-bottom:20px}
.rp-title{
  font-size:11.5px;font-weight:700;color:var(--white);
  margin-bottom:10px;
  display:flex;align-items:center;gap:6px;
}
.rp-title::after{
  content:'';flex:1;height:1px;background:rgba(77,179,216,0.15);
}
.util-link{
  display:flex;align-items:center;gap:10px;
  padding:9px 10px;border-radius:var(--radius-sm);
  color:rgba(255,255,255,0.65);font-size:12.5px;
  cursor:pointer;transition:all 0.15s;margin-bottom:2px;
  text-decoration:none;
}
.util-link:hover{
  background:rgba(77,179,216,0.10);
  color:var(--white);
}
.util-icon{
  width:30px;height:30px;border-radius:6px;
  display:flex;align-items:center;justify-content:center;
  font-size:14px;flex-shrink:0;
}
.ui-blue{background:rgba(17,99,168,0.25)}
.ui-green{background:rgba(29,164,98,0.2)}
.ui-orange{background:rgba(212,129,26,0.25)}

/* ── LOG BOX ── */
.log-box{
  background:rgba(2,10,18,0.8);
  border:1px solid rgba(77,179,216,0.12);
  border-radius:var(--radius-sm);
  padding:10px 12px;max-height:200px;overflow-y:auto;
  font-family:'DM Mono',monospace;font-size:10.5px;
  color:rgba(77,179,216,0.8);line-height:1.6;
}
.log-box .log-err{color:#ef6c66}
.log-box .log-ok{color:#2de08a}

/* ── TOAST ── */
.toast-wrap{
  position:fixed;top:70px;right:20px;z-index:9999;
  display:flex;flex-direction:column;gap:8px;
  pointer-events:none;
}
.toast{
  background:rgba(6,29,51,0.95);
  border:1px solid rgba(77,179,216,0.25);
  backdrop-filter:blur(12px);
  border-radius:var(--radius-sm);
  padding:10px 16px;font-size:12px;
  color:var(--white);
  box-shadow:0 4px 20px rgba(0,0,0,0.4);
  display:flex;align-items:center;gap:8px;
  animation:slide-in 0.25s ease;pointer-events:all;
}
.toast.success{border-color:rgba(29,164,98,0.4)}
.toast.error{border-color:rgba(192,48,42,0.4)}
@keyframes slide-in{from{opacity:0;transform:translateX(20px)}to{opacity:1;transform:none}}

/* ── SCROLLBAR ── */
::-webkit-scrollbar{width:5px;height:5px}
::-webkit-scrollbar-track{background:transparent}
::-webkit-scrollbar-thumb{background:rgba(77,179,216,0.25);border-radius:3px}


</style>
</head>
<body>

<!-- BG -->
<div class="bg-canvas">
  <svg class="dna-strand" viewBox="0 0 100 600" fill="none">
    <path d="M20 0 Q80 30 20 60 Q80 90 20 120 Q80 150 20 180 Q80 210 20 240 Q80 270 20 300 Q80 330 20 360 Q80 390 20 420 Q80 450 20 480 Q80 510 20 540 Q80 570 20 600" stroke="white" stroke-width="1.5"/>
    <path d="M80 0 Q20 30 80 60 Q20 90 80 120 Q20 150 80 180 Q20 210 80 240 Q20 270 80 300 Q20 330 80 360 Q20 390 80 420 Q20 450 80 480 Q20 510 80 540 Q20 570 80 600" stroke="white" stroke-width="1.5"/>
  </svg>
</div>

<!-- TOAST -->
<div class="toast-wrap" id="toastWrap"></div>

<!-- SHELL -->
<div class="shell">

  <!-- TOPBAR -->
  <header class="topbar">
    <div class="logo-icon">
      <svg class="logo-helix" viewBox="0 0 40 40" fill="none" xmlns="http://www.w3.org/2000/svg">
        <path d="M8 4 Q32 10 8 20 Q32 30 8 36" stroke="#4db3d8" stroke-width="2.2" fill="none"/>
        <path d="M32 4 Q8 10 32 20 Q8 30 32 36" stroke="#1a8fc1" stroke-width="2.2" fill="none"/>
        <line x1="11" y1="9" x2="29" y2="11" stroke="rgba(77,179,216,0.4)" stroke-width="1.2"/>
        <line x1="8" y1="20" x2="32" y2="20" stroke="rgba(77,179,216,0.6)" stroke-width="1.4"/>
        <line x1="11" y1="31" x2="29" y2="29" stroke="rgba(77,179,216,0.4)" stroke-width="1.2"/>
      </svg>
      <div>
        <div class="logo-text">MITO<span>CLIN</span></div>
      </div>
    </div>
    <span class="logo-sub">Mitochondrial Genome Analysis Platform</span>
    <div class="topbar-right">
      <div style="display:flex;align-items:center;gap:6px">
        <div class="status-dot" id="serverDot"></div>
        <span class="status-label" id="serverLabel">Connecting…</span>
      </div>
    </div>
  </header>

  <!-- SIDEBAR -->
  <nav class="sidebar">
    <div style="margin-bottom:16px"></div>
    <div class="nav-section">Navigation</div>
    <div class="nav-item active" data-panel="dashboard">
      <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="3" y="3" width="7" height="7" rx="1"/><rect x="14" y="3" width="7" height="7" rx="1"/><rect x="3" y="14" width="7" height="7" rx="1"/><rect x="14" y="14" width="7" height="7" rx="1"/></svg>
      Dashboard
    </div>
    <div class="nav-item" data-panel="upload">
      <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="17 8 12 3 7 8"/><line x1="12" y1="3" x2="12" y2="15"/></svg>
      Upload &amp; Run
    </div>
    <div class="nav-item" data-panel="status">
      <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="10"/><polyline points="12 6 12 12 16 14"/></svg>
      Pipeline Status
    </div>
    <div class="nav-item" data-panel="results">
      <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><line x1="18" y1="20" x2="18" y2="10"/><line x1="12" y1="20" x2="12" y2="4"/><line x1="6" y1="20" x2="6" y2="14"/></svg>
      Results
    </div>
    <div class="nav-divider"></div>
    <div class="nav-section">Tools</div>
    <div class="nav-item" data-panel="report">
      <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/><line x1="16" y1="13" x2="8" y2="13"/><line x1="16" y1="17" x2="8" y2="17"/></svg>
      Generate Report
    </div>
  </nav>

  <!-- MAIN -->
  <main class="main">

    <!-- ═══ DASHBOARD ═══ -->
    <div id="panel-dashboard" class="panel active">
      <div class="card-header" style="padding:0 0 14px;border:none">
        <div>
          <h2 style="font-size:18px;color:#fff">Welcome to MITOCLIN</h2>
          <p style="margin-top:3px">Automated mitochondrial genome variant analysis pipeline</p>
        </div>
      </div>

      <div class="stat-grid">
        <div class="stat-card">
          <div class="stat-icon">
            <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><path d="M17 21v-2a4 4 0 00-4-4H5a4 4 0 00-4 4v2"/><circle cx="9" cy="7" r="4"/><path d="M23 21v-2a4 4 0 00-3-3.87"/><path d="M16 3.13a4 4 0 010 7.75"/></svg>
          </div>
          <div class="stat-val" id="dash-total-jobs">0</div>
          <div class="stat-lbl">Total Jobs</div>
        </div>
        <div class="stat-card">
          <div class="stat-icon">
            <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#2de08a" stroke-width="2"><polyline points="20 6 9 17 4 12"/></svg>
          </div>
          <div class="stat-val" id="dash-completed">0</div>
          <div class="stat-lbl">Completed</div>
        </div>
        <div class="stat-card">
          <div class="stat-icon">
            <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><circle cx="12" cy="12" r="10"/><polyline points="12 6 12 12 16 14"/></svg>
          </div>
          <div class="stat-val" id="dash-running">0</div>
          <div class="stat-lbl">Running</div>
        </div>
      </div>

      <div class="card">
        <div class="card-header">
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg>
          <div><h2>Pipeline Overview</h2><p>Analysis workflow stages</p></div>
        </div>
        <div class="card-body">
          <div style="display:flex;gap:0;align-items:center;flex-wrap:wrap">
            <div id="pipeline-step" style="display:none">
              Loading pipeline info…
            </div>
          </div>
          <div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;row-gap:8px">
            {steps}
          </div>
        </div>
      </div>

      <div class="card">
        <div class="card-header">
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><rect x="3" y="3" width="18" height="18" rx="2"/><path d="M3 9h18M9 21V9"/></svg>
          <div><h2>Recent Jobs</h2></div>
        </div>
        <div class="card-body" style="padding:0">
          <table class="job-table">
            <thead><tr>
              <th>Sample</th><th>Status</th><th>Progress</th>
              <th>Started</th><th>Action</th>
            </tr></thead>
            <tbody id="dash-job-rows">
              <tr><td colspan="5" style="text-align:center;color:var(--muted);padding:20px">
                No jobs yet. Submit your first analysis.
              </td></tr>
            </tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- ═══ UPLOAD & RUN ═══ -->
    <div id="panel-upload" class="panel">
      <div class="card">
        <div class="card-header">
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="17 8 12 3 7 8"/><line x1="12" y1="3" x2="12" y2="15"/></svg>
          <div>
            <h2>New Analysis</h2>
            <p>Illumina mtDNA sequencing — FASTQ upload &amp; automated variant calling</p>
          </div>
        </div>
        <div class="card-body">

          <!-- Sample info row -->
          <div class="form-grid form-grid-2" style="margin-bottom:12px">
            <div>
              <label class="field-label">Sample Name *</label>
              <input id="f-sample" class="form-input" type="text"
                     placeholder="e.g. S105 (letters, numbers, _ or -)"/>
            </div>
            <div>
              <label class="field-label">Case / Patient ID</label>
              <input id="f-caseid" class="form-input" type="text"
                     placeholder="e.g. CASE-001"/>
            </div>
          </div>

          <!-- Email + SE/PE -->
          <div class="form-grid form-grid-2" style="margin-bottom:14px">
            <div>
              <label class="field-label">
                Notification Email
                <span style="font-weight:400;color:var(--muted);font-size:10px">
                  &nbsp;— report link sent on completion (valid 7 days)
                </span>
              </label>
              <input id="f-email" class="form-input" type="email"
                     placeholder="your@email.com (optional)"/>
            </div>
            <div>
              <label class="field-label">Sequencing Mode</label>
              <div style="display:flex;gap:8px;margin-top:2px">
                <label style="display:flex;align-items:center;gap:8px;cursor:pointer;
                       padding:9px 12px;border-radius:var(--radius-sm);
                       border:2px solid var(--blue3);flex:1;
                       background:rgba(77,179,216,0.08)" id="lbl-pe">
                  <input type="radio" name="seq_mode" id="mode-pe" value="PE"
                         onchange="toggleSeqMode()" checked
                         style="accent-color:var(--blue3);width:14px;height:14px"/>
                  <div>
                    <div style="font-size:12px;font-weight:600;color:#fff">Paired-End</div>
                    <div style="font-size:9.5px;color:var(--muted)">R1 + R2 files</div>
                  </div>
                </label>
                <label style="display:flex;align-items:center;gap:8px;cursor:pointer;
                       padding:9px 12px;border-radius:var(--radius-sm);
                       border:2px solid rgba(77,179,216,0.2);flex:1" id="lbl-se">
                  <input type="radio" name="seq_mode" id="mode-se" value="SE"
                         onchange="toggleSeqMode()"
                         style="accent-color:var(--blue3);width:14px;height:14px"/>
                  <div>
                    <div style="font-size:12px;font-weight:600;color:#fff">Single-End</div>
                    <div style="font-size:9.5px;color:var(--muted)">R1 only</div>
                  </div>
                </label>
              </div>
            </div>
          </div>

          <!-- FASTQ upload zones -->
          <div class="sec-divider"><span>Illumina FASTQ Files</span></div>
          <div style="display:flex;gap:12px;margin-bottom:6px">
            <div style="flex:1">
              <label class="field-label">R1 FASTQ *</label>
              <div class="file-zone" id="zone-r1">
                <input type="file" id="f-r1" accept=".fastq,.gz,.fq,.fastq.gz,.fq.gz"/>
                <svg class="file-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/><line x1="12" y1="18" x2="12" y2="12"/><line x1="9" y1="15" x2="15" y2="15"/></svg>
                <p>Click or drag R1 file</p>
                <div class="fname" id="r1-name"></div>
              </div>
            </div>
            <div style="flex:1" id="r2-zone-wrap">
              <label class="field-label">
                R2 FASTQ
                <span id="r2-req-label" style="color:#ef6c66;font-size:10px">&nbsp;*</span>
                <span id="r2-se-label" style="display:none;font-weight:400;color:var(--muted);font-size:10px">&nbsp;(not needed for SE)</span>
              </label>
              <div class="file-zone" id="zone-r2">
                <input type="file" id="f-r2" accept=".fastq,.gz,.fq,.fastq.gz,.fq.gz"/>
                <svg class="file-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/><line x1="12" y1="18" x2="12" y2="12"/><line x1="9" y1="15" x2="15" y2="15"/></svg>
                <p id="r2-zone-text">Click or drag R2 file</p>
                <div class="fname" id="r2-name"></div>
              </div>
            </div>
          </div>
          <div style="font-size:10.5px;color:var(--muted);margin-bottom:14px">
            &#9432; Platform: <b style="color:rgba(255,255,255,0.7)">Illumina only</b>.
            Accepted: .fastq.gz&nbsp;/&nbsp;.fq.gz
          </div>

          <!-- Pipeline config -->
          <div class="sec-divider"><span>Pipeline Configuration</span></div>
          <div class="form-grid form-grid-2" style="margin-bottom:14px">
            <div>
              <label class="field-label">Reference Genome</label>
              <select id="f-ref" class="form-input form-select">
                <option value="">Loading references…</option>
              </select>
              <div style="font-size:10px;color:var(--muted);margin-top:3px">
                rCRS (NC_012920.1) recommended for mtDNA analysis
              </div>
            </div>
            <div>
              <label class="field-label">Variant Caller</label>
              <div class="form-input" style="opacity:0.65;cursor:not-allowed;display:flex;align-items:center;gap:6px">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><polyline points="20 6 9 17 4 12"/></svg>
                GATK Mutect2 (mtDNA mode)
              </div>
            </div>
          </div>

          <!-- Advanced params -->
          <div class="sec-divider"><span>Advanced Parameters</span></div>
          <div class="form-grid form-grid-3">
            <div>
              <label class="field-label">Heteroplasmy Threshold (%)</label>
              <input id="f-het" class="form-input" type="number" value="5" min="1" max="50"/>
              <div style="font-size:10px;color:var(--muted);margin-top:3px">Min allele fraction (default 5%)</div>
            </div>
            <div>
              <label class="field-label">Base Quality (Phred)</label>
              <input id="f-bq" class="form-input" type="number" value="30" min="10" max="40"/>
              <div style="font-size:10px;color:var(--muted);margin-top:3px">Trim quality cutoff (default 30)</div>
            </div>
            <div>
              <label class="field-label">Min Read Length (bp)</label>
              <input id="f-mq" class="form-input" type="number" value="30" min="15" max="100"/>
              <div style="font-size:10px;color:var(--muted);margin-top:3px">After trimming (default 30 bp)</div>
            </div>
          </div>

          <div style="margin-top:20px">
            <button class="btn-primary" id="btn-submit" onclick="submitJob()">
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polygon points="5 3 19 12 5 21 5 3"/></svg>
              Run mtDNA Analysis
            </button>
            <div class="status-row">
              <span class="status-badge sb-ready" id="submit-status">
                <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="3"><polyline points="20 6 9 17 4 12"/></svg>
                Ready to run
              </span>
            </div>
          </div>
        </div>
      </div>
    </div>

    <!-- ═══ PIPELINE STATUS ═══ -->
    <div id="panel-status" class="panel">
      <div class="card">
        <div class="card-header">
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><circle cx="12" cy="12" r="10"/><polyline points="12 6 12 12 16 14"/></svg>
          <div><h2>Active Jobs</h2><p>Real-time pipeline progress</p></div>
          <button class="btn-secondary btn-sm" onclick="refreshJobs()" style="margin-left:auto">
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="23 4 23 10 17 10"/><path d="M20.49 15a9 9 0 11-2.12-9.36L23 10"/></svg>
            Refresh
          </button>
        </div>
        <div class="card-body" id="jobs-container">
          <div style="text-align:center;color:var(--muted);padding:30px">
            No jobs submitted yet
          </div>
        </div>
      </div>
    </div>

    <!-- ═══ RESULTS ═══ -->
    <div id="panel-results" class="panel">
      <div class="card">
        <div class="card-header">
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><line x1="18" y1="20" x2="18" y2="10"/><line x1="12" y1="20" x2="12" y2="4"/><line x1="6" y1="20" x2="6" y2="14"/></svg>
          <div><h2>Completed Results</h2><p>Download HTML and PDF reports</p></div>
        </div>
        <div class="card-body" id="results-container">
          <div style="text-align:center;color:var(--muted);padding:30px">
            No completed jobs yet
          </div>
        </div>
      </div>
    </div>

    <!-- ═══ GENERATE REPORT ═══ -->
    <div id="panel-report" class="panel">
      <div class="card">
        <div class="card-header">
          <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/></svg>
          <div>
            <h2>Generate Report from Existing Data</h2>
            <p>Upload pipeline output files to generate HTML + PDF reports instantly</p>
          </div>
        </div>
        <div class="card-body">
          <!-- Mode toggle: CSV+COV  vs  VCF -->
          <div style="display:flex;gap:8px;margin-bottom:16px">
            <button id="btn-mode-csv" onclick="setReportMode('csv')"
              style="flex:1;padding:9px 12px;border-radius:var(--radius-sm);
                     border:2px solid var(--blue3);background:rgba(77,179,216,0.12);
                     color:#fff;font-size:12px;font-weight:600;cursor:pointer">
              CSV + Coverage
            </button>
            <button id="btn-mode-vcf" onclick="setReportMode('vcf')"
              style="flex:1;padding:9px 12px;border-radius:var(--radius-sm);
                     border:2px solid rgba(77,179,216,0.2);background:transparent;
                     color:rgba(255,255,255,0.6);font-size:12px;font-weight:600;cursor:pointer">
              VCF File
            </button>
          </div>

          <!-- Info note (CSV mode) -->
          <div id="rg-info-csv" style="background:rgba(77,179,216,0.08);border:1px solid rgba(77,179,216,0.2);
               border-radius:var(--radius-sm);padding:10px 14px;margin-bottom:16px;
               font-size:12px;color:rgba(255,255,255,0.7);display:flex;gap:10px;align-items:flex-start">
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2" style="flex-shrink:0;margin-top:1px"><circle cx="12" cy="12" r="10"/><line x1="12" y1="8" x2="12" y2="12"/><line x1="12" y1="16" x2="12.01" y2="16"/></svg>
            <span>Upload <b>S105_final_report.csv</b> from
            <code style="background:rgba(0,0,0,0.3);padding:1px 5px;border-radius:3px;font-size:11px">output_S105/10_final_report_data/</code>
            and <b>S105_coverage.txt</b> from
            <code style="background:rgba(0,0,0,0.3);padding:1px 5px;border-radius:3px;font-size:11px">output_S105/06_coverage/</code></span>
          </div>

          <!-- Info note (VCF mode) — hidden until VCF button clicked -->
          <div id="rg-info-vcf" style="display:none;background:rgba(77,179,216,0.08);border:1px solid rgba(77,179,216,0.2);
               border-radius:var(--radius-sm);padding:10px 14px;margin-bottom:16px;
               font-size:12px;color:rgba(255,255,255,0.7);gap:10px;align-items:flex-start">
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2" style="flex-shrink:0;margin-top:1px"><circle cx="12" cy="12" r="10"/><line x1="12" y1="8" x2="12" y2="12"/><line x1="12" y1="16" x2="12.01" y2="16"/></svg>
            <span>Upload your <b>filtered VCF file</b>
            (e.g. <code style="background:rgba(0,0,0,0.3);padding:1px 5px;border-radius:3px;font-size:11px">S105_AF5_snps.vcf</code>
            from <code style="background:rgba(0,0,0,0.3);padding:1px 5px;border-radius:3px;font-size:11px">output_S105/04_variant_calling/</code>).
            <b>Full MITOMAP + MITOMASTER annotation runs automatically</b> —
            same as the FASTQ pipeline output.</span>
          </div>

          <div class="form-grid form-grid-2" style="margin-bottom:14px">
            <div>
              <label class="field-label">Sample Name *</label>
              <input id="rg-sample" class="form-input" type="text" placeholder="e.g. S105"/>
            </div>
            <div>
              <label class="field-label">Case ID</label>
              <input id="rg-caseid" class="form-input" type="text" placeholder="e.g. CASE-001"/>
            </div>
          </div>

          <!-- CSV + COV upload zones -->
          <div id="rg-csv-section">
            <div class="sec-divider"><span>Upload Output Files</span></div>
            <div class="form-grid form-grid-2" style="margin-bottom:14px">
              <div>
                <label class="field-label">Final Report CSV *
                  <span style="font-weight:400;color:var(--muted);font-size:10px">
                    (S105_final_report.csv)
                  </span>
                </label>
                <div class="file-zone" id="zone-csv">
                  <input type="file" id="rg-csv-file" accept=".csv,.tsv,.txt"/>
                  <svg class="file-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/><line x1="16" y1="13" x2="8" y2="13"/><line x1="16" y1="17" x2="8" y2="17"/><polyline points="10 9 9 9 8 9"/></svg>
                  <p>Click to upload CSV report</p>
                  <div class="fname" id="csv-fname"></div>
                </div>
              </div>
              <div>
                <label class="field-label">Coverage File *
                  <span style="font-weight:400;color:var(--muted);font-size:10px">
                    (S105_coverage.txt)
                  </span>
                </label>
                <div class="file-zone" id="zone-cov">
                  <input type="file" id="rg-cov-file" accept=".txt,.tsv,.bed"/>
                  <svg class="file-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg>
                  <p>Click to upload coverage file</p>
                  <div class="fname" id="cov-fname"></div>
                </div>
              </div>
            </div>
          </div>

          <!-- VCF upload zone -->
          <div id="rg-vcf-section" style="display:none">
            <div class="sec-divider"><span>Upload VCF File</span></div>
            <div style="margin-bottom:14px">
              <label class="field-label">VCF File *
                <span style="font-weight:400;color:var(--muted);font-size:10px">
                  (GATK Mutect2 or compatible .vcf)
                </span>
              </label>
              <div class="file-zone" id="zone-vcf">
                <input type="file" id="rg-vcf-file" accept=".vcf,.vcf.gz"/>
                <svg class="file-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/><polyline points="8 13 12 17 16 13"/><line x1="12" y1="17" x2="12" y2="9"/></svg>
                <p>Click to upload VCF file</p>
                <div class="fname" id="vcf-fname"></div>
              </div>
            </div>
          </div>

          <button class="btn-primary" id="btn-rg" onclick="generateReport()" style="max-width:340px">
            <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/><line x1="16" y1="13" x2="8" y2="13"/></svg>
            Generate HTML &amp; PDF Reports
          </button>
          <div id="rg-status" style="margin-top:12px"></div>
        </div>
      </div>
    </div>

  </main>

  <!-- RIGHT PANEL -->
  <aside class="right-panel">
    <div class="rp-section">
      <div class="rp-title">Platform Utilities</div>
      <a class="util-link" href="#" onclick="showPanel('upload');return false">
        <span class="util-icon ui-blue">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="17 8 12 3 7 8"/><line x1="12" y1="3" x2="12" y2="15"/></svg>
        </span>
        New Analysis
      </a>
      <a class="util-link" href="#" onclick="showPanel('report');return false">
        <span class="util-icon ui-green">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#2de08a" stroke-width="2"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/></svg>
        </span>
        Generate Report
      </a>
      <a class="util-link" href="https://gatk.broadinstitute.org" target="_blank">
        <span class="util-icon ui-orange">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#d4811a" stroke-width="2"><circle cx="12" cy="12" r="10"/><line x1="2" y1="12" x2="22" y2="12"/><path d="M12 2a15.3 15.3 0 010 20"/></svg>
        </span>
        GATK Documentation
      </a>
      <a class="util-link" href="https://www.mitomap.org" target="_blank">
        <span class="util-icon ui-blue">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><path d="M9 19c-5 1.5-5-2.5-7-3m14 6v-3.87a3.37 3.37 0 00-.94-2.61c3.14-.35 6.44-1.54 6.44-7A5.44 5.44 0 0020 4.77 5.07 5.07 0 0019.91 1S18.73.65 16 2.48a13.38 13.38 0 00-7 0C6.27.65 5.09 1 5.09 1A5.07 5.07 0 005 4.77a5.44 5.44 0 00-1.5 3.78c0 5.42 3.3 6.61 6.44 7A3.37 3.37 0 009 18.13V22"/></svg>
        </span>
        MITOMAP Database
      </a>
    </div>

    <div class="rp-section">
      <div class="rp-title">Pipeline Steps</div>
      <div style="display:flex;flex-direction:column;gap:5px" id="rp-steps">
        <div id="rp-step-fastqc"   class="rp-step">
          <div style="display:flex;justify-content:space-between;font-size:11px;color:rgba(255,255,255,0.55);padding:5px 0;border-bottom:1px solid rgba(77,179,216,0.07)">
            <span>01 FastQC</span><span style="color:var(--muted)">QC</span></div>
        </div>
        <div id="rp-step-trim"     class="rp-step">
          <div style="display:flex;justify-content:space-between;font-size:11px;color:rgba(255,255,255,0.55);padding:5px 0;border-bottom:1px solid rgba(77,179,216,0.07)">
            <span>02 Trim Galore</span><span style="color:var(--muted)">Trim</span></div>
        </div>
        <div id="rp-step-bwa"      class="rp-step">
          <div style="display:flex;justify-content:space-between;font-size:11px;color:rgba(255,255,255,0.55);padding:5px 0;border-bottom:1px solid rgba(77,179,216,0.07)">
            <span>03 BWA MEM</span><span style="color:var(--muted)">Align</span></div>
        </div>
        <div id="rp-step-dedup"    class="rp-step">
          <div style="display:flex;justify-content:space-between;font-size:11px;color:rgba(255,255,255,0.55);padding:5px 0;border-bottom:1px solid rgba(77,179,216,0.07)">
            <span>04 MarkDups</span><span style="color:var(--muted)">Dedup</span></div>
        </div>
        <div id="rp-step-mutect2"  class="rp-step">
          <div style="display:flex;justify-content:space-between;font-size:11px;color:rgba(255,255,255,0.55);padding:5px 0;border-bottom:1px solid rgba(77,179,216,0.07)">
            <span>05 Mutect2</span><span style="color:var(--muted)">Call</span></div>
        </div>
        <div id="rp-step-haplo"    class="rp-step">
          <div style="display:flex;justify-content:space-between;font-size:11px;color:rgba(255,255,255,0.55);padding:5px 0;border-bottom:1px solid rgba(77,179,216,0.07)">
            <span>06 HaploCheck</span><span style="color:var(--muted)">Haplo</span></div>
        </div>
        <div id="rp-step-mito"     class="rp-step">
          <div style="display:flex;justify-content:space-between;font-size:11px;color:rgba(255,255,255,0.55);padding:5px 0">
            <span>07 MITOMASTER</span><span style="color:var(--muted)">Annotate</span></div>
        </div>
      </div>
    </div>

    <div class="rp-section">
      <div class="rp-title">About MITOCLIN</div>
      <div style="font-size:11.5px;color:rgba(255,255,255,0.55);line-height:1.7;padding:6px 4px">
        Automated mtDNA variant calling pipeline for Illumina sequencing data.<br/><br/>
        <b style="color:rgba(255,255,255,0.8)">Team:</b><br/>
        Dr. Ajay Kumar Mahato<br/>
        Ms. Aastha Yadav<br/>
        <span style="font-size:10px;color:rgba(255,255,255,0.35)">CDFD, Hyderabad</span>
      </div>
      <div style="margin-top:8px;padding:4px">
        <a href="https://github.com/lgi/mitoclin" target="_blank"
           style="color:var(--blue3);font-size:11px;text-decoration:none">
          &#128279; github.com/lgi/mitoclin
        </a>
      </div>
      <div style="font-size:10px;color:rgba(255,255,255,0.25);margin-top:6px;padding:4px;line-height:1.5">
        &copy; 2026 CDFD &amp; MITOCLIN Team.<br/>
        Reports retained 7 days on server.<br/>
        For research use only.
      </div>
    </div>

  </aside>
</div>

<script>
// ── STATE
const API = '';
let pollTimer = null;

// ── PANEL NAVIGATION
function showPanel(name) {
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));
  const p = document.getElementById('panel-'+name);
  if(p) p.classList.add('active');
  document.querySelectorAll(`.nav-item[data-panel="${name}"]`)
    .forEach(n => n.classList.add('active'));
  if(name === 'status' || name === 'results') refreshJobs();
  if(name === 'dashboard') refreshDashboard();
}

document.querySelectorAll('.nav-item[data-panel]').forEach(el => {
  el.addEventListener('click', () => showPanel(el.dataset.panel));
});

// ── FILE INPUTS
['r1','r2'].forEach(id => {
  const inp = document.getElementById('f-'+id);
  const nm  = document.getElementById(id+'-name');
  inp.addEventListener('change', () => {
    if(inp.files[0]) {
      nm.textContent = inp.files[0].name;
      document.getElementById('zone-'+id).style.borderColor = 'var(--blue3)';
    }
  });
});

// ── HEALTH CHECK
async function checkHealth() {
  try {
    const r = await fetch(API+'/api/health');
    const d = await r.json();
    const ok = d.status === 'ok';
    document.getElementById('serverDot').style.background = ok ? 'var(--success)' : '#ef6c66';
    document.getElementById('serverDot').style.boxShadow  = ok
      ? '0 0 6px var(--success)' : '0 0 6px #ef6c66';
    document.getElementById('serverLabel').textContent = ok ? 'Server online' : 'Server error';
    if(d.email_configured){
      document.getElementById('serverLabel').textContent = 'Server online • Email ready';
    }
  } catch(e) {
    document.getElementById('serverDot').style.background = '#ef6c66';
    document.getElementById('serverLabel').textContent = 'Offline';
  }
}

// ── LOAD REFERENCES
async function loadRefs() {
  try {
    const r = await fetch(API+'/api/references');
    const d = await r.json();
    const sel = document.getElementById('f-ref');
    sel.innerHTML = '';
    if(!d.references || d.references.length===0) {
      sel.innerHTML = '<option value="">No references found</option>';
      return;
    }
    d.references.forEach(ref => {
      const o = document.createElement('option');
      o.value = ref; o.textContent = ref;
      if(ref.toLowerCase().includes('rcrs') || ref.toLowerCase().includes('mt')) o.selected=true;
      sel.appendChild(o);
    });
  } catch(e) {
    document.getElementById('f-ref').innerHTML = '<option value="">Error loading</option>';
  }
}

// ── TOAST
function toast(msg, type='info') {
  const w = document.getElementById('toastWrap');
  const t = document.createElement('div');
  t.className = `toast ${type}`;
  const icon = type==='success'
    ? '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#2de08a" stroke-width="2.5"><polyline points="20 6 9 17 4 12"/></svg>'
    : type==='error'
    ? '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#ef6c66" stroke-width="2.5"><circle cx="12" cy="12" r="10"/><line x1="15" y1="9" x2="9" y2="15"/><line x1="9" y1="9" x2="15" y2="15"/></svg>'
    : '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><circle cx="12" cy="12" r="10"/><line x1="12" y1="8" x2="12" y2="12"/><line x1="12" y1="16" x2="12.01" y2="16"/></svg>';
  t.innerHTML = icon + msg;
  w.appendChild(t);
  setTimeout(()=>{ t.style.opacity='0'; t.style.transition='opacity 0.3s';
    setTimeout(()=>t.remove(), 300); }, 4000);
}

// ── SE/PE TOGGLE
function toggleSeqMode() {
  const isPE = document.getElementById('mode-pe').checked;
  const r2wrap = document.getElementById('r2-zone-wrap');
  const r2req  = document.getElementById('r2-req-label');
  const r2se   = document.getElementById('r2-se-label');
  const lblPe  = document.getElementById('lbl-pe');
  const lblSe  = document.getElementById('lbl-se');
  if(isPE) {
    r2wrap.style.opacity = '1';
    r2wrap.style.pointerEvents = 'auto';
    r2req.style.display = 'inline';
    r2se.style.display  = 'none';
    lblPe.style.border  = '2px solid var(--blue3)';
    lblSe.style.border  = '2px solid rgba(77,179,216,0.2)';
  } else {
    r2wrap.style.opacity = '0.4';
    r2wrap.style.pointerEvents = 'none';
    r2req.style.display = 'none';
    r2se.style.display  = 'inline';
    lblPe.style.border  = '2px solid rgba(77,179,216,0.2)';
    lblSe.style.border  = '2px solid var(--blue3)';
  }
}

// ── SUBMIT JOB
async function submitJob() {
  const sample  = document.getElementById('f-sample').value.trim();
  const r1      = document.getElementById('f-r1').files[0];
  const r2      = document.getElementById('f-r2').files[0];
  const ref     = document.getElementById('f-ref').value;
  const caseid  = document.getElementById('f-caseid').value.trim() || 'N/A';
  const email   = document.getElementById('f-email').value.trim();
  const het     = document.getElementById('f-het').value;
  const bq      = document.getElementById('f-bq').value;
  const mq      = document.getElementById('f-mq').value;
  const isPE    = document.getElementById('mode-pe').checked;
  const seqMode = isPE ? 'PE' : 'SE';

  if(!sample){ toast('Sample name is required','error'); return; }
  if(!/^[A-Za-z0-9_\-]+$/.test(sample)){ toast('Sample name: letters, numbers, _ or - only','error'); return; }
  if(!r1){ toast('R1 FASTQ file is required','error'); return; }
  if(isPE && !r2){ toast('R2 FASTQ file is required for Paired-End mode','error'); return; }

  const btn = document.getElementById('btn-submit');
  btn.disabled = true;
  btn.innerHTML = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" style="animation:spin 1s linear infinite"><polyline points="23 4 23 10 17 10"/><path d="M20.49 15a9 9 0 11-2.12-9.36L23 10"/></svg> Uploading…';

  const sb = document.getElementById('submit-status');
  sb.className = 'status-badge sb-running';
  sb.innerHTML = '⏳ Uploading &amp; submitting…';

  const fd = new FormData();
  fd.append('sample_name',   sample);
  fd.append('r1',            r1);
  if(isPE && r2) fd.append('r2', r2);
  fd.append('reference',     ref);
  fd.append('case_id',       caseid);
  fd.append('platform',      'Illumina');
  fd.append('seq_mode',      seqMode);
  fd.append('email',         email);
  fd.append('het_threshold', het);
  fd.append('base_quality',  bq);
  fd.append('map_quality',   mq);

  try {
    const res  = await fetch(API+'/api/submit', {method:'POST', body:fd});
    const data = await res.json();
    if(res.ok) {
      const emailMsg = email ? ` — confirmation sent to ${email}` : '';
      toast(`Job submitted! ID: ${data.job_id.slice(0,8)}…${emailMsg}`, 'success');
      sb.className = 'status-badge sb-done';
      sb.innerHTML = `✓ Submitted (${seqMode}) — see Pipeline Status`;
      btn.disabled = false;
      btn.innerHTML = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polygon points="5 3 19 12 5 21 5 3"/></svg> Run mtDNA Analysis';
      startPolling();
      setTimeout(()=>showPanel('status'), 1500);
    } else {
      toast(data.error || 'Submission failed','error');
      btn.disabled = false;
      btn.innerHTML = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polygon points="5 3 19 12 5 21 5 3"/></svg> Run mtDNA Analysis';
      sb.className = 'status-badge sb-failed';
      sb.textContent = data.error || 'Submission failed';
    }
  } catch(e) {
    toast('Network error: '+e.message, 'error');
    btn.disabled = false;
    btn.innerHTML = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><polygon points="5 3 19 12 5 21 5 3"/></svg> Run mtDNA Analysis';
    sb.className = 'status-badge sb-failed';
    sb.textContent = 'Network error';
  }
}

// ── STATUS BADGE HTML
function statusBadge(status, progress) {
  if(status==='running') return `<span class="status-badge sb-running">● Running ${progress}%</span>`;
  if(status==='completed') return `<span class="status-badge sb-done">✓ Completed</span>`;
  if(status==='failed') return `<span class="status-badge sb-failed">✗ Failed</span>`;
  if(status==='queued') return `<span class="status-badge" style="background:rgba(255,255,255,0.06);color:rgba(255,255,255,0.5);border:1px solid rgba(255,255,255,0.1)">⏳ Queued</span>`;
  return `<span class="status-badge">${status}</span>`;
}

// ── REFRESH JOBS
async function refreshJobs() {
  try {
    const r = await fetch(API+'/api/jobs');
    const d = await r.json();
    const jobs = d.jobs || [];

    // Dashboard
    document.getElementById('dash-total-jobs').textContent = jobs.length;
    document.getElementById('dash-completed').textContent  = jobs.filter(j=>j.status==='completed').length;
    document.getElementById('dash-running').textContent    = jobs.filter(j=>j.status==='running').length;

    // Dashboard job rows
    const tr = document.getElementById('dash-job-rows');
    if(jobs.length===0) {
      tr.innerHTML='<tr><td colspan="5" style="text-align:center;color:var(--muted);padding:20px">No jobs yet</td></tr>';
    } else {
      tr.innerHTML = jobs.slice(0,6).map(j=>`
        <tr>
          <td><b style="color:#fff">${j.sample}</b></td>
          <td>${statusBadge(j.status,j.progress)}</td>
          <td style="min-width:100px">
            <div class="progress-wrap"><div class="progress-bar" style="width:${j.progress}%"></div></div>
            <span style="font-size:10px;color:var(--muted)">${j.progress}%</span>
          </td>
          <td style="color:var(--muted);font-size:11px">${j.created_at.slice(0,16).replace('T',' ')}</td>
          <td>
            ${j.status==='completed'
              ? `<button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','html')">HTML</button>
                 <button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','pdf')" style="margin-left:4px">PDF</button>
                 ${j.vcf ? `<button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','vcf')" style="margin-left:4px;border-color:rgba(77,179,216,0.35)">VCF</button>` : ''}`
              : `<button class="btn-secondary btn-sm" onclick="viewLogs('${j.job_id}')">Logs</button>`}
          </td>
        </tr>`).join('');
    }

    // Status panel
    const jc = document.getElementById('jobs-container');
    if(jobs.length===0) {
      jc.innerHTML='<div style="text-align:center;color:var(--muted);padding:30px">No jobs submitted yet</div>';
    } else {
      jc.innerHTML = jobs.map(j=>`
        <div style="border:1px solid rgba(77,179,216,0.12);border-radius:8px;padding:14px;margin-bottom:10px;background:rgba(6,29,51,0.3)">
          <div style="display:flex;justify-content:space-between;align-items:flex-start;margin-bottom:8px">
            <div>
              <span style="font-size:14px;font-weight:600;color:#fff">${j.sample}</span>
              <span style="font-size:10px;color:var(--muted);margin-left:8px;font-family:'DM Mono',monospace">${j.job_id.slice(0,12)}…</span>
            </div>
            <div>${statusBadge(j.status,j.progress)}</div>
          </div>
          <div class="progress-wrap" style="margin-bottom:6px">
            <div class="progress-bar" style="width:${j.progress}%"></div>
          </div>
          <div style="display:flex;justify-content:space-between;align-items:center">
            <span style="font-size:11px;color:var(--muted)">${j.message||''}</span>
            <div style="display:flex;gap:6px">
              ${j.status==='completed'
                ? `<button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','html')">
                    <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
                    HTML</button>
                   <button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','pdf')">
                    <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
                    PDF</button>
                   ${j.vcf ? `<button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','vcf')" style="border-color:rgba(77,179,216,0.35)">
                    <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
                    VCF</button>` : ''}`
                : `<button class="btn-secondary btn-sm" onclick="viewLogs('${j.job_id}')">View Logs</button>`}
            </div>
          </div>
        </div>`).join('');
    }

    // Results panel
    const rc = document.getElementById('results-container');
    const done = jobs.filter(j=>j.status==='completed');
    if(done.length===0) {
      rc.innerHTML='<div style="text-align:center;color:var(--muted);padding:30px">No completed jobs yet</div>';
    } else {
      rc.innerHTML = `<table class="job-table">
        <thead><tr><th>Sample</th><th>Job ID</th><th>Completed</th><th>Download</th></tr></thead>
        <tbody>${done.map(j=>`
          <tr>
            <td><b style="color:#fff">${j.sample}</b></td>
            <td style="font-family:'DM Mono',monospace;font-size:11px;color:var(--muted)">${j.job_id.slice(0,16)}…</td>
            <td style="font-size:11px;color:var(--muted)">${(j.completed_at||j.updated_at||'').slice(0,16).replace('T',' ')}</td>
            <td>
              <button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','html')" style="margin-right:4px">
                HTML Report
              </button>
              <button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','pdf')" style="margin-right:4px">
                PDF Report
              </button>
              ${j.vcf ? `<button class="btn-secondary btn-sm" onclick="downloadReport('${j.job_id}','vcf')" style="border-color:rgba(77,179,216,0.35)">
                VCF File
              </button>` : ''}
            </td>
          </tr>`).join('')}
        </tbody></table>`;
    }
  } catch(e) {
    console.error('refreshJobs:', e);
  }
}

function downloadReport(jid, fmt) {
  window.open(API+'/api/download/'+jid+'/'+fmt);
}

async function viewLogs(jid) {
  const r = await fetch(API+'/api/logs/'+jid+'?tail=50');
  const d = await r.json();
  const logs = (d.logs||[]).join('\n') || 'No logs yet.';
  alert(logs);
}

// ── REPORT MODE TOGGLE (CSV+COV vs VCF)
let _reportMode = 'csv';
function setReportMode(mode) {
  _reportMode = mode;
  const isCsv = mode === 'csv';
  document.getElementById('rg-csv-section').style.display = isCsv ? '' : 'none';
  document.getElementById('rg-vcf-section').style.display = isCsv ? 'none' : '';
  document.getElementById('rg-info-csv').style.display    = isCsv ? 'flex' : 'none';
  document.getElementById('rg-info-vcf').style.display    = isCsv ? 'none' : 'flex';
  const bCsv = document.getElementById('btn-mode-csv');
  const bVcf = document.getElementById('btn-mode-vcf');
  if(isCsv){
    bCsv.style.border = '2px solid var(--blue3)';
    bCsv.style.background = 'rgba(77,179,216,0.12)';
    bCsv.style.color = '#fff';
    bVcf.style.border = '2px solid rgba(77,179,216,0.2)';
    bVcf.style.background = 'transparent';
    bVcf.style.color = 'rgba(255,255,255,0.6)';
  } else {
    bVcf.style.border = '2px solid var(--blue3)';
    bVcf.style.background = 'rgba(77,179,216,0.12)';
    bVcf.style.color = '#fff';
    bCsv.style.border = '2px solid rgba(77,179,216,0.2)';
    bCsv.style.background = 'transparent';
    bCsv.style.color = 'rgba(255,255,255,0.6)';
  }
}

// ── FILE INPUT HANDLERS FOR GENERATE REPORT
document.addEventListener('DOMContentLoaded', () => {
  // ── CSV file input
  const csvInp = document.getElementById('rg-csv-file');
  if (csvInp) csvInp.addEventListener('change', () => {
    const file = csvInp.files[0];
    if (!file) return;
    const kb = (file.size/1024).toFixed(1);
    document.getElementById('csv-fname').textContent = file.name + ' (' + kb + ' KB)';
    document.getElementById('zone-csv').style.borderColor = 'var(--success)';
    document.getElementById('zone-csv').style.background  = 'rgba(29,164,98,0.06)';
    // Auto-fill sample name from filename
    const sn = document.getElementById('rg-sample');
    const patterns = [/^(.+?)_final_report/, /^(.+?)_report/, /^(.+?)\.csv$/];
    for (const pat of patterns) {
      const m = file.name.match(pat);
      if (m) { sn.value = m[1]; break; }
    }
  });

  // ── Coverage file input
  const covInp = document.getElementById('rg-cov-file');
  if (covInp) covInp.addEventListener('change', () => {
    const file = covInp.files[0];
    if (!file) return;
    const kb = (file.size/1024).toFixed(1);
    document.getElementById('cov-fname').textContent = file.name + ' (' + kb + ' KB)';
    document.getElementById('zone-cov').style.borderColor = 'var(--success)';
    document.getElementById('zone-cov').style.background  = 'rgba(29,164,98,0.06)';
  });

  // ── VCF file input
  const vcfInp = document.getElementById('rg-vcf-file');
  if (vcfInp) vcfInp.addEventListener('change', () => {
    const file = vcfInp.files[0];
    if (!file) return;
    const kb = (file.size/1024).toFixed(1);
    document.getElementById('vcf-fname').textContent = file.name + ' (' + kb + ' KB)';
    document.getElementById('zone-vcf').style.borderColor = 'var(--blue3)';
    document.getElementById('zone-vcf').style.background  = 'rgba(77,179,216,0.08)';
    // Auto-fill sample name from VCF filename
    // Handles: S105_AF5_snps.vcf, S105_filtered.vcf, S105.vcf
    const sn = document.getElementById('rg-sample');
    const vcfPatterns = [
      /^(.+?)_AF5_snps/i,
      /^(.+?)_AF5/i,
      /^(.+?)_filtered/i,
      /^(.+?)_snps/i,
      /^(.+?)\.vcf/i,
    ];
    for (const pat of vcfPatterns) {
      const m = file.name.match(pat);
      if (m) { sn.value = m[1]; break; }
    }
  });

  // Initialise mode
  setReportMode('csv');
});

// ── GENERATE REPORT (CSV+COV or VCF mode)
async function generateReport() {
  const sample = document.getElementById('rg-sample').value.trim();
  const caseid = document.getElementById('rg-caseid').value.trim() || 'N/A';
  const isCsv  = (_reportMode === 'csv');

  if (!sample) { toast('Sample name is required', 'error'); return; }

  const fd = new FormData();
  fd.append('sample',   sample);
  fd.append('case_id',  caseid);
  fd.append('platform', 'NGS');

  if (isCsv) {
    const csvFile = document.getElementById('rg-csv-file').files[0];
    const covFile = document.getElementById('rg-cov-file').files[0];
    if (!csvFile) { toast('Please upload the final_report.csv file', 'error'); return; }
    if (!covFile) { toast('Please upload the coverage.txt file', 'error'); return; }
    fd.append('csv_file', csvFile);
    fd.append('cov_file', covFile);
  } else {
    const vcfFile = document.getElementById('rg-vcf-file').files[0];
    if (!vcfFile) { toast('Please upload a VCF file', 'error'); return; }
    fd.append('vcf_file', vcfFile);
  }

  const btn = document.getElementById('btn-rg');
  btn.disabled = true;
  const spinSvg = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" style="animation:spin 1s linear infinite"><polyline points="23 4 23 10 17 10"/><path d="M20.49 15a9 9 0 11-2.12-9.36L23 10"/></svg>';
  const genMsg  = isCsv
    ? 'Uploading &amp; generating reports…'
    : 'Uploading VCF → running MITOMAP + MITOMASTER annotation… (may take ~60s)';
  btn.innerHTML = spinSvg + ' Generating…';

  const st = document.getElementById('rg-status');
  st.innerHTML = `<span class="status-badge sb-running">⏳ ${genMsg}</span>`;

  try {
    const r = await fetch(API+'/api/report-only', { method:'POST', body:fd });
    const d = await r.json();
    if (r.ok) {
    const vcfBtn = (d.has_vcf || !isCsv)
        ? `<button class="btn-secondary" onclick="downloadReport('${d.job_id}','vcf')" style="gap:6px;border-color:rgba(77,179,216,0.35)">
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="#4db3d8" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
            ${sample}_annotated.vcf
          </button>`
        : '';
      st.innerHTML = `
        <div style="display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin-bottom:8px">
          <span class="status-badge sb-done">✓ Reports ready — ${d.variants} variants</span>
        </div>
        <div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center">
          <button class="btn-secondary" onclick="downloadReport('${d.job_id}','html')" style="gap:6px">
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
            ${sample}_mitoclin_report.html
          </button>
          <button class="btn-secondary" onclick="downloadReport('${d.job_id}','pdf')" style="gap:6px">
            <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/></svg>
            ${sample}_mitoclin_report.pdf
          </button>
          ${vcfBtn}
          <button class="btn-secondary" onclick="resetReportForm()" style="gap:6px;border-color:rgba(77,179,216,0.15)">
            ↩ New Report
          </button>
        </div>`;
      toast('Reports generated successfully!', 'success');
      refreshJobs();
      ['csv','cov','vcf'].forEach(id => {
        const z = document.getElementById('zone-'+id);
        if (z) { z.style.borderColor=''; z.style.background=''; }
      });
    } else {
      const msg = d.error || 'Report generation failed';
      st.innerHTML = `<span class="status-badge sb-failed">✗ ${msg}</span>`;
      toast(msg, 'error');
    }
  } catch(e) {
    st.innerHTML = `<span class="status-badge sb-failed">✗ ${e.message}</span>`;
    toast('Network error: ' + e.message, 'error');
  } finally {
    btn.disabled = false;
    btn.innerHTML = '<svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M14 2H6a2 2 0 00-2 2v16a2 2 0 002 2h12a2 2 0 002-2V8z"/><polyline points="14 2 14 8 20 8"/><line x1="16" y1="13" x2="8" y2="13"/></svg> Generate HTML &amp; PDF Reports';
  }
}

// ── RESET REPORT FORM for next sample
function resetReportForm() {
  document.getElementById('rg-sample').value = '';
  document.getElementById('rg-caseid').value = '';
  ['rg-csv-file','rg-cov-file','rg-vcf-file'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.value = '';
  });
  ['csv-fname','cov-fname','vcf-fname'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.textContent = '';
  });
  document.getElementById('rg-status').innerHTML = '';
  ['csv','cov','vcf'].forEach(id => {
    const z = document.getElementById('zone-'+id);
    if (z) { z.style.borderColor=''; z.style.background=''; }
  });
  setReportMode('csv');
}

// ── DASHBOARD REFRESH
function refreshDashboard() { refreshJobs(); }

// ── AUTO-POLL running jobs
function startPolling() {
  if(pollTimer) return;
  pollTimer = setInterval(async ()=>{
    const r = await fetch(API+'/api/jobs').then(r=>r.json()).catch(()=>({jobs:[]}));
    const running = (r.jobs||[]).some(j=>j.status==='running'||j.status==='queued');
    refreshJobs();
    if(!running){ clearInterval(pollTimer); pollTimer=null; }
  }, 3500);
}

// ── CSS SPIN KEYFRAME
const style = document.createElement('style');
style.textContent='@keyframes spin{from{transform:rotate(0deg)}to{transform:rotate(360deg)}}';
document.head.appendChild(style);

// ── INIT
checkHealth();
loadRefs();
refreshDashboard();
setInterval(checkHealth, 15000);
setInterval(()=>{ refreshJobs(); }, 8000);
</script>
</body>
</html>
"""

# Add the pipeline steps HTML into the template
PIPELINE_STEPS = "".join([
    f'<div style="display:flex;align-items:center;gap:0">'
    f'<div style="background:rgba(17,99,168,0.2);border:1px solid rgba(77,179,216,0.25);'
    f'border-radius:6px;padding:7px 12px;font-size:11px;color:rgba(255,255,255,0.7);white-space:nowrap">'
    f'{s}</div>'
    f'{"<div style=&quot;width:18px;height:1px;background:rgba(77,179,216,0.3)&quot;></div>" if i<6 else ""}'
    f'</div>'
    for i,s in enumerate(["FastQC","Trim Galore","BWA MEM","MarkDups",
                           "Mutect2","HaploCheck","MITOMASTER"])
])
FRONTEND_HTML = FRONTEND_HTML.replace("{steps}", PIPELINE_STEPS)

@app.route("/")
def index():
    return FRONTEND_HTML, 200, {"Content-Type": "text/html; charset=utf-8"}

# ══════════════════════════════════════════════════════════════
# CLI — REPORT-ONLY MODE
# ══════════════════════════════════════════════════════════════

def cli_report_only(args):
    print(f"[MITOCLIN] Generating reports for sample: {args.sample}")
    html_p, pdf_p, n, vcf_out = generate_reports(
        args.sample, args.csv, args.cov,
        args.case_id, args.platform,
        args.outdir or os.path.join(os.path.dirname(args.csv), "..", "reports")
    )
    print(f"  HTML : {html_p}")
    print(f"  PDF  : {pdf_p}")
    if vcf_out:
        print(f"  VCF  : {vcf_out}")
    print(f"  Variants: {n}")
    print("[MITOCLIN] Done.")

# ══════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════
@app.route("/upload", methods=["POST"])
def upload():
    file = request.files.get("file")
    sample = request.form.get("sample", "Sample1")
    case_id = request.form.get("case_id", "Case1")
    platform = request.form.get("platform", "Illumina")

    job_id = str(uuid.uuid4())
    outdir = os.path.join(REPORTS_DIR, job_id)
    os.makedirs(outdir, exist_ok=True)

    if not file:
        return jsonify({"error": "No file uploaded"}), 400

    filename = file.filename
    filepath = os.path.join(UPLOAD_DIR, f"{job_id}_{filename}")
    file.save(filepath)

    JOBS[job_id] = {
        "status": "running",
        "sample": sample,
        "outdir": outdir
    }

    try:
        # ===============================
        # VCF MODE
        # ===============================
        if filename.lower().endswith(".vcf"):
            html_path, pdf_path, nvar, vcf_out = generate_reports(
                sample=sample,
                csv_path=None,
                cov_path=None,
                case_id=case_id,
                platform=platform,
                outdir=outdir,
                vcf_path=filepath
            )
            final_vcf = vcf_out if (vcf_out and os.path.isfile(vcf_out)) else filepath
            JOBS[job_id].update({
                "status": "done",
                "html": html_path,
                "pdf": pdf_path,
                "vcf": final_vcf,
                "variants": nvar
            })

        else:
            # ===============================
            # CSV PIPELINE MODE
            # ===============================
            csv_path = filepath
            cov_path = None

            html_path, pdf_path, nvar, _vcf_out = generate_reports(
                sample, csv_path, cov_path,
                case_id, platform, outdir
            )

            JOBS[job_id].update({
                "status": "done",
                "html": html_path,
                "pdf": pdf_path,
                "variants": nvar
            })
    
    except Exception as e:
        JOBS[job_id]["status"] = "error"
        JOBS[job_id]["error"] = str(e)
        return jsonify({"error": str(e)}), 500

    return jsonify({"job_id": job_id})
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="MITOCLIN — Mitochondrial Genome Analysis Platform"
    )
    parser.add_argument("--report-only", action="store_true",
                        help="CLI mode: generate reports only (no web server)")
    parser.add_argument("--sample",   help="Sample name (CLI mode)")
    parser.add_argument("--csv",      help="Path to final_report.csv (CLI mode)")
    parser.add_argument("--cov",      help="Path to coverage.txt (CLI mode)")
    parser.add_argument("--case-id",  default="N/A")
    parser.add_argument("--platform", default="NGS")
    parser.add_argument("--outdir",   default="")
    parser.add_argument("--host",     default="0.0.0.0")
    parser.add_argument("--port",     type=int, default=5001)
    parser.add_argument("--debug",    action="store_true")
    args = parser.parse_args()

    if args.report_only:
        if not args.sample or not args.csv or not args.cov:
            print("ERROR: --sample, --csv and --cov are required in --report-only mode")
            sys.exit(1)
        cli_report_only(args)
    else:
        print("╔══════════════════════════════════════════════╗")
        print("║  MITOCLIN — Mitochondrial Genome Analysis    ║")
        print(f"║  Open: http://localhost:{args.port:<19}║")
        print("╚══════════════════════════════════════════════╝")
        print(f"  Pipeline script : {PIPELINE_SH}")
        print(f"  Reference dir   : {REFERENCE_DIR}")
        print(f"  Uploads dir     : {UPLOAD_DIR}")
        print()
        app.run(host=args.host, port=args.port,
                debug=args.debug, threaded=True)
