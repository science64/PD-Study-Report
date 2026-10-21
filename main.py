"""
pd-study-report

Finds Proteome Discoverer files (.pdResult, .msf, .pdStudy) and writes a
self-contained HTML report (optionally also a PDF) with the study design,
raw files, search / validation / quantification settings, identification
numbers and software versions - ready to send to collaborators.

All files are opened read-only, so nothing is changed.

Usage:
    python main.py                          # files in the folder of this script
    python main.py D:/Studies/MyStudy       # files in another folder
    python main.py my.pdResult --pdf        # one file, HTML + PDF
    python main.py --help                   # all options
"""

import argparse
import glob
import html
import os
import re
import shutil
import sqlite3
import struct
import subprocess
import sys
import tempfile
import webbrowser
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

__version__ = "1.0.0"

# ---------------- settings ----------------
MAX_VALUE_LENGTH = 300   # long values (e.g. XML) are shortened in the report
MAX_LIST = 50            # max number of version hits listed for unknown files
LIST_ALL_TABLES = False  # True = also list every table name in each file
COUNT_PEPTIDES_AND_PSMS = True  # False = skip the slow counts (same as --fast)
ANONYMIZE = False        # True = hide user/computer names and folder paths
PATTERNS = ["*.pdResult", "*.pdresult", "*.msf", "*.MSF",
            "*.pdStudy", "*.pdstudy"]
# ------------------------------------------

HIGH = 3  # PD confidence enum: 1 = Low, 2 = Medium, 3 = High
CONFIDENCE_NAMES = {1: "Low", 2: "Medium", 3: "High"}
DATE_FORMATS = ["%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y %H:%M:%S", "%d.%m.%Y %H:%M:%S",
                "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S"]

# (report section, node category/name regex, parameter name regex)
KEY_SETTINGS = [
    ("Database search", r"sequence database search",
     r"protein database|enzyme name|max\. missed cleavage sites|min\. peptide length|"
     r"max\. peptide length|precursor mass tolerance|fragment mass tolerance"),
    ("PSM validation", r"psm validation",
     r"target fdr.*|validation based on|maximum delta cn"),
    ("Peptide and protein level", r"peptide validation|protein validation|"
     r"peptide and protein filter|protein grouping|protein marker",
     r"validation mode|target fdr.*|peptide confidence at least|minimum peptide length|"
     r"minimum number of peptide sequences|apply strict parsimony principle|protein database"),
    ("Quantification", r"quantif",
     r"integration tolerance|integration method|mass analyzer|ms order|activation type|"
     r"reporter abundance based on|apply quan value corrections|co-isolation threshold|"
     r"average reporter s/n threshold|sps mass matches.*|peptides to use|"
     r"reject quan results with missing channels|normalization mode|scaling mode|"
     r"ratio calculation|imputation mode|hypothesis test|precursor abundance based on"),
]

# processing-log lines worth showing in the report
LOG_PATTERNS = [r"^Sequence Database: ", r"\bversion\b", r"^Sent \d+ spectra from \d+ files",
                r"Total search time", r"confident target/decoy", r"^Stored \d+ (decoy )?PSMs",
                r"^Found \d+ (decoy )?(peptides|protein groups)", r"^Start quantitation with method"]

try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:  # e.g. running inside Jupyter
    SCRIPT_DIR = os.getcwd()


# ======================= small helpers =======================

class Html(str):
    """Text that is already HTML and must not be escaped again."""


def esc(value):
    if isinstance(value, Html):
        return value
    return html.escape("" if value is None else str(value))


def short(value, limit=MAX_VALUE_LENGTH):
    if isinstance(value, bytes):
        return f"<binary data, {len(value)} bytes>"
    text = " ".join(str(value).split())
    if len(text) > limit:
        text = text[:limit] + " ..."
    return text


def basename_any(path):
    """File name from a Windows or POSIX path (PD stores Windows paths)."""
    return re.split(r"[\\/]", str(path))[-1]


def show_path(path):
    if not path:
        return ""
    return basename_any(path) if ANONYMIZE else str(path)


def show_name(name):
    if not name:
        return ""
    return "(hidden)" if ANONYMIZE else str(name)


def scrub(text):
    """Remove folder paths from free text when anonymizing."""
    if not ANONYMIZE or not text:
        return text
    return re.sub(r"(?:[A-Za-z]:)?(?:[\\/][^\\/\s]+)+[\\/]([^\\/\s]+)", r"\1", text)


def natural_key(text):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(text or ""))]


def fmt_int(n):
    if n is None:
        return "–"
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


def fmt_compact(n):
    if n is None:
        return "–"
    for limit, suffix in ((1e12, "T"), (1e9, "B"), (1e6, "M"), (1e3, "K")):
        if abs(n) >= limit:
            return f"{n / limit:.1f} {suffix}"
    return f"{n:.0f}"


def fmt_size(n):
    try:
        n = float(re.match(r"\s*([\d.]+)", str(n)).group(1))
    except (AttributeError, TypeError, ValueError):
        return ""
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def fmt_date(text):
    if not text:
        return ""
    raw = str(text).strip()
    core = re.sub(r"(:\d\d)\.\d+", r"\1", raw).rstrip("Z")
    for fmt in DATE_FORMATS:
        try:
            stamp = datetime.strptime(core, fmt).strftime("%Y-%m-%d %H:%M")
            return stamp + (" UTC" if raw.endswith("Z") else "")
        except ValueError:
            pass
    return raw


def fmt_duration(delta):
    minutes = int(delta.total_seconds() // 60)
    days, minutes = divmod(minutes, 1440)
    hours, minutes = divmod(minutes, 60)
    parts = [f"{days} d" if days else "", f"{hours} h" if hours else "", f"{minutes} min"]
    return " ".join(p for p in parts if p)


def ticks_to_datetime(ticks):
    """.NET DateTime ticks (used in WorkflowMessages.Time) to datetime."""
    try:
        return datetime(1, 1, 1) + timedelta(microseconds=int(ticks) // 10)
    except (TypeError, ValueError, OverflowError):
        return None


def version_number(text):
    match = re.search(r"\d+(?:\.\d+){1,3}", str(text or ""))
    return match.group(0) if match else ""


def decode_text(raw):
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        try:
            return raw.decode("utf-16")
        except UnicodeDecodeError:
            return None
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        return None


def parse_xml(text):
    if isinstance(text, bytes):
        text = decode_text(text)
    if not text or not text.lstrip("\ufeff \r\n\t").startswith("<"):
        return None
    text = re.sub(r"^\s*<\?xml[^>]*\?>", "", text.lstrip("\ufeff"))
    try:
        return ET.fromstring(text)
    except (ET.ParseError, ValueError):
        return None


def find_files(folder):
    found = {}
    for pattern in PATTERNS:
        for path in glob.glob(os.path.join(folder, pattern)):
            found[os.path.normcase(path)] = path  # removes duplicates on Windows
    return sorted(found.values())


def collect_inputs(paths):
    files = {}
    for path in paths:
        if os.path.isdir(path):
            for f in find_files(path):
                files[os.path.normcase(os.path.abspath(f))] = f
        elif os.path.isfile(path):
            files[os.path.normcase(os.path.abspath(path))] = path
        else:
            print(f"Not found, skipped: {path}")
    return sorted(files.values(), key=lambda p: natural_key(os.path.basename(p)))


# ======================= SQLite (.pdResult / .msf) =======================

def is_sqlite(path):
    with open(path, "rb") as f:
        return f.read(16) == b"SQLite format 3\x00"


def connect_read_only(path):
    uri = Path(os.path.abspath(path)).as_uri() + "?mode=ro"
    try:
        return sqlite3.connect(uri, uri=True)
    except sqlite3.OperationalError:
        return sqlite3.connect(path)


def get_tables(con):
    rows = con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    ).fetchall()
    return [r[0] for r in rows]


def get_columns(con, table):
    return [r[1] for r in con.execute(f'PRAGMA table_info("{table}")')]


def fetch_dicts(con, sql, params=()):
    cur = con.execute(sql, params)
    names = [d[0] for d in cur.description]
    return [dict(zip(names, row)) for row in cur.fetchall()]


def read_versions(con, tables):
    """PD version that created the file and the latest version that upgraded it."""
    if "SchemaInfo" not in tables:
        return {}
    cols = get_columns(con, "SchemaInfo")
    if "SoftwareVersion" not in cols:
        return {}
    order = '"Date"' if "Date" in cols else "rowid"
    history = fetch_dicts(con, f"SELECT * FROM SchemaInfo ORDER BY {order}")
    versions = [r["SoftwareVersion"] for r in history if r.get("SoftwareVersion")]
    if not versions:
        return {"history": history}
    latest_date = next((r.get("Date") for r in reversed(history)
                        if r.get("SoftwareVersion") == versions[-1]), "")
    return {"created": versions[0], "latest": versions[-1],
            "latest_date": fmt_date(latest_date), "history": history}


def parse_workflow(source):
    """Workflow name, template and all node parameters from workflow XML."""
    root = source if isinstance(source, ET.Element) else parse_xml(source)
    if root is None:
        return None
    wf = root if root.tag == "Workflow" else root.find(".//Workflow")
    if wf is None:
        return None
    nodes = []
    for node in wf.iter("WorkflowNode"):
        params = []
        for p in node.iter("ProcessingNodeParameter"):
            if p.get("IsHidden") == "True":
                continue
            value = p.get("DisplayValue")
            if value is None:
                value = p.text or ""
            params.append({"group": re.sub(r"^\d+\.\s*", "", p.get("Category", "")),
                           "name": p.get("FriendlyName") or p.get("Name", ""),
                           "value": value.strip(),
                           "advanced": p.get("IsAdvanced") == "True"})
        nodes.append({"number": node.get("ProcessingNodeNumber", ""),
                      "parent": node.get("ParentProcessingNodeNumber", ""),
                      "name": node.get("FriendlyName") or node.get("ProcessingNodeName", ""),
                      "category": node.get("Category", ""),
                      "version": f'{node.get("MainVersion", "?")}.{node.get("MinorVersion", "?")}',
                      "params": params})
    return {"name": wf.get("Name", ""), "template": wf.get("TemplateName", ""),
            "description": wf.get("Description", ""), "type": wf.get("WorkflowType", ""),
            "user": wf.get("CreatedFromUser", ""), "nodes": nodes}


def read_workflows(con, tables):
    if "Workflows" not in tables:
        return []
    durations = {}
    if "WorkflowMessages" in tables and "Time" in get_columns(con, "WorkflowMessages"):
        for wid, first, last in con.execute(
                "SELECT WorkflowID, MIN(Time), MAX(Time) FROM WorkflowMessages GROUP BY WorkflowID"):
            durations[wid] = (ticks_to_datetime(first), ticks_to_datetime(last))
    workflows = []
    for row in fetch_dicts(con, "SELECT * FROM Workflows"):
        wf = parse_workflow(row.get("WorkflowXML")) or {"nodes": []}
        start, end = durations.get(row.get("WorkflowID"), (None, None))
        wf.update({"type": row.get("WorkflowType") or wf.get("type", ""),
                   "name": row.get("WorkflowName") or wf.get("name", ""),
                   "description": row.get("WorkflowDescription") or wf.get("description", ""),
                   "level": row.get("Level") or 0,
                   "started": fmt_date(row.get("WorkflowStartDate")),
                   "user": row.get("User") or wf.get("user", ""),
                   "machine": row.get("MachineName", ""),
                   "software": version_number(row.get("SoftwareVersion")),
                   "ran": (f"{start:%Y-%m-%d %H:%M} → {end:%Y-%m-%d %H:%M} ({fmt_duration(end - start)})"
                           if start and end else "")})
        workflows.append(wf)
    workflows.sort(key=lambda w: -int(w["level"]))  # processing first, then consensus
    return workflows


def read_input_files(con, tables):
    if "WorkflowInputFiles" not in tables:
        return []
    files = []
    for r in fetch_dicts(con, "SELECT * FROM WorkflowInputFiles"):
        if r.get("ChildWorkflowID") not in (0, None) or "msf" in str(r.get("FileType", "")).lower():
            continue  # .msf files produced by a processing step, not spectrum files
        files.append({"fraction": r.get("StudyFileID") or "",
                      "name": basename_any(r.get("FileName") or ""),
                      "path": show_path(r.get("FileName")),
                      "acquired": fmt_date(r.get("CreationDate")),
                      "size": r.get("FileSize"),
                      "instrument": r.get("InstrumentDescription") or "",
                      "instrument_software": r.get("SoftwareRevision") or "",
                      "rt_range": r.get("RetentionTimeRange") or ""})
    files.sort(key=lambda f: natural_key(f["fraction"] or f["name"]))
    return files


def read_messages(con, tables):
    if "WorkflowMessages" not in tables:
        return [], []
    cols = get_columns(con, "WorkflowMessages")
    if "Message" not in cols:
        return [], []
    order = "WorkflowID, MessageID" if "MessageID" in cols else "rowid"
    highlights, warnings = {}, {}
    for r in fetch_dicts(con, f"SELECT * FROM WorkflowMessages ORDER BY {order}"):
        text = scrub(" ".join(str(r.get("Message") or "").split()))
        node = re.sub(r"^\(\d+\):", "", str(r.get("ProcessingNodeName") or ""))
        if (r.get("MessageKind") or 0) >= 2:
            warnings[(node, text)] = None
        elif any(re.search(p, text, re.I) for p in LOG_PATTERNS):
            highlights[(node, text)] = None
    return list(highlights), list(warnings)


def confidence_value(value):
    if isinstance(value, bytes):
        return int.from_bytes(value[:4], "little") if len(value) >= 4 else None
    return value


def decode_abundances(blob):
    """PD stores per-channel values as 8-byte doubles followed by a 1-byte 'has value' flag."""
    if not isinstance(blob, bytes) or len(blob) % 9:
        return None
    values = []
    for i in range(0, len(blob), 9):
        number = struct.unpack("<d", blob[i:i + 8])[0]
        values.append(number if blob[i + 8] else None)
    return values


def count_results(con, tables, detailed):
    counts, per_file = {}, {}

    def one(sql):
        try:
            return con.execute(sql).fetchone()
        except sqlite3.Error:
            return None

    def total_and_filtered(table):
        cols = get_columns(con, table)
        filtered = 'SUM("ExcludedBy"=-1)' if "ExcludedBy" in cols else "NULL"
        return one(f'SELECT COUNT(*), {filtered} FROM "{table}"') or (None, None)

    if "TargetProteinGroups" in tables:
        counts["protein_groups"] = total_and_filtered("TargetProteinGroups")

    channel_sums = None
    if "TargetProteins" in tables:
        cols = get_columns(con, "TargetProteins")
        counts["proteins"] = (one('SELECT COUNT(*) FROM "TargetProteins"') or [None])[0]
        if "IsMasterProtein" in cols:
            wanted = [c for c in ("ProteinFDRConfidence", "Contaminant", "ExcludedBy", "Abundances")
                      if c in cols]
            select = ", ".join(f'"{c}"' for c in wanted) or "1"
            by_conf, master, contaminants, sums = {}, 0, 0, None
            try:
                rows = con.execute(f'SELECT {select} FROM "TargetProteins" WHERE "IsMasterProtein"=0')
                for row in rows:
                    row = dict(zip(wanted, row))
                    master += 1
                    conf = confidence_value(row.get("ProteinFDRConfidence"))
                    by_conf[conf] = by_conf.get(conf, 0) + 1
                    if row.get("Contaminant"):
                        contaminants += 1
                        continue
                    if row.get("ExcludedBy", -1) != -1:
                        continue
                    values = decode_abundances(row.get("Abundances"))
                    if values:
                        if sums is None:
                            sums = [0.0] * len(values)
                        if len(values) == len(sums):
                            sums = [s + (v or 0.0) for s, v in zip(sums, values)]
            except sqlite3.Error:
                pass
            counts.update({"master": master, "master_by_conf": by_conf,
                           "master_contaminants": contaminants})
            channel_sums = sums if sums and len(sums) > 1 and any(sums) else None

    if "MSnSpectrumInfo" in tables:
        counts["spectra"] = (one('SELECT COUNT(*) FROM "MSnSpectrumInfo"') or [None])[0]

    if detailed and "TargetPeptideGroups" in tables:
        print("  counting peptide groups ...", flush=True)
        cols = get_columns(con, "TargetPeptideGroups")
        high = f'SUM("Confidence"={HIGH})' if "Confidence" in cols else "NULL"
        filtered = 'SUM("ExcludedBy"=-1)' if "ExcludedBy" in cols else "NULL"
        counts["peptides"] = one(f'SELECT COUNT(*), {high}, {filtered} FROM "TargetPeptideGroups"')

    if detailed and "TargetPsms" in tables:
        print("  counting PSMs per raw file (can take a minute for large files) ...", flush=True)
        cols = get_columns(con, "TargetPsms")
        conf_col = next((c for c in ("MatchConfidence", "Confidence") if c in cols), None)
        high = f'SUM("{conf_col}"={HIGH})' if conf_col else "NULL"
        filtered = 'SUM("ExcludedBy"=-1)' if "ExcludedBy" in cols else "NULL"
        group = '"SpectrumFileName"' if "SpectrumFileName" in cols else "''"
        totals = [0, 0, 0]
        try:
            for name, n, n_high, n_filtered in con.execute(
                    f'SELECT {group}, COUNT(*), {high}, {filtered} FROM "TargetPsms" GROUP BY 1'):
                per_file.setdefault(str(name).lower(), {}).update(
                    {"psms": n, "psms_high": n_high, "psms_filtered": n_filtered})
                totals = [t + (v or 0) for t, v in zip(totals, (n, n_high, n_filtered))]
            counts["psms"] = (totals[0], totals[1] if conf_col else None,
                              totals[2] if "ExcludedBy" in cols else None)
        except sqlite3.Error:
            pass

    if detailed and "MSnSpectrumInfo" in tables and "SpectrumFileName" in get_columns(con, "MSnSpectrumInfo"):
        try:
            for name, n in con.execute(
                    'SELECT "SpectrumFileName", COUNT(*) FROM "MSnSpectrumInfo" GROUP BY 1'):
                per_file.setdefault(str(name).lower(), {})["spectra"] = n
        except sqlite3.Error:
            pass
    return counts, per_file, channel_sums


def inspect_sqlite(path, detailed):
    con = connect_read_only(path)
    try:
        tables = get_tables(con)
        info = {"kind": "sqlite", "tables": tables}
        info["versions"] = read_versions(con, tables)
        info["workflows"] = read_workflows(con, tables)
        info["raw_files"] = read_input_files(con, tables)
        info["log"], info["warnings"] = read_messages(con, tables)
        info["study"] = None
        if "AnalysisDefinition" in tables:
            row = con.execute("SELECT AnalysisDefinitionXML FROM AnalysisDefinition").fetchone()
            root = parse_xml(row[0]) if row else None
            if root is not None and root.find("StudyDefinition") is not None:
                info["study"] = parse_study_definition(root.find("StudyDefinition"), root)
        info["counts"], info["per_file"], info["channel_sums"] = count_results(con, tables, detailed)
        if not info["versions"].get("created"):
            version = next((w["software"] for w in info["workflows"] if w.get("software")), "")
            info["versions"]["created"] = version
        return info
    finally:
        con.close()


# ======================= .pdStudy (XML) =======================

def read_reporter_masses(sd):
    masses = {}
    for qm in sd.iter("QuantitationMethod"):
        method = parse_xml(qm.text)
        if method is None:
            continue
        for part in method.iter("MethodPart"):
            params = {p.get("name"): (p.text or "").strip() for p in part.findall("Parameter")}
            if params.get("TagName"):
                masses[params["TagName"]] = params.get("MonoisotopicMZ", "")
    return masses


def parse_study_definition(sd, extra=None):
    """Quan method(s), file sets and samples from a <StudyDefinition> element."""
    factor_parent = sd.find("Factors")
    factors = {}
    for f in (list(factor_parent) if factor_parent is not None else []):
        factors[f.get("Id") or f.get("Name") or f.tag] = f.get("Name") or f.tag

    methods, method_names, channel_names = [], {}, {}
    for qm in sd.findall("QuanMethods/QuanMethod"):
        channels = sorted(qm.iter("QuanChannel"), key=lambda c: natural_key(c.get("Position")))
        for c in channels:
            channel_names[(qm.get("Id"), c.get("Id"))] = c.get("Name", "")
        method_names[qm.get("Id")] = qm.get("Name", "")
        methods.append({"name": qm.get("Name", ""),
                        "description": (qm.get("Description") or "").replace("\ufffd", ""),
                        "kind": qm.get("MethodKind", ""),
                        "channels": [c.get("Name", "") for c in channels]})

    file_sets = []
    for fs in sd.findall("FileSets/FileSet"):
        files = [{"fraction": f.get("FractionNumber", ""), "name": basename_any(f.get("FileName", "")),
                  "path": show_path(f.get("FileName")), "time": fmt_date(f.get("FileTime")),
                  "size": f.get("FileSize")} for f in fs.iter("File")]
        files.sort(key=lambda f: natural_key(f["fraction"] or f["name"]))
        file_sets.append({"id": fs.get("Id"), "name": fs.get("Name", ""),
                          "label": f'F{fs.get("SampleFileSetNumber")}' if fs.get("SampleFileSetNumber") else "",
                          "method": method_names.get(fs.get("QuanMethodId"), ""), "files": files})
    set_labels = {fs["id"]: fs["label"] or fs["name"] for fs in file_sets}

    sample_info = {}
    if extra is not None:
        for si in extra.iter("SampleInfo"):
            sample_info[si.get("SampleId")] = {
                "label": si.get("SampleNumberDisplayString", ""),
                "variables": [(v.get("Name", ""), v.get("Value", ""), v.get("Kind", ""))
                              for v in si.iter("StudyVariable")]}

    samples, factor_columns = [], {}
    for s in sd.findall("Samples/Sample"):
        q = s.find("QuanChannelInformation")
        channel = channel_names.get((q.get("QuanMethodId"), q.get("QuanChannelId")), "") if q is not None else ""
        values = {}
        fv_parent = s.find("FactorValues")
        for fv in (list(fv_parent) if fv_parent is not None else []):
            key = next((v for k, v in fv.attrib.items() if "factor" in k.lower()), fv.get("Id"))
            value = fv.get("Value") or (fv.text or "").strip()
            if value:
                values[factors.get(key, fv.get("Name") or key or fv.tag)] = value
        extra_info = sample_info.get(s.get("Id"), {})
        for name, value, kind in extra_info.get("variables", []):
            if kind not in ("QuanChannel", "SampleType", "File") and value:
                values[name] = value
        for name in values:
            factor_columns[name] = None
        samples.append({"number": s.get("SampleNumber", ""),
                        "label": extra_info.get("label") or f'S{s.get("SampleNumber", "")}',
                        "name": s.get("Name", ""), "type": s.get("SampleType", ""),
                        "channel": channel, "file_set": set_labels.get(s.get("FileSetId"), ""),
                        "factors": values})
    samples.sort(key=lambda s: natural_key(s["number"]))

    grouping, ratios = [], []
    if extra is not None:
        grouping = list(dict.fromkeys(
            v.get("Name", "") for v in extra.findall(".//StudyVariablesForGrouping/StudyVariable")))
        for r in extra.findall(".//QuanRatios/*"):
            ratios.append(" ".join(f"{k}={v}" for k, v in r.attrib.items()) or r.tag)

    return {"name": sd.get("Name", ""), "description": sd.get("Description", ""),
            "created": fmt_date(sd.get("CreationDate")), "changed": fmt_date(sd.get("LastChangeDate")),
            "directory": show_path(sd.get("Directory")),
            "factors": list(factors.values()) or list(factor_columns),
            "factor_columns": list(factor_columns), "methods": methods,
            "reporter_mz": read_reporter_masses(sd), "file_sets": file_sets,
            "samples": samples, "grouping": grouping, "ratios": ratios}


def parse_analyses(root):
    analyses = []
    for ar in root.findall("AnalysisResults/AnalysisResult"):
        analysis = ar.find("Analysis")
        if analysis is None:
            continue
        steps, workflows = [], []
        for step in analysis.iter("SequenceStep"):
            wi = step.find("WorkflowInfo")
            if wi is None:
                continue
            steps.append({"type": wi.findtext("WorkflowType", ""), "workflow": wi.findtext("WorkflowName", ""),
                          "file": basename_any(step.findtext("FileName", "")),
                          "state": step.findtext("WorkflowExecutionState", "")})
            wf = parse_workflow(wi.find("WorkflowDefinition/Workflow")) \
                if wi.find("WorkflowDefinition/Workflow") is not None else None
            if wf:
                wf["type"] = wf["type"] or wi.findtext("WorkflowType", "")
                workflows.append(wf)
        steps.sort(key=lambda s: s["type"] != "Processing")
        workflows.sort(key=lambda w: w["type"] != "Processing")
        analyses.append({"number": ar.get("ResultNumber", ""),
                         "created": fmt_date(analysis.get("CreationDate")),
                         "description": analysis.get("Description", ""),
                         "steps": steps, "workflows": workflows})
    analyses.sort(key=lambda a: natural_key(a["number"]))
    return analyses


def read_study_root(path):
    with open(path, "rb") as f:
        data = f.read()
    if data[:2] == b"PK":
        with zipfile.ZipFile(path) as z:
            for name in z.namelist():
                root = parse_xml(z.read(name))
                if root is not None and root.tag == "Study":
                    return root
        return None
    return parse_xml(data)


def inspect_study_file(path):
    root = read_study_root(path)
    if root is None or root.find("StudyDefinition") is None:
        return inspect_other_file(path)
    return {"kind": "study", "study": parse_study_definition(root.find("StudyDefinition")),
            "analyses": parse_analyses(root)}


def inspect_other_file(path):
    """Best-effort search for version entries in files of unknown format."""
    with open(path, "rb") as f:
        data = f.read()
    entries = []
    if data[:2] == b"PK":
        with zipfile.ZipFile(path) as z:
            entries = [(name, z.read(name)) for name in z.namelist()]
    else:
        entries.append((os.path.basename(path), data))
    hits = {}
    for name, raw in entries:
        text = decode_text(raw)
        if text is None:
            continue
        found = re.findall(r'[\w:.-]*version[\w:.-]*\s*=\s*"[^"]*"', text, re.I)
        found += [f"{tag} = {val}" for tag, val in
                  re.findall(r"<([\w:.-]*version[\w:.-]*)>([^<]*)</", text, re.I)]
        found = list(dict.fromkeys(found))
        if found:
            hits[name] = found[:MAX_LIST]
    return {"kind": "other", "hits": hits}


# ======================= settings summary & methods text =======================

def search_nodes(nodes):
    return [n for n in nodes if re.search(r"sequence database search", n["category"], re.I)]


def pick(nodes, node_re, param_re):
    for node in nodes:
        if re.search(node_re, f'{node["category"]} {node["name"]}', re.I):
            for p in node["params"]:
                if p["value"] and re.fullmatch(param_re, p["name"], re.I):
                    return p["value"]
    return ""


def collect_modifications(nodes):
    static, dynamic = {}, {}
    for node in search_nodes(nodes):
        for p in node["params"]:
            value, where = p["value"], f'{p["group"]} {p["name"]}'.lower()
            if not re.search(r"/\s*[+-]?\d", value) or "modification" not in where + " terminus":
                continue
            if "protein terminus" in where:
                value = value.replace("(N-Terminus)", "(protein N-terminus)")
                value = value.replace("(C-Terminus)", "(protein C-terminus)")
            (static if "static" in where else dynamic)[value] = None
    return list(static), list(dynamic)


def collect_key_settings(workflows):
    sections = {}
    for title, node_re, param_re in KEY_SETTINGS:
        rows, seen = [], set()
        for wf in workflows:
            for node in wf["nodes"]:
                if not re.search(node_re, f'{node["category"]} {node["name"]}', re.I):
                    continue
                for p in node["params"]:
                    key = (p["name"], p["value"], node["name"])
                    if p["value"] and re.fullmatch(param_re, p["name"], re.I) and key not in seen:
                        seen.add(key)
                        rows.append(key)
        if title == "Database search":
            nodes = [n for wf in workflows for n in wf["nodes"]]
            engines = list(dict.fromkeys(n["name"] for n in search_nodes(nodes)))
            static, dynamic = collect_modifications(nodes)
            if engines:
                rows.insert(0, ("Search engine", ", ".join(engines), ""))
            if static:
                rows.append(("Static modifications", Html("<br>".join(esc(m) for m in static)), ""))
            if dynamic:
                rows.append(("Dynamic modifications", Html("<br>".join(esc(m) for m in dynamic)), ""))
        if rows:
            sections[title] = rows
    return sections


def methods_summary(rep):
    """A draft methods paragraph built only from values stored in the file."""
    nodes = [n for wf in rep.get("workflows", []) for n in wf["nodes"]]
    sentences = []
    raw = rep.get("raw_files", [])
    instruments = sorted({f["instrument"] for f in raw if f.get("instrument")})
    version = rep.get("versions", {}).get("created", "")
    if raw or version:
        text = f"{len(raw)} raw file{'s' if len(raw) != 1 else ''}" if raw else "Raw files"
        if instruments:
            text += f" acquired on {' / '.join(instruments)}"
        text += " were processed with Proteome Discoverer" + (f" {version}" if version else "") + "."
        sentences.append(text)

    engines = list(dict.fromkeys(n["name"] for n in search_nodes(nodes)))
    if engines:
        S = r"sequence database search"
        db, enzyme = pick(nodes, S, r"protein database"), pick(nodes, S, r"enzyme name")
        missed = pick(nodes, S, r"max\. missed cleavage sites")
        min_len, max_len = pick(nodes, S, r"min\. peptide length"), pick(nodes, S, r"max\. peptide length")
        text = f"Spectra were searched with {' and '.join(engines)}"
        text += f" against {db}" if db else ""
        if enzyme:
            text += f", using {enzyme}" + (f" with up to {missed} missed cleavages" if missed else "")
        if min_len and max_len:
            text += f" and a peptide length of {min_len}–{max_len} amino acids"
        sentences.append(text + ".")
        precursor = pick(nodes, S, r"precursor mass tolerance")
        fragment = pick(nodes, S, r"fragment mass tolerance")
        if precursor and fragment:
            sentences.append(f"Precursor and fragment mass tolerances were {precursor} and {fragment}.")
        static, dynamic = collect_modifications(nodes)
        if static:
            sentences.append("Static modifications: " + "; ".join(static) + ".")
        if dynamic:
            sentences.append("Dynamic modifications: " + "; ".join(dynamic) + ".")

    for node in nodes:
        if re.search(r"psm validation", node["category"], re.I):
            strict = pick([node], ".", r"target fdr \(strict\)")
            relaxed = pick([node], ".", r"target fdr \(relaxed\)")
            based = pick([node], ".", r"validation based on")
            if strict:
                text = f"PSMs were validated with {node['name']} (strict / relaxed target FDR {strict} / {relaxed}"
                sentences.append(text + (f", based on {based})." if based else")."))

    peptide_fdr = pick(nodes, r"peptide validation", r"target fdr \(strict\) for peptides")
    protein_fdr = pick(nodes, r"protein validation", r"target fdr \(strict\)")
    if peptide_fdr or protein_fdr:
        parts = [f"{peptide_fdr} at peptide level" if peptide_fdr else "",
                 f"{protein_fdr} at protein level" if protein_fdr else ""]
        sentences.append("The strict target FDR was " + " and ".join(p for p in parts if p) + ".")
    confidence = pick(nodes, r"peptide and protein filter", r"peptide confidence at least")
    min_peptides = pick(nodes, r"peptide and protein filter", r"minimum number of peptide sequences")
    if confidence:
        text = f"Only peptides with at least {confidence.lower()} confidence were kept"
        text += f", and proteins needed at least {min_peptides} peptide sequence(s)" if min_peptides else ""
        sentences.append(text + ".")
    if pick(nodes, r"protein grouping", r"apply strict parsimony principle") == "True":
        sentences.append("Protein groups were formed with the strict parsimony principle.")

    study = rep.get("study") or {}
    method = next((m["name"] for m in study.get("methods", []) if m.get("name")), "")
    Q = r"quantif"
    tolerance = pick(nodes, Q, r"integration tolerance")
    co_isolation = pick(nodes, Q, r"co-isolation threshold")
    sn = pick(nodes, Q, r"average reporter s/n threshold")
    if tolerance or co_isolation:
        text = "Reporter ions" + (f" ({method})" if method else "") + " were quantified"
        text += f" with an integration tolerance of {tolerance}" if tolerance else ""
        details = [f"co-isolation threshold {co_isolation}%" if co_isolation else "",
                   f"average reporter S/N threshold {sn}" if sn else "",
                   f"peptides used: {pick(nodes, Q, r'peptides to use')}" if pick(nodes, Q, r"peptides to use") else ""]
        details = [d for d in details if d]
        sentences.append(text + (" (" + ", ".join(details) + ")" if details else "") + ".")
    normalization = pick(nodes, Q, r"normalization mode")
    scaling = pick(nodes, Q, r"scaling mode")
    if normalization or scaling:
        parts = [f"normalization: {normalization}" if normalization else "",
                 f"scaling: {scaling}" if scaling else ""]
        sentences.append("Abundances – " + "; ".join(p for p in parts if p) + ".")
    return " ".join(sentences)


# ======================= HTML rendering =======================

CSS = """
:root {
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink-2: #52514e; --muted: #77756f;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,.10);
  --accent: #2a78d6; --accent-soft: #e9f1fb; --warn-ink: #7a4d00; --warn-bg: #fff6e0;
}
@media screen and (prefers-color-scheme: dark) {
  :root {
    --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink-2: #c3c2b7; --muted: #9a988f;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,.10);
    --accent: #3987e5; --accent-soft: #17263a; --warn-ink: #fab219; --warn-bg: #2a2110;
  }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
  font: 14px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; }
main { max-width: 1060px; margin: 0 auto; padding: 40px 24px 64px; }
h1 { font-size: 28px; line-height: 1.2; margin: 4px 0 8px; letter-spacing: -.01em; }
h2 { font-size: 21px; margin: 0 0 4px; display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }
h3 { font-size: 16px; margin: 32px 0 10px; padding-bottom: 6px; border-bottom: 1px solid var(--grid); }
h4 { font-size: 13px; margin: 18px 0 6px; color: var(--ink-2); text-transform: uppercase; letter-spacing: .04em; }
p { margin: 6px 0; }
a { color: var(--accent); }
.eyebrow { color: var(--accent); font-weight: 600; font-size: 13px; letter-spacing: .06em; text-transform: uppercase; }
.muted { color: var(--muted); }
.small { font-size: 12px; }
.lede { color: var(--ink-2); margin: 0; }
.card { background: var(--surface); border: 1px solid var(--border); border-radius: 12px; padding: 24px 28px; margin: 24px 0; }
.badge { font-size: 11px; font-weight: 600; letter-spacing: .05em; text-transform: uppercase;
  color: var(--accent); background: var(--accent-soft); border-radius: 999px; padding: 2px 10px; }
.facts { display: grid; grid-template-columns: max-content 1fr; gap: 4px 20px; margin: 12px 0 0; }
.facts dt { color: var(--muted); }
.facts dd { margin: 0; overflow-wrap: anywhere; }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 12px; margin: 20px 0 4px; }
.tile { border: 1px solid var(--border); border-radius: 10px; padding: 12px 14px; background: var(--page); }
.tile .label { color: var(--ink-2); font-size: 12px; }
.tile .value { font-size: 22px; font-weight: 600; margin-top: 2px; }
.tile .note { color: var(--muted); font-size: 11px; }
.table-wrap { overflow-x: auto; margin: 8px 0; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--grid); vertical-align: top; }
th { color: var(--ink-2); font-weight: 600; font-size: 12px; border-bottom: 1px solid var(--axis); }
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
td.nowrap { white-space: nowrap; }
tr.group th { background: var(--page); color: var(--ink); font-size: 12px; text-transform: uppercase; letter-spacing: .04em; }
tbody tr:hover td { background: var(--accent-soft); }
td.barcell { width: 22%; min-width: 90px; vertical-align: middle; }
.bar { height: 12px; }
.bar span { display: block; height: 100%; background: var(--accent); border-radius: 0 4px 4px 0; }
.adv { font-size: 10px; color: var(--muted); border: 1px solid var(--grid); border-radius: 4px; padding: 0 4px; margin-left: 4px; }
.methods { background: var(--page); border-left: 3px solid var(--accent); padding: 12px 16px; border-radius: 0 8px 8px 0; }
.warn { background: var(--warn-bg); color: var(--warn-ink); border-radius: 8px; padding: 8px 12px; margin: 6px 0; }
ul.log { margin: 6px 0; padding-left: 18px; }
ul.log li { margin: 2px 0; }
.tree ul { list-style: none; margin: 0; padding-left: 18px; border-left: 1px solid var(--grid); }
.tree > ul { padding-left: 0; border-left: 0; }
.tree li { margin: 3px 0; }
.tree .node { font-weight: 600; }
.wf { border: 1px solid var(--border); border-radius: 10px; padding: 14px 18px; margin: 12px 0; background: var(--page); }
details { margin: 8px 0; }
summary { cursor: pointer; font-weight: 600; }
.toc ol { margin: 6px 0; padding-left: 20px; }
footer { color: var(--muted); font-size: 12px; margin-top: 32px; }
@page { size: A4; margin: 14mm 12mm 16mm;
  @bottom-right { content: "Page " counter(page) " of " counter(pages); font-size: 9px; color: #77756f; } }
@media print {
  body { background: #fff; font-size: 11px; }
  main { max-width: none; padding: 0; }
  .card { border: 0; padding: 0; margin: 0 0 12px; background: #fff; }
  .card.file { break-before: page; }
  .tiles { grid-template-columns: repeat(3, 1fr); }
  h3, h4 { break-after: avoid; }
  tr, .tile, .wf { break-inside: avoid; }
  .table-wrap { overflow: visible; }
  table { font-size: 9px; }
  th { font-size: 8.5px; }
  th, td { padding: 3px 5px; }
  td.barcell { width: 14%; min-width: 50px; }
  h1 { font-size: 22px; }
  h2 { font-size: 16px; }
  h3 { font-size: 13px; margin-top: 20px; }
  .tile .value { font-size: 17px; }
  .bar span, .badge, .methods, .warn, .tile, tr.group th { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
  tbody tr:hover td { background: none; }
  summary { list-style: none; }
  .toc { display: none; }
}
@media (max-width: 640px) {
  main { padding: 20px 16px 40px; }
  .card { padding: 16px; }
  .facts { grid-template-columns: 1fr; }
  .facts dt { margin-top: 6px; }
}
"""


def table(headers, rows, num=(), titles=None, bars=(), nowrap=()):
    out = ['<div class="table-wrap"><table><thead><tr>']
    for i, h in enumerate(headers):
        out.append(f'<th class="num">{esc(h)}</th>' if i in num else f"<th>{esc(h)}</th>")
    out.append("</tr></thead><tbody>")
    for r, row in enumerate(rows):
        title = f' title="{esc(titles[r])}"' if titles else ""
        cells = []
        for i, cell in enumerate(row):
            cls = ("num" if i in num else "barcell" if i in bars else "nowrap" if i in nowrap else "")
            cells.append(f'<td class="{cls}">{esc(cell)}</td>' if cls else f"<td>{esc(cell)}</td>")
        out.append(f"<tr{title}>{''.join(cells)}</tr>")
    out.append("</tbody></table></div>")
    return Html("".join(out))


def bar(value, maximum, label):
    pct = 100 * value / maximum if value and maximum else 0
    return Html(f'<div class="bar" role="img" aria-label="{esc(label)}"><span style="width:{pct:.1f}%"></span></div>')


def facts(items):
    rows = "".join(f"<dt>{esc(k)}</dt><dd>{esc(v)}</dd>" for k, v in items if v not in (None, ""))
    return Html(f'<dl class="facts">{rows}</dl>') if rows else Html("")


def tiles(items):
    cells = "".join(f'<div class="tile"><div class="label">{esc(label)}</div>'
                    f'<div class="value">{esc(value)}</div><div class="note">{esc(note)}</div></div>'
                    for label, value, note in items if value not in (None, "", "–"))
    return Html(f'<div class="tiles">{cells}</div>') if cells else Html("")


class Section:
    """Collects sub-headings of one file section for the table of contents."""

    def __init__(self, prefix):
        self.prefix, self.toc = prefix, []

    def h3(self, title):
        anchor = f"{self.prefix}-{len(self.toc) + 1}"
        self.toc.append((anchor, title))
        return f'<h3 id="{anchor}">{esc(title)}</h3>'


def render_study_design(sec, study, channel_sums=None):
    out = [sec.h3("Study design")]
    methods = study.get("methods", [])
    out.append(facts([
        ("Study", study.get("name")),
        ("Description", study.get("description")),
        ("Created", study.get("created")),
        ("Quantification", "; ".join(f'{m["name"]} – {m["kind"]}' for m in methods)),
        ("Method description", "; ".join(m["description"] for m in methods if m["description"])),
        ("Channels", "; ".join(", ".join(m["channels"]) for m in methods if m["channels"])),
        ("Study factors", ", ".join(study.get("factors", [])) or "none defined"),
        ("Grouping (for ratios)", ", ".join(study.get("grouping", []))),
        ("Quan ratios", "; ".join(study.get("ratios", [])) or ("none defined" if study.get("grouping") else "")),
    ]))
    samples = study.get("samples", [])
    if samples:
        factor_cols = study.get("factor_columns", [])
        many_sets = len({s["file_set"] for s in samples}) > 1
        mz = study.get("reporter_mz", {})
        loading = channel_sums if channel_sums and len(channel_sums) == len(samples) else None
        headers = ["Sample", "Name", "Channel"] + (["Reporter m/z"] if mz else []) + ["Type"]
        headers += factor_cols + (["File set"] if many_sets else [])
        num = {headers.index("Reporter m/z")} if mz else set()
        if loading:
            headers += ["Abundance (sum)", "vs. mean", "Channel loading"]
            num |= {len(headers) - 3, len(headers) - 2}
        rows, titles = [], []
        mean = sum(loading) / len(loading) if loading else None
        for i, s in enumerate(samples):
            row = [s["label"], s["name"], s["channel"]] + ([mz.get(s["channel"], "")] if mz else [])
            row += [s["type"]] + [s["factors"].get(f, "") for f in factor_cols]
            row += [s["file_set"]] if many_sets else []
            title = f'{s["label"]} {s["channel"]}'
            if loading:
                row += [fmt_compact(loading[i]), f"{100 * loading[i] / mean:.0f} %" if mean else "",
                        bar(loading[i], max(loading), f'{s["channel"]}: {fmt_compact(loading[i])}')]
                title += f": summed abundance {fmt_compact(loading[i])} ({100 * loading[i] / mean:.0f} % of mean)" if mean else ""
            rows.append(row)
            titles.append(title)
        out.append("<h4>Samples</h4>")
        out.append(table(headers, rows, num=num, titles=titles, bars={len(headers) - 1} if loading else ()))
        if loading:
            out.append('<p class="muted small">Abundance (sum): raw reporter abundances summed over all master '
                       'proteins that pass the PD filters (contaminants excluded), before normalization. '
                       'Bars are relative to the largest channel.</p>')
    return "".join(out)


def render_file_sets(study):
    out = []
    for fs in study.get("file_sets", []):
        out.append(f'<h4>File set {esc(fs["label"])} – {esc(fs["name"])}</h4>')
        out.append(f'<p class="muted small">{len(fs["files"])} files'
                   + (f' · quan method {esc(fs["method"])}' if fs["method"] else "") + "</p>")
        rows = [[f["fraction"], f["path"], f["time"], fmt_size(f["size"])] for f in fs["files"]]
        out.append(table(["Fraction", "File", "File time", "Size"], rows, num={3}, nowrap={2}))
    return "".join(out)


def render_results(sec, rep):
    c = rep.get("counts", {})
    if not c:
        return ""
    out = [sec.h3("Identification results")]
    groups = c.get("protein_groups") or (None, None)
    peptides = c.get("peptides") or (None, None, None)
    psms = c.get("psms") or (None, None, None)
    rows = []
    if c.get("spectra") is not None:
        rows.append(["MS/MS spectra in result", fmt_int(c["spectra"]), "–", "–"])
    if psms[0] is not None:
        rows.append(["PSMs", fmt_int(psms[0]), fmt_int(psms[1]), fmt_int(psms[2])])
    if peptides[0] is not None:
        rows.append(["Peptide groups", fmt_int(peptides[0]), fmt_int(peptides[1]), fmt_int(peptides[2])])
    if c.get("proteins") is not None:
        rows.append(["Proteins (incl. non-master)", fmt_int(c["proteins"]), "–", "–"])
    if c.get("master") is not None:
        rows.append(["Master proteins", fmt_int(c["master"]), fmt_int(c["master_by_conf"].get(HIGH)), "–"])
    if groups[0] is not None:
        rows.append(["Protein groups", fmt_int(groups[0]), "–", fmt_int(groups[1])])
    out.append(table(["Level", "Total", "High confidence", "Passing PD filters"], rows, num={1, 2, 3}))
    if c.get("master_by_conf"):
        parts = [f"{CONFIDENCE_NAMES.get(k, 'not set' if k is None else k)} {fmt_int(v)}"
                 for k, v in sorted(c["master_by_conf"].items(), key=lambda kv: -(kv[0] or 0))]
        out.append(f'<p class="muted small">Master proteins by protein FDR confidence: {esc(" · ".join(parts))}'
                   f' · flagged as contaminant: {fmt_int(c.get("master_contaminants"))}.</p>')
    if not rep.get("detailed"):
        out.append('<p class="muted small">Peptide and PSM counts were skipped (--fast).</p>')
    return "".join(out)


def render_raw_files(sec, rep):
    files = rep.get("raw_files", [])
    if not files:
        return ""
    per_file = rep.get("per_file", {})
    out = [sec.h3("Raw files")]
    instruments = sorted({f"{f['instrument']} (software {f['instrument_software']})" if f["instrument_software"]
                          else f["instrument"] for f in files if f["instrument"]})
    dates = sorted(f["acquired"] for f in files if f["acquired"])
    total = sum(float(f["size"] or 0) for f in files)
    rt = [float(x) for f in files for x in re.findall(r"\d+(?:\.\d+)?", str(f["rt_range"]))]
    folders = sorted({str(f["path"])[:-len(f["name"])].rstrip("\\/") for f in files
                      if f["name"] and str(f["path"]).endswith(f["name"])} - {""})
    out.append(facts([("Number of files", str(len(files))), ("Instrument", "; ".join(instruments)),
                      ("Acquired", f"{dates[0]} → {dates[-1]}" if dates else ""),
                      ("Total size", fmt_size(total) if total else ""),
                      ("Retention time range", f"{min(rt):.1f}–{max(rt):.1f} min" if rt else ""),
                      ("Folder", "; ".join(folders))]))
    has_counts = any(per_file.get(f["name"].lower()) for f in files)
    if has_counts:
        stats = [per_file.get(f["name"].lower(), {}) for f in files]
        top = max((s.get("psms_filtered") or 0 for s in stats), default=0)
        headers = ["Fraction", "Raw file", "Acquired", "Size", "MS/MS spectra", "Filtered PSMs", "ID rate", "Filtered PSMs"]
        rows, titles = [], []
        for f, s in zip(files, stats):
            spectra, good = s.get("spectra"), s.get("psms_filtered")
            rate = f"{100 * good / spectra:.1f} %" if spectra and good is not None else "–"
            rows.append([f["fraction"], f["name"], f["acquired"], fmt_size(f["size"]), fmt_int(spectra),
                         fmt_int(good), rate, bar(good, top, f"{f['fraction']}: {fmt_int(good)} PSMs")])
            titles.append(f"{f['fraction']} {f['name']}: {fmt_int(good)} filtered PSMs of "
                          f"{fmt_int(s.get('psms'))} PSMs, {fmt_int(spectra)} MS/MS spectra")
        out.append(table(headers, rows, num={3, 4, 5, 6}, titles=titles, bars={7}, nowrap={0, 2}))
        out.append('<p class="muted small">Filtered PSMs: PSMs that pass the PD result filters. '
                   'ID rate: filtered PSMs / MS/MS spectra stored in the result.</p>')
    else:
        rows = [[f["fraction"], f["name"], f["acquired"], fmt_size(f["size"]), f["instrument"]] for f in files]
        out.append(table(["Fraction", "Raw file", "Acquired", "Size", "Instrument"], rows, num={3}, nowrap={0, 2}))
    return "".join(out)


def render_key_settings(sec, workflows, title="Search, validation and quantification settings"):
    sections = collect_key_settings(workflows)
    if not sections:
        return ""
    out = [sec.h3(title) if sec else f"<h4>{esc(title)}</h4>",
           '<div class="table-wrap"><table><thead><tr><th>Setting</th><th>Value</th><th>Node</th></tr></thead><tbody>']
    for name, rows in sections.items():
        out.append(f'<tr class="group"><th colspan="3">{esc(name)}</th></tr>')
        for label, value, node in rows:
            out.append(f'<tr><td>{esc(label)}</td><td>{esc(value)}</td><td class="muted nowrap">{esc(node)}</td></tr>')
    out.append("</tbody></table></div>")
    return "".join(out)


def render_node_tree(nodes):
    numbers = {n["number"] for n in nodes}
    children = {}
    for n in nodes:
        parents = [p.strip() for p in re.split(r"[;,]", n["parent"] or "") if p.strip() in numbers]
        children.setdefault(parents[0] if parents else None, []).append(n)
    seen = set()

    def walk(key):
        items = []
        for n in children.get(key, []):
            if n["number"] in seen:
                continue
            seen.add(n["number"])
            items.append(f'<li><span class="node">{esc(n["name"])}</span> '
                         f'<span class="muted small">{esc(n["category"])} · v{esc(n["version"])}</span>'
                         f'{walk(n["number"])}</li>')
        return f"<ul>{''.join(items)}</ul>" if items else ""

    return Html(f'<div class="tree">{walk(None)}</div>')


def render_workflows(sec, workflows):
    if not workflows:
        return ""
    out = [sec.h3("Workflows")]
    for wf in workflows:
        out.append('<div class="wf">')
        out.append(f'<p><span class="badge">{esc(wf.get("type") or "workflow")}</span> '
                   f'<strong>{esc(wf.get("name"))}</strong></p>')
        out.append(facts([("Template", wf.get("template")), ("Description", wf.get("description")),
                          ("Started", wf.get("started")), ("Ran (log)", wf.get("ran")),
                          ("User", show_name(wf.get("user"))), ("Computer", show_name(wf.get("machine"))),
                          ("PD version", wf.get("software"))]))
        if wf["nodes"]:
            out.append("<h4>Processing nodes</h4>")
            out.append(render_node_tree(wf["nodes"]))
        out.append("</div>")
    return "".join(out)


def render_log(sec, rep):
    log, warnings = rep.get("log", []), rep.get("warnings", [])
    if not log and not warnings:
        return ""
    out = [sec.h3("Processing log")]
    for node, text in warnings:
        out.append(f'<div class="warn"><strong>Warning</strong> · {esc(node)}: {esc(text)}</div>')
    if log:
        out.append('<ul class="log">' + "".join(
            f'<li><span class="muted">{esc(node)}:</span> {esc(text)}</li>' for node, text in log) + "</ul>")
    return "".join(out)


def render_versions(sec, rep):
    versions = rep.get("versions", {})
    workflows = rep.get("workflows", [])
    out = [sec.h3("Software versions")]
    upgraded = ""
    if versions.get("latest") and versions.get("latest") != versions.get("created"):
        upgraded = f'{versions["latest"]} ({versions.get("latest_date", "")})'
    out.append(facts([("Created with Proteome Discoverer", versions.get("created")),
                      ("File format last upgraded by PD", upgraded)]))
    if upgraded:
        out.append('<p class="muted small">The file was opened in a newer Proteome Discoverer version, which '
                   'upgraded its format. The processing itself was done with the version it was created with.</p>')
    rows = [[wf.get("type", ""), n["name"], n["version"]] for wf in workflows for n in wf["nodes"]]
    if rows:
        out.append(table(["Workflow", "Node", "Node version"], rows))
    history = versions.get("history") or []
    if history:
        rows = [[r.get("Kind", ""), r.get("Version", ""), fmt_date(r.get("Date")), r.get("SoftwareVersion", ""),
                 short(r.get("Comment", ""), 140)] for r in history]
        out.append(f"<details><summary>File schema history ({len(rows)} entries)</summary>")
        out.append(table(["Kind", "Schema", "Date", "PD version", "Comment"], rows))
        out.append("</details>")
    return "".join(out)


def render_all_parameters(sec, workflows, expand, title="Appendix: all workflow parameters"):
    if not any(wf["nodes"] for wf in workflows):
        return ""
    out = [sec.h3(title),
           '<p class="muted small">Every visible parameter of every node, as stored in the file. '
           'Empty parameters are left out; <span class="adv">adv.</span> marks advanced parameters.</p>']
    for wf in workflows:
        out.append(f'<h4>{esc(wf.get("type", ""))} workflow – {esc(wf.get("name", ""))}</h4>')
        for n in wf["nodes"]:
            params = [p for p in n["params"] if p["value"]]
            count = f"{len(params)} parameter{'s' if len(params) != 1 else ''}" if params else "no parameters"
            out.append(f'<details{" open" if expand and params else ""}><summary>{esc(n["name"])} '
                       f'<span class="muted small">({count})</span></summary>')
            if params:
                rows = [[p["group"], Html(esc(p["name"]) + ('<span class="adv">adv.</span>' if p["advanced"] else "")),
                         scrub(short(p["value"]))] for p in params]
                out.append(table(["Group", "Parameter", "Value"], rows))
            out.append("</details>")
    return "".join(out)


def file_header(rep, sec, kind_label):
    return (f'<section class="card file" id="{sec.prefix}"><h2><span class="badge">{esc(kind_label)}</span>'
            f'{esc(rep["name"])}</h2>'
            f'<p class="lede">{esc(fmt_size(rep["size"]))} · last modified {esc(rep["modified"])}</p>')


def render_result_file(rep, sec, expand):
    kind = "pdResult" if rep["name"].lower().endswith(".pdresult") else (
        "msf" if rep["name"].lower().endswith(".msf") else "SQLite")
    out = [file_header(rep, sec, kind)]
    c = rep.get("counts", {})
    psms = c.get("psms") or (None, None, None)
    peptides = c.get("peptides") or (None, None, None)
    groups = c.get("protein_groups") or (None, None)
    out.append(tiles([("Raw files", fmt_int(len(rep.get("raw_files", []))) if rep.get("raw_files") else "", ""),
                      ("MS/MS spectra", fmt_int(c.get("spectra")), "stored in result"),
                      ("PSMs", fmt_int(psms[2]), "passing PD filters"),
                      ("Peptide groups", fmt_int(peptides[2]), "passing PD filters"),
                      ("Protein groups", fmt_int(groups[1]), "passing PD filters"),
                      ("Master proteins", fmt_int((c.get("master_by_conf") or {}).get(HIGH)),
                       "high protein FDR confidence")]))
    versions = rep.get("versions", {})
    processing = next((w for w in rep.get("workflows", []) if w.get("type") == "Processing"), {})
    out.append(facts([("Study", (rep.get("study") or {}).get("name")),
                      ("Created with", f'Proteome Discoverer {versions["created"]}' if versions.get("created") else ""),
                      ("Processed", processing.get("started")),
                      ("Processed by", show_name(processing.get("user"))),
                      ("Computer", show_name(processing.get("machine")))]))
    if rep.get("study"):
        out.append(render_study_design(sec, rep["study"], rep.get("channel_sums")))
    summary = methods_summary(rep)
    if summary:
        out.append(sec.h3("Methods summary (draft)"))
        out.append(f'<div class="methods">{esc(summary)}</div>')
        out.append('<p class="muted small">Written automatically from the parameters stored in the file. '
                   'Please check it against the settings below before using it in a manuscript.</p>')
    out.append(render_results(sec, rep))
    out.append(render_raw_files(sec, rep))
    out.append(render_key_settings(sec, rep.get("workflows", [])))
    out.append(render_workflows(sec, rep.get("workflows", [])))
    out.append(render_log(sec, rep))
    out.append(render_versions(sec, rep))
    out.append(render_all_parameters(sec, rep.get("workflows", []), expand))
    if LIST_ALL_TABLES:
        out.append(sec.h3(f'All tables ({len(rep.get("tables", []))})'))
        out.append(f'<p class="small">{esc(", ".join(rep.get("tables", [])))}</p>')
    out.append("</section>")
    return "".join(out)


def render_study_file(rep, sec, expand, shown_elsewhere):
    study = rep["study"]
    out = [file_header(rep, sec, "pdStudy")]
    out.append(facts([("Study", study.get("name")), ("Created", study.get("created")),
                      ("Last changed", study.get("changed")), ("Study folder", study.get("directory")),
                      ("Description", study.get("description")),
                      ("Samples", str(len(study.get("samples", [])))),
                      ("Input files", str(sum(len(fs["files"]) for fs in study.get("file_sets", []))))]))
    if shown_elsewhere:
        out.append('<p class="muted small">Samples, channels, raw files and workflow settings of this study are '
                   'shown in the result file section of this report.</p>')
    else:
        out.append(render_study_design(sec, study))
        if study.get("file_sets"):
            out.append(sec.h3("Input files"))
            out.append(render_file_sets(study))
    analyses = rep.get("analyses", [])
    if analyses:
        out.append(sec.h3("Analyses in this study"))
        rows = []
        for a in analyses:
            for s in a["steps"]:
                rows.append([a["number"], a["created"], s["type"], s["workflow"], s["file"], s["state"]])
        out.append(table(["#", "Created", "Workflow type", "Workflow", "Output file", "State"], rows,
                         nowrap={1}))
        latest = analyses[-1]
        if latest["workflows"] and not shown_elsewhere:
            out.append(render_key_settings(sec, latest["workflows"],
                                           f'Settings of the latest analysis (#{latest["number"]})'))
            out.append(render_all_parameters(sec, latest["workflows"], expand,
                                             f'Appendix: all parameters of analysis #{latest["number"]}'))
    out.append("</section>")
    return "".join(out)


def render_other_file(rep, sec):
    out = [file_header(rep, sec, rep["kind"] if rep["kind"] == "error" else "file")]
    if rep.get("error"):
        out.append(f'<div class="warn">Could not read file: {esc(rep["error"])}</div>')
    hits = rep.get("hits") or {}
    if not hits and not rep.get("error"):
        out.append("<p>No readable version information found in this file.</p>")
    for name, found in hits.items():
        out.append(f"<h4>Version entries in {esc(name)}</h4>")
        out.append('<ul class="log">' + "".join(f"<li>{esc(short(h))}</li>" for h in found) + "</ul>")
    out.append("</section>")
    return "".join(out)


def render_report(reports, title, expand=False):
    result_studies = {(r.get("study") or {}).get("name") for r in reports if r["kind"] == "sqlite"}
    sections_html, toc = [], []
    for i, rep in enumerate(reports, 1):
        sec = Section(f"f{i}")
        if rep["kind"] == "sqlite":
            body = render_result_file(rep, sec, expand)
        elif rep["kind"] == "study":
            body = render_study_file(rep, sec, expand, rep["study"].get("name") in result_studies)
        else:
            body = render_other_file(rep, sec)
        sections_html.append(body)
        toc.append((sec.prefix, rep["name"], sec.toc))

    overview_rows = []
    for rep in reports:
        study = (rep.get("study") or {}).get("name", "")
        version = (rep.get("versions") or {}).get("created", "")
        kind = {"sqlite": "Result / processing file", "study": "Study file"}.get(rep["kind"], "Other")
        overview_rows.append([rep["name"], kind, study, version, fmt_size(rep["size"])])
    toc_html = "".join(
        f'<li><a href="#{anchor}">{esc(name)}</a><ol>'
        + "".join(f'<li><a href="#{a}">{esc(t)}</a></li>' for a, t in subs) + "</ol></li>"
        for anchor, name, subs in toc)
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="generator" content="pd-study-report {__version__}">
<title>{esc(title)}</title>
<style>{CSS}</style>
</head>
<body>
<main>
<header class="card">
  <div class="eyebrow">Proteome Discoverer study report</div>
  <h1>{esc(title)}</h1>
  <p class="lede">Generated {generated} · {len(reports)} file{'s' if len(reports) != 1 else ''}
  {'· anonymized' if ANONYMIZE else ''}</p>
  {table(["File", "Type", "Study", "Created with PD", "Size"], overview_rows, num={4})}
  <nav class="toc"><h4>Contents</h4><ol>{toc_html}</ol></nav>
</header>
{''.join(sections_html)}
<footer>
  Generated with pd-study-report {__version__}. All files were opened read-only. Numbers reflect the state
  stored in the files (PD filters as applied when the result was processed).
</footer>
</main>
</body>
</html>
"""


# ======================= PDF =======================

def find_browser():
    for name in ("msedge", "chrome", "google-chrome", "google-chrome-stable",
                 "chromium", "chromium-browser", "microsoft-edge"):
        found = shutil.which(name)
        if found:
            return found
    candidates = []
    for base in (os.environ.get("PROGRAMFILES(X86)"), os.environ.get("PROGRAMFILES"),
                 os.environ.get("LOCALAPPDATA")):
        if base:
            candidates += [os.path.join(base, "Microsoft", "Edge", "Application", "msedge.exe"),
                           os.path.join(base, "Google", "Chrome", "Application", "chrome.exe")]
    candidates += ["/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
                   "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
                   "/Applications/Chromium.app/Contents/MacOS/Chromium"]
    return next((c for c in candidates if os.path.isfile(c)), None)


def write_pdf(html_text, pdf_path):
    """Print the report to PDF with a headless Edge / Chrome / Chromium."""
    browser = find_browser()
    if not browser:
        print("PDF: no Edge, Chrome or Chromium found. Open the HTML report and use Print -> Save as PDF.")
        return False
    work = tempfile.mkdtemp(prefix="pd_report_")
    try:
        html_file = os.path.join(work, "report.html")
        with open(html_file, "w", encoding="utf-8") as f:
            f.write(html_text)
        if os.path.exists(pdf_path):
            os.remove(pdf_path)
        cmd = [browser, "--headless", "--disable-gpu", "--no-first-run", "--disable-extensions",
               f"--user-data-dir={os.path.join(work, 'profile')}", "--no-pdf-header-footer",
               "--print-to-pdf-no-header", f"--print-to-pdf={os.path.abspath(pdf_path)}",
               Path(html_file).as_uri()]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=180)
    except (OSError, subprocess.SubprocessError) as err:
        print(f"PDF: printing failed ({err}).")
    finally:
        shutil.rmtree(work, ignore_errors=True)
    if os.path.isfile(pdf_path) and os.path.getsize(pdf_path) > 0:
        return True
    print("PDF: the browser did not create a PDF. Open the HTML report and use Print -> Save as PDF.")
    return False


# ======================= main =======================

def inspect_file(path, detailed):
    rep = {"path": path, "name": os.path.basename(path), "size": os.path.getsize(path),
           "modified": datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M"),
           "detailed": detailed}
    try:
        if is_sqlite(path):
            rep.update(inspect_sqlite(path, detailed))
        elif path.lower().endswith(".pdstudy"):
            rep.update(inspect_study_file(path))
        else:
            rep.update(inspect_other_file(path))
    except Exception as err:
        rep.update({"kind": "error", "error": str(err)})
    return rep


def safe_filename(text):
    return re.sub(r"[^\w.-]+", "_", text).strip("_") or "report"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create an HTML/PDF report from Proteome Discoverer .pdResult, .msf and .pdStudy files.")
    parser.add_argument("paths", nargs="*", default=[SCRIPT_DIR],
                        help="files or folders to read (default: the folder of this script)")
    parser.add_argument("-o", "--output", help="path of the HTML report (default: next to the input files)")
    parser.add_argument("--pdf", action="store_true", help="also write a PDF (needs Edge, Chrome or Chromium)")
    parser.add_argument("--anonymize", action="store_true",
                        help="hide user names, computer names and folder paths")
    parser.add_argument("--fast", action="store_true",
                        help="skip the slow peptide and PSM counts (useful for very large files)")
    parser.add_argument("--title", help="report title (default: the study name)")
    parser.add_argument("--open", action="store_true", help="open the report when it is ready")
    parser.add_argument("--version", action="version", version=f"pd-study-report {__version__}")
    args, _unknown = parser.parse_known_args()  # tolerate extra args, e.g. from Jupyter
    return args


def main():
    global ANONYMIZE
    args = parse_args()
    ANONYMIZE = ANONYMIZE or args.anonymize
    detailed = COUNT_PEPTIDES_AND_PSMS and not args.fast

    files = collect_inputs(args.paths)
    if not files:
        print("No .pdResult, .msf or .pdStudy files found in: " + ", ".join(args.paths))
        return 1

    reports = []
    for path in files:
        size_mb = os.path.getsize(path) / 1024 ** 2
        print(f"Reading {os.path.basename(path)} ({size_mb:,.1f} MB) ...", flush=True)
        reports.append(inspect_file(path, detailed))
    reports.sort(key=lambda r: {"sqlite": 0, "study": 1}.get(r["kind"], 2))  # results first
    studies = list(dict.fromkeys(s for s in ((r.get("study") or {}).get("name") for r in reports) if s))
    title = args.title or (", ".join(studies) if studies else os.path.splitext(reports[0]["name"])[0])
    if args.output:
        html_path = args.output if args.output.lower().endswith((".html", ".htm")) else args.output + ".html"
    else:
        html_path = os.path.join(os.path.dirname(os.path.abspath(files[0])),
                                 f"PD_report_{safe_filename(title)}.html")

    try:
        with open(html_path, "w", encoding="utf-8") as f:
            f.write(render_report(reports, title))
    except OSError as err:  # e.g. read-only data folder
        print(f"Cannot write {html_path} ({err}); writing to the current folder instead.")
        html_path = os.path.join(os.getcwd(), os.path.basename(html_path))
        with open(html_path, "w", encoding="utf-8") as f:
            f.write(render_report(reports, title))

    print("\n" + "=" * 70)
    print("SUMMARY")
    for rep in reports:
        version = (rep.get("versions") or {}).get("created") or {"study": "study file"}.get(rep["kind"], rep["kind"])
        print(f"  {rep['name']}: {version}")
    print(f"\nHTML report: {os.path.abspath(html_path)}")
    if args.pdf:
        pdf_path = os.path.splitext(html_path)[0] + ".pdf"
        if write_pdf(render_report(reports, title, expand=True), pdf_path):
            print(f"PDF report:  {os.path.abspath(pdf_path)}")
    if args.open:
        webbrowser.open(Path(os.path.abspath(html_path)).as_uri())
    return 0


if __name__ == "__main__":
    sys.exit(main())
