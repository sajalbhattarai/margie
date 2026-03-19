#!/usr/bin/env python3
"""Shared concordance utilities for model building and scoring."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Dict, List, Optional
import math


def first_existing_column(columns: List[str], candidates: List[str]) -> Optional[str]:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def normalized_string(value) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value).strip()


def extract_gene_set(row, tool_column_map: Dict[str, List[str]], columns: List[str]) -> Dict[str, str]:
    """Return per-tool annotation identifiers for one gene row."""
    gene_set: Dict[str, str] = {}

    for tool, candidates in tool_column_map.items():
        column = first_existing_column(columns, candidates)
        if column is None:
            continue

        value = normalized_string(row.get(column, ""))
        if not value:
            continue

        # RAST IDs are feature-specific and can dominate signatures; use generalized tag.
        if tool == "RAST":
            value = "RAST_ID_GENERIC"

        gene_set[tool] = value

    return gene_set


def get_operon_info(row, columns: List[str]) -> Dict[str, object]:
    """Extract operon context using whichever OPERON columns exist in the table."""
    id_col = first_existing_column(columns, ["OPERON_operon_id", "OPERON_id", "operon_id", "Operon_ID"])
    size_col = first_existing_column(columns, ["OPERON_operon_size", "OPERON_size", "operon_size"])
    pos_col = first_existing_column(columns, ["OPERON_operon_position", "OPERON_position", "operon_position"])

    operon_id = normalized_string(row.get(id_col, "")) if id_col else ""
    if not operon_id:
        return {
            "has_operon": False,
            "operon_id": "",
            "operon_size": 0,
            "position_in_operon": 0,
        }

    size_raw = normalized_string(row.get(size_col, "0")) if size_col else "0"
    pos_raw = normalized_string(row.get(pos_col, "0")) if pos_col else "0"

    try:
        operon_size = int(float(size_raw))
    except ValueError:
        operon_size = 0

    try:
        position = int(float(pos_raw))
    except ValueError:
        position = 0

    return {
        "has_operon": True,
        "operon_id": operon_id,
        "operon_size": operon_size,
        "position_in_operon": position,
    }


def build_signature(gene_set: Dict[str, str]) -> str:
    if not gene_set:
        return ""
    return " ## ".join(f"{tool}|{gene_set[tool]}" for tool in sorted(gene_set.keys()))


def now_utc_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def calculate_concordance_score(
    signature: str,
    pattern_metadata: Dict[str, dict],
    gene_set: Dict[str, str],
    operon_info: Dict[str, object],
    total_genomes: int,
) -> int:
    """Compute concordance score in range 0-100."""
    if not signature or signature not in pattern_metadata:
        num_tools = len(gene_set)
        if num_tools >= 5:
            return 25
        if num_tools >= 3:
            return 15
        if num_tools >= 1:
            return 5
        return 0

    pattern = pattern_metadata[signature]
    total_count = int(pattern.get("total_count", 0))
    num_genomes = int(pattern.get("num_genomes", 0))

    # Frequency component (0-40)
    if total_count >= 1000:
        frequency_score = 40
    elif total_count >= 500:
        frequency_score = 35
    elif total_count >= 100:
        frequency_score = 30
    elif total_count >= 50:
        frequency_score = 25
    elif total_count >= 20:
        frequency_score = 20
    elif total_count >= 10:
        frequency_score = 15
    elif total_count >= 5:
        frequency_score = 10
    else:
        frequency_score = 5

    # Breadth component (0-30)
    denom = max(1, total_genomes)
    breadth_pct = (num_genomes / denom) * 100.0
    if breadth_pct >= 50:
        breadth_score = 30
    elif breadth_pct >= 25:
        breadth_score = 25
    elif breadth_pct >= 10:
        breadth_score = 20
    elif breadth_pct >= 5:
        breadth_score = 15
    elif breadth_pct >= 2:
        breadth_score = 10
    else:
        breadth_score = 5

    # Operon consistency (0-15)
    operon_freq = float(pattern.get("operon_frequency", 0.0))
    gene_in_operon = bool(operon_info.get("has_operon", False))
    if total_count >= 5:
        if operon_freq >= 0.7 and gene_in_operon:
            operon_score = 15
        elif operon_freq <= 0.3 and not gene_in_operon:
            operon_score = 15
        elif 0.3 < operon_freq < 0.7:
            operon_score = 10
        else:
            operon_score = 5
    else:
        operon_score = 10

    # Tool agreement (0-15)
    tools = len(gene_set)
    if tools >= 8:
        tool_score = 15
    elif tools >= 6:
        tool_score = 12
    elif tools >= 4:
        tool_score = 10
    elif tools >= 2:
        tool_score = 7
    elif tools == 1:
        tool_score = 3
    else:
        tool_score = 0

    return min(100, frequency_score + breadth_score + operon_score + tool_score)
