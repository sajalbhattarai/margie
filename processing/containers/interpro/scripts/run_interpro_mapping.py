#!/usr/bin/env python3
"""Map Pfam and TIGRFAM signatures to InterPro IDs per gene."""

import argparse
import csv
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple


def log(message: str) -> None:
    """Log with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}]  {message}", flush=True)


def split_values(value: str) -> Iterable[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text:
        return []
    for token in re.split(r"[;,|]", text):
        token = token.strip()
        if token:
            yield token


def load_mapping(path: Path) -> Dict[Tuple[str, str], Set[str]]:
    mapping: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            db_name = (row.get("signature_db") or "").strip()
            signature = (row.get("signature_acc") or "").strip()
            ipr = (row.get("ipr_acc") or "").strip()
            if not db_name or not signature or not ipr:
                continue
            mapping[(db_name, signature)].add(ipr)
    return mapping


def load_signatures(path: Path, feature_col: str, signature_col: str, description_col: str = None) -> Tuple[Dict[str, Set[str]], Dict[str, Dict[str, str]]]:
    """
    Load signatures and descriptions from annotation file.
    
    Returns:
        Tuple of (signatures dict, descriptions dict)
        - signatures: {feature_id: set of signature names}
        - descriptions: {feature_id: {signature_name: description}}
    """
    signatures: Dict[str, Set[str]] = defaultdict(set)
    descriptions: Dict[str, Dict[str, str]] = defaultdict(dict)
    
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            feature_id = (row.get(feature_col) or "").strip()
            if not feature_id:
                continue
            
            # Get signature name(s)
            signature_list = list(split_values(row.get(signature_col, "")))
            for signature in signature_list:
                signatures[feature_id].add(signature)
            
            # Get description if column specified
            if description_col and len(signature_list) > 0:
                desc = (row.get(description_col) or "").strip()
                # For each signature in this row, associate the description
                # (In domain outputs, each row is one domain with one description)
                for signature in signature_list:
                    descriptions[feature_id][signature] = desc
    
    return signatures, descriptions


def mapped_iprs(mapping: Dict[Tuple[str, str], Set[str]], db_name: str, signatures: Set[str]) -> Set[str]:
    iprs: Set[str] = set()
    for signature in signatures:
        iprs.update(mapping.get((db_name, signature), set()))
    return iprs


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate InterPro per-gene mapping")
    parser.add_argument("--organism", required=True)
    parser.add_argument("--pfam", required=True)
    parser.add_argument("--tigrfam", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    pfam_path = Path(args.pfam)
    tigrfam_path = Path(args.tigrfam)
    mapping_path = Path(args.mapping)
    output_path = Path(args.output)

    log("InterPro Signature Mapping")
    log("=" * 60)
    log(f"Organism:     {args.organism}")
    log(f"Pfam input:   {pfam_path}")
    log(f"TIGRfam input: {tigrfam_path}")
    log(f"Mapping DB:   {mapping_path}")
    log(f"Output:       {output_path}")
    log("")

    log("Loading signature to InterPro mapping...")
    mapping = load_mapping(mapping_path)
    log(f"  Loaded {len(mapping)} signature mappings")
    log("")

    log("Loading Pfam signatures...")
    pfam_signatures, pfam_descriptions = load_signatures(pfam_path, "feature_id", "PFAM_accession", "PFAM_description")
    log(f"  Loaded Pfam signatures for {len(pfam_signatures)} features")
    log("")

    log("Loading TIGRfam signatures...")
    tigrfam_signatures, tigrfam_descriptions = load_signatures(tigrfam_path, "feature_id", "TIGRFAM_id", "TIGRFAM_description")
    log(f"  Loaded TIGRfam signatures for {len(tigrfam_signatures)} features")
    log("")

    genes = sorted(set(pfam_signatures.keys()) | set(tigrfam_signatures.keys()))
    log(f"Processing {len(genes)} total features with Pfam/TIGRfam annotations...")
    log("")

    # Track statistics
    pfam_only = 0
    tigrfam_only = 0
    both = 0
    agreement = 0
    conflict = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "organism",
                "feature_id",
                "INTERPRO_ipr_pfam_ids",
                "INTERPRO_ipr_tigrfam_ids",
                "INTERPRO_ipr_overlap_count",
                "INTERPRO_ipr_support_flags",
                "INTERPRO_pfam_signatures",
                "INTERPRO_tigrfam_signatures",
                "INTERPRO_pfam_descriptions",
                "INTERPRO_tigrfam_descriptions",
            ]
        )

        for feature_id in genes:
            pfam_vals = pfam_signatures.get(feature_id, set())
            tigrfam_vals = tigrfam_signatures.get(feature_id, set())

            # Track statistics
            if pfam_vals and tigrfam_vals:
                both += 1
            elif pfam_vals:
                pfam_only += 1
            elif tigrfam_vals:
                tigrfam_only += 1

            ipr_pfam = mapped_iprs(mapping, "Pfam", pfam_vals)
            ipr_tigrfam = mapped_iprs(mapping, "TIGRFAM", tigrfam_vals)
            overlap = ipr_pfam & ipr_tigrfam

            flags = []
            if overlap:
                flags.append("IPR_PFAM_TIGR_AGREE")
                agreement += 1
            elif ipr_pfam and ipr_tigrfam:
                flags.append("IPR_PFAM_TIGR_CONFLICT")
                conflict += 1
            elif ipr_pfam:
                flags.append("IPR_PFAM_ONLY")
            elif ipr_tigrfam:
                flags.append("IPR_TIGRFAM_ONLY")
            
            # Build description strings (signature: description pairs)
            pfam_desc_list = []
            for sig in sorted(pfam_vals):
                desc = pfam_descriptions.get(feature_id, {}).get(sig, "")
                if desc:
                    pfam_desc_list.append(f"{sig}: {desc}")
                else:
                    pfam_desc_list.append(sig)
            
            tigrfam_desc_list = []
            for sig in sorted(tigrfam_vals):
                desc = tigrfam_descriptions.get(feature_id, {}).get(sig, "")
                if desc:
                    tigrfam_desc_list.append(f"{sig}: {desc}")
                else:
                    tigrfam_desc_list.append(sig)

            writer.writerow(
                [
                    args.organism,
                    feature_id,
                    ";".join(sorted(ipr_pfam)),
                    ";".join(sorted(ipr_tigrfam)),
                    len(overlap),
                    ";".join(flags),
                    ";".join(sorted(pfam_vals)),
                    ";".join(sorted(tigrfam_vals)),
                    "; ".join(pfam_desc_list),
                    "; ".join(tigrfam_desc_list),
                ]
            )

    log("Statistics:")
    log(f"  Total features:            {len(genes)}")
    log(f"  Features with Pfam only:   {pfam_only}")
    log(f"  Features with TIGRfam only: {tigrfam_only}")
    log(f"  Features with both:        {both}")
    log(f"  InterPro agreement:        {agreement}")
    log(f"  InterPro conflict:         {conflict}")
    log("")
    log(f"✓ Wrote {len(genes)} gene mappings to {output_path}")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
