#!/usr/bin/env python3
"""Audit acquisition-source overlap from version-pinned provider metadata.

Inputs remain in the private audit workspace. Published outputs contain the
crosswalk and evidence, not redistributed raw expression or image data.
"""
import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlsplit

import pandas as pd


def normalized_url(value):
    u = urlsplit(str(value).strip())
    return re.sub(r"^www\.", "", u.netloc.lower()) + u.path.lower().rstrip("/")


def normalized_label(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def valid_id(value):
    return isinstance(value, str) and len(value.strip()) > 10


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    root = args.root
    inp = root / "inputs"
    out = root / "outputs"
    out.mkdir(exist_ok=True)
    release = pd.read_csv(inp / "release_samples.csv").fillna("")
    inventory = pd.read_csv(inp / "sample_inventory.csv").fillna("")
    hest = pd.read_csv(inp / "HEST_v1_1_0.csv").fillna("")
    stimage = pd.read_csv(inp / "STimage_meta_all_gene02122025.csv").fillna("")
    native = json.loads((inp / "hest_metadata_all.json").read_text())
    st_verification = json.loads((root / "stimage_coordinate_verification.json").read_text())
    image_verification = json.loads((root / "ffmb028_image_verification.json").read_text())
    assert image_verification["same_image_confirmed"]
    hest_revision = (inp / "hest_native_revision.txt").read_text().strip()
    st_revision = json.loads((inp / "stimage_hub_api.json").read_text())["sha"]
    assert len(release) == release.sample_id.nunique() == 124
    assert int(release.gene_positive_observations.sum()) == 38472192
    assert len(hest) == hest.id.nunique() == 1229
    assert sha256(inp / "HEST_v1_1_0.csv") == "61367f7991591550672cd214aa9d70263a65aff751a722ccb513dcab439110c5"
    assert Counter(hest.st_technology) == {"Visium": 602, "Spatial Transcriptomics": 552, "Visium HD": 10, "Xenium": 65}
    assert len(stimage) == stimage.slide.nunique() == 1149

    # Full source inventory determines whether a page represents one acquisition;
    # using only the processed subset would incorrectly collapse multi-sample pages.
    page_acquisitions = defaultdict(set)
    for row in inventory.to_dict("records"):
        for url in row["source_page_aliases"].split(";"):
            page_acquisitions[normalized_url(url)].add(row["sample_id"])
    inv_by_resource = {row["official_resource_id"]: row for row in inventory.to_dict("records")}
    assert len(inv_by_resource) == len(inventory)
    candidates, matched = [], []
    relevant_hest = hest[hest.st_technology.str.contains("HD|Xenium")]
    assert len(relevant_hest) == 75 and set(relevant_hest.id) == set(native)
    for rr in release.to_dict("records"):
        ii = inv_by_resource[rr["official_resource_id"]]
        aliases = set(ii["resource_aliases"].split(";"))
        urls = set(map(normalized_url, ii["source_page_aliases"].split(";")))
        for hh in relevant_hest.to_dict("records"):
            if (rr["platform"] == "HD") != ("HD" in hh["st_technology"]):
                continue
            mm = native[hh["id"]]
            source_url = normalized_url(hh["download_page_link1"])
            shared_url = source_url in urls
            equal_dims = (int(hh["fullres_px_width"]), int(hh["fullres_px_height"])) == (
                int(rr["image_width_px"]), int(rr["image_height_px"]))
            reasons, rules = [], []
            if mm.get("Sample ID") in aliases:
                reasons.append("native_sample_id")
                rules.append("native_resource_identity")
            roi1, roi2 = mm.get("roi_uuid"), rr["native_roi_uuid"]
            roi_conflict = valid_id(roi1) and valid_id(roi2) and roi1 != roi2
            if valid_id(roi1) and roi1 == roi2:
                reasons.append("native_roi_uuid")
                rules.append("native_roi_identity")
            if valid_id(mm.get("experiment_uuid")) and mm["experiment_uuid"] == rr["native_experiment_uuid"]:
                reasons.append("same_experiment_not_sufficient")
            if shared_url:
                reasons.append("source_url_or_frozen_official_alias")
                if normalized_label(hh["subseries"]) and normalized_label(hh["subseries"]) == normalized_label(ii["official_sample_label"]):
                    reasons.append("explicit_section_or_panel_label")
                    rules.append("source_and_section_identity")
                if len(page_acquisitions[source_url]) == 1 and equal_dims:
                    rules.append("single_acquisition_source_and_image_geometry")
            if equal_dims:
                reasons.append("exact_full_resolution_image_dimensions")
            if rr["sample_id"] == "FFMB028" and hh["id"] == "TENX159":
                assert hh["dataset_title"] == ii["source_title"] and equal_dims
                reasons.append("exact_provider_title_and_verified_image_content")
                rules.append("provider_title_and_verified_image_content")
            if roi_conflict:
                # Different ROIs are different acquisitions even in the same run.
                rules = []
                reasons.append("different_native_roi_uuid")
            if not reasons:
                continue
            decision = "confirmed_overlap" if rules else "not_same_acquisition_or_insufficient_identity"
            record = {
                "sample_id": rr["sample_id"], "external_resource": "HEST",
                "external_id": hh["id"], "decision": decision,
                "decision_rules": ";".join(rules), "evidence": ";".join(reasons),
                "source_url": rr["source_url"], "external_source_url": hh["download_page_link1"],
                "our_resource_id": rr["source_resource_id"], "our_acquisition_identity": rr["acquisition_identity"],
                "our_section_label": ii["official_sample_label"], "external_section_label": hh["subseries"],
                "our_roi_uuid": roi2, "external_roi_uuid": roi1 if isinstance(roi1, str) else "",
                "external_metadata_url": "https://huggingface.co/datasets/MahmoodLab/hest/blob/" + hest_revision + "/metadata/" + hh["id"] + ".json",
                "notes": "",
            }
            if rr["sample_id"] == "MBDC037" and hh["id"] == "TENX137":
                record["notes"] = "Formic-acid percentage label differs (8% vs 10%); native ROI and image dimensions agree."
            if rr["sample_id"] == "HD-043" and hh["id"] == "TENX155":
                record["notes"] = "Additional same-dimension CRC P2 image candidate; not needed to establish overlap, already confirmed via TENX128. No donor/identity inference from dimensions alone."
            candidates.append(record)
            if rules:
                matched.append(record)

    # Full-barcode coordinates and original provider-report fingerprints identify
    # all four slides. Two current Space Ranger reanalyses changed coordinates;
    # retain the measured differences rather than claiming byte-identical output.
    st_hd_ids = set(stimage.loc[stimage.tech == "VisiumHD", "slide"])
    assert st_hd_ids == {x["stimage_slide"] for x in st_verification}
    assert not stimage.tech.str.contains("Xenium", case=False).any()
    for ver in st_verification:
        assert ver["identity_confirmed"]
        assert ver["matched_barcodes"] == ver["stimage_rows"]
        assert ver["max_abs_coordinate_error_px"] < 25
        rr = release.loc[release.sample_id == ver["sample_id"]].iloc[0]
        row = {
            "sample_id": ver["sample_id"], "external_resource": "STimage",
            "external_id": ver["stimage_slide"], "decision": "confirmed_overlap",
            "decision_rules": ver["identity_method"],
            "evidence": ver["identity_evidence"],
            "source_url": rr.source_url,
            "external_source_url": "https://huggingface.co/datasets/jiawennnn/STimage-1K4M/blob/" + st_revision + "/VisiumHD/coord/" + ver["stimage_slide"] + "_coord.csv",
            "our_resource_id": rr.source_resource_id, "our_acquisition_identity": rr.acquisition_identity,
            "our_section_label": "", "external_section_label": "", "our_roi_uuid": "", "external_roi_uuid": "",
            "external_metadata_url": "https://github.com/JiawenChenn/STimage-1K4M/blob/56b6fc36d313e161a031fb86a15dd9e293fde9f2/meta/meta_all_gene02122025.csv",
            "notes": "Acquisition overlap, not identical segmentation, bin size, image modality, or prediction targets.",
        }
        matched.append(row)
        candidates.append(row)

    pairs = pd.DataFrame(matched).sort_values(["sample_id", "external_resource", "external_id"])
    assert not pairs.duplicated(["sample_id", "external_resource", "external_id"]).any()
    # An external acquisition record must not be assigned to two frozen acquisitions.
    assert pairs.groupby(["external_resource", "external_id"]).sample_id.nunique().max() == 1
    pairs.to_csv(out / "matched_record_pairs.csv", index=False)
    pd.DataFrame(candidates).sort_values(["sample_id", "external_resource", "external_id"]).to_csv(out / "candidate_decisions.csv", index=False)
    crosswalk = []
    for rr in release.to_dict("records"):
        ps = pairs[pairs.sample_id == rr["sample_id"]]
        hs = sorted(ps.loc[ps.external_resource == "HEST", "external_id"])
        ss = sorted(ps.loc[ps.external_resource == "STimage", "external_id"])
        category = "both" if hs and ss else "HEST_only" if hs else "STimage_only" if ss else "no_catalogue_match"
        crosswalk.append({k: rr[k] for k in ["sample_id", "platform", "species", "organ", "segmentation_kind", "stain", "gene_positive_observations", "source_resource_id", "source_url", "acquisition_identity"]})
        crosswalk[-1].update({"hest_ids": ";".join(hs), "stimage_ids": ";".join(ss), "overlap_group": category, "hest_status": "confirmed" if hs else "no_match_in_audited_catalogue", "stimage_status": "confirmed" if ss else "no_match_in_audited_catalogue", "release_version": "2026-10-28", "audit_date": "2026-11-02"})
    crosswalk = pd.DataFrame(crosswalk)
    crosswalk.to_csv(out / "acquisition_crosswalk.csv", index=False)
    hest_publication = crosswalk[["sample_id", "platform", "source_resource_id", "gene_positive_observations", "hest_ids", "hest_status"]]
    hest_publication.to_csv(out / "hest_publication_crosswalk.csv", index=False)
    directly_matched = hest_publication[hest_publication.hest_ids != ""]
    directly_unmatched = hest_publication[hest_publication.hest_ids == ""]
    pd.DataFrame([{
        "reference": "HEST 2024 publication",
        "catalogue_records": len(hest),
        "matched_HD": int((directly_matched.platform == "HD").sum()),
        "matched_Xenium": int((directly_matched.platform == "Xenium").sum()),
        "matched_acquisitions": len(directly_matched),
        "matched_our_observations": int(directly_matched.gene_positive_observations.sum()),
        "unmatched_acquisitions": len(directly_unmatched),
        "unmatched_our_observations": int(directly_unmatched.gene_positive_observations.sum()),
    }]).to_csv(out / "hest_publication_summary.csv", index=False)
    summary = []
    for category in ["HEST_only", "STimage_only", "both", "no_catalogue_match"]:
        group = crosswalk[crosswalk.overlap_group == category]
        summary.append({"overlap_group": category, "HD_acquisitions": int((group.platform == "HD").sum()), "Xenium_acquisitions": int((group.platform == "Xenium").sum()), "acquisitions": len(group), "gene_positive_observations_in_our_release": int(group.gene_positive_observations.sum())})
    pd.DataFrame(summary).to_csv(out / "overlap_summary.csv", index=False)
    input_names = ["release_samples.csv", "sample_inventory.csv", "sample_versions.csv", "HEST_v1_1_0.csv", "STimage_meta_all_gene02122025.csv", "hest_metadata_all.json"]
    versions = {
        "audit_date": "2026-11-02", "release_version": "2026-10-28",
        "unit": "molecular acquisition / capture area / Xenium ROI; not donor, organ, study or processed cell",
        "HEST": {
            "reference": "NeurIPS 2024 publication cohort, Figure 1a",
            "catalogue_version": "v1.1.0",
            "metadata_records": len(hest),
            "relevant_assay_records": len(relevant_hest),
            "metadata_url": "https://github.com/mahmoodlab/HEST/blob/56b242670f4bd27e30fd2438a2aafc414cae1460/assets/HEST_v1_1_0.csv",
            "paper_url": "https://proceedings.neurips.cc/paper_files/paper/2024/hash/60a899cc31f763be0bde781a75e04458-Abstract-Datasets_and_Benchmarks_Track.html",
            "native_evidence_revision": hest_revision,
        },
        "STimage": {"metadata_version": "meta_all_gene02122025.csv", "metadata_records": 1149, "github_revision": "56b6fc36d313e161a031fb86a15dd9e293fde9f2", "hub_revision": st_revision, "HD_records": len(st_hd_ids), "Xenium_records": 0, "coordinate_file_sha256": {f.name: sha256(f) for f in sorted((inp / "stimage_coords").glob("*.csv"))}},
        "input_sha256": {name: sha256(inp / name) for name in input_names},
        "negative_status": "No acquisition match found in the specified catalogues. Not a claim of first-ever data generation, donor independence, or absence from every external resource.",
        "observation_count_scope": "Counts are our frozen gene-positive cell/nucleus observations associated with shared-source acquisitions; they are not counts of identical cells or target vectors present in HEST/STimage.",
    }
    versions.update(json.loads((inp / "release_index_receipt.json").read_text()))
    (out / "source_versions.json").write_text(json.dumps(versions, indent=2) + "\n")
    no_match_ids = set(crosswalk.loc[crosswalk.overlap_group == "no_catalogue_match", "sample_id"])
    unresolved = [c for c in candidates if c["sample_id"] in no_match_ids and "different_native_roi_uuid" not in c["evidence"]]
    assert not unresolved, unresolved
    validation = {"acquisitions": len(crosswalk), "gene_positive_observations": int(crosswalk.gene_positive_observations.sum()), "hest_matched_acquisitions": int((crosswalk.hest_ids != "").sum()), "stimage_matched_acquisitions": int((crosswalk.stimage_ids != "").sum()), "matched_either_acquisitions": int((crosswalk.overlap_group != "no_catalogue_match").sum()), "hest_matching_catalogue_records": int((pairs.external_resource == "HEST").sum()), "unresolved_candidates_among_no_catalogue_match": len(unresolved), "stimage_coordinate_checks": st_verification, "extra_hest_image_check": image_verification, "summary": summary}
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps({k: v for k, v in validation.items() if k not in ["stimage_coordinate_checks", "extra_hest_image_check"]}, indent=2))


if __name__ == "__main__":
    main()
