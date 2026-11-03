Source-acquisition overlap with HEST-1k and STimage-1K4M

The 2026-10-28 freeze contains 124 acquisitions and 38,472,192 gene-positive
cell or nucleus records. We compare it with the 1,229 profiles reported in
Figure 1a of the 2024 HEST-1k paper and the corrected 1,149-entry STimage
catalogue. The HEST authors' September 2024 manifest has 602 Visium, 552
Spatial Transcriptomics, 10 Visium HD and 65 Xenium profiles.

HEST shares 42 acquisitions (6 HD and 36 Xenium), accounting for 13,883,605
frozen records. The other 82 acquisitions account for 24,588,587 records.
STimage shares four HD acquisitions. Three occur in both catalogues and
HD-028 is the STimage-only match. The union covers 43 of 124 acquisitions.

                         HD  Xenium  Total   Our gene-positive records
HEST only                 3      36     39                  13,317,654
STimage only              1       0      1                     483,460
Both                      3       0      3                     565,951
No catalogue match       46      35     81                  24,105,127

The unit is a source acquisition, capture area or Xenium ROI. Native sample
and ROI identifiers, official page aliases, section and panel labels, and
image geometry identify matching records. Original image content corroborates
the FFMB028/TENX159 match. The four STimage matches are HD-011, HD-022,
HD-028 and HD-043. Barcode coordinates and provider reports identify the
source acquisitions across processing versions. Counts of cells and nuclei
refer to the processed records in our freeze.

acquisition_crosswalk.csv and overlap_summary.csv contain the joint audit.
hest_publication_crosswalk.csv and hest_publication_summary.csv contain the
direct comparison with the HEST paper. matched_record_pairs.csv and
candidate_decisions.csv give the supporting identifiers and decisions.
source_versions.json pins the reference manifests and input checksums.
validation.json and independent_validation.json record the checks.

Run python build_overlap_audit.py --root AUDIT_WORKSPACE to regenerate the
tables under outputs/. Supply the versioned inputs listed in source_versions.json
under inputs/, the aggregated native metadata for the 75 eligible HEST profiles,
hest_native_revision.txt, the STimage Hub revision and coordinate inputs,
release_index_receipt.json, and the two verified evidence files in the workspace
root (stimage_coordinate_verification.json and ffmb028_image_verification.json).
