# Main overview figure

[`assets/dataset.png`](../assets/dataset.png) is the author's supplied **sMMC-40M
Structure** figure, preserved byte-for-byte at its original 5232 × 2686 resolution.
It describes the intended connection between morphology, molecular measurements,
and context. The following details distinguish the illustration from the
**2026-10-28 frozen implementation**:

- The figure's scale panel retains source label variants such as Ovarian/Ovary,
  Colon/Colorectal, and Kidney/kidney. The current release uses 23 harmonized
  reporting categories. Its species panel displays human/mouse only, whereas the
  release includes seven source species. Use the catalog and exact summary
  rather than these illustrative plot labels.
- The technology percentages displayed in the illustration are not the frozen
  release proportions. The current gene-positive counts are HD **13,660,756**
  and Xenium **24,811,436**. HD includes **2,041,153 nuclei**.
- The diagram shows the earlier CellViT/HoVer-Net outline. Current HD H&E
  segmentation uses **VitaminP** and 2 µm full-bin-area aggregation; HD IF is
  nuclear. Xenium uses native cell measurements and boundaries.
- The diagram's 10 µm and 200 µm contexts are conceptual scales. The reported
  CONCH benchmark uses **150 native-pixel** crops followed by a 448-pixel encoder
  transform. Physical field of view therefore varies with image pixel size.
- Generated descriptions, communication graphs, and illustrated metadata fields
  describe the framework. This frozen release does not claim complete generated
  captions, validated cell-cell communication predictions, or verified disease
  states for every cell. See the actual field schema and annotation evidence.

Exact release values are in [`release_summary.json`](../data/release_summary.json).
The figure has not been regenerated or edited as part of this code publication.
