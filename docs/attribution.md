# Attribution, citation, and usage

This repository accompanies the sMMC-40M resource, developed from sMMC-22M. The
established repository URL and some data-side filenames retain earlier names.
The reviewed manuscript and discussion are on
[OpenReview](https://openreview.net/forum?id=Hk85ZGSGbQ). Use the final published
paper's author list and bibliographic record when available; no DOI, author list,
or final proceedings citation is invented here.

When reporting experiments, identify **data version 2026-10-28**, the protocol
track, input image encoder, target definition, and sample IDs. Historical Xenium
and current HD/Xenium extension scores have different protocols.

The image encoder is [CONCH](https://github.com/mahmoodlab/CONCH). Its official
checkpoint is obtained from [MahmoodLab/CONCH](https://huggingface.co/MahmoodLab/CONCH)
under the upstream access/usage terms. The immutable revision and file SHA256
are recorded in [`conch_v1_checkpoint.json`](../benchmark/models/conch_v1_checkpoint.json).
No checkpoint, tokenizer, or third-party CONCH implementation is redistributed
in this Git repository. The complete frozen data-side benchmark package retains
the encoder source snapshot used in the experiments.

Original acquisition providers, publications, and source URLs are retained in
the [catalog](../data/sample_catalog.csv). Upstream datasets, models, images,
software, and source publications retain their respective terms. This repository
does not assert a blanket license over those third-party assets. A separate
formal data-distribution and repository-code licensing statement has not yet
been supplied; no new license is assigned by this publication.

The overview figure was supplied by the resource author. If reusing figures,
data, or models, consult the corresponding authors/providers and applicable
terms. Please report reproducibility issues with the exact data version and
protocol rather than assuming that all acquisitions share one gene panel,
imaging branch, or donor-independent evaluation.
