# AGENTS Guide for CRISP

## 1. What This Repository Does
CRISP predicts post-perturbation single-cell gene expression from:
- a control-state cell representation (typically `scGPT` embedding),
- a drug condition (drug identity + dose),
- optional covariates.

Primary use cases in this repo:
- train perturbation response models on AnnData (`.h5ad`) datasets,
- evaluate IID and OOD perturbation generalization,
- run zero-shot predictions on unseen cell contexts,
- run downstream drug-screening style GSEA analysis on predictions.

## 2. High-Level Repository Organization
The repo is small and organized by workflow stage:
- Core library (`CRISP/`): data object construction, model, training loop, evaluation, embeddings utilities.
- Experiment configs (`experiments/configs/`): dataset/training key mappings and paths.
- Experiment launch scripts (`experiments/`): shell loops for multi-seed/multi-split runs.
- Data preprocessing assets (`data/`): notebooks documenting raw-to-AnnData pipeline and embedding generation.
- Tutorials (`tutorials/`): practical training, prediction, and drug-screening walkthroughs.

Important runtime entrypoints:
- `CRISP/train_script.py`: CLI training entrypoint.
- `CRISP/trainer.py`: orchestration class (`Trainer`) for dataset/model/train/eval/predict.
- `CRISP/data.py`: `Dataset`/`SubDataset` used by training/evaluation.
- `CRISP/model.py`: core model (`PertAE`) and loss composition.
- `CRISP/eval.py`: grouped perturbation evaluation metrics.

## 3. Environment and Run Basics
Install:
```bash
pip install -r requirement.txt
pip install -e .
```

Typical training command:
```bash
python CRISP/train_script.py \
  --config experiments/configs/nips.yaml \
  --split split \
  --seed 1327 \
  --savedir experiments/results/nips_split_1327 \
  --MMD 0.1
```

`--split` is the **column name** in `adata.obs` that stores `train/test/ood` labels (e.g. `split`, `split2`, `split3`).

## 4. Data Contract (What the Model Expects)
The training path expects an AnnData with standardized fields.

### 4.1 Required `adata` components
- `adata.X`: normalized/log1p gene expression matrix (cells x genes).
- `adata.obsm[FM_key]`: precomputed foundation-model cell embedding (default `X_scGPT`).
- `adata.uns[degs_key]`: dict of DE gene lists per perturbation group.

### 4.2 Required `adata.obs` fields (names are config-driven)
- perturbation name (e.g. `condition`)
- dose scalar (e.g. `dose_val`)
- cell type (e.g. `cell_type`)
- SMILES (e.g. `SMILES`)
- control indicator (`1` for control; e.g. `neg_control` or `control`)
- split label column (values: `train`, `test`, `ood`)
- perturbation category key for grouped eval (e.g. `cov_drug_name` or `cov_drug_dose_name`)
- paired-control covariate key (e.g. `type_donor` or `cell_type`)

### 4.3 Drug embedding table
A parquet dataframe indexed by canonical SMILES is expected for initializing drug embeddings.
- Training uses it to map training-set drugs -> fixed chemical vectors.
- Inference can also use it for unseen drugs via SMILES lookup.

## 5. End-to-End Data Pipeline (Raw Data -> Trainable Inputs)
Canonical flow documented in `data/nips_data.ipynb` and `data/sciplex_dataprocess.ipynb`:

1. Load raw count and metadata tables.
2. Build AnnData.
3. Standardize perturbation fields:
- normalize/clean drug names,
- derive normalized dose (`dose_val`),
- build composite keys such as `cov_drug_name`, `cov_drug_dose_name`, `drug_dose_name`,
- build control indicator (`neg_control` or `control`).
4. Filter sparse perturbation groups (example notebooks drop groups with <5 samples).
5. Compute DE genes per perturbation group using `rank_genes_groups_by_cov(...)` and store in `adata.uns`.
6. Canonicalize SMILES strings.
7. Precompute FM embeddings and store in `adata.obsm[...]` via `CRISP/scFM.py::calc_gpt`.
8. Build one or more split columns (`split`, `split2`, `split3`, etc.) with `train/test/ood` assignment.
9. Persist `.h5ad` and optional parquet drug-embedding table.

## 6. Training Task Pipeline (Exact Data Flow)
CLI entry (`CRISP/train_script.py`) performs:

1. Load YAML config.
2. Override config with CLI args (`split`, `seed`, `savedir`, optional paths/coefficients).
3. `Trainer.init_dataset(...)`:
- create `Dataset` object from AnnData,
- compute tensors: genes, FM embeddings, drug indices, doses, DE masks, group indices,
- compute paired-control FM embedding per sample (`get_paired_mean`),
- sample negative examples (`sample_neg`) for contrastive loss,
- create subsets: training/test control/test treated/(optional ood control/ood treated).
4. `Trainer.init_drug_embedding(...)`:
- load pretrained chemical embeddings (parquet),
- create frozen torch embedding aligned to training drugs.
5. `Trainer.init_model(...)`:
- construct `PertAE` with dataset dimensions and hyperparameters.
6. `Trainer.load_train()`:
- create DataLoader over training subset using custom collate.
7. `Trainer.train(...)` loop:
- for each batch, call `PertAE.iter_update(...)`,
- periodically run grouped evaluation on test/(ood) sets,
- save model and evaluation artifacts at stop conditions.

## 7. Model Consumption Pipeline (Inside `PertAE`)
Per batch, model consumes:
- target genes (`genes`) for supervised reconstruction,
- paired control FM embedding (`cell_embeddings`),
- drug index + dose (`drugs_idx`, `dosages`) or precomputed drug vector,
- optional covariates,
- DE mask + negative pair batch for auxiliary losses.

Forward logic (`predict`):
1. Encode control FM embedding -> latent Gaussian params (`mu`, `logvar`).
2. Reparameterize to sampled basal latent.
3. Encode drug embedding and scale it with learned dose-response (`dosers`).
4. Concatenate basal latent + drug latent (+ optional covariate latents).
5. Decode to predicted perturbed gene expression.

Training loss composition (`iter_update`):
- reconstruction = weighted MSE + DE-focused AFMSE,
- MMD between predicted and true expression distributions,
- cell-type classifier loss,
- cosine-embedding negative-pair loss,
- KL divergence regularization.

## 8. Evaluation Task Pipeline
`CRISP/eval.py::evaluate(...)`:
1. Iterate perturbation categories in treated set.
2. Skip tiny groups and control-like perturbations.
3. For each group, pick matching control cells (same cell type).
4. Generate predictions with fixed perturbation condition applied to control pool.
5. Compare mean predicted vs mean true profiles.
6. Aggregate metrics (global + per-group), including:
- `r2score`, `r2score_de`,
- `pearson`, `pearson_de`,
- `mse`, `mse_de`,
- `pearson_delta`, `pearson_delta_de`,
- `sinkhorn_de`.

## 9. Inference/Prediction Task Pipeline
`Trainer.get_prediction(...)` supports two perturbation encodings:

1. **Seen drug by name**:
- provide `drug_name`, `dose`, and `ref_drug_dict` (`drug_name -> training index`).

2. **SMILES-based (can support unseen drugs)**:
- provide `smile`, `dose`, and `smile_df` (SMILES-indexed embedding dataframe).

Common flow:
- input: control AnnData with FM embedding in `.obsm[FM_emb]`,
- build repeated perturbation tensors for all control cells,
- run model prediction,
- return predicted expression AnnData + latent outputs (or tensors).

## 10. Artifacts Produced by Training
In `save_dir`, training writes:
- `model.pt`: tuple of model state dict, init args, and history.
- `config.yaml`: resolved run config.
- `log.txt`: training/evaluation logs.
- `eval_stats.pkl`: aggregated eval metrics.
- `eval_stats_all.pkl`: per-group eval metrics.
- `pred_mean.pkl`: per-group mean true/pred/control profiles.

## 11. Practical Caveats / Implementation Notes
- Code is GPU-centric in multiple places (`custom_collate` and evaluation paths hardcode CUDA moves).
- `Dataset.subset(...)` uses set intersections; ordering is not explicitly stabilized.
- The optional generic covariate path in `SubDataset` is less exercised than default configs.
- No formal test suite is included; validation is mostly notebook/script driven.

## 12. Suggested Workflow for Future Agents
When making changes:
1. Start from config + `train_script.py` to understand active key mappings.
2. Verify AnnData schema against Section 4 before debugging model behavior.
3. Trace issues through `Trainer.init_dataset -> Dataset/SubDataset -> PertAE.predict/iter_update`.
4. Re-run a short training/eval smoke run on one split/seed before broader experiments.

