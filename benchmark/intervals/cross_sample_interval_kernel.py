"""Fixed-prediction spatial reweighting; no training/generalization CI claim."""
from __future__ import annotations

import numpy as np
from scipy import sparse

METRICS = ('pearson_gene_macro', 'spearman_gene_macro', 'pearson_cell_macro',
           'f1_gene_macro', 'r2_gene_macro', 'rmse_log1p')


def finite_mean(values, axis=0):
    good = np.isfinite(values)
    count = good.sum(axis=axis)
    return np.divide(np.where(good, values, 0).sum(axis=axis), count,
                     out=np.full(np.shape(count), np.nan), where=count > 0)


def exact_spearman(y, p, block, weights, batch=64):
    """Midranks of literally repeated observations, including prediction ties."""
    reps, nb = weights.shape
    sums, counts = np.zeros(reps), np.zeros(reps, dtype=int)
    eligibility = np.zeros((reps, y.shape[1]), dtype=bool)
    for gene in range(y.shape[1]):
        _, yi = np.unique(y[:, gene], return_inverse=True)
        order = np.argsort(p[:, gene], kind='stable')
        ps = p[order, gene]
        starts = np.r_[0, np.flatnonzero(ps[1:] != ps[:-1]) + 1]
        if yi.max() == 0 or len(starts) == 1:
            continue
        yc = np.zeros((yi.max() + 1, nb), dtype=float)
        np.add.at(yc, (yi, block), 1)
        yorder, border = yi[order], block[order]
        tied = len(starts) != len(p)
        if tied:
            pg = np.cumsum(np.r_[True, ps[1:] != ps[:-1]]) - 1
        for lo in range(0, reps, batch):
            hi = min(reps, lo + batch)
            w = weights[lo:hi]
            ymass = w @ yc.T
            yr = (np.cumsum(ymass, axis=1) - .5 * ymass)[:, yorder]
            pw = w[:, border].astype(float)
            if tied:
                mass = np.add.reduceat(pw, starts, axis=1)
                pr = (np.cumsum(mass, axis=1) - .5 * mass)[:, pg]
            else:
                pr = np.cumsum(pw, axis=1) - .5 * pw
            center = .5 * pw.sum(axis=1, keepdims=True)
            yr -= center
            pr -= center
            cov = np.einsum('ij,ij,ij->i', pw, yr, pr)
            vy = np.einsum('ij,ij,ij->i', pw, yr, yr)
            vp = np.einsum('ij,ij,ij->i', pw, pr, pr)
            den = np.sqrt(vy * vp)
            good = den > 0
            sums[lo:hi] += np.divide(cov, den, out=np.zeros(hi-lo), where=good)
            counts[lo:hi] += good
            eligibility[lo:hi, gene] = good
    values = np.divide(sums, counts, out=np.full(reps, np.nan), where=counts > 0)
    return values, eligibility


def metrics(y, p, raw, block, weights, with_spearman=False):
    """First weight row must be all ones; return original point plus draws/QC.

    The exploratory macro functional matches the original metric (finite-gene
    mean). Changes in its contributing gene set are separately diagnosed and
    can suppress the published interval even if that macro remains finite.
    """
    y, p, weights = y.astype(float), p.astype(float), weights.astype(float)
    assert weights.ndim == 2 and np.all(weights[0] == 1)
    assert np.isfinite(y).all() and np.isfinite(p).all() and np.all(weights >= 0)
    nb, ng = weights.shape[1], y.shape[1]
    reducer = sparse.csr_matrix((np.ones(len(block)), (block, np.arange(len(block)))),
                                shape=(nb, len(block)))
    yc, pc = y-y.mean(axis=1, keepdims=True), p-p.mean(axis=1, keepdims=True)
    den = np.sqrt((yc*yc).sum(axis=1) * (pc*pc).sum(axis=1))
    rc = np.divide((yc*pc).sum(axis=1), den,
                   out=np.full(len(y), np.nan), where=den > 0)
    truth, called = raw > 0, p > np.log1p(.5)
    const_y, const_p = np.ptp(y, axis=0) == 0, np.ptp(p, axis=0) == 0
    ym, pm = y-y.mean(axis=0), p-p.mean(axis=0)
    ym[:, const_y], pm[:, const_p] = 0, 0
    features = np.concatenate((ym, ym*ym, pm, pm*pm, ym*pm, (y-p)**2,
        (truth & called).astype(float), (~truth & called).astype(float),
        (truth & ~called).astype(float), np.ones((len(y), 1)),
        np.where(np.isfinite(rc), rc, 0)[:, None], np.isfinite(rc)[:, None]), axis=1)
    block_stats = np.asarray(reducer @ features)
    stat = weights @ block_stats
    sy, sy2, sp, sp2, cross, sse, tp, fp, fn = [stat[:, k*ng:(k+1)*ng] for k in range(9)]
    n = stat[:, 9*ng:9*ng+1]
    safe_n = np.where(n > 0, n, 1)
    vy, vp = sy2-sy*sy/safe_n, sp2-sp*sp/safe_n
    for v, second, constant in [(vy, sy2, const_y), (vp, sp2, const_p)]:
        v[v <= 256*np.finfo(float).eps*np.maximum(second, np.finfo(float).tiny)] = 0
        v[:, constant] = 0
    gd = np.sqrt(vy*vp)
    corr = np.divide(cross-sy*sp/safe_n, gd,
                     out=np.full_like(cross, np.nan), where=gd > 0)
    corr = np.clip(corr, -1, 1)
    fd = 2*tp+fp+fn
    f1 = np.divide(2*tp, fd, out=np.full_like(tp, np.nan), where=fd > 0)
    r2 = 1-np.divide(sse, vy, out=np.full_like(sse, np.nan), where=vy > 0)
    cell = np.divide(stat[:, 9*ng+1], stat[:, 9*ng+2],
                     out=np.full(len(weights), np.nan), where=stat[:, 9*ng+2] > 0)
    rmse = np.sqrt(sse.sum(axis=1)/(safe_n[:, 0]*ng))
    out = np.column_stack((finite_mean(corr, 1), np.full(len(weights), np.nan),
                           cell, finite_mean(f1, 1), finite_mean(r2, 1), rmse))
    eligible = {0: np.isfinite(corr), 3: np.isfinite(f1), 4: np.isfinite(r2)}
    if with_spearman:
        out[:, 1], eligible[1] = exact_spearman(y, p, block, weights)
    out[n[:, 0] == 0] = np.nan
    occupancy = np.bincount(block, minlength=nb).astype(float)
    ess = float(occupancy.sum()**2/(occupancy*occupancy).sum())
    occupied_count = int(np.count_nonzero(occupancy))
    variance_mass = block_stats[:, ng:2*ng]
    variance_ess = np.divide(variance_mass.sum(axis=0)**2,
        (variance_mass*variance_mass).sum(axis=0),
        out=np.zeros(ng), where=(variance_mass*variance_mass).sum(axis=0) > 0)
    quality = []
    for j, metric in enumerate(METRICS):
        q = dict(metric=metric, occupied_blocks=occupied_count, effective_blocks=ess,
                 finite_fraction=float(np.isfinite(out[1:, j]).mean()) if len(out) > 1 else float(np.isfinite(out[0, j])),
                 eligibility_changed_fraction=0., status='ok', reasons=[])
        if j == 1 and not with_spearman:
            q.update(status='not_computed', reasons=['exact_spearman_not_requested'])
        elif not np.isfinite(out[0, j]):
            q.update(status='undefined', reasons=['original_metric_undefined'])
        else:
            if j in eligible:
                e = eligible[j]
                q['original_gene_eligibility'] = e[0].tolist()
                q['gene_finite_fraction'] = e[1:].mean(axis=0).tolist() if len(e) > 1 else e[0].astype(float).tolist()
                q['eligibility_changed_fraction'] = float(np.any(e[1:] != e[0], axis=1).mean()) if len(e) > 1 else 0.
            if j == 4:
                q['gene_variance_support_ess'] = variance_ess.tolist()
                q['min_gene_variance_support_ess'] = float(variance_ess[eligible[4][0]].min())
            if occupied_count < 10 or ess < 10:
                q['reasons'].append('fewer_than_10_spatial_blocks')
            if q['finite_fraction'] < .99:
                q['reasons'].append('fewer_than_99pct_finite_replicates')
            if q['eligibility_changed_fraction'] > .01:
                q['reasons'].append('gene_eligibility_changes_in_over_1pct_replicates')
            if j == 4 and q['min_gene_variance_support_ess'] < 10:
                q['reasons'].append('r2_variance_supported_by_fewer_than_10_effective_blocks')
            if q['reasons']:
                q['status'] = 'unsupported'
            elif occupied_count < 20 or ess < 20 or (j == 4 and q['min_gene_variance_support_ess'] < 20):
                q.update(status='limited', reasons=['fewer_than_20_effective_support_blocks'])
        quality.append(q)
    return out, quality
