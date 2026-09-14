"""Connect the unchanged drawing functions to package-local RET-01 data."""
from pathlib import Path
import pandas as pd
from s1_helpers import PROTEIN_ORDER, CONDITIONS, CATEGORIES

ROOT = Path(__file__).resolve().parents[2]


def retention_tables():
    frame = pd.read_csv(ROOT/'02_postprocessing/results/trajectory_summaries/retention_manifest.csv')
    conditions = pd.read_csv(ROOT/'02_postprocessing/results/aggregate_statistics/retention_conditions.csv')
    if len(frame) != 126 or len(conditions) != 42 or not conditions.trajectory_n.eq(3).all():
        raise ValueError('Invalid retention population')
    return frame, conditions


def s1_tables():
    frame, conditions = retention_tables()
    counts = pd.read_csv(ROOT/'02_postprocessing/results/aggregate_statistics/retention_qc_counts.csv').set_index('protein_short').reindex(PROTEIN_ORDER)
    counts.columns = CATEGORIES
    retained = conditions.pivot(index='protein_short', columns='condition', values='retained_ready_n').reindex(index=PROTEIN_ORDER, columns=CONDITIONS).astype(int)
    if not counts.sum(axis=1).eq(21).all():
        raise ValueError('Invalid QC category population')
    return frame, counts, retained
