"""
Load a Sei diffs HDF5 file into a pandas DataFrame with rsIDs as the row index
and Sei target names as column labels.

Example:
    python view_diffs_h5.py \\
        /work_dir/results/<gwas>/sei_outputs_main_gwas/chromatin-profiles-hdf5/main_gwas_snps_correct_ref_final_diffs.h5 \\
        --row_labels /work_dir/results/<gwas>/sei_outputs_main_gwas/chromatin-profiles-hdf5/main_gwas_snps_correct_ref_final_row_labels.txt \\
        --target_names /work_dir/sei-framework/model/target.names \\
        --out diffs_with_labels.csv

If --row_labels / --target_names are omitted, the script tries the conventional
sibling path (*_row_labels.txt) and the default Sei target.names location.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def load_diffs(h5_path: Path) -> np.ndarray:
    with h5py.File(h5_path, "r") as f:
        return f["data"][:]


def load_row_labels(row_labels_path: Path) -> pd.Series:
    df = pd.read_csv(row_labels_path, sep="\t")
    return df["name"].astype(str)


def load_target_names(target_names_path: Path) -> pd.DataFrame:
    return pd.read_csv(
        target_names_path,
        sep="|",
        header=None,
        names=["Tissue", "Feature", "ProjectID", "Else"],
    ).apply(lambda c: c.str.strip() if c.dtype == "object" else c)


def build_dataframe(h5_path: Path, row_labels_path: Path, target_names_path: Path) -> pd.DataFrame:
    diffs = load_diffs(h5_path)
    rsids = load_row_labels(row_labels_path)
    targets = load_target_names(target_names_path)

    if diffs.shape[0] != len(rsids):
        raise ValueError(f"Row count mismatch: diffs={diffs.shape[0]}, row_labels={len(rsids)}")
    if diffs.shape[1] != len(targets):
        raise ValueError(f"Column count mismatch: diffs={diffs.shape[1]}, target_names={len(targets)}")

    columns = (
        targets["Tissue"].fillna("") + " | " +
        targets["Feature"].fillna("") + " | " +
        targets["ProjectID"].fillna("")
    )

    return pd.DataFrame(diffs, index=pd.Index(rsids, name="rsID"), columns=columns.values)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("h5", type=Path, help="Path to *_diffs.h5")
    p.add_argument("--row_labels", type=Path, default=None,
                   help="Path to *_row_labels.txt (default: sibling of h5)")
    p.add_argument("--target_names", type=Path, default=None,
                   help="Path to Sei target.names (default: target_names from paths.yaml).")
    p.add_argument("--paths_yaml", default=None,
                   help="Path to paths.yaml (default: <repo>/paths.yaml).")
    p.add_argument("--out", type=Path, default=None,
                   help="Optional output CSV path. If omitted, prints a preview only.")
    return p.parse_args()


def main():
    args = parse_args()

    row_labels = args.row_labels or Path(str(args.h5).replace("_diffs.h5", "_row_labels.txt"))

    target_names = args.target_names
    if target_names is None:
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from paths import load_paths
        target_names = Path(load_paths(args.paths_yaml)['target_names'])

    df = build_dataframe(args.h5, row_labels, target_names)
    print(f"Loaded DataFrame: {df.shape[0]} SNPs x {df.shape[1]} Sei profiles")
    print(df.iloc[:5, :3])

    if args.out:
        df.to_csv(args.out)
        print(f"Saved to {args.out}")


if __name__ == "__main__":
    main()
