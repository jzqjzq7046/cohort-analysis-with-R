"""Read-only SAS ingestion; selected clinical columns only, no plotting."""
import argparse
import json
from pathlib import Path


def export(source: Path, out: Path, fields: list[str]) -> None:
    from pandas.io.sas.sas7bdat import SAS7BDATReader
    if source.resolve() == out.resolve():
        raise ValueError("Source and output must differ")
    if out.exists():
        raise FileExistsError(out)
    with SAS7BDATReader(source, encoding="infer") as reader:
        expected = (reader.row_count, reader.column_count)
        encoding = reader.encoding
        available = {c.lower(): c for c in reader.column_names}
        if "tutknom18" not in available:
            raise ValueError("Clinical participant ID is absent")
        df = reader.read()
    if df.shape != expected:
        raise ValueError(f"Incomplete read: {df.shape} versus header {expected}")
    selected = [available[f.lower()] for f in fields if f.lower() in available]
    out.parent.mkdir(parents=True, exist_ok=True)
    df.loc[:, selected].to_csv(out, index=False, encoding="utf-8-sig", date_format="%Y-%m-%d")
    receipt = dict(source=str(source), rows=len(df), source_columns=len(df.columns), encoding=encoding,
      exported_columns=selected, unavailable_fields=[f for f in fields if f.lower() not in available])
    out.with_suffix(".metadata.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(f"Clinical input read: {len(df)} records; {len(selected)} selected columns exported")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--fields", required=True)
    args = parser.parse_args()
    export(args.source, args.out, args.fields.split(","))
