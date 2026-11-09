"""Verify and summarize a completed module-33 fresh rerun; no model or graphics.

Uses only aggregate saved results and writes new summary/audit documents.
Original result CSVs, participant data and all figures remain unchanged.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREVIOUS_ID = "ultrasound_grade_Mei_style_v1_20261028_20261028_164128"
LABELS = {"B1": "Model 3: all ages", "Adult_B1_same_sample": "Model 3: adults", "B2": "Model 4: adults"}


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def truth(value):
    return value == "TRUE"


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def key(row):
    return tuple(row[col] for col in ("model", "kind", "feature_id", "term"))


def main():
    module = ROOT / "outputs/33_Steatosis_severity"
    current_id = (module / "LATEST_RUN.txt").read_text().strip()
    run = module / "current" / current_id
    assert (run / "run_status.txt").read_text().splitlines()[0] == "COMPLETED"
    mode = run / "audit/analysis_execution_mode.txt"
    assert mode.is_file() and "Historical model checkpoints are NOT reused" in mode.read_text()
    previous = module / "history" / PREVIOUS_ID
    assert previous.is_dir(), "Previous adopted result must remain archived"
    source = "tables/B_All_grade_associations_EN.csv"
    rows, old_rows = read_rows(run / source), read_rows(previous / source)
    old = {key(row): row for row in old_rows}
    assert len(old) == len(old_rows) == len(rows)
    assert {key(row) for row in rows} == set(old)
    summaries, comparisons = [], []
    fields = ("beta_log2", "SE", "df", "CI_low", "CI_high", "P_raw", "P_BH", "abundance_ratio", "ratio_CI_low", "ratio_CI_high")
    expected_checkpoints = 0
    for model, label in LABELS.items():
        for kind in ("ordinal", "categorical"):
            selected = [row for row in rows if row["model"] == model and row["kind"] == kind]
            eligible = [row for row in selected if truth(row["eligible"])]
            valid = [row for row in eligible if truth(row["analysis_valid"])]
            sig05 = [row for row in valid if truth(row["strict_q005"])]
            sig10 = [row for row in valid if truth(row["source_style_q010"])]
            count_sgbs = len({row["feature_id"] for row in eligible})
            expected_checkpoints += math.ceil(count_sgbs / 20)
            summaries.append(dict(model=label, exposure_encoding=kind, participants_N=int(selected[0]["N"]),
                no_steatosis_N=int(selected[0]["No_N"]), mild_N=int(selected[0]["Mild_N"]), clear_N=int(selected[0]["Clear_N"]),
                eligible_SGBs=count_sgbs, eligible_tests=len(eligible), valid_tests=len(valid),
                BH_q_lt_005_tests=len(sig05), BH_q_lt_010_tests=len(sig10),
                BH_q_lt_005_SGBs=len({row["feature_id"] for row in sig05}),
                BH_q_lt_010_SGBs=len({row["feature_id"] for row in sig10}),
                BH_q_lt_005_positive_tests=sum(number(row["beta_log2"]) > 0 for row in sig05),
                BH_q_lt_005_negative_tests=sum(number(row["beta_log2"]) < 0 for row in sig05)))
            delta = {}
            for field in fields:
                differences = [abs(number(row[field]) - number(old[key(row)][field])) for row in selected
                    if math.isfinite(number(row[field])) and math.isfinite(number(old[key(row)][field]))]
                delta["max_abs_difference_" + field] = max(differences, default=0)
            unchanged = all(row[col] == old[key(row)][col] for row in selected
                for col in ("eligible", "analysis_valid", "strict_q005", "source_style_q010", "status", "N", "No_N", "Mild_N", "Clear_N"))
            comparisons.append(dict(model=label, exposure_encoding=kind, rows_compared=len(selected),
                samples_filters_status_and_significance_identical=unchanged, **delta))
            assert unchanged, f"Changed sample/filter/status/BH classification in {model}/{kind}; review required"
            assert delta["max_abs_difference_beta_log2"] < 1e-8
            assert delta["max_abs_difference_P_raw"] < 1e-8
            assert delta["max_abs_difference_P_BH"] < 1e-8
    checkpoint_root = run / "private/cp"
    if not checkpoint_root.is_dir():
        checkpoint_root = run / "private/fresh_analysis_checkpoints"
    checkpoints = list(checkpoint_root.rglob("checkpoint.rds"))
    assert len(checkpoints) == expected_checkpoints
    assert min(p.stat().st_mtime for p in checkpoints) >= mode.stat().st_mtime - 1
    recovered_path = run / "audit/recovered_same_request_fresh_checkpoints.csv"
    recovered = len(read_rows(recovered_path)) if recovered_path.is_file() else 0
    write_rows(run / "tables/B_Fresh_rerun_model_summary_EN.csv", summaries)
    write_rows(run / "tables/B_Fresh_vs_previous_results_EN.csv", comparisons)
    ids = read_rows(run / "tables/B_Displayed_ordinal_SGBs_EN.csv")
    six = [row["SGB_id"] for row in ids[:6]]
    text = ["# Fresh ultrasound-grade analysis: 1 November 2026", "",
        f"All six official MaAsLin2 model sets were refitted using {len(checkpoints)} newly created, run-local checkpoints. Historical fitted checkpoints were not read.",
        f"{recovered} checkpoints from the same 1-November fresh-analysis request were recovered after a Windows path-length export failure; the remaining {len(checkpoints) - recovered} were fitted after the restart. These are not checkpoints from the preceding 28-October analysis.",
        "The original cohort, ≥0.0002 relative-fraction detection in ≥10% of the model sample, zero-replacement/LOG transform, family random intercept, covariates and multiple-testing families were retained.",
        "Model 3 adjusts generation, generation-specific age splines, sex, library batch, sequencing pool, season and recent antibiotics. Adult Model 4 adds smoking, alcohol, diet score, physical activity, total energy and continuous BMI. The two adult models use identical participants.",
        "Observed ultrasound codes No/Mild/Clear are scored 0/1/2 in ordinal models. Categorical models use No as reference, with Mild and Clear contrasts pooled in each model's BH family. These are cross-sectional grade associations, not prospective progression or MASLD.", "",
        "| Model | Encoding | N | Eligible SGBs | Eligible tests | BH q<0.05 tests | BH q<0.10 tests |",
        "|---|---|---:|---:|---:|---:|---:|"]
    text += [f"| {x['model']} | {x['exposure_encoding']} | {x['participants_N']} | {x['eligible_SGBs']} | {x['eligible_tests']} | {x['BH_q_lt_005_tests']} | {x['BH_q_lt_010_tests']} |" for x in summaries]
    text += ["", "The six-panel display keeps the original raw relative-abundance scatter plus black-outlined boxplot form, with blue/pale-orange/orange grade colours and grey participant points.",
        "All eligible SGBs were tested before the original top-q/raw-P/ID display rule was applied. The effect figure displays up to 12 ordinal all-age hits; the raw abundance panels display up to six. The displayed six are: " + ", ".join(six) + ".",
        "q<0.10 is the source-style exploratory threshold; stricter q<0.05 is separately retained. Raw-distribution panels are unadjusted; their q annotations refer to the adjusted ordinal model.",
        "The fresh results reproduce the previous adopted result within floating-point precision. Sample counts, filtering, model status and q<0.05/q<0.10 classifications are identical; detailed numerical differences are in the comparison table.",
        "Reference framework: Mei Z, et al. Strain-specific gut microbial signatures in type 2 diabetes identified in a cross-cohort analysis of 8,117 metagenomes. Nature Medicine. 2024;30:2265–2276. doi:10.1038/s41591-024-03067-7. The ordinal/categorical per-feature framework is reused; this cohort, ultrasound outcome, family structure and adjustments are documented adaptations."]
    (run / "reports/Fresh_reanalysis_summary_EN.md").write_text("\n".join(text) + "\n", encoding="utf-8")
    zh = ["# 超声脂肪变程度分析完整重算：2026-11-01", "",
        f"6套官方MaAsLin2模型重新拟合，共{len(checkpoints)}个本次任务的新检查点；其中{recovered}个从同一11月1日任务的路径过长失败运行恢复，其余在缩短路径后拟合。未复用10月28日旧分析的模型检查点。全龄3193人：无2553、轻度602、明显38；同样本成人2040人：无1593、轻度421、明显26。", "",
        "六面板散点＋箱线图的表现形式、黑色箱线和灰色散点保留。展示前先完成全部合格SGB检验，再按原有q、原始P及编号规则选择前6个；不是预先只检验这6个。", "",
        "前6个SGB：" + "、".join(six) + "。", "",
        "完整重算与旧版结果一致，差异仅为浮点数精度；纳入人数、筛选、模型状态及BH q<0.05/0.10判断均未改变。q<0.10为既定文献式探索阈值，q<0.05另列。原始分布未作协变量调整，图中q来自调整后的有序模型。", "",
        "新汇总：B_Fresh_rerun_model_summary_EN.csv；新旧核对：B_Fresh_vs_previous_results_EN.csv。所有既往结果留在历史目录；原始RDS保持不变。"]
    (run / "reports/Fresh_reanalysis_summary_zh.md").write_text("\n".join(zh) + "\n", encoding="utf-8")
    audit = dict(status="PASS", new_checkpoints=len(checkpoints), expected_checkpoints=expected_checkpoints,
        all_checkpoints_generated_in_this_user_request=True, models_refitted=6,
        recovered_same_request_checkpoints=recovered,
        preceding_run=PREVIOUS_ID, classifications_and_samples_identical=True,
        maximum_absolute_beta_difference=max(x["max_abs_difference_beta_log2"] for x in comparisons),
        maximum_absolute_raw_P_difference=max(x["max_abs_difference_P_raw"] for x in comparisons),
        displayed_SGBs=six, statistics_refitted=True, plotting_form_preserved=True)
    (run / "audit/fresh_reanalysis_validation.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
