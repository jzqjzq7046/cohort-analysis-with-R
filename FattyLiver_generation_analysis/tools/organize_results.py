"""Publish one study-ordered viewer from adopted outputs. Never fits models.

The first migration archives every previous viewing file in a verified ZIP.
Authoritative outputs, raw data, RDS, numerical results and figures are read-only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import re
import shutil
import sys
import tempfile
import time
import zipfile
from collections import Counter
from datetime import datetime, timezone, timedelta
from pathlib import Path
from urllib.parse import quote, unquote, urlsplit

GALLERY = "Results_表格与图片"
LOGS = "99_来源与整理记录"
ARCHIVE = "90_历史版本"
CURRENT = "01_当前分析"
SUPPORT = "02_数据质控与处理"
KEY = "Results_重点表格"
ALLOWED = {".csv", ".tsv", ".xlsx", ".xls", ".xlsm", ".ods", ".html", ".htm", ".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".svg", ".eps", ".webp", ".gif", ".bmp", ".md", ".txt", ".docx", ".rtf", ".pptx"}
PRIVATE = {"private", "audit", "history", "_running", "_checkpoints"}
STYLE = "body{font:15px Arial,'Microsoft YaHei',sans-serif;max-width:1150px;margin:32px auto;padding:0 24px;color:#243746;line-height:1.7}a{color:#23658e}h1{font-size:25px}h2{font-size:20px;margin-top:30px}.card{border:1px solid #dbe3e9;padding:15px 20px;margin:16px 0;border-radius:6px}small,.muted{color:#647784}table{width:100%;border-collapse:collapse}td,th{border-bottom:1px solid #e1e6ea;text-align:left;padding:8px;word-break:break-word}th{background:#f3f6f8}img{max-width:100%;max-height:510px;object-fit:contain}nav{display:flex;gap:20px;flex-wrap:wrap}input{padding:9px;width:95%;border:1px solid #aebdc6}details{margin:16px 0}code{word-break:break-all}"


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def within(path: Path, parent: Path) -> Path:
    resolved, base = path.resolve(), parent.resolve()
    if not resolved.is_relative_to(base) or resolved == base:
        raise RuntimeError(f"Unsafe target outside its intended parent: {path}")
    return resolved


def files_in(root: Path):
    if not root.exists():
        return
    for path in sorted(root.rglob("*")):
        if path.is_symlink() or (hasattr(path, "is_junction") and path.is_junction()):
            raise RuntimeError(f"Linked directories/files are not followed: {path}")
        if path.is_file():
            within(path, root)
            yield path


def write(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def write_csv(path: Path, rows: list[dict], fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def doc(title: str, body: str) -> str:
    return f'<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>{html.escape(title)}</title><style>{STYLE}</style></head><body>{body}</body></html>'


def link(destination: Path, page: Path, label: str) -> str:
    relative = os.path.relpath(destination, page.parent).replace("\\", "/")
    return f'<a href="{quote(relative, safe="/.:_-~")}">{html.escape(label)}</a>'


def category(relative: Path):
    if relative.parts[0] == "tables" or relative.suffix.lower() in {".csv", ".tsv", ".xlsx", ".xls", ".ods"}:
        return "Tables_表格"
    if relative.parts[0] == "figures" or relative.suffix.lower() in {".png", ".pdf", ".svg", ".tiff", ".tif"}:
        return "Figures_图片"
    return "Reports_报告"


def collect(root: Path, cfg: dict):
    records, modules = [], []
    for module in sorted((root / "outputs").iterdir()):
        if not module.is_dir() or not re.match(r"^\d+_", module.name):
            continue
        number = int(module.name.split("_", 1)[0])
        selected = module.name in cfg["current"] or module.name in cfg["support"] or number > 30
        if not selected:
            continue
        runs = list((module / "current").iterdir()) if (module / "current").exists() else []
        runs = [p for p in runs if p.is_dir()]
        if not runs:
            continue
        if len(runs) != 1:
            raise RuntimeError(f"Expected exactly one adopted run: {module}")
        run = runs[0]
        if not (run / "run_status.txt").exists() or (run / "run_status.txt").read_text(encoding="utf-8-sig").splitlines()[0] != "COMPLETED":
            raise RuntimeError(f"Adopted run is not complete: {run}")
        if module.name in cfg["support"]:
            folder = Path(SUPPORT) / cfg["support"][module.name]
            title = cfg["support"][module.name]
            note = "质控及预处理来源记录；不是新增统计分析。"
        else:
            meta = cfg["current"].get(module.name, {"folder": f"{number:02d}_{module.name}", "title": module.name, "note": "后续新增的成功分析模块。"})
            folder, title, note = Path(CURRENT) / meta["folder"], meta["title"], meta["note"]
        modules.append({"module": module.name, "run_id": run.name, "folder": folder.as_posix(), "title": title, "note": note})
        for source in files_in(run):
            relative = source.relative_to(run)
            if PRIVATE.intersection(relative.parts) or source.suffix.lower() not in ALLOWED or source.name in {"run_status.txt", "STATUS", "CURRENT", "CURRENT.txt"} or source.name.startswith("~$"):
                continue
            # The polished clinical Table 1/2 live only in the key-table folder.
            if relative.parts[:2] == ("tables", "Manuscript_Table1_Table2_20261024"):
                viewing = Path(KEY) / "tables/English/Table1_Table2" / Path(*relative.parts[2:])
            else:
                viewing = folder / relative
            records.append({"module": module.name, "run_id": run.name, "category": category(relative), "file": source.name, "viewing_path": viewing.as_posix(), "source_path": source.as_posix(), "source_md5": md5(source), "copied_md5": "", "bytes": source.stat().st_size})
    paths = [r["viewing_path"].casefold() for r in records]
    if len(paths) != len(set(paths)):
        raise RuntimeError("Colliding public result names")
    if len({r["source_path"] for r in records}) != len(records):
        raise RuntimeError("A source file would be copied more than once")
    return records, modules


def key_table(record):
    name = record["file"]
    if name.endswith(".xlsx") or name.startswith(("Main_Table", "Table1_", "Table2_", "Table3_", "Table4_", "Table5_", "A_Genus_model_summary", "A_Adjusted_genus_Shannon", "A_Genus_community_tests", "B_Global_and_taxon_summary", "B_BH_supported", "C_External_health_set_associations", "C_External_set_coverage", "Robustness_model_summary", "Effect_comparison_summary", "Taxon_concordance_summary", "Family_pair_sample_sizes")):
        return name.endswith((".csv", ".xlsx", ".html", ".pdf"))
    return "Clinical_association" in name and name.endswith((".csv", ".html", ".pdf"))


def make_pages(stage: Path, records: list[dict], modules: list[dict]):
    root_page = stage / "index.html"
    current_page = stage / CURRENT / "index.html"
    support_page = stage / SUPPORT / "index.html"
    key_page = stage / KEY / "index.html"
    top = f'<h1>脂肪肝研究 · 当前结果</h1><p>按研究顺序浏览；只纳入已成功完成的当前方案。旧版结果单独归档。</p><nav>{link(current_page,root_page,"当前分析")}{link(key_page,root_page,"重点表格")}{link(support_page,root_page,"数据质控与处理")}{link(stage/ARCHIVE/"index.html",root_page,"历史版本")}{link(stage/"file_manifest.csv",root_page,"来源清单")}</nav>'
    cards = []
    for module in modules:
        folder = stage / module["folder"]
        page = folder / "index.html"
        subset = [r for r in records if r["module"] == module["module"]]
        body = f'<nav>{link(root_page,page,"返回总目录")}{link(key_page,page,"重点表格")}</nav><h1>{html.escape(module["title"])}</h1><p>{html.escape(module["note"])}</p><small>来源：{html.escape(module["module"])} / {html.escape(module["run_id"])}</small>'
        # Only PNG previews appear in the reading flow; all export formats stay together.
        previews = [r for r in subset if r["category"] == "Figures_图片" and r["file"].lower().endswith(".png")]
        for record in previews:
            image = stage / record["viewing_path"]
            src = quote(os.path.relpath(image, page.parent).replace("\\", "/"), safe="/._-~")
            siblings = [r for r in subset if Path(r["viewing_path"]).with_suffix("") == Path(record["viewing_path"]).with_suffix("")]
            downloads = " · ".join(link(stage/r["viewing_path"],page,Path(r["file"]).suffix[1:].upper()) for r in siblings)
            body += f'<div class="card"><b>{html.escape(Path(record["file"]).stem)}</b><p>{downloads}</p><img loading="lazy" src="{src}" alt="{html.escape(record["file"])}"></div>'
        body += '<input type="search" placeholder="搜索本分析的表格、图注或文件" oninput="document.querySelectorAll(\'tr.result\').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(this.value.toLowerCase()))">'
        for cat in ("Tables_表格", "Reports_报告", "Figures_图片"):
            items = [r for r in subset if r["category"] == cat]
            body += f'<details {"open" if cat != "Figures_图片" else ""}><summary>{cat}（{len(items)}）</summary><table>'
            for record in items:
                body += f'<tr class="result"><td>{link(stage/record["viewing_path"],page,record["file"])}</td></tr>'
            body += '</table></details>'
        write(page, doc(module["title"], body))
        card = f'<div class="card"><h2>{link(page,root_page,module["title"])}</h2><p>{html.escape(module["note"])}</p><small>{len(subset)}个文件 · {html.escape(module["module"])}</small></div>'
        if module["folder"].startswith(CURRENT):
            cards.append(card)
    top += '<p><b>阅读顺序：</b>基础队列与临床 → 个人脂肪肝相关菌 → 单菌跨代相似性 → 未患病年轻成员 → 家系内患病比较 → 外部健康菌集合。</p>' + "".join(cards)
    top += '<p class="muted">超声肝脏脂肪变；不进行MASLD分类。仅要求相关两代有可用成员。不同分析样本互有重叠，人数不可相加。这里是查看副本；原始结果保留在outputs。</p>'
    write(root_page, doc("脂肪肝研究 · 当前结果", top))
    for prefix, page, title in ((CURRENT,current_page,"当前分析"),(SUPPORT,support_page,"数据质控与处理")):
        body = f'<nav>{link(root_page,page,"返回总目录")}</nav><h1>{title}</h1>'
        for module in modules:
            if module["folder"].startswith(prefix):
                body += f'<div class="card">{link(stage/module["folder"]/"index.html",page,module["title"])}<p>{html.escape(module["note"])}</p></div>'
        write(page, doc(title, body))
    body = f'<nav>{link(root_page,key_page,"返回总目录")}</nav><h1>重点表格 · Current English tables</h1><p>Table 1/2采用当前两代可用方案，临床Table 1为5,660人。其余表格链接至当前分析中的唯一查看副本，避免多份文件更新不同步。</p><table><tr><th>分析</th><th>Table / workbook</th></tr>'
    selected = [r for r in records if key_table(r) and not r["module"] in {m["module"] for m in modules if m["folder"].startswith(SUPPORT)}]
    for record in selected:
        body += f'<tr><td>{html.escape(record["module"])}</td><td>{link(stage/record["viewing_path"],key_page,record["file"])}</td></tr>'
    write(key_page, doc("重点表格", body+'</table>'))
    english = stage/KEY/"tables/English/index.html"
    write(english, doc("Current English tables",f'<h1>Current English tables</h1>{link(key_page,english,"Open current key tables")}'))
    write_csv(stage/KEY/"tables/English/Current_key_table_locations_EN.csv",[{k:r[k] for k in ("module","file","viewing_path","source_path","source_md5")} for r in selected], fields=["module","file","viewing_path","source_path","source_md5"])
    write(stage/KEY/"CURRENT.txt","Current two-generation protocol; Table 1 clinical N=5660.\nOpen index.html. Historical three-generation-only workbooks are in 90_历史版本.\n")


def preserve_previous(target: Path, stage: Path, stamp: str):
    marker = target/LOGS/"managed_files.json"
    if marker.exists():
        managed = set(json.loads(marker.read_text(encoding="utf-8"))["paths"])
        kept = []
        for source in files_in(target):
            relative = source.relative_to(target)
            if source.name.startswith("~$"):
                continue
            if relative.parts[0] == ARCHIVE or relative.as_posix() not in managed:
                destination = stage / relative
                if destination.exists():
                    if md5(source) == md5(destination):
                        continue
                    destination = stage/"98_用户补充文件"/stamp/relative
                destination.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(source,destination)
                if md5(source) != md5(destination):
                    raise RuntimeError(f"Preserved user/archive copy differs: {relative}")
                kept.append(relative.as_posix())
        return {"mode":"refresh","preserved_files":len(kept)}
    if not target.exists():
        return {"mode":"new","archived_files":0}
    archive = stage/ARCHIVE/f"整理前完整目录_{stamp}.zip"
    archive.parent.mkdir(parents=True,exist_ok=True)
    inventory = []
    with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as bundle:
        for source in files_in(target):
            relative=source.relative_to(target).as_posix()
            if source.name.startswith("~$"):
                continue
            inventory.append({"previous_path":relative,"md5":md5(source),"bytes":source.stat().st_size,"archive":archive.name})
            bundle.write(source,relative)
    with zipfile.ZipFile(archive) as bundle:
        if bundle.testzip() is not None:
            raise RuntimeError("Historical ZIP CRC check failed")
        for record in inventory:
            if hashlib.md5(bundle.read(record["previous_path"])).hexdigest() != record["md5"]:
                raise RuntimeError(f"Archive hash mismatch: {record['previous_path']}")
    write_csv(stage/ARCHIVE/f"整理前文件清单_{stamp}.csv",inventory)
    return {"mode":"migration","archived_files":len(inventory),"archive_md5":md5(archive),"archive":archive.name}


def validate_links(stage: Path, pages: list[Path]):
    checked=0
    for page in pages:
        text=page.read_text(encoding="utf-8")
        for raw in re.findall(r'(?:href|src)="([^"]+)"',text):
            url=urlsplit(html.unescape(raw))
            if url.scheme or not url.path:
                continue
            destination=(page.parent/unquote(url.path)).resolve()
            if not destination.is_relative_to(stage.resolve()) or not destination.is_file():
                raise RuntimeError(f"Broken generated link: {page} -> {raw}")
            checked+=1
    return checked


def publish_in_place(target: Path, stage: Path, root: Path, preservation: dict):
    """Windows Explorer may hold the root directory name open.

    Publish the verified archive first, then payloads, retire only preserved old
    viewing files, and publish the index/manifest last. Never touches outputs.
    """
    within(target,root)
    previous={p.relative_to(target).as_posix():md5(p) for p in files_in(target) if not p.name.startswith("~$")}
    staged={p.relative_to(stage).as_posix():p for p in files_in(stage)}
    archived={}
    if preservation["mode"]=="migration":
        inventories=list((stage/ARCHIVE).glob("整理前文件清单_*.csv"))
        for inventory in inventories:
            with inventory.open(encoding="utf-8-sig",newline="") as stream:
                archived.update({r["previous_path"]:r["md5"] for r in csv.DictReader(stream)})
    def copy(relative):
        source=staged[relative]
        destination=target/relative
        within(destination,target)
        destination.parent.mkdir(parents=True,exist_ok=True)
        if destination.exists() and md5(destination)==md5(source):
            return
        # Different temporary basename avoids overwriting an active Office lock.
        temporary=destination.with_name(destination.name+".result_publish_tmp")
        shutil.copy2(source,temporary)
        if md5(temporary)!=md5(source):
            raise RuntimeError(f"In-place copy differs: {relative}")
        os.replace(temporary,destination)
    archive_paths=[p for p in staged if p.startswith(ARCHIVE+"/")]
    for relative in archive_paths:
        copy(relative)
    controls={"index.html","file_manifest.csv","README_zh.md","README_en.md",".managed_results_gallery",f"{LOGS}/managed_files.json",f"{LOGS}/validation.json"}
    for relative in staged:
        if relative not in archive_paths and relative not in controls:
            copy(relative)
    retained=[]
    for relative,old_hash in previous.items():
        if relative in staged or relative.startswith(ARCHIVE+"/"):
            continue
        destination=target/relative
        within(destination,target)
        if preservation["mode"]=="migration" and archived.get(relative)!=old_hash:
            retained.append(relative)
            continue
        try:
            if md5(destination)!=old_hash:
                retained.append(relative)
                continue
            destination.unlink()
        except PermissionError:
            retained.append(relative)
    # Only empty old directories are removed; no recursive delete in the viewer.
    for directory in sorted([p for p in target.rglob("*") if p.is_dir()],key=lambda p:len(p.parts),reverse=True):
        within(directory,target)
        try:
            directory.rmdir()
        except OSError:
            pass
    validation=json.loads((stage/LOGS/"validation.json").read_text(encoding="utf-8"))
    validation.update({"publication":"in_place_root_directory_open","retained_open_or_changed_files":retained})
    write(stage/LOGS/"validation.json",json.dumps(validation,ensure_ascii=False,indent=2))
    for relative in sorted(controls-{"index.html"}):
        copy(relative)
    copy("index.html")
    for relative,source in staged.items():
        if md5(target/relative)!=md5(source):
            raise RuntimeError(f"Published file differs: {relative}")
    return validation


def build(root: Path, check=False):
    root=root.resolve()
    cfg=json.loads((root/"config/results_viewer.json").read_text(encoding="utf-8"))
    records,modules=collect(root,cfg)
    target=root/GALLERY
    if check:
        old=json.loads((target/LOGS/"validation.json").read_text(encoding="utf-8"))
        for record in records:
            if md5(target/record["viewing_path"]) != record["source_md5"]:
                raise RuntimeError(f"Viewing copy differs: {record['viewing_path']}")
        pages=[target/"index.html",target/CURRENT/"index.html",target/SUPPORT/"index.html",target/KEY/"index.html",target/KEY/"tables/English/index.html",target/ARCHIVE/"index.html"]+[target/m["folder"]/"index.html" for m in modules]
        links=validate_links(target,pages)
        print(json.dumps({"validation":"PASS","source_and_copy_hashes":len(records),"modules":len(modules),"generated_local_links":links,"previous_build":old},ensure_ascii=False,indent=2))
        return
    if target.exists() and not (target/".managed_results_gallery").exists():
        raise RuntimeError("Refusing to replace an unmanaged result folder")
    stamp=datetime.now(timezone(timedelta(hours=8))).strftime("%Y%m%d_%H%M%S")+f"_{os.getpid()}"
    stage=Path(tempfile.mkdtemp(prefix=".result_gallery_build_",dir=root))
    backup=root/f".result_gallery_previous_{stamp}"
    try:
        for record in records:
            source=Path(record["source_path"])
            destination=stage/record["viewing_path"]
            destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(source,destination)
            record["copied_md5"]=md5(destination)
            if record["copied_md5"] != record["source_md5"]:
                raise RuntimeError(f"Copy checksum mismatch: {source}")
        preservation=preserve_previous(target,stage,stamp)
        make_pages(stage,records,modules)
        archive_page=stage/ARCHIVE/"index.html"
        archive_files=[p for p in (stage/ARCHIVE).iterdir() if p.suffix.lower() in {".zip",".csv"}]
        body=f'<nav>{link(stage/"index.html",archive_page,"返回当前结果")}</nav><h1>历史版本</h1><p>整理前的所有查看文件均保存在完整ZIP中，并逐文件核对MD5。包括旧三代限制、旧调整集、旧绘图版本及原来重复的查看副本。这里的数值不应与当前方案混用。原始outputs未移动。</p><ul>'
        body+=''.join(f'<li>{link(p,archive_page,p.name)}</li>' for p in sorted(archive_files))+'</ul>'
        write(archive_page,doc("历史版本",body))
        zh="""# 当前结果查看目录

双击 index.html。按研究顺序阅读，PNG用于预览，PDF / TIFF / SVG保留在同一分析的figures文件夹。

- 01_当前分析：24–30号当前方案，每个分析只保留一个查看位置；文件名和数值不变。
- 02_数据质控与处理：协变量资料核查、菌群QC和预处理记录。
- Results_重点表格：当前英文重点表格；Table 1/2单列，其余链接到唯一查看副本。
- 90_历史版本：整理前全部查看文件的完整ZIP及逐文件MD5清单。旧三代限制和旧模型与当前结果分开。
- 99_来源与整理记录：本次迁移、校验、当前来源模块和托管文件清单。

原始RDS、outputs权威结果、模型和数值均未修改。本次只整理文件，不运行统计模型。
不汇入private、audit、逐人数据和RDS。不要在托管副本上长期保存修改；用户补充文件会保留。
当前临床Table 1为5,660人；各菌群模块按自己的实际纳入样本解释。24号AST和25号LOG结果有不同方法，不能混成一套结果。
更新入口仍是 scripts/refresh_results_gallery.R。各后续脚本自动使用本结构，避免再创建重复专题目录。
"""
        write(stage/"README_zh.md",zh)
        write(stage/"README_en.md","# Current study results\n\nOpen index.html. Results are ordered by the research workflow, with one viewing location per source file.\nCurrent analyses (24–30), QC/preprocessing, current key tables, historical ZIP archives and provenance are separated.\nOriginal RDS, authoritative outputs, statistics and figures are unchanged. No statistical model is executed.\nThe complete previous viewer is retained in a checksum-verified ZIP. Current clinical Table 1: N=5,660. Sample sizes differ across analyses.\nUse scripts/refresh_results_gallery.R to refresh existing results only.\n")
        write(stage/".managed_results_gallery","schema_version=2; study-ordered viewing copies; authoritative analyses remain in outputs\n")
        write_csv(stage/"file_manifest.csv",records)
        write_csv(stage/LOGS/"Current_modules_EN.csv",modules)
        write(stage/LOGS/"migration.json",json.dumps(preservation,ensure_ascii=False,indent=2))
        generated=[stage/"index.html",stage/CURRENT/"index.html",stage/SUPPORT/"index.html",stage/KEY/"index.html",stage/KEY/"tables/English/index.html",archive_page]+[stage/m["folder"]/"index.html" for m in modules]
        checked=validate_links(stage,generated)
        # Re-read every source after copying, so publication cannot hide source changes.
        if any(md5(Path(r["source_path"])) != r["source_md5"] for r in records):
            raise RuntimeError("A source result changed during organization")
        validation={"status":"PASS","built_at":datetime.now(timezone(timedelta(hours=8))).isoformat(),"files":len(records),"modules":len(modules),"generated_local_links":checked,"current_modules":sum(m["folder"].startswith(CURRENT) for m in modules),"preservation":preservation,"source_results_unchanged":True,"no_statistics_rerun":True}
        write(stage/LOGS/"validation.json",json.dumps(validation,ensure_ascii=False,indent=2))
        owned=[p.relative_to(stage).as_posix() for p in files_in(stage) if not p.relative_to(stage).parts[0] in {ARCHIVE,"98_用户补充文件"}]
        owned.append(f"{LOGS}/managed_files.json")
        write(stage/LOGS/"managed_files.json",json.dumps({"schema_version":2,"paths":sorted(set(owned))},ensure_ascii=False,indent=2))
        in_place=False
        if target.exists():
            within(target,root);within(backup,root)
            try:
                target.rename(backup)
            except PermissionError:
                validation=publish_in_place(target,stage,root,preservation)
                in_place=True
        if not in_place:
            try:
                stage.rename(target)
            except Exception:
                if backup.exists():
                    backup.rename(target)
                raise
        if backup.exists():
            # This is only the verified, archived viewing directory, never outputs/data.
            within(backup,root)
            if not backup.name.startswith(".result_gallery_previous_"):
                raise RuntimeError("Unexpected cleanup target")
            try:
                shutil.rmtree(backup)
            except PermissionError:
                print(f"Verified previous viewer remains recoverable at {backup}")
        print(json.dumps(validation,ensure_ascii=False,indent=2))
        print(f"Result viewer updated: {target/'index.html'}")
    finally:
        if stage.exists():
            within(stage,root)
            if stage.name.startswith(".result_gallery_build_"):
                shutil.rmtree(stage)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,required=True)
    parser.add_argument("--check",action="store_true")
    options=parser.parse_args()
    build(options.root,options.check)
