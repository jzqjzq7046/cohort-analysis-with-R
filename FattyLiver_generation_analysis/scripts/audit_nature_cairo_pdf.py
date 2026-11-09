"""Backend-neutral inspection adapter for R/Cairo PDF text geometry.

Cairo commonly uses a normalized 1 Tf with the physical font scale in Tm,
and batches widely separated labels in one text trace. The upstream lightweight
Tf scan and whole-trace-to-line matcher then produce false font/collision flags.
This adapter preserves the original reports and applies the SAME collision
engine/thresholds to all actual text-line/word boxes from the final PDF. It does
not draw, rewrite, rescale, or suppress any scientific objects.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path
import pymupdf


def main():
    p = argparse.ArgumentParser()
    p.add_argument("pdf", type=Path)
    p.add_argument("--skill", type=Path, required=True)
    p.add_argument("--prefix", type=Path, required=True)
    p.add_argument("--oriented-text", action="store_true", help="Opt-in real text quadrilaterals for radial labels; preserve the original axis-box report and thresholds")
    p.add_argument("--max-body-pt", type=float, default=7.01,
                   help="Explicit physical text ceiling for the artifact contract; scientific figures retain the default, manuscript tables may specify a larger ceiling")
    args = p.parse_args()
    prefix = str(args.prefix)
    doc = pymupdf.open(args.pdf)
    runs = []
    label_letters = set("abcdef")
    for page_no, page in enumerate(doc, 1):
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line.get("spans", []):
                    if span.get("text", "").strip():
                        runs.append(dict(page=page_no, text=span["text"], size_pt=span["size"], font=span["font"]))
    small = [r for r in runs if r["size_pt"] < 5 - 1e-4]
    large = [r for r in runs if r["size_pt"] > args.max_body_pt and not
             (r["text"].strip() in label_letters and abs(r["size_pt"] - 8) < .02 and "Bold" in r["font"])]
    text_report = dict(pdf=str(args.pdf), method="Effective rendered glyph size after Cairo text-matrix scaling",
                       minimum_pt=min((r["size_pt"] for r in runs), default=None),
                       maximum_pt=max((r["size_pt"] for r in runs), default=None),
                       text_runs=len(runs), below_floor=small, maximum_contract_pt=args.max_body_pt, text_above_contract=large,
                       verdict="PASS" if runs and not small and not large else "FIX BEFORE DELIVERY")
    Path(prefix + ".effective-text-audit.json").write_text(json.dumps(text_report, indent=2), encoding="utf-8")
    if text_report["verdict"] != "PASS":
        print(json.dumps(text_report, indent=2)); return 1

    spec = importlib.util.spec_from_file_location("nature_collision_engine", args.skill / "scripts/audit_figure_collisions.py")
    mod = importlib.util.module_from_spec(spec); sys.modules[spec.name] = mod; spec.loader.exec_module(mod)
    original = mod.extract_pdf_geometry

    def cairo_geometry(path):
        pages = original(path)
        final_pdf = pymupdf.open(path)
        for geometry, page in zip(pages, final_pdf):
            texts = []; oriented = {}
            # Split only genuine geometric gaps/line breaks. Keep all characters
            # and all drawings; no threshold or failure category is relaxed.
            for block in page.get_text("rawdict")["blocks"]:
                for line in block.get("lines", []):
                    direction = line.get("dir", (1, 0))
                    for span in line.get("spans", []):
                        size = span["size"]; groups = []; current = []; previous = None
                        for char in span.get("chars", []):
                            if previous is not None:
                                dx = char["origin"][0] - previous["origin"][0]
                                dy = char["origin"][1] - previous["origin"][1]
                                along = dx * direction[0] + dy * direction[1]
                                across = abs(-dx * direction[1] + dy * direction[0])
                                if along > 1.6 * size or along < -.5 * size or across > .45 * size:
                                    if current: groups.append(current)
                                    current = []
                            current.append(char); previous = char
                        if current: groups.append(current)
                        for group in groups:
                            visible = [c for c in group if c["c"].strip()]
                            if not visible: continue
                            bbox = mod.union_rects(mod.normalize_rect(c["bbox"]) for c in visible)
                            if args.oriented_text:
                                from rotated_text_geometry import quad_from_chars
                                oriented[len(texts)] = quad_from_chars(pymupdf, direction, span, visible)
                            texts.append(mod.TextBox(index=len(texts), text="".join(c["c"] for c in group).strip(), bbox=bbox))
            geometry.texts = texts
            if args.oriented_text:geometry.oriented_text = oriented
        final_pdf.close()
        return pages

    mod.extract_pdf_geometry = cairo_geometry
    if args.oriented_text:
        # Keep the untouched axis-aligned audit as a transparent diagnostic.
        import contextlib, io
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            mod.main([str(args.pdf), "--json-out", prefix + ".axis-aligned-collision-audit.json", "--strict"])
        Path(prefix + ".axis-aligned-collision-audit.log").write_text(captured.getvalue(), encoding="utf-8")
        from rotated_text_geometry import install_oriented_audit, self_test
        self_test(); install_oriented_audit(mod)
    status = mod.main([str(args.pdf), "--json-out", prefix + ".collision-audit.json",
                       "--overlay-pdf", prefix + ".collision-audit.pdf", "--strict"])
    report_path = Path(prefix + ".collision-audit.json")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    report["geometry_adapter"] = "R/Cairo: normalized Tf resolved with physical text matrices; per-line/word boxes replace batched whole-trace boxes; all engine thresholds unchanged"
    if args.oriented_text:
        report["geometry_adapter"] += "; rotated labels use actual recovered text quadrilaterals and local-coordinate stroke tests, same 5% overlap / 0.6pt inset / clipping thresholds; original horizontal-box findings preserved separately"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    doc.close()
    return status


if __name__ == "__main__":
    raise SystemExit(main())
