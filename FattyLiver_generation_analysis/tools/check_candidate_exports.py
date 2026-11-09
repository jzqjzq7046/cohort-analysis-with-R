"""Read-only figure metadata and export QA; never redraws scientific content."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / '.figure_qa_python'))
import pymupdf as fitz
from PIL import Image


def check(run: Path) -> dict:
    records = []
    figures = sorted((run / "figures").glob("*.svg"))
    assert len(figures) == 3
    for svg in figures:
        root = ET.parse(svg).getroot()
        width = float(root.attrib["width"].removesuffix("pt"))
        height = float(root.attrib["height"].removesuffix("pt"))
        assert any(node.tag.endswith("}text") for node in root.iter()), svg
        stem = svg.stem
        with fitz.open(svg.with_suffix(".pdf")) as doc:
            assert len(doc) == 1 and doc[0].get_text().strip()
            assert abs(doc[0].rect.width - width) < 1
            assert abs(doc[0].rect.height - height) < 1
        for extension in (".png", ".tiff"):
            path = svg.with_suffix(extension)
            with Image.open(path) as im:
                dpi = im.info.get("dpi")
                assert dpi and all(abs(float(value) - 600) < 0.1 for value in dpi)
                assert abs(im.width - width / 72 * 600) < 2
                assert abs(im.height - height / 72 * 600) < 2
                records.append(dict(figure=stem, extension=extension, pixels=im.size,
                                    dpi=tuple(float(value) for value in dpi)))
        audit_path = run / "audit" / f"{stem}.collision-audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        assert audit.get("verdict") == "PASS", audit_path
    hashes = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted((run / "figures").iterdir()) if path.is_file()
    }
    assert len(hashes) == 12
    return dict(status="PASS", figures=3, exports=12, raster_metadata=records,
                editable_vector_text=True, strict_collision_audit="PASS",
                sha256=hashes, original_images_modified=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    result = check(args.run.resolve())
    (args.run / "audit" / "final_export_verification.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("status", "figures", "exports")}))
