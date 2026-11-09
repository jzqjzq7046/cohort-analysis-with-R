"""Inventory input filenames and ZIP directories without extracting private data."""
import json
import os
import re
import zipfile
from pathlib import Path

root = Path('E:/YFS')
project = root / 'workspace/FattyLiver_generation_analysis'
skip = {'.git', 'node_modules', '.r-library', '.r-library-sgb', '.r-library-4.5',
        '.figure_qa_python', '.python_packages', 'tools', '.venv', 'renv', 'code_history', '90_历史版本'}
raw = re.compile(r'\.(?:fastq|fq)(?:\.gz|\.bz2|\.zst)?$|\.(?:sam|bam|cram)(?:\.bz2|\.gz)?$', re.I)
strain = re.compile(r'(?:strainphlan|instrain|strain[_-](?:distance|sharing|profile)|sample.*\.markers|\.pkl$|\.markers\.(?:bz2|gz)$|snv[_-])', re.I)
matches, archives, errors = [], [], []
for folder, dirs, files in os.walk(root):
    dirs[:] = [d for d in dirs if d not in skip]
    for name in files:
        path = Path(folder) / name
        if raw.search(name) or strain.search(name):
            matches.append({'path': str(path), 'kind': 'raw_reads' if raw.search(name) else 'candidate_strain_input'})
        if path.suffix.lower() == '.zip' and not path.is_relative_to(project / 'Results_表格与图片'):
            try:
                with zipfile.ZipFile(path) as z:
                    names = z.namelist()
                    hits = [n for n in names if raw.search(n) or strain.search(n)]
                    archives.append({'path': str(path), 'entries': len(names), 'candidate_entries': hits})
            except (OSError, zipfile.BadZipFile) as exc:
                errors.append({'path': str(path), 'error': str(exc)})
report = {'scope': str(root), 'skipped_software_or_history_directories': sorted(skip),
          'filesystem_candidate_inputs': matches, 'zip_inventory': archives, 'errors': errors,
          'raw_or_marker_inputs_found': bool(matches or any(a['candidate_entries'] for a in archives))}
out = project / 'inventory/ADE_20261101'
out.mkdir(parents=True, exist_ok=True)
(out / 'sequence_input_inventory.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps(report, ensure_ascii=False, indent=2))
