"""Collect adopted MI tables without replacing the existing primary analyses."""
from pathlib import Path
from datetime import datetime
import csv, hashlib, html, shutil

ROOT = Path(__file__).resolve().parents[1]

def adopted(module):
    runs = list((ROOT/'outputs'/module/'current').iterdir())
    runs = [p for p in runs if p.is_dir()]
    if len(runs) != 1 or (runs[0]/'run_status.txt').read_text(encoding='utf-8-sig').splitlines()[0] != 'COMPLETED':
        raise RuntimeError(f'No completed adopted run: {module}')
    return runs[0]

def main():
    r12 = adopted('12_Adult_covariate_MI')
    r13 = adopted('13_Clinical_MI')
    with (r12/'audit/final_integrity_checks.csv').open(encoding='utf-8-sig', newline='') as f:
        if any(r['passed'].upper() != 'TRUE' for r in csv.DictReader(f)):
            raise RuntimeError('MI integrity failure')
    base = ROOT/'Results_表格与图片/Results_重点表格'
    folder = base/('adult_MI_'+datetime.now().strftime('%Y%m%d_%H%M%S'))
    folder.mkdir(parents=True, exist_ok=False)
    files = [
        (r13/'tables/pooled_adult_models.csv', '01_MI成人合并结果.csv'),
        (r13/'tables/complete_case_vs_MI.csv', '02_完整资料与MI对照.csv'),
        (r12/'tables/adult_missingness.csv', '03_成人原始缺失情况.csv'),
        (r12/'tables/imputation_summary.csv', '04_插补概况.csv'),
        (r12/'tables/observed_imputed_distribution.csv', '05_观测与插补分布.csv'),
        (r12/'tables/distribution_by_generation.csv', '06_分代观测与插补分布.csv'),
        (r12/'tables/complete_vs_missing_adult_characteristics.csv', '07_完整与缺失资料者特征.csv'),
        (r13/'reports/results_summary_zh.md', '结果摘要_zh.md'),
        (ROOT/'config/imputation_design_20261020.md', '方法与文献_zh.md'),
    ]
    rows = []
    for source, name in files:
        dest = folder/name
        shutil.copy2(source, dest)
        before = hashlib.md5(source.read_bytes()).hexdigest()
        after = hashlib.md5(dest.read_bytes()).hexdigest()
        if before != after:
            raise RuntimeError('Copy checksum mismatch')
        rows.append(dict(file=name, source=str(source), md5=after))
    with (folder/'file_manifest.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['file','source','md5']);writer.writeheader();writer.writerows(rows)
    links='\n'.join(f'<li><a href="{html.escape(r["file"])}">{html.escape(r["file"])}</a></li>' for r in rows)
    with (r13/'tables/pooled_adult_models.csv').open(encoding='utf-8-sig',newline='') as f:
        estimates=list(csv.DictReader(f))
    order={'G0-G1':0,'G1-G2':1,'G0-G2':2}
    labels={'R2_adult_core_same_sample':'成人核心参照','R3_adult_lifestyle':'＋生活方式','R4_adult_lifestyle_BMI':'＋生活方式及BMI'}
    estimates.sort(key=lambda r:(order[r['comparison']],r['model']))
    table=[]
    for r in estimates:
        interval=(f"{float(r['OR']):.2f}（{float(r['CI_low']):.2f}–{float(r['CI_high']):.2f}）" if r['status']=='POOLED' else '不可可靠合并')
        p=f"{float(r['p_raw']):.4f}" if r['status']=='POOLED' else '—'
        vals=[r['comparison'],labels[r['model']],r['n']+'/'+r['cases'],interval,p,r['m_estimable']+'/'+r['m_total']]
        table.append('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in vals)+'</tr>')
    intro='<p>50套家系多层MICE数据、40轮，最大Rhat=1.009。原始RDS和全龄主结果保留；仅插补成人缺失协变量。G0–G1恢复为319人、88例；涉及G2的成人比较恢复为224人、12例。</p>'
    interpretation='<p>原先6个不可估模型中，核心和生活方式层的4个现可合并。G1–G2的BMI层29/50套分离，不能只合并21套成功结果；G0–G2的BMI层50/50套分离。G1–G2成人生活方式OR约1.92，区间较宽；不能据此证明无关联，也不能确认独立于BMI。</p><p>核心调整：年轻年龄/性别、年长观察人数/平均年龄/女性比例、调查地区/城乡；恒定项按原规则移除。生活方式增加饮酒、吸烟、饮食评分、体力活动及log10能量；BMI层进一步增加BMI。全部采用原家系GEE与稳健方差，Rubin合并log OR，不平均OR。</p><p>全龄主模型与成人模型样本不同，OR变化不能全部归因于生活方式。MI依赖MAR及插补模型假设，不能增加病例或保证解决分离。</p>'
    header='<tr>'+''.join('<th>'+v+'</th>' for v in ['比较','模型','人数/病例','合并OR（95%CI）','原始P','可估套数'])+'</tr>'
    page='<!doctype html><meta charset="utf-8"><title>成人MI重点表格</title><style>body{font:15px/1.7 Arial,"Microsoft YaHei";max-width:1200px;margin:30px auto;padding:16px;color:#243746}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:8px}th{background:#eef2f7}</style><h1>成人协变量多重插补与重分析</h1>'+intro+'<table>'+header+''.join(table)+'</table>'+interpretation+'<h2>下载表格与方法</h2><ul>'+links+'</ul>'
    (folder/'index.html').write_text(page,encoding='utf-8')
    (r13/'reports/results_readable_zh.html').write_text(page.split('<h2>下载表格与方法</h2>')[0],encoding='utf-8')
    index=base/'index.html';backup=base/'index_before_MI_20261020.html'
    if not backup.exists():shutil.copy2(index,backup)
    index.write_text('<!doctype html><meta charset="utf-8"><title>脂肪肝重点表格</title><h1>脂肪肝重点表格</h1><ul><li><a href="index_before_MI_20261020.html">原主要结果及补充表格</a></li><li><a href="'+folder.name+'/index.html">新增：成人协变量多重插补与重分析</a></li></ul>',encoding='utf-8')
    from consolidate_key_table_folder import rebuild_index
    rebuild_index()
    print('MI key tables:', folder)
    print('Verified copies:',len(rows))

if __name__ == '__main__':main()
