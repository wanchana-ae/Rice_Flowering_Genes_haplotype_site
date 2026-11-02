import pandas as pd,numpy as np,matplotlib,html,sys; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.stdout.reconfigure(encoding='utf-8')
R='C:/Users/wanch/OneDrive/Work/01 Projects/Rice_Flowering_Genes/results/haplotype_2026-11-02/'
S=pd.read_csv(R+'hap232_samples.csv',index_col=0); ps=pd.read_csv(R+'hap232_per_sample.csv',index_col=0); sm=pd.read_csv(R+'hap232_summary.csv').set_index('gene')
ids={'OsLHY':'Os08g0157600','OsFKF1':'Os11g0547000','Ghd7':'Os07g0261200','Hd2_OsPRR37':'Os07g0695100','Ehd1':'Os10g0463400','Hd3a':'Os06g0157700','RFT1':'Os06g0157500','OsMADS14':'Os03g0752800','OsMADS15':'Os07g0108900','osFTL4':'Os09g0513500','OsATG8a':'Os07g0512200','OsATG8b':'Os04g0624000','OsATG8c':'Os08g0191600','OsATG8e':'Os11g0100100','OsCOL10':'Os03g0711100','OsZHD1':'Os09g0466400','DTH8':'Os08g0174500','Ehd2':'Os10g0419200','OsHAPL1':'Os05g0494100','FIE2':'Os08g0137100','SDG711':'Os06g0275500','SDG718':'Os03g0307800','OsDof12':'Os03g0169600','OsDof4':'Os01g0277500','bZIP71':'Os09g0306400','LHD3':'Os03g0752100','TAC1':'Os09g0529300','OsLF':'Os05g0541400'}
zone={'OsLHY':1,'OsFKF1':1,'OsCOL10':1,'OsZHD1':1,'Ghd7':2,'DTH8':2,'Hd2_OsPRR37':2,'Ehd1':3,'Ehd2':3,'OsHAPL1':3,'FIE2':3,'SDG711':3,'SDG718':3,'Hd3a':4,'RFT1':4,'OsMADS14':5,'OsMADS15':5,'OsDof12':5,'OsDof4':5,'bZIP71':5,'LHD3':5,'osFTL4':5,'TAC1':5,'OsLF':5,'OsATG8a':5,'OsATG8b':5,'OsATG8c':5,'OsATG8e':5}
cols={'PS':'#c0392b','NPS':'#2471a3'}; rng=np.random.default_rng(3); cards=[]; rows=[]
disp=lambda g:g.replace('_','/')
for g in sorted(ids,key=lambda x:(zone[x],x)):
    if g not in ps: continue
    h=ps[g].dropna(); vc=h.value_counts(); keep=vc[vc>=5].index.tolist()
    if len(keep)<1: continue
    hh=h[h.isin(keep)]
    t=S.loc[hh.index].assign(hap=hh.map({k:f'H{i+1}' for i,k in enumerate(keep)}))
    k=len(keep); fig,ax=plt.subplots(1,2,figsize=(max(7,1.1*k+5),4.6),sharey=True)
    for a,(hd,tt) in zip(ax,[('P','Hd1 = P (no deletion)'),('NP','Hd1 = NP (deletion)')]):
        for i in range(k):
            s=t[(t.hd1==hd)&(t.hap==f'H{i+1}')]
            if len(s): a.boxplot(s.dDTF,positions=[i],widths=.55,showfliers=False,medianprops=dict(color='k',lw=2))
            for gr in ['PS','NPS']:
                q=s[s.grp==gr]; a.scatter(i+rng.uniform(-.18,.18,len(q)),q.dDTF,c=cols[gr],s=22,edgecolor='k',lw=.3,zorder=3,label=gr if i==0 else None)
            a.text(i,150,f'n={len(s)}\n{int((s.grp=="PS").sum())}/{int((s.grp=="NPS").sum())}',ha='center',va='bottom',fontsize=7.5)
        a.axhline(34.5,ls='--',c='gray',lw=.8); a.set_xticks(range(k)); a.set_xticklabels([f'H{i+1}' for i in range(k)]); a.set_title(tt,fontsize=10); a.set_ylim(-85,185); a.grid(axis='y',alpha=.25); a.set_xlim(-.6,k-.4)
    ax[0].set_ylabel('ΔDTF (days)'); ax[0].legend(loc='lower left',fontsize=7,title='PS / NPS',title_fontsize=7)
    fig.suptitle(f'{disp(g)} ({ids[g]}) — ΔDTF by haplotype; labels: n and PS/NPS',fontsize=10)
    fig.savefig(f'site/img/{g}.png',dpi=105,bbox_inches='tight'); plt.close(fig)
    r=sm.loc[g]; rows.append((g,r))
    ht=t.groupby('hap').agg(n=('y','size'),PS=('y','sum'),HD=('hd1NP','sum'),med=('dDTF','median')).reindex([f'H{i+1}' for i in range(k)])
    hr=''.join(f'<tr><td>{i}</td><td>{int(x.n)}</td><td>{int(x.PS)}</td><td>{int(x.n-x.PS)}</td><td>{int(x.HD)}</td><td>{x.med:.0f}</td></tr>' for i,x in ht.iterrows())
    lk=''
    if g in ('RFT1','Hd3a'):
        other='Hd3a' if g=='RFT1' else 'RFT1'
        lk=f'<p class="note">† Linked with {other} (~12 kb apart): the all-reference haplotype is carried by exactly the same accessions in both genes, so this panel is the same signal as {other}, not an independent result.</p>'
    cards.append(f'<section class="card" id="{g}"><h3>{html.escape(disp(g))} <span class="m">{ids[g]} · zone {zone[g]} · {int(r.n_sites_used)} sites</span></h3>{lk}<img src="img/{g}.png" alt="{g} boxplot" loading="lazy"><table class="ht"><tr><th>hap</th><th>n</th><th>PS</th><th>NPS</th><th>Hd1 NP</th><th>median ΔDTF</th></tr>{hr}</table></section>')
def f(x): return '—' if pd.isna(x) else (f'{x:.2g}' if x<0.01 else f'{x:.3f}')
tr=''
for g,r in rows:
    flag=''
    if pd.notna(r.maxabs_coef) and r.maxabs_coef>20: flag=' <span class="w" title="large coefficient: a small haplotype is almost one group; logistic LRT unreliable">sep.</span>'
    cls=' class="hi"' if g in ('RFT1','Hd3a') else ''
    nt=int(r.n_in_test) if pd.notna(r.n_in_test) else '—'; nh=int(r.n_hap_ge5) if pd.notna(r.n_hap_ge5) else '—'
    dag=' †' if g in ('RFT1','Hd3a') else ''
    tr+=f'<tr{cls}><td><a href="#{g}">{html.escape(disp(g))}</a>{dag}</td><td>{ids[g]}</td><td>{zone[g]}</td><td>{int(r.n_sites_used)}</td><td>{nt}</td><td>{nh}</td><td>{f(r.logit_LRT_p)}{flag}</td><td>{f(r.logit_LRT_FDR)}</td><td>{f(r.OLS_dDTF_p)}</td><td>{f(r.OLS_dDTF_FDR)}</td></tr>'
SNPX=R+'SNP_haplotype_tables.xlsx'
def _tbl(df,hapcols=()):
    h='<tr>'+''.join(f'<th>{html.escape(str(c))}</th>' for c in df.columns)+'</tr>'
    b=''
    for _,x in df.iterrows():
        fx=' class="fx"' if str(x.get('impact',''))in('MODERATE','HIGH') else ''
        b+=f'<tr{fx}>'+''.join((f'<td class="{"r" if v=="REF" else "a"}">{v}</td>' if c in hapcols and isinstance(v,str) else f'<td>{"" if pd.isna(v) else html.escape(str(v))}</td>') for c,v in x.items())+'</tr>'
    return '<div class="sc"><table>'+h+b+'</table></div>'
def snp_section():
    out='<h2>SNP / indel haplotype tables</h2><div class="card"><p>Alleles at each site (gene ±2 kb, sites used in the analysis) for each haplotype with n≥5. <b>REF</b> = reference-genome (Nipponbare) allele (blue), <b>ALT</b> = alternative (red); bold rows change the protein (MODERATE/HIGH, snpEff). H1 = most frequent haplotype. Accessions heterozygous or missing at a used site are excluded. Download: <a href="data/SNP_haplotype_tables.xlsx">SNP_haplotype_tables.xlsx</a>.</p>'
    for g,title in [('RFT1','RFT1 (Os06g0157500)'),('Hd3a','Hd3a (Os06g0157700)')]:
        hs=pd.read_excel(SNPX,g+'_haplotypes'); sn=pd.read_excel(SNPX,g+'_SNPs')
        hc=[c for c in sn.columns if c.startswith('H') and c[1:].isdigit()]
        out+=f'<details><summary>{title} — {len(sn)} sites, {len(hs)} haplotypes</summary>'+_tbl(hs.rename(columns={'median_dDTF':'median ΔDTF'}))+_tbl(sn[['pos','REF','ALT','region','effect','impact','hgvsp','MAF_in_analysis_set']+hc],hc)+'</details>'
    out+='<details><summary>Hd1 deletions (chr6)</summary>'+_tbl(pd.read_excel(SNPX,'Hd1_sites'))+_tbl(pd.read_excel(SNPX,'Hd1_groups').rename(columns={'median_dDTF':'median ΔDTF'}))+'<p>Hd1 = NP if the accession carries deletion A or B; no accession carries both.</p></details>'
    return out+'<p class="note">Sites shown are those passing the analysis filters (missing ≤10%, heterozygous ≤5%, MAF ≥0.05 in the 221 accessions); association only, no causal claim.</p></div>'
css='''<style>:root{--bg:#F6F8F3;--surface:#fff;--border:#DCE3D4;--text:#17201A;--muted:#566257;--accent:#2F7A46;--hero1:#123A22;--hero2:#2F7A46;--code:#EAEFE4;--warn:#B45309}
@media (prefers-color-scheme:dark){:root{--bg:#0F1611;--surface:#172019;--border:#2A3A2E;--text:#E8EFE6;--muted:#9FB0A2;--accent:#5DBB7A;--hero1:#0B2214;--hero2:#1E5A30;--code:#223326;--warn:#F59E0B}}
*{box-sizing:border-box}body{margin:0;font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;background:var(--bg);color:var(--text);line-height:1.6}
header{background:linear-gradient(135deg,var(--hero1),var(--hero2));color:#fff;padding:36px 20px}header .in,main{max-width:1100px;margin:0 auto;padding:0 20px}header h1{margin:0 0 6px;font-size:1.7rem}header p{margin:0;opacity:.85}
main{padding-top:20px;padding-bottom:50px}.card{background:var(--surface);border:1px solid var(--border);border-radius:12px;padding:16px 18px;margin:16px 0}h2{margin:28px 0 8px}h3{margin:0 0 8px}.m{font-weight:400;font-size:.82rem;color:var(--muted);margin-left:8px}
img{max-width:100%;height:auto;display:block;margin:6px 0}table{border-collapse:collapse;width:100%;font-size:.86rem}th,td{border-bottom:1px solid var(--border);padding:5px 8px;text-align:left}th{color:var(--muted)}.ht{width:auto;margin-top:6px}
.sc{overflow-x:auto}tr.hi td{background:rgba(47,122,70,.12)}.w{color:var(--warn);font-size:.75rem;font-weight:700}.note{border-left:4px solid var(--warn);padding-left:12px}a{color:var(--accent)}details{margin:10px 0}summary{cursor:pointer;font-weight:600}td.r{background:rgba(36,113,163,.14)}td.a{background:rgba(192,57,43,.14)}tr.fx td:first-child{font-weight:700}.sc table td,.sc table th{white-space:nowrap}</style>'''
body=('<div class="card note"><b>Preliminary / not validated.</b> Exploratory results from one panel (221 accessions with ΔDTF-based groups), computed and checked only by the analysis team. Associations are not evidence of causation: the strongest signal (RFT1) is partly entangled with population structure, and there is no external validation yet.</div>'
'<h2>Method in brief</h2><div class="card"><ul>'
'<li><b>Samples:</b> accessions classified from the flowering-date sheet: photosensitive (PS, ΔDTF&gt;35 d, n=150) or non-photosensitive (NPS, ΔDTF&lt;34 d, n=71); 34–35 d and missing ΔDTF excluded; one accession with conflicting rows dropped, and two accessions with implausible negative ΔDTF (long-day flowering date earlier than the short-day date, which looks implausible; values awaiting confirmation by the data owner) excluded. Accession names/IDs are not shown.</li>'
'<li><b>Hd1 call:</b> NP if any non-reference allele at the three Hd1 deletion/indel sites (chr6:9338004, 9338005, 9338220), otherwise P.</li>'
'<li><b>Haplotype per gene:</b> string of homozygous genotypes at sites within the gene ±2 kb that pass missing≤10%, heterozygosity≤5%, MAF≥0.05; samples heterozygous or missing at a used site are left out for that gene; haplotypes with n≥5 are shown (H1 = most frequent). Genotypes are unphased.</li>'
'<li><b>Tests (per gene):</b> logistic (PS ~ PC1–5 + Hd1 + haplotype) and OLS (ΔDTF ~ PC1–5 + Hd1 + haplotype), likelihood-ratio / F-test for the haplotype term; Benjamini–Hochberg FDR across the 28 genes. "sep." marks genes with a very large coefficient (a small haplotype almost entirely one group), where the logistic test is unreliable.</li></ul></div>'
'<h2>Summary</h2><div class="card"><p>The <b>chr6 ~2.93 Mb locus (RFT1 / Hd3a)</b> is the clearest signal: RFT1 passes FDR&lt;0.05 in both tests, and <b>RFT1 and Hd3a cannot be separated</b> (see note below; Hd3a also passes both here, but it is the same signal). The all-reference haplotype of this locus is enriched in non-photosensitive accessions, also within Hd1 = P (5 of 14 NPS vs 2 of 125 for the main haplotype). Matching on ancestry was poor, so this stays a lead, not a conclusion. Other genes passing a test: <b>LHD3</b> passes both, but its logistic result shows separation (unreliable) and its ΔDTF result appeared only after the two implausible accessions were excluded (FDR 0.14 → 0.025), so it is unconfirmed; <b>Hd2/OsPRR37</b> passes only the logistic test (FDR 0.016); Ehd1, OsMADS15, OsATG8e, DTH8, Ehd2, SDG711, SDG718, OsDof4 and OsLF pass only the logistic test and most show separation ("sep."), so they are not supported. <b>The ΔDTF (OLS) test is sensitive to extreme values</b>: excluding the two accessions changed which genes pass it. In the earlier analysis restricted to photosensitive Hd1 = P accessions (n=144), no gene passed FDR.</p>'
'<p class="note"><b>† RFT1 and Hd3a are one signal here.</b> The two genes lie about 12 kb apart on chr6 (RFT1 2,926,823–2,928,474; Hd3a 2,940,004–2,942,452). Among the 155 accessions with a haplotype for both, the all-reference haplotype of RFT1 (n=34) and of Hd3a (n=34) are carried by exactly the same accessions, with no discordant accession. This data set therefore cannot say whether the association comes from RFT1, Hd3a, or the locus as a whole, and a separate "Hd3a" result in the table is not independent evidence.</p>'
'<div class="sc"><table><tr><th>Gene</th><th>RAP-DB</th><th>Zone</th><th>Sites</th><th>n tested</th><th>Hap n≥5</th><th>Logistic p</th><th>FDR</th><th>ΔDTF p</th><th>FDR</th></tr>'+tr+'</table></div></div>'
'<h2>Hd1 combined with the RFT1/Hd3a locus</h2><div class="card"><img src="img/hd1_rft1_6combos.png" alt="ΔDTF by combined Hd1 and RFT1/Hd3a haplotype" loading="lazy">'
'<p>Hd1 genotype (P = no deletion; deletion A = chr6:9338220 CAAGA&gt;C; deletion B = 9338004 CTT&gt;C / 9338005) combined with the RFT1/Hd3a haplotype (ref = all-reference, alt = main haplotype; the two genes have identical ref haplotypes, so they are shown as one). Accession counts (PS / NPS) and median ΔDTF:</p>'
'<div class="sc"><table><tr><th>Hd1</th><th>RFT1/Hd3a</th><th>n</th><th>PS</th><th>NPS</th><th>median ΔDTF</th></tr>'
'<tr><td>P</td><td>alt</td><td>125</td><td>123</td><td>2</td><td>99</td></tr><tr><td>P</td><td>ref</td><td>14</td><td>9</td><td>5</td><td>49</td></tr>'
'<tr><td>NP, deletion A</td><td>alt</td><td>24</td><td>3</td><td>21</td><td>25.5</td></tr><tr><td>NP, deletion A</td><td>ref</td><td>15</td><td>0</td><td>15</td><td>17</td></tr>'
'<tr><td>NP, deletion B</td><td>alt</td><td>6</td><td>2</td><td>4</td><td>32</td></tr><tr><td>NP, deletion B</td><td>ref</td><td>13</td><td>0</td><td>13</td><td>15</td></tr></table></div>'
'<p>Within both Hd1 = P (Mann–Whitney p = 0.0006; PC1–5-adjusted OLS p = 1.4e-09) and Hd1 = NP (p = 7.7e-06; adjusted p = 0.0013), the ref haplotype has lower ΔDTF than alt. The overall difference among the six groups (Kruskal–Wallis p = 6.7e-24) is mostly the known Hd1 P vs NP difference. p-values are per comparison and uncorrected for multiple testing.</p>'
'<p class="note">Association, not causation: ref carriers are mostly Hd1 NP (28 of 42, vs 30 of 155 for alt) and could not be matched well on ancestry; groups are small (n = 6–15 in four of the six); no external validation.</p></div>'
+snp_section()+
'<h2>Per-gene boxplots</h2><p style="color:var(--muted)">Each point is one accession; red = PS, blue = NPS; dashed line = 34.5 d (group cut-off). Labels above boxes: n and PS/NPS counts. Tables: PS and NPS counts, number with Hd1 NP, median ΔDTF.</p>'+''.join(cards))
page='<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Flowering-time haplotypes — Results</title>'+css+'</head><body><header><div class="in"><h1>Haplotype analysis of 28 flowering-pathway genes</h1><p>ΔDTF and photosensitivity by haplotype, controlling for Hd1 — preliminary results (draft, 2026-11-02)</p></div></header><main>'+body+'</main></body></html>'
open('site/index.html','w',encoding='utf-8').write(page)
print(len(cards),'cards')
