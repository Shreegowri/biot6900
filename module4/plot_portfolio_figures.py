#!/usr/bin/env python3
"""Plot BIOT 6900 PDAC portfolio figures from existing Lab 2–4 results.
The file is kept in module 4 as it is the final module for part I.
No model fitting, network downloads or new target ranking. Paths supplied on
command line are resolved against --repo-root, never the working directory.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
import platform
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Patch
from matplotlib.collections import LineCollection

PALETTE={'RNA':'#2f6f97','Protein':'#d18a35','Genomic':'#63936b'}
AF_COLORS=['#ff7d45','#ffdb13','#65cbf3','#0053d6']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def intersections(sets):
    a,b,c=sets
    return {'100':len(a-b-c),'010':len(b-a-c),'110':len((a&b)-c),'001':len(c-a-b),'101':len((a&c)-b),'011':len((b&c)-a),'111':len(a&b&c)}

def venn(ax, sets, title):
    # Equal circles deliberately do not imply area proportionality.
    centers=[(-.52,.28),(.52,.28),(0,-.58)]
    for (name,col),xy in zip(PALETTE.items(),centers):ax.add_patch(Circle(xy,.94,facecolor=col,edgecolor=col,alpha=.21,lw=1.3))
    positions={'100':(-.95,.48),'010':(.95,.48),'110':(0,.66),'001':(0,-1.13),'101':(-.49,-.43),'011':(.49,-.43),'111':(0,-.02)}
    for k,n in intersections(sets).items():ax.text(*positions[k],f'{n:,}',ha='center',va='center',fontsize=11,fontweight='bold')
    for x,y,n in [(-1.0,1.37,'RNA'),(1.,1.37,'Protein'),(0.,-1.78,'Genomic')]:
        ix=['RNA','Protein','Genomic'].index(n);ax.text(x,y,f'{n}\n{len(sets[ix]):,}',ha='center',va='center',fontsize=10,color=PALETTE[n])
    ax.set(xlim=(-1.65,1.65),ylim=(-2.,1.6),aspect='equal');ax.axis('off');ax.set_title(title,loc='left',fontsize=12,fontweight='bold',pad=12)

def chain_plot(ax,ca,pocket,title):
    xyz=ca[['x','y','z']].to_numpy();xyz=xyz-xyz.mean(axis=0);_,_,v=np.linalg.svd(xyz,full_matrices=False);xy=xyz@v[:2].T
    colors=np.array(AF_COLORS)[np.digitize(ca.plddt.to_numpy(),[50,70,90])]
    seg=np.stack([xy[:-1],xy[1:]],axis=1);ax.add_collection(LineCollection(seg,colors=colors[:-1],linewidths=1.1,alpha=.95))
    keep=ca.residue.isin(pocket);ax.scatter(xy[keep,0],xy[keep,1],s=21,facecolors='none',edgecolors='#8b1697',linewidths=.8,zorder=4)
    ax.autoscale();ax.set_aspect('equal',adjustable='datalim');ax.axis('off');ax.margins(.12);ax.set_title(title,loc='left',fontsize=12,fontweight='bold')


def load_table(path, required):
    if not path.is_file():
        raise FileNotFoundError(str(path))
    df = pd.read_csv(path, sep='\t' if path.suffix.lower() == '.tsv' else ',')
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f'{path.name}: missing columns {sorted(missing)}')
    if 'gene' in df:
        if df.gene.isna().any():
            raise ValueError(f'{path.name}: missing gene identifiers')
        df['gene'] = df.gene.astype(str).str.strip().str.upper()
        if df.gene.eq('').any() or df.gene.duplicated().any():
            raise ValueError(f'{path.name}: empty or duplicate gene identifiers')
    return df


def numbers(df, names, filename):
    for name in names:
        converted = pd.to_numeric(df[name], errors='coerce')
        if (df[name].notna() & converted.isna()).any():
            raise ValueError(f'{filename}: nonnumeric {name}')
        df[name] = converted
    return df


def save_figure(fig, out, stem, outputs):
    for extension in ['png', 'pdf']:
        path = out / f'{stem}.{extension}'
        fig.savefig(path, dpi=240)
        outputs.append(path.name)
    plt.close(fig)


def scatter(ax, df, title='B  RNA and protein direction'):
    finite = np.isfinite(df.rna_effect) & np.isfinite(df.prot_effect)
    good = df.loc[finite]
    if len(good) < 2:
        raise ValueError('At least two finite RNA/protein effect pairs are needed')
    concordant = np.sign(good.rna_effect) == np.sign(good.prot_effect)
    ax.scatter(good.rna_effect, good.prot_effect,
               c=np.where(concordant, '#587f91', '#d48c66'), s=8, alpha=.45, rasterized=True)
    ax.axhline(0, color='#999999', lw=.7); ax.axvline(0, color='#999999', lw=.7)
    for gene, offset in [('ERBB4', (14, -13)), ('BRSK2', (-50, -15)),
                         ('KRAS', (12, -16)), ('COL11A1', (10, 9))]:
        row = good.loc[good.gene.eq(gene)]
        if row.empty:
            continue
        r = row.iloc[0]
        ax.scatter(r.rna_effect, r.prot_effect, s=31, color='#20252b')
        ax.annotate(gene, (r.rna_effect, r.prot_effect), xytext=offset,
                    textcoords='offset points', fontsize=9,
                    arrowprops={'arrowstyle': '-', 'color': '#777', 'lw': .6})
    rho = good.rna_effect.corr(good.prot_effect, method='spearman')
    ax.set(xlabel='RNA mean difference on log2 abundance scale',
           ylabel='Protein mean difference on normalized scale')
    ax.set_title(title, loc='left', fontsize=12, fontweight='bold')
    ax.text(.02, .98, f'n = {len(good):,} evaluable genes\nSpearman rho = {rho:.3f}\n'
            f'Same direction = {concordant.mean():.1%}', transform=ax.transAxes,
            va='top', fontsize=9, bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .85})
    return {'joined_genes': len(df), 'evaluable_pairs': len(good),
            'same_direction': int(concordant.sum()), 'opposite_direction': int((~concordant).sum()),
            'unassessable_pairs': int((~finite).sum()),
            'same_direction_fraction': float(concordant.mean()),
            'spearman_rho': float(rho) if np.isfinite(rho) else None}


def plot_heatmap(df, out, outputs):
    top = df.loc[np.isfinite(df.score)].sort_values('rank', kind='stable').head(15)
    if top.empty:
        raise ValueError('No finite ranked Lab 2 scores')
    if not np.isfinite(top[['rna_effect', 'prot_effect', 'association']]).all().all():
        raise ValueError('Top-ranked genes have missing or infinite heatmap values')
    if not top.association.between(0, 1).all():
        raise ValueError('association must be a mutation frequency between zero and one')
    fig, axs = plt.subplots(1, 3, figsize=(8.7, 4.5), sharey=True)
    for j, (col, title, cmap) in enumerate([
        ('rna_effect', 'RNA mean difference\nlog2 abundance', 'RdBu_r'),
        ('prot_effect', 'Protein mean difference\nnormalized abundance', 'RdBu_r'),
        ('association', 'Mutation frequency\npercent of observed cases', 'YlGnBu')]):
        values = top[col].to_numpy() * (100 if col == 'association' else 1)
        limit = max(float(np.max(np.abs(values))), 1e-12)
        im = axs[j].imshow(values[:, None], cmap=cmap, aspect='auto',
                           vmin=0 if col == 'association' else -limit, vmax=limit)
        axs[j].set_xticks([]); axs[j].set_title(title, fontsize=10, pad=10)
        axs[j].set_yticks(np.arange(len(top)), [f'{int(r)}  {g}' for r, g in zip(top['rank'], top.gene)])
        axs[j].tick_params(axis='y', length=0, labelsize=9)
        for i, value in enumerate(values):
            axs[j].text(0, i, f'{value:.2f}' if col == 'association' else f'{value:+.2f}',
                        ha='center', va='center', fontsize=9,
                        color='white' if abs(value) / limit > .65 else '#111111')
        for spine in axs[j].spines.values():
            spine.set_visible(False)
        fig.colorbar(im, ax=axs[j], orientation='horizontal', fraction=.045,
                     pad=.04, aspect=25).ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=.18, right=.99, top=.85, bottom=.08, wspace=.18)
    save_figure(fig, out, 'figure2_omics_heatmap', outputs)
    return top.gene.tolist()


def read_ca(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith('ENDMDL'):
            break
        if line.startswith('ATOM  ') and line[12:16].strip() == 'CA' and line[16] in (' ', 'A'):
            if line[26].strip():
                raise ValueError('Insertion codes are unsupported for canonical AlphaFold residue numbering')
            rows.append((line[21], int(line[22:26]), float(line[30:38]),
                         float(line[38:46]), float(line[46:54]), float(line[60:66])))
    ca = pd.DataFrame(rows, columns=['chain', 'residue', 'x', 'y', 'z', 'plddt'])
    if ca.empty or ca.chain.nunique() != 1 or ca.residue.duplicated().any():
        raise ValueError(f'{path.name}: expected one AlphaFold chain with unique C-alpha residues')
    if ca.residue.tolist() != list(range(1, len(ca) + 1)):
        raise ValueError(f'{path.name}: canonical residues must run consecutively from 1 for PAE mapping')
    if not np.isfinite(ca[['x','y','z','plddt']]).all().all() or not ca.plddt.between(0,100).all():
        raise ValueError(f'{path.name}: invalid coordinates or confidence values')
    return ca


def check_structure(gene, acc, pdb, pae_file, pocket_file):
    ca = read_ca(pdb)
    js = json.loads(pae_file.read_text())
    js = js[0] if isinstance(js, list) else js
    pae = np.asarray(js['predicted_aligned_error'], dtype=float)
    if pae.shape != (len(ca), len(ca)) or not np.isfinite(pae).all() or (pae < 0).any():
        raise ValueError(f'{gene}: PAE shape/values do not match structure')
    pk = load_table(pocket_file, ['pocket','Druggability Score','lining_residues','mean_pLDDT','PAE_within','verdict'])
    pk = numbers(pk, ['pocket','Druggability Score','mean_pLDDT','PAE_within'], pocket_file.name)
    if pk.empty or pk.pocket.duplicated().any() or not np.isfinite(pk[['pocket','Druggability Score','mean_pLDDT','PAE_within']]).all().all():
        raise ValueError(f'{gene}: invalid pocket table')
    linings = {}
    for _, row in pk.iterrows():
        residues = ast.literal_eval(row.lining_residues)
        if not isinstance(residues, list) or not residues or any(type(r) is not int for r in residues):
            raise ValueError(f'{gene}: lining_residues must be a nonempty list of integers')
        if len(set(residues)) != len(residues) or min(residues) < 1 or max(residues) > len(ca):
            raise ValueError(f'{gene}: invalid pocket residue mapping')
        ii = np.asarray(residues) - 1
        mean_conf = float(ca.iloc[ii].plddt.mean()); mean_pae = float(pae[np.ix_(ii, ii)].mean())
        if not np.isclose(mean_conf, row.mean_pLDDT, rtol=1e-6, atol=1e-6) or not np.isclose(mean_pae, row.PAE_within, rtol=1e-6, atol=1e-6):
            raise ValueError(f'{gene} pocket {row.pocket}: table confidence disagrees with PDB/PAE')
        verdict = 'reject' if mean_conf < 70 else ('suspect' if mean_pae > 10 else 'keep')
        if verdict != row.verdict:
            raise ValueError(f'{gene} pocket {row.pocket}: confidence verdict disagrees')
        linings[int(row.pocket)] = residues
    retained = pk.loc[pk.verdict.eq('keep')].sort_values(['Druggability Score','pocket'], ascending=[False,True])
    best = retained.iloc[0] if not retained.empty else None
    return ca, ([] if best is None else linings[int(best.pocket)]), {
        'uniprot': acc, 'residues': len(ca), 'mean_plddt': float(ca.plddt.mean()),
        'percent_below50': float(100 * ca.plddt.lt(50).mean()), 'pockets_checked': len(pk),
        'pockets_retained': len(retained), 'best_pocket': None if best is None else int(best.pocket),
        'best_pocket_score': None if best is None else float(best['Druggability Score']),
        'best_pocket_mean_plddt': None if best is None else float(best.mean_pLDDT),
        'best_pocket_mean_pae': None if best is None else float(best.PAE_within)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    here = Path(__file__).resolve().parent
    ap.add_argument('--repo-root', type=Path, default=here.parent if here.name == 'portfolio' else here,
                    help='biot6900 root; default parent of portfolio/, otherwise the script folder')
    ap.add_argument('--targets', default='../data/pdac/targets_pdac.csv', help='Lab 2 CSV')
    ap.add_argument('--week3', default='../data/pdac/targets_pdac_w3.csv', help='Lab 3 CSV')
    ap.add_argument('--rna', default='../data/pdac/pdac_transcriptomics.tsv')
    ap.add_argument('--protein', default='../data/pdac/pdac_proteomics.tsv')
    ap.add_argument('--genomics', default='../data/pdac/pdac_genomics.tsv')
    ap.add_argument('--structures-dir', default='structures')
    ap.add_argument('--pockets-dir', default='results')
    for gene, acc in [('erbb4','Q15303'),('brsk2','Q8IWQ3')]:
        ap.add_argument(f'--{gene}-pdb', help=f'Override AF-{acc}.pdb path')
        ap.add_argument(f'--{gene}-pae', help=f'Override AF-{acc}-pae.json path')
        ap.add_argument(f'--{gene}-pockets', help=f'Override {gene.upper()}_all_pockets.csv path')
    ap.add_argument('--output-dir', default='portfolio/figures')
    ap.add_argument('--figures', nargs='+', choices=['1','1b','2','3','4'], default=['1','2','3','4'],
                    help='Default all four; use 1b 2 3 if only the two target CSVs are available')
    ap.add_argument('--overwrite', action='store_true', help='Replace managed figure files and plotting_status.json')
    args = ap.parse_args(argv)
    root = args.repo_root.expanduser().resolve()
    def resolve(value):
        p = Path(value).expanduser()
        return p.resolve() if p.is_absolute() else (root / p).resolve()
    out = resolve(args.output_dir)
    paths = {k: resolve(getattr(args,k)) for k in ['targets','week3','rna','protein','genomics']}
    requested = list(dict.fromkeys(args.figures))
    stems = {'1':['figure1_integration_venn','figure1b_rna_protein'], '1b':['figure1b_rna_protein'],
             '2':['figure2_omics_heatmap'], '3':['figure3_druggability'], '4':['figure4_structures']}
    managed = {out / f'{s}.{ext}' for f in requested for s in stems[f] for ext in ['png','pdf']}
    managed.add(out / 'plotting_status.json')
    existing = sorted(p for p in managed if p.exists())
    if existing and not args.overwrite:
        ap.error(f'Output already exists: {existing[0]}. Choose another --output-dir or use --overwrite.')
    out.mkdir(parents=True, exist_ok=True)
    if args.overwrite:
        for p in existing:
            p.unlink()
    report = {'scope':'Plots from existing Lab 2–4 outputs; no new model training or structure prediction',
              'repo_root':str(root), 'requested_figures':requested, 'outputs':[], 'figures':{},
              'input_sha256':{}, 'environment':{'python':platform.python_version(),
              'numpy':np.__version__, 'pandas':pd.__version__, 'scipy':scipy.__version__, 'matplotlib':matplotlib.__version__}}
    def used(path):
        report['input_sha256'][str(path)] = sha(path)
    def lab2():
        used(paths['targets'])
        df = load_table(paths['targets'], ['gene','rna_effect','prot_effect','association','score','rank'])
        df = numbers(df, ['rna_effect','prot_effect','association','score','rank'], paths['targets'].name)
        if not df['rank'].gt(0).all() or not np.isfinite(df['rank']).all() or df['rank'].duplicated().any() or not (df['rank'] % 1 == 0).all():
            raise ValueError('Lab 2 rank must be a unique positive integer')
        return df
    for number in requested:
        try:
            if number in ['1','1b']:
                df = lab2()
                fig, ax = plt.subplots(figsize=(6.5,5.4))
                report['concordance'] = scatter(ax, df)
                fig.tight_layout()
                save_figure(fig, out, 'figure1b_rna_protein', report['outputs'])
                if number == '1':
                    layers = []
                    for key in ['rna','protein','genomics']:
                        t = load_table(paths[key], ['gene']); used(paths[key]); layers.append(t)
                    sets = [set(t.gene) for t in layers]
                    if set.intersection(*sets) != set(df.gene):
                        raise ValueError('Full-layer intersection does not match Lab 2 target genes; check source versions')
                    # Verify effects and mutation frequency against the unmerged layer tables.
                    for t, source_col, target_col in [(layers[0],'effect','rna_effect'),(layers[1],'effect','prot_effect'),(layers[2],'association','association')]:
                        if source_col not in t:
                            raise ValueError(f'Full layer lacks {source_col}')
                        a = pd.to_numeric(t.set_index('gene').loc[df.gene,source_col], errors='raise').to_numpy()
                        if not np.allclose(a,df[target_col].to_numpy(),equal_nan=True,rtol=1e-7,atol=1e-10):
                            raise ValueError(f'Full layer and target table disagree for {target_col}')
                    report['venn'] = {'definition':'Gene presence in complete source tables, not significance',
                                      'set_sizes':dict(zip(['RNA','Protein','Genomic'],map(len,sets))),
                                      'regions':intersections(sets)}
                    fig, axs = plt.subplots(1,2,figsize=(10.8,4.6),gridspec_kw={'width_ratios':[1,1.25]})
                    venn(axs[0],sets,'A  Cross-layer gene availability'); scatter(axs[1],df)
                    fig.tight_layout();save_figure(fig,out,'figure1_integration_venn',report['outputs'])
            elif number == '2':
                report['heatmap_genes'] = plot_heatmap(lab2(),out,report['outputs'])
            elif number == '3':
                w3 = load_table(paths['week3'], ['gene','p_sm','p_ab','rank']);used(paths['week3'])
                w3 = numbers(w3,['p_sm','p_ab','rank'],paths['week3'].name)
                if not np.isfinite(w3[['p_sm','p_ab','rank']]).all().all() or not w3.p_sm.between(0,1).all() or not w3.p_ab.between(0,1).all():
                    raise ValueError('Invalid Lab 3 scores/ranks')
                if w3['rank'].duplicated().any() or not w3['rank'].gt(0).all() or not (w3['rank'] % 1 == 0).all():
                    raise ValueError('Lab 3 rank must be a unique positive integer')
                if paths['targets'].is_file():
                    w2 = lab2().set_index('gene')
                    if not set(w3.gene).issubset(w2.index):raise ValueError('Lab 3 genes absent from Lab 2 table')
                    if 'week2_score' in w3 and not np.allclose(w3.week2_score,w2.loc[w3.gene,'score'],rtol=1e-8,atol=1e-10):
                        raise ValueError('Lab 3 week2_score disagrees with Lab 2 input')
                w3 = w3.sort_values('rank',kind='stable')
                fig,ax=plt.subplots(figsize=(8.6,4.9));y=np.arange(len(w3))
                ax.hlines(y,w3.p_sm,w3.p_ab,color='#bcc9ce',lw=1.2)
                ax.scatter(w3.p_sm,y,s=35,color='#207f92',label='Small-molecule precedence')
                ax.scatter(w3.p_ab,y,s=30,color='#c76535',marker='D',label='Antibody precedence')
                ax.set_yticks(y,w3.gene);ax.invert_yaxis()
                ax.set(xlabel='Uncalibrated model score',xlim=(-.02,max(.85,float(w3[['p_sm','p_ab']].max().max())+.04)))
                ax.legend(loc='lower right',frameon=False,fontsize=9);ax.grid(axis='x',alpha=.15)
                fig.tight_layout();save_figure(fig,out,'figure3_druggability',report['outputs'])
                report['druggability_order']=w3.gene.tolist()
            else:
                structures={}
                for gene,acc in [('ERBB4','Q15303'),('BRSK2','Q8IWQ3')]:
                    g=gene.lower();sd=resolve(args.structures_dir);pd_=resolve(args.pockets_dir)
                    pdb=resolve(getattr(args,g+'_pdb')) if getattr(args,g+'_pdb') else sd/f'AF-{acc}.pdb'
                    pae=resolve(getattr(args,g+'_pae')) if getattr(args,g+'_pae') else sd/f'AF-{acc}-pae.json'
                    pk=resolve(getattr(args,g+'_pockets')) if getattr(args,g+'_pockets') else pd_/f'{gene}_all_pockets.csv'
                    for p in [pdb,pae,pk]:
                        if not p.is_file():raise FileNotFoundError(str(p))
                        used(p)
                    structures[gene]=check_structure(gene,acc,pdb,pae,pk)
                report['structures']={g:v[2] for g,v in structures.items()}
                fig=plt.figure(figsize=(10.6,7.3));gs=fig.add_gridspec(2,2,height_ratios=[2.4,1.0],hspace=.29,wspace=.24)
                for j,(gene,(ca,lining,st)) in enumerate(structures.items()):
                    ax=fig.add_subplot(gs[0,j]);chain_plot(ax,ca,lining,f'{chr(65+j)}  {gene} predicted model')
                    label=f'Purple rings: pocket {st["best_pocket"]} lining C-alpha atoms' if st['best_pocket'] is not None else 'No pocket passes the confidence screen'
                    ax.text(.02,-.03,f'Mean pLDDT {st["mean_plddt"]:.1f}; {st["percent_below50"]:.1f}% below 50\n'+label,
                            transform=ax.transAxes,fontsize=9,va='bottom',bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
                    ax=fig.add_subplot(gs[1,j]);colors=np.array(AF_COLORS)[np.digitize(ca.plddt,[50,70,90])]
                    ax.scatter(ca.residue,ca.plddt,s=2.4,c=colors);ax.axhline(70,color='#777',ls='--',lw=.7)
                    ax.set(xlabel='Canonical residue number',ylabel='pLDDT',ylim=(0,103));ax.grid(axis='y',alpha=.12)
                handles=[Patch(facecolor=c,label=t) for c,t in zip(AF_COLORS,['<50','50–<70','70–<90','≥90'])]
                fig.legend(handles=handles,title='Local confidence pLDDT',loc='upper center',ncol=4,frameon=False,bbox_to_anchor=(.5,.985),fontsize=9)
                fig.subplots_adjust(top=.89,bottom=.09);save_figure(fig,out,'figure4_structures',report['outputs'])
            report['figures'][number]={'status':'generated'}
            print(f'Figure {number}: generated')
        except (FileNotFoundError,ValueError,KeyError,IndexError,SyntaxError,TypeError) as error:
            plt.close('all')
            state='missing_input' if isinstance(error,FileNotFoundError) else 'invalid_input'
            report['figures'][number]={'status':state,'reason':str(error)}
            print(f'Figure {number}: {state}: {error}',file=sys.stderr)
    incomplete=any(x['status']!='generated' for x in report['figures'].values())
    report['status']='incomplete' if incomplete else 'complete'
    report['outputs']=list(dict.fromkeys(report['outputs']))
    (out/'plotting_status.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(f'{report["status"].upper()}: {out}')
    if incomplete:print('Requested figures are missing. See plotting_status.json; no missing evidence was invented.',file=sys.stderr)
    return 2 if incomplete else 0


if __name__ == '__main__':
    sys.exit(main())
