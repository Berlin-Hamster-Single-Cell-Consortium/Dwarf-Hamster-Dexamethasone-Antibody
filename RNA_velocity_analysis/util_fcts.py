import matplotlib.pyplot as pl
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import scvelo as scv
from scipy.sparse import issparse
import matplotlib.gridspec as gridspec
from scipy.stats import gaussian_kde, spearmanr, pearsonr
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
import seaborn as sns
import re
import os
import gzip
import mygene
from csv import Sniffer

signatures_path_= os.path.join(os.path.dirname(os.path.realpath(__file__)), 'metadata/')

def get_genefamily_percentage(adata, key='MT-', start=True, name='mito'):
    keys = key if isinstance(key, list) else [key, '____ignore____']
    if start:
        family_genes = np.logical_or(*[adata.var_names.str.startswith(k) for k in keys])
    else:
        family_genes = np.logical_or(*[adata.var_names.str.endswith(k) for k in keys])
    if issparse(adata.X):
        adata.obs['percent_'+name] = np.sum(
            adata[:, family_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    else:
        adata.obs['percent_'+name] = np.sum(
            adata[:, family_genes].X, axis=1) / np.sum(adata.X, axis=1)

def get_mito_percentage(adata, species='human'):
    key = 'MT-' if species == 'human' else 'mt-'
    get_genefamily_percentage(adata, key=key, start=True, name='mito')

def get_ribo_percentage(adata, species='human'):
    key = specify_genes(['RPS', 'RPL'], species=species)
    get_genefamily_percentage(adata, key=key, start=True, name='ribo')

def get_hemo_percentage(adata, species='human'):
    key = specify_genes(['HBA', 'HBB'], species=species)
    get_genefamily_percentage(adata, key=key, start=True, name='hemo')

def score_cell_cycle(adata, signatures_path=signatures_path_, species='human'):
    adatas = adata if isinstance(adata, list) else [adata]
    for i in range(len(adatas)):
        adata = adatas[i]
        # score cell cycle
        # cc score with genes from Kowalczyk, Monika S., et al. “Single-Cell RNA-Seq Reveals Changes in Cell Cycle and Differentiation Programs upon Aging of Hematopoietic Stem Cells.” Genome Research, vol. 25, no. 12, 2015, pp. 1860–72, doi:10.1101/gr.192237.115.
        cell_cycle_genes = [x.strip() for x in open(signatures_path+'/regev_lab_cell_cycle_genes.txt')]
        cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
        # Split into 2 lists
        s_genes = cell_cycle_genes[:43]
        g2m_genes = cell_cycle_genes[43:]

        # score
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
        adatas[i] = adata
    return adatas[0] if len(adatas)==1 else adatas

def score_smillie_str_epi_imm(adata, signatures_path=signatures_path_, species='human'):
    tab=pd.read_excel(signatures_path+'/colonoid_cancer_uhlitz_markers_revised.xlsx', skiprows=1, index_col=0)
    score_genes(adata, np.array(tab.index[tab['Epithelial']==1].values, dtype='str'), score_name='epi_score', species=species)
    score_genes(adata, np.array(tab.index[tab['Stromal']==1].values, dtype='str'), score_name='str_score', species=species)
    score_genes(adata, np.array(tab.index[tab['Immune']==1].values, dtype='str'), score_name='imm_score', species=species)


def score_tumor_immune_cells(adata, signatures_path=signatures_path_, species='human'):
    # ImSigGenes immune tumor signatures
    tab=pd.read_excel(signatures_path+'/ImSigGenes_immunetumor.xlsx', skiprows=2, index_col=1)
    annot = dict()
    for ct in pd.unique(tab.Signature):
        annot[ct] = tab[tab.Signature==ct].index.values
    for ct in annot.keys():
        score_genes(adata, annot[ct], score_name=ct, species=species)

def calc_qc_scvelo(adata):
    adatas = adata if isinstance(adata, list) else [adata]
    for adata in adatas:
        # obs qc
        adata.obs['ucounts'] = rsum(adata.layers['unspliced'], axis=1)
        adata.obs['scounts'] = rsum(adata.layers['spliced'], axis=1)
        adata.obs['ufeatures'] = rsum(adata.layers['unspliced']>0, axis=1)
        adata.obs['sfeatures'] = rsum(adata.layers['spliced']>0, axis=1)
        # var qc
        adata.var['ucounts'] = rsum(adata.layers['unspliced'], axis=0)
        adata.var['scounts'] = rsum(adata.layers['spliced'], axis=0)
        adata.var['ucells'] = rsum(adata.layers['unspliced']>0, axis=0)
        adata.var['scells'] = rsum(adata.layers['spliced']>0, axis=0)

def calc_qc(adata, extended_genesets=False, species='detect'):
    adatas = adata if isinstance(adata, list) else [adata]
    for adata in adatas:
        # qc counts
        adata.obs['ncounts'] = rsum(adata.X, axis=1)
        adata.obs['ngenes'] = rsum(adata.X>0, axis=1)
        adata.var['ncounts'] = rsum(adata.X, axis=0)
        adata.var['ncells'] = rsum(adata.X>0, axis=0)

        species = detect_organism(adata) if species == 'detect' else species

        # gene modules
        # mitochondrial genes
        get_mito_percentage(adata, species)
        # ribosomal genes
        get_ribo_percentage(adata, species)
        # hemoglobin genes
        get_hemo_percentage(adata, species)

        if extended_genesets:
            if species is not 'human':
                raise ValueError(species,' species is not known. Pls do not use extended_genesets=True.')
            # interferon genes, immune response
            get_genefamily_percentage(adata, key='IFIT', start=True, name='ifit')
            # Cell adhesion molecules genes
            get_genefamily_percentage(adata, key='CAM', start=False, name='cam')
            # HLA genes encode MHC I and MHC II
            get_genefamily_percentage(adata, key='HLA-', start=True, name='hla')  # genome specific sometimes!!!
            # S100 genes, saw them often in organoids
            get_genefamily_percentage(adata, key='S100', start=True, name='s100')
            # FOX genes, TFs
            get_genefamily_percentage(adata, key='FOX', start=True, name='fox')
            # Heat shock protein genes
            get_genefamily_percentage(adata, key='HSP', start=True, name='heatshock')
            # ABC transporter genes, can lead to multi-drug resistance in cancer
            get_genefamily_percentage(adata, key='ABC', start=True, name='abc')

def specify_genes(genes, species='human'):
    genes = genes if isinstance(genes, list) else list(genes) if isinstance(genes, np.ndarray) else [genes]
    if species is 'human':
        return [x.upper() for x in genes]
    elif species is 'mouse':
        return [x.capitalize() for x in genes]
    else:
        raise ValueError('Species '+species+' not known.')

def score_genes(adata, gene_list, score_name, species='human', **kwargs):
    gene_list_ = specify_genes(gene_list, species=species)
    sc.tl.score_genes(adata, gene_list_, score_name=score_name)

def score_hallmarks(adata, subset='organoid', signatures_path=signatures_path_, species='human'):
    sc.settings.verbosity = 0
    # subset can be a list of hallmarks, 'organoid' (), 'CRC' (~18) or 'all' (50 scores)
    tab = pd.read_csv(signatures_path + 'h.all.v6.2.symbols.gmt', sep='\t', index_col=0, header=None).drop(1, axis=1).T
    hallsigs={hallmark : tab[hallmark][~pd.isna(tab[hallmark])].values for hallmark in tab.columns}
    if isinstance(subset, list):
        selection = subset
    elif subset == 'organoid':
        selection = ['HALLMARK_DNA_REPAIR', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'HALLMARK_APOPTOSIS']
    elif subset == 'CRC':  # TODO this list is bugged, some entries do not exist
        selection = ['HALLMARK_DNA_REPAIR', 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'HALLMARK_APOPTOSIS',
        'HALLMARK_NOTCH_SIGNALING', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_HYPOXIA', 'HALLMARK_TGF_BETA_SIGNALING',
        'HALLMARK_MITOTIC_SPINDLE', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_PROTEIN_SECRETION'
        'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
        'HALLMARK_P53_PATHWAY', 'HALLMARK_ANGIOGENESIS', 'HALLMARK_KRAS_SIGNALING_UP', 'HALLMARK_KRAS_SIGNALING_DN',
        'HALLMARK_GLYCOLYSIS']
    elif subset == 'all':
        selection = hallsigs.keys()
    else:
        raise ValueError('Please select a valid subset of hallmark to use. You can also choose "all".')
    for hm in selection:
        score_genes(adata, hallsigs[hm], score_name=hm, species=species)

def lin_corr_adata(adata, x, keys, method='spearman'):
    """Linearly correlates features (genes/obs_keys) of adata with a given array.
    Computes pearson linear correlation (r and p value) for each selected feature
    with the given values in x.
    ----------
    adata: An adata object.
    x: numeric numpy array
        For example x = adata.obsm['X_diffmap'][:,1] or another gene's
        expression.
    keys: Either a list of genes or a list of adata.obs.columns.
    method: Either 'spearman' or 'pearson'.
    Returns
    -------
    df: A pandas DataFrame
        The dataframe has genes as index and columns pearson_r and pearson_p. It
        is sorted by correlation coefficient (pearson_r).
    """
    # input keys may be list or str, make list
    keys = [keys] if isinstance(keys, str) else keys
    # select correlation method
    if method == 'spearman':
        correlate = spearmanr
    elif method == 'pearsonr':
        correlate = pearsonr
    else:
        raise ValueError(f'Method {method} not valid (pearson or spearman only).')
    # feature set
    if all(np.isin(keys, adata.obs.columns)):
        feature_type = 'obs_keys'
        Y = adata.obs[keys].values
    elif any(np.isin(keys, adata.var_names)):
        feature_type = 'genes'
        Y = adata.X.A if issparse(adata.X) else adata.X
    else:
        raise ValueError('Keys must be list of genes or adata.obs keys.')
    # linearly correlated
    lincors = []
    for i, key in enumerate(keys):
        y = Y[:, i]
        r, p = correlate(x, y)
        lincors.append([key, r, p])
    # format result as pandas.DataFrame
    df = pd.DataFrame(lincors, columns=[feature_type, f'{method}_r', f'{method}_p']).set_index(feature_type)
    df = df.sort_values(f'{method}_r', ascending=False)  # sort by correlation
    return df

def kde_trajectory(adata, key, groupby, velocity=False, rug_keys=[], component=1,
                   figsize=[15,5], n_convolve=10, ax=None, show=True, n_eval=200,
                   range_percs=[0,100], linewidth=4, rug_alpha=0.1,
                   n=19, ylim=30, scale=8):
    X = adata.obsm['X_'+key] if 'X_'+key in adata.obsm.keys() else adata.obsm[key] if key in adata.obsm.keys() else adata.obs[key]
    X = X[:, component] if len(X.shape)>1 else X
    if velocity:
        V = adata.obsm['velocity_'+key] if 'velocity_'+key in adata.obsm.keys() else adata.obsm[key] if key in adata.obsm.keys() else adata.obs[key]
        V = V[:, component] if len(V.shape)>1 else V
    ax = pl.figure(figsize=figsize).gca() if ax is None else ax

    xmin = np.percentile(X, range_percs[0])
    xmax = np.percentile(X, range_percs[1])
    ev = np.linspace(xmin, xmax, n_eval)

    # plot density per group
    for i, cond in enumerate(np.sort(pd.unique(adata.obs[groupby]))):
        mask = adata.obs[groupby] == cond
        kernel = gaussian_kde(X[mask])
        ax.plot(ev, kernel(ev), label=cond, linewidth=linewidth, color=adata.uns[groupby+'_colors'][i])

        if velocity:
            # arrow projections
            edges = np.linspace(ev[0], ev[-1], n)
            bins = [(edges[k]<X) & (X<edges[k+1]) for k in range(n-1)]
            in_bin = bins[2]
            xs = np.array([np.mean(X[mask & in_bin]) for in_bin in bins])
            ys = np.array([np.mean(kernel(X[mask & in_bin])) for in_bin in bins])
            vs = np.array([np.mean(V[mask & in_bin]) if y>ylim else 0 for y, in_bin in zip(ys, bins)])
            vs = vs / np.max(np.abs(vs))
            # pl.plot(X[mask & in_bin], kernel(X[mask & in_bin]), label=cond, linewidth=linewidth, color='red')
            # pl.quiver(X[mask & in_bin], kernel(X[mask & in_bin]), V[mask & in_bin], 0)
            ix = np.abs(vs) > 0
            ax.quiver(xs[ix], ys[ix], vs[ix], 0 , zorder=100, scale_units='width', scale=scale, color=adata.uns[groupby+'_colors'][i])

    # plot categorical annotations as rug
    rug_y = ax.get_ylim()[1]/10
    rug_keys = rug_keys if isinstance(rug_keys, list) else [rug_keys]
    for i, rug in enumerate(rug_keys):
        for j, cond in enumerate(np.sort(pd.unique(adata.obs[rug]))):
            mask = (adata.obs[rug] == cond) & (X>xmin) & (X<xmax)
            plot = ax.plot(X[mask], np.zeros(np.sum(mask)) - rug_y * (i+1), '|', color=adata.uns[rug+'_colors'][j], ms=10, alpha=rug_alpha)

    ax.set_xticks([])
    ax.set_xlabel(key + ' component '+str(component))
    ax.set_ylabel(f'Cell density (KDE) by {groupby}')
    ax.set_yticks(ax.get_yticks()[ax.get_yticks()>=0])
    ax.axhline(y=0, c='k')
    ax.legend()
    if show:
        pl.show()
    else:
        return

def diffusion_analysis_(adata, groupby, species='human', component=1, corr_cutoff=0.1, figsize=[10,8], range_percs=[3,97], velocity_mode=None, show=True):
    """Performs a diffusion analysis on adata for a specific diffusion component.
    velocity_mode may be None, 'on density', 'average' or , 'single'
    ----------
    adata: An adata object.
    Returns
    -------
    None
    """
    ckey = 'DC'+str(component)
    add_velocity_subplot = velocity_mode!=None and velocity_mode!='on density'

    # set layout
    fig = pl.figure(constrained_layout=True, figsize=figsize)
    widths = [1, 1, 1]
    n_rows = 3 + add_velocity_subplot
    heights = [1] * n_rows
    spec = fig.add_gridspec(ncols=3, nrows=n_rows, width_ratios=widths,
                              height_ratios=heights)

    ax0 = fig.add_subplot(spec[0, :])
    kde_trajectory(adata, key='diffmap', groupby=groupby, range_percs=range_percs, ax=ax0,
                   show=False, component=component,
                   velocity=velocity_mode=='on density'
                  )
    ax0.set_xlabel('diffusion pseudotime')
    ax0.set_ylabel('cell density')

    def add_annotation(row, keys, fig, df, name):
        n_top=8
        ax_0 = fig.add_subplot(spec[row, 0])
        ax_1 = fig.add_subplot(spec[row, 1], sharey=ax_0)
        ax_2 = fig.add_subplot(spec[row, 2], sharey=ax_0)
        ax_0.set_axis_off()
        ax_2.set_axis_off()

        # Arrows
        ax_0.annotate('', xy=(.4, 1), xytext=(.6, 1),
                    arrowprops=dict(facecolor='black', shrink=0.05), rotation=90)
        ax_2.annotate('', xy=(.6, 1), xytext=(.4, 1),
                    arrowprops=dict(facecolor='black', shrink=0.05), rotation=90)

        # Texts
        neg_df = df['spearman_r'][df['spearman_r']<-corr_cutoff].iloc[::-1][:n_top]
        pos_df = df['spearman_r'][df['spearman_r']>corr_cutoff][:n_top]
        for i, hallmark in enumerate(neg_df.index):
            ax_0.text(0.5, .8 - i/len(hallmarks), hallmark.replace('HALLMARK_','').replace('_',' '), ha='center', va='center')
        for i, hallmark in enumerate(pos_df.index):
            ax_2.text(0.5, .8 - i/len(hallmarks), hallmark.replace('HALLMARK_','').replace('_',' '), ha='center', va='center')

        # Barplot
        ax_1.barh([.8- i/10 for i in range(len(neg_df))], neg_df.values, align='center', height=0.08, color='tab:blue')
        ax_1.barh([.8- i/10 for i in range(len(pos_df))], pos_df.values, align='center', height=0.08, color='tab:red')
        ax_1.spines['right'].set_visible(False)
        ax_1.spines['left'].set_visible(False)
        ax_1.spines['top'].set_visible(False)
        ax_1.set_yticks([])
        m = np.max(np.abs(df['spearman_r']))
        ax_1.set_xlim([-m,m])
        ax_1.set_xlabel(f'correlation between diffusion axis \n and {name} expression \n (spearman R)')
        ax_1.set_ylim([0,1])

    ### Pathways
    # aggregate hallmarks
    dfs = []
    hallmarks = ['HALLMARK_ANGIOGENESIS', 'HALLMARK_APOPTOSIS', 'HALLMARK_COAGULATION', 'HALLMARK_COMPLEMENT',
    'HALLMARK_IL2_STAT5_SIGNALING', 'HALLMARK_INFLAMMATORY_RESPONSE',
    'HALLMARK_INTERFERON_ALPHA_RESPONSE', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING',
    'HALLMARK_TGF_BETA_SIGNALING', 'HALLMARK_XENOBIOTIC_METABOLISM']
    if not all(np.isin(hallmarks, adata.obs.keys())): score_hallmarks(adata, species=species, subset=hallmarks)
    df_hallmarks = lin_corr_adata(adata, adata.obsm['X_diffmap'][:, component], hallmarks)
    df_hallmarks = df_hallmarks[~pd.isna(df_hallmarks.spearman_r)]
    add_annotation(-2, hallmarks, fig, df_hallmarks, 'signature score')

    ### Genes
    df_genes = lin_corr_adata(adata, adata.obsm['X_diffmap'][:, component], adata.var_names)
    df_genes = df_genes[~pd.isna(df_genes.spearman_r)]
    add_annotation(-1, hallmarks, fig, df_genes, 'gene')

    ### velocities
    if add_velocity_subplot:
        ax1 = fig.add_subplot(spec[1, :], sharex=ax0)
        groups = list(adata.obs[groupby].cat.categories)
        colors = adata.uns[f'{groupby}_colors']
        x = adata.obsm['X_diffmap'][:, component]
        v = adata.obsm['velocity_diffmap'][:, component]
        mask0 = (x>np.percentile(x, range_percs[0])) & (x<np.percentile(x, range_percs[1]))
        if velocity_mode=='single':
            for i, group in enumerate(groups):
                mask = (adata.obs[groupby] == group) & mask0
                ax1.quiver(x[mask], np.random.uniform(1-i, -i, x.shape)[mask], v[mask], np.zeros_like(v)[mask], color=colors[i], scale=0.4, edgecolor='k', linewidth = .5)
            ax1.set_ylabel(f'RNA velocity\nby {groupby}')
        else:
            from scipy.interpolate import interp1d
            n_evals = 10
            xint=np.linspace(np.percentile(x, range_percs[0]), np.percentile(x, range_percs[1]), n_evals)
            for i, group in enumerate(groups):
                mask = (adata.obs[groupby] == group) & mask0
                f = interp1d(x[mask], v[mask])
                x_int = xint[(xint >= np.min(x[mask])) & (xint <= np.max(x[mask]))]
                v_int = f(x_int)
                # Normalize
                v_absmax = np.max(np.abs(v_int))
                x_segment = (x_int[1] - x_int[0]) / (n_evals/5)
                v_int = v_int * x_segment / v_absmax
                ax1.quiver(x_int, i * np.ones_like(x_int), v_int, np.zeros_like(v_int),
                           headwidth=4, color=colors[i], edgecolor='k', linewidth = .5, angles='xy', scale_units='xy', scale=1)
            ax1.set_ylim(-1, len(groups))
            ax1.set_ylabel(f'Average RNA velocity\nby {groupby}')
        ax1.set_yticks([])
    # pl.suptitle('Neutrophil cell density on diffusion pseudotime')
    if show: pl.show()

def identify_barcode_overlap(df1, df2, key1, key2, reg1='[ACGT]+-', reg2='[ACGT]+-', kick=-1, plot=True):
    # clear index
    x1 = np.array([re.findall(reg1, txt)[0][:kick] for txt in df1.index])
    x2 = np.array([re.findall(reg2, txt)[0][:kick] for txt in df2.index])

    # count co-occurences of barcodes by key categories
    c1 = pd.unique(df1[key1])
    c2 = pd.unique(df2[key2])
    Z = np.zeros((len(c1), len(c2)))
    for i, ci in enumerate(c1):
        for j, cj in enumerate(c2):
            Z[i,j] = np.sum(np.isin(x1[df1[key1]==ci], x2[df2[key2]==cj]))
    X = pd.DataFrame(Z, index=c1, columns=c2)
    if plot: sns.heatmap(X, annot=False), pl.show()
    return X

def get_subfolders(d, full_path=True):
    prefix=d if full_path else ''
    return [os.path.join(prefix, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]

def get_files(d, full_path=True):
    prefix=d if full_path else ''
    return [os.path.join(prefix, f) for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

def force_merge(df1, df2):
    # pd.concat([adata.obs, tab], axis=0, ignore_index=True) is not working as one would think
    # see https://stackoverflow.com/questions/32801806/pandas-concat-ignore-index-doesnt-work
    df = pd.DataFrame(np.concatenate([df1, df2], axis=1),
                             index=df1.index,
                             columns=list(df1.columns) + list(df2.columns))
    df = df.loc[:,~df.columns.duplicated()]  # remove duplicate columns
    return df

def peek(f):
    opener = open if f.split('.')[-1]!='gz' else lambda x: gzip.open(x, 'rb')
    # peek into file to find out length and separator
    file_length = sum(1 for line in opener(f))

    for line in opener(f):
        try:
            first_line = line.decode()
        except (UnicodeDecodeError, AttributeError):
            first_line = line
        break
    sniffer = Sniffer()
    dialect = sniffer.sniff(first_line)
    separator=dialect.delimiter
    return file_length, separator

def gene_symbols_to_entrezid(gene_list, species='human', verbose=False):
    mg = mygene.MyGeneInfo()
    out = mg.querymany(gene_list, scopes='symbol', fields='entrezgene', species=species, verbose=verbose)
    df = pd.DataFrame([[o['query'], o['_id']] if '_id' in o.keys() else [o['query'], None] for o in out], columns=['gene_symbol', 'entrez_id']).set_index('gene_symbol')
    return df

# NOTE: you need to run beforehand:
# from goatools.base import download_go_basic_obo, download_ncbi_associations
# obo_fname = download_go_basic_obo(goa_path+'go-basic.obo')
# fin_gene2go = download_ncbi_associations(goa_path+'gene2go')
# also see:
# https://github.com/tanghaibao/goatools/blob/main/notebooks/goea_nbt3102.ipynb
goa_path = '/fast/work/users/peidlis_c/utils/goa/'  # replace with your path

def GOEA(gene_list, species='human', namespaces=['BP'], sig_alpha=0.05, verbose=False):
    """Performs GO enrichment analysis with goatools.
    Based on https://github.com/tanghaibao/goatools/blob/main/notebooks/goea_nbt3102.ipynb.
    Note that you must ensure goa_path is filled (see jnb link above) by running:
        from goatools.base import download_go_basic_obo, download_ncbi_associations
        obo_fname = download_go_basic_obo(goa_path+'go-basic.obo')
        fin_gene2go = download_ncbi_associations(goa_path+'gene2go')
    ----------
    gene_list: `list` of `str`
        A list of gene symbols.
    species: `str`, either `'human'` or `'mouse'` (default: `'human'`)
        The species the genes came from.
    namespaces: `list` of `str` (default: `['BP']`)
        A `list` of strings from `['BP', 'CC', 'MF']`.
        BP: Biological Process (larger processes, e.g. immun)
        CC: Cellular Component (location in the cell)
        MF: Molecular Function (small process)
        See http://geneontology.org/docs/ontology-documentation/.
    sig_alpha: `float` in `[0,1]` (default: `0.05`)
        Significance cut-off for multiple testing corrected p values.
    Returns
    -------
    df: A pandas DataFrame
        The dataframe has GO term ids as index and is sorted by multiple testing
        corrected p values. Gene ids of the respective GO term are found in the
        column 'study_items'.
    """
    largs = {} if verbose else {'prt': None}

    # species handling
    taxids_dics = {'human' : 9606, 'mouse': 10090}
    if species not in taxids_dics.keys():
        raise ValueError('Species ', species, ' not known...')
    else:
        taxid = taxids_dics[species]

    # convert gene symbols to entrez ids
    df = gene_symbols_to_entrezid(gene_list, species=species, verbose=verbose)
    genes = [x for x in df.entrez_id if x!=None]
    geneids_study = [int(x) for x in genes if x.isdigit()]

    # read GO-relevant databases
    obodag = GODag(goa_path+"go-basic.obo", **largs)
    objanno = Gene2GoReader(goa_path+'gene2go', taxids=[taxid], **largs)
    ns2assoc = objanno.get_ns2assc()

    # define background
    if species == 'mouse':
        from metadata.goatools_bg_genes.genes_ncbi_10090_proteincoding import GENEID2NT
    elif species == 'human':
        from metadata.goatools_bg_genes.genes_ncbi_9606_proteincoding import GENEID2NT

    # initialize GOEA object
    goeaobj = GOEnrichmentStudyNS(
        GENEID2NT.keys(), # List of protein-coding genes
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = sig_alpha, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method
    # run analysis
    goea_results_all = goeaobj.run_study(geneids_study, **largs)
    goea_results = [r for r in goea_results_all if r.p_fdr_bh < sig_alpha]

    if len(goea_results) > 0:
        df = pd.DataFrame([o.__dict__ for o in goea_results])
        df = df.set_index('GO').drop(['kws', 'method_flds'], axis=1)
        df = df[np.isin(df.NS, namespaces)]
    else:
        df = None
        print('Results are empty. Check if gene set or species wrong. Or set only_significant=False.')
    return df
