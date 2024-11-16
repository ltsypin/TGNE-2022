import os
import pandas as pd
import numpy as np
import warnings
import scipy.spatial
import scipy.cluster.hierarchy
from copy import deepcopy

import bokeh
from bokeh import events
from bokeh.events import Reset
from bokeh.plotting import show as show_interactive
from bokeh.plotting import output_file, output_notebook
from bokeh.layouts import column, row
from bokeh.models import ResetTool, TabPanel, Tabs, Circle, Div, ColumnDataSource, CustomJS, TextInput, LassoSelectTool, Select, MultiSelect, ColorBar, Legend, LegendItem, Spinner
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn, Button, HTMLTemplateFormatter
from bokeh.events import SelectionGeometry
from bokeh.transform import linear_cmap, jitter
from bokeh.core.enums import ResetPolicy

import umap
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.manifold import MDS

from .clustering_utils import compute_pairwise_distance_matrix

# bokeh_ui_utils

file_dir = os.path.dirname(os.path.abspath(__file__))
main_dir = '../TGNE/'

# The two functions below are taken and adapted from the UMAP package
def _get_umap_embedding(umap_object):
    if hasattr(umap_object, "embedding_"):
        return umap_object.embedding_
    elif hasattr(umap_object, "embedding"):
        return umap_object.embedding
    else:
        raise ValueError("Could not find embedding attribute of umap_object")
        
def plot_enrichment(enrich_column_data_source, plot_sizing_mode='stretch_both'):
    
    # pdb.set_trace()
    
    # y_range = FactorRange(factors=[str(y) for y in enrich_df['module'].unique()])
    
    # grouped = enrich_df.groupby('module')
    
    hover = [
        ('module', '@module'),
        ('term', '@term'),
        ('info', '@info'),
        ('fold-change', '@fold_change'),
        ('bonferroni', '@bonferroni')
    ]
    
    p = bokeh.plotting.figure(
        # y_range=y_range,
        title='Functional term enrichment in modules',
        x_axis_label='Fold change',
        y_axis_label='Module',
        x_axis_type='log',
        tooltips=hover,
        # background_fill_color='black'
        sizing_mode=plot_sizing_mode,
        output_backend="webgl",
        tools=["ypan", "box_zoom", "ywheel_zoom", "undo", "reset", "save"]
    )
    
    # cds = bokeh.models.ColumnDataSource(enrich_df)
    # print(enrich_df.head())
    
    renderer = p.circle(y=jitter('module', width=0.4), x='fold_change', source=enrich_column_data_source, alpha='alpha', size='size', color='color', line_color='line_color')
    # p.xaxis.major_label_orientation = 45
    p.ygrid.minor_grid_line_color = 'navy'
    p.ygrid.minor_grid_line_alpha = 0.1
    # p.xgrid.band_fill_alpha = 0.1
    # p.xgrid.band_fill_color = "navy"
    ticker = []
    for m in enrich_column_data_source.data['module']:
        if m not in ticker:
            ticker.append(m)
    p.yaxis.ticker = ticker
    p.y_range.flipped = True
    p.xaxis.major_label_text_font_size = '12pt'
    p.yaxis.major_label_text_font_size = '8pt'
    p.yaxis.axis_label_text_font_size = '8pt'
    p.xaxis.axis_label_text_font_size = '12pt'

    selected_circle = Circle(fill_alpha=0.3, 
                            #  size=7, 
                             line_color='black', fill_color='color')
    nonselected_circle = Circle(fill_alpha=0.05, 
                                # size=1, 
                                line_color=None, fill_color='color')

    renderer.selection_glyph = selected_circle
    renderer.nonselection_glyph = nonselected_circle
    
    return p

def heatmap(column_data_source, ls_color_palette, r_low, r_high, x_axis_factors, y_axis_factors, s_z="normalized_expression", index_name='TTHERM_ID', col_name='phase', plot_sizing_mode='stretch_both', hover_module=True):
    # adapted from https://gitlab.com/biotransistor/bokehheat/-/blob/master/bokehheat/heat.py
    """
    input:
        df_matrx: a dataframe in same xy orientation as the final heatmap.
          the index should cary the y axis label.
          the column should cary the x axis label.
          the matrix as such should only cary the z axis values.

        ls_color_palette: a list color strings to specify the color spectrum.
            this variable is compatible with the ordinary bokeh palettes:
            https://bokeh.pydata.org/en/latest/docs/reference/palettes.html

        r_low: quantitative minimum value. the dataset can contain lower values,
            but for color labeling they will be mapped to this minimum value.
            e.g.: -8.

        r_high: quantitative maximum value. the dataset can contain lower values,
            but for color labeling they will be mapped to this maximum value.
            e.g.: 8.

        s_z: string. label that specifies what the values in the matrix actually
            are. e.g.: 'gene expression [log2]'

    output:
        p: bokeh plot object.

    description:
        this function will return a bokeh based interactive heatmap plot.
        the color are representing the z value.
    """
    # index as string
#     df_matrix.index = df_matrix.index.astype(str)
#     df_matrix.columns = df_matrix.columns.astype(str)

#     # handle y and x axis name
#     if (df_matrix.index.name == None):
#         df_matrix.index.name = "y_axis"
#     if (df_matrix.columns.name == None):
#         df_matrix.columns.name = "x_axis"
    # pdb.set_trace()
    s_y = index_name
    
    # df_matrix.columns.name = 'phase'
    s_x = col_name
    
    
    # print(df_matrix.head())
    
    # melt dataframe
    # df_tidy = df_matrix.reset_index().melt(
    #     id_vars=[df_matrix.index.name],
    #     value_name=s_z
    # )
    # print(df_tidy.head())
    # color declaration
    d_zcolormapper = linear_cmap(
        field_name=s_z,
        palette=ls_color_palette,
        low=r_low,
        high=r_high
    )
    # tooltip declaration
    lt_tooltip = [
        (s_y, f"@{s_y}"),
        (s_x, f"@{s_x}"),
        (s_z, f"@{s_z}"),
    ]

    if hover_module:
        lt_tooltip.append(('module', f'@module'))

    # generate figure
    o_colorbar = ColorBar(color_mapper=d_zcolormapper['transform'])
    p = bokeh.plotting.figure(
        y_range=y_axis_factors,
        x_range=x_axis_factors,
        tools = "box_zoom,hover,pan,reset,ywheel_zoom,save",  # have to be set hardcoded
        active_drag = "box_zoom",  # have to be set hardcoded
        tooltips=lt_tooltip,
        title='Whole genome normalized expression',
        toolbar_location='right',
        sizing_mode=plot_sizing_mode,
        output_backend="webgl"
    )
    
    p.rect(
        source=column_data_source,
        x=s_x,
        y=s_y,
        color=d_zcolormapper,
        width=1,
        height=1,
        fill_alpha='fill_alpha',
        line_alpha='line_alpha',
        # line_color='white',
        # nonselection_fill_alpha=0.01,
        # nonselection_line_alpha=0.01,
        # nonselection_line_color="white"
    )
    p.add_layout(o_colorbar, place='left')
    # p.yaxis.major_label_orientation = "horizontal"
    p.xaxis.major_label_orientation = 45
    # p.yaxis.major_label_text_font_size = '0pt'
    p.yaxis.visible = False
    p.xaxis.major_label_text_font_size = '12pt'

    # out
    return(p)
        
def interactive(
    embedding_df,
    enrich_df,
    num_genes,
    x,
    # mean_expression_df,
    title=None,
    labels=None,
    values=None,
    hover_data=None,
    theme=None,
    cmap="Blues",
    color_key=None,
    color_key_cmap="Spectral",
    background="white",
#     width=800,
#     height=800,
    point_size=None,
    radius=None, # My contribution
#     subset_points=None,
    interactive_text_search=False,
    interactive_text_search_columns=[
        'TTHERM_ID', 
        'Description', 
        'TGD2021_description',
        'InterPro_description',
        'module', 
        'common_name', 
        'Preferred_name'
    ],
    interactive_text_search_columns2=[
        # 'COG_category',
        'EC',
        'GOs',
        'PFAMs',
        'KEGG_ko',
        'InterPro',
        'KEGG_Pathway',
        'KEGG_Module',
        'KEGG_Reaction',
        'KEGG_rclass',
        # 'BRITE',
        'KEGG_TC',
        # 'CAZy',
        # 'BiGG_Reaction'
    ],
    interactive_text_search_alpha_contrast=0.999,
    alpha=None,
    expr_min = 0,
    expr_max = 1,
    plot_sizing_mode='stretch_both',
    table_sizing_mode = 'stretch_both',
    search_sizing_mode = 'stretch_both',
    avg_df=None,
    avg_radius=None,
):
    """Create an interactive bokeh plot of a UMAP embedding.
    While static plots are useful, sometimes a plot that
    supports interactive zooming, and hover tooltips for
    individual points is much more desireable. This function
    provides a simple interface for creating such plots. The
    result is a bokeh plot that will be displayed in a notebook.
    Note that more complex tooltips etc. will require custom
    code -- this is merely meant to provide fast and easy
    access to interactive plotting.
    Parameters
    ----------
    embedding_df: pandas DataFrame
        A expression dataframe with columns x and y, which are the
        2D embedding of a model (e.g., UMAP or pyMDE) on the expression data, and all the
        annotations, geometric means of expression, etc.
    x: list
        The categories for the x-axes of the heatmap and expression profiles
    labels: array, shape (n_samples,) (optional, default None)
        An array of labels (assumed integer or categorical),
        one for each data sample.
        This will be used for coloring the points in
        the plot according to their label. Note that
        this option is mutually exclusive to the ``values``
        option.
    values: array, shape (n_samples,) (optional, default None)
        An array of values (assumed float or continuous),
        one for each sample.
        This will be used for coloring the points in
        the plot according to a colorscale associated
        to the total range of values. Note that this
        option is mutually exclusive to the ``labels``
        option.
    hover_data: DataFrame, shape (n_samples, n_tooltip_features)
    (optional, default None)
        A dataframe of tooltip data. Each column of the dataframe
        should be a Series of length ``n_samples`` providing a value
        for each data point. Column names will be used for
        identifying information within the tooltip.
    theme: string (optional, default None)
        A color theme to use for plotting. A small set of
        predefined themes are provided which have relatively
        good aesthetics. Available themes are:
           * 'blue'
           * 'red'
           * 'green'
           * 'inferno'
           * 'fire'
           * 'viridis'
           * 'darkblue'
           * 'darkred'
           * 'darkgreen'
    cmap: string (optional, default 'Blues')
        The name of a matplotlib colormap to use for coloring
        or shading points. If no labels or values are passed
        this will be used for shading points according to
        density (largely only of relevance for very large
        datasets). If values are passed this will be used for
        shading according the value. Note that if theme
        is passed then this value will be overridden by the
        corresponding option of the theme.
    color_key: dict or array, shape (n_categories) (optional, default None)
        A way to assign colors to categoricals. This can either be
        an explicit dict mapping labels to colors (as strings of form
        '#RRGGBB'), or an array like object providing one color for
        each distinct category being provided in ``labels``. Either
        way this mapping will be used to color points according to
        the label. Note that if theme
        is passed then this value will be overridden by the
        corresponding option of the theme.
    color_key_cmap: string (optional, default 'Spectral')
        The name of a matplotlib colormap to use for categorical coloring.
        If an explicit ``color_key`` is not given a color mapping for
        categories can be generated from the label list and selecting
        a matching list of colors from the given colormap. Note
        that if theme
        is passed then this value will be overridden by the
        corresponding option of the theme.
    background: string (optional, default 'white')
        The color of the background. Usually this will be either
        'white' or 'black', but any color name will work. Ideally
        one wants to match this appropriately to the colors being
        used for points etc. This is one of the things that themes
        handle for you. Note that if theme
        is passed then this value will be overridden by the
        corresponding option of the theme.
    width: int (optional, default 800)
        The desired width of the plot in pixels.
    height: int (optional, default 800)
        The desired height of the plot in pixels
    point_size: int (optional, default None)
        The size of each point marker
    radius: int (optional, default None)
        The radius of each point marker (adjusts the point size while zooming)
    subset_points: array, shape (n_samples,) (optional, default None)
        A way to select a subset of points based on an array of boolean
        values.
    interactive_text_search: bool (optional, default False)
        Whether to include a text search widget above the interactive plot
    interactive_text_search_columns: list (optional, default None)
        Columns of data source to search. Searches labels and hover_data by default.
    interactive_text_search_alpha_contrast: float (optional, default 0.95)
        Alpha value for points matching text search. Alpha value for points
        not matching text search will be 1 - interactive_text_search_alpha_contrast
    alpha: float (optional, default: None)
        The alpha blending value, between 0 (transparent) and 1 (opaque).
    Returns
    -------
    """
    import textwrap
    
    if 'TTHERM_IDs' in list(embedding_df.columns):
        interactive_text_search_columns=[
        'TTHERM_ID', 
        'Description', 
        'TGD2021_description', 
        'module', 
        'common_name', 
        'Preferred_name'
    ]

    if theme is not None:
        cmap = _themes[theme]["cmap"]
        color_key_cmap = _themes[theme]["color_key_cmap"]
        background = _themes[theme]["background"]

    if labels is not None and values is not None:
        raise ValueError(
            "Conflicting options; only one of labels or values should be set"
        )
        
    if alpha is not None:
        if not 0.0 <= alpha <= 1.0:
            raise ValueError("Alpha must be between 0 and 1 inclusive")

    if point_size is None and radius is None:
        point_size = 100.0 / np.sqrt(points.shape[0])
        
    data = embedding_df
    # data = data.set_index('TTHERM_ID')
    # pdb.set_trace()
    if radius is not None:
        data['radius'] = radius

    if labels is not None:
        data["label"] = labels

        if color_key is None:
            unique_labels = np.unique(labels)
            num_labels = unique_labels.shape[0]
            color_key = _to_hex(
                plt.get_cmap(color_key_cmap)(np.linspace(0, 1, num_labels))
            )

        if isinstance(color_key, dict):
            data["color"] = pd.Series(labels).map(color_key)
        else:
            # print('here')
            unique_labels = np.unique(labels)
            if len(color_key) < unique_labels.shape[0]:
                # raise ValueError(
                #     "Color key must have enough colors for the number of labels"
                # )
                
                print('Color key has fewer colors than labels. Making all green')
                data['color'] = ['green'] * len(labels)
                avg_df["color"] = ['green'] * avg_df.shape[0]
            else:
                new_color_key = {k: color_key[i] for i, k in enumerate(unique_labels)}
                data["color"] = pd.Series(labels).map(new_color_key)
                avg_df["color"] = pd.Series(avg_df["label"].values).map(new_color_key)

        colors = "color"

    elif values is not None:
        data["value"] = values
        palette = _to_hex(plt.get_cmap(cmap)(np.linspace(0, 1, 256)))
        colors = btr.linear_cmap(
            "value", palette, low=np.min(values), high=np.max(values)
        )

    else:
        colors = matplotlib.colors.rgb2hex(plt.get_cmap(cmap)(0.5))

    # print(data['color'].unique())
    # print(colors)

    if hover_data is not None:
        tooltip_dict = {}
        for col_name in hover_data:
            data[col_name] = hover_data[col_name]
            tooltip_dict[col_name] = "@{" + col_name + "}"
        tooltips = list(tooltip_dict.items())
    else:
        tooltips = None

    if alpha is not None:
        data["alpha"] = alpha
    else:
        alpha = 1
        data["alpha"] = alpha
    
    # print(list(hover_data['module'].values)[0:10])
    # print(list(hover_data['module'].values)[len(list(hover_data['module'].values))-11:len(list(hover_data['module'].values))-1])

    # print('RADIUS:', radius)

    data_source = bokeh.plotting.ColumnDataSource(data)
    data_source.data['module'] = hover_data['module']
    data_source.data['ID'] = hover_data['ID']
    # data_source.data['YF_ID'] = hover_data['YF_ID']
    data_source.data['radius'] = np.ones_like(hover_data['ID']) * radius
    data_source.data['alpha'] = np.ones_like(hover_data['ID']) * alpha
    
    # print(data_source.data['ID'][:5])

    x_pad = (data['x'].max() - data['x'].min()) * 0.05
    y_pad = (data['y'].max() - data['y'].min()) * 0.05

    plot = bokeh.plotting.figure(
        tooltips=tooltips,
        tools="tap,lasso_select,box_select,pan,wheel_zoom,box_zoom,reset,save",
        background_fill_color=background,
        title="",
        sizing_mode=plot_sizing_mode,
        output_backend="webgl",
        x_range=(data['x'].min() - x_pad, data['x'].max() + x_pad),
        y_range=(data['y'].min() - y_pad, data['y'].max() + y_pad),
    )

    ix_start = data['x'].min() - x_pad
    ix_end = data['x'].max() + x_pad
    iy_start = data['y'].min() - y_pad
    iy_end = data['y'].max() + y_pad

    if point_size is not None:

        renderer = plot.circle(
            x="x",
            y="y",
            source=data_source,
            color=colors,
            size=point_size,
            alpha="alpha",
            line_color='black'
        )

    elif radius is not None:
        renderer = plot.circle(
            x="x",
            y="y",
            source=data_source,
            color=colors,
            radius=radius,
            alpha="alpha",
            line_color='black'
        )

    plot.grid.visible = False
    plot.axis.visible = False

    plot.reset_policy = 'event_only'

    reset_callback = CustomJS(args=dict(plot=plot, ix_start=ix_start, ix_end=ix_end, iy_start=iy_start, iy_end=iy_end), code="""
    console.log('Custom reset callback triggered!');
    plot.x_range.start = ix_start;
    plot.x_range.end = ix_end;
    plot.y_range.start = iy_start;
    plot.y_range.end = iy_end;
    """)

    plot.js_on_event(Reset, reset_callback)

    selected_circle = Circle(fill_alpha=1, radius=radius, line_color='black', fill_color='color')
    nonselected_circle = Circle(fill_alpha=0.05, radius=radius/20, line_color=None, fill_color='color')

    renderer.selection_glyph = selected_circle
    renderer.nonselection_glyph = nonselected_circle


    avg_tooltip_dict = {}
    for col_name in avg_df:
        if col_name not in ['x', 'y', 'color']:
            avg_tooltip_dict[col_name if col_name != 'label' else 'module'] = "@{" + col_name + "}"
    avg_tooltips = list(avg_tooltip_dict.items())


    avg_data_source = bokeh.plotting.ColumnDataSource(avg_df)
    avg_data_source.data['radius'] = np.ones_like(avg_df['label']) * avg_radius
    avg_data_source.data['alpha'] = np.ones_like(avg_df['label']) * alpha
    avg_data_source.data['line_color'] = avg_df.shape[0] * ['black']

    avg_x_pad = (avg_df['x'].max() - avg_df['x'].min()) * 0.05
    avg_y_pad = (avg_df['y'].max() - avg_df['y'].min()) * 0.05

    plot_avg = bokeh.plotting.figure(
    tooltips=avg_tooltips,
    tools="tap,lasso_select,box_select,pan,wheel_zoom,box_zoom,reset,save",
    background_fill_color=background,
    title="",
    sizing_mode=plot_sizing_mode,
    output_backend="webgl",
    x_range=(avg_df['x'].min() - avg_x_pad, avg_df['x'].max() + avg_x_pad),
    y_range=(avg_df['y'].min() - avg_y_pad, avg_df['y'].max() + avg_y_pad),
    )

    avg_ix_start = avg_df['x'].min() - avg_x_pad
    avg_ix_end = avg_df['x'].max() + avg_x_pad
    avg_iy_start = avg_df['y'].min() - avg_y_pad
    avg_iy_end = avg_df['y'].max() + avg_y_pad

    avg_renderer = plot_avg.circle(
        x="x",
        y="y",
        source=avg_data_source,
        color=colors,
        radius=radius,
        alpha="alpha",
        line_color='black'
    )

    plot_avg.reset_policy = 'event_only'

    avg_reset_callback = CustomJS(args=dict(plot=plot_avg, ix_start=avg_ix_start, ix_end=avg_ix_end, iy_start=avg_iy_start, iy_end=avg_iy_end), code="""
    console.log('Custom avg reset callback triggered!');
    plot.x_range.start = ix_start;
    plot.x_range.end = ix_end;
    plot.y_range.start = iy_start;
    plot.y_range.end = iy_end;
    """)

    plot_avg.js_on_event(Reset, avg_reset_callback)

    plot_avg.grid.visible = False
    plot_avg.axis.visible = False

    avg_selected_circle = Circle(fill_alpha=1, radius=radius, line_color='black', fill_color='color')
    avg_nonselected_circle = Circle(fill_alpha=0.05, radius=radius/20, line_color=None, fill_color='color')

    avg_renderer.selection_glyph = avg_selected_circle
    avg_renderer.nonselection_glyph = avg_nonselected_circle


    
    x_heatmap_profile = x

    # rna_seq_phase_dict = {
    #         '000min': '(000min) G1',
    #         '030min': '(030min) G1',
    #         '060min': '(060min) S',
    #         '090min': '(090min) S',
    #         '120min': '(120min) G2',
    #         '150min': '(150min) M',
    #         '180min': '(180min) M',
    #         '210min': '(210min) G1',
    #         '240min': '(240min) S',
    # }

    # for idx in range(len(x_heatmap_profile)):
    #     if x_heatmap_profile[idx] in rna_seq_phase_dict:
    #         x_heatmap_profile[idx] = rna_seq_phase_dict[x_heatmap_profile[idx]]

    # for FIXME add the rna seq phases
    
    # if normalized:
    hm_min = expr_min
    hm_max = expr_max
        
    # else:
    #     hm_min = 2
    #     hm_max = 16
    
    # For companion heatmap plot
    ttherm_ids = embedding_df['TTHERM_ID'].values
    hm_df = embedding_df.loc[:, ['TTHERM_ID'] + x]
    hm_df['module'] = hover_data.loc[:,:]['module'].values
    hm_df_tidy = hm_df.melt(id_vars=['TTHERM_ID', 'module'], var_name='phase', value_name='normalized_expression')
    hm_cds = bokeh.plotting.ColumnDataSource(hm_df_tidy)
    hm_cds.data['fill_alpha'] = [0.7]*len(hm_df_tidy)
    hm_cds.data['line_alpha'] = [0.7]*len(hm_df_tidy)
    # hm_cds.data['y_axis'] = ttherm_ids

    
    hm = heatmap(hm_cds, bokeh.palettes.Inferno256, hm_min, hm_max, x_heatmap_profile, ttherm_ids, plot_sizing_mode=plot_sizing_mode)
    
    

    # For companion expression plot

    expr_source = bokeh.plotting.ColumnDataSource(dict(
        TTHERM_ID=['blah'],
        module=['blah'],
        ID=['blah'], 
        expr_xs=[['Ll']], 
        expr_ys=[[0]],
        alpha=[0],
        color=['black'],
        line_dash = ['solid'],
        ))
    
    # print(expr_source.data)
    
    # if normalized:
    y_axis_label = 'Mean normalized expression'
    y_range = (expr_min - 0.01, expr_max + 0.01)
        
    # else:
    #     y_axis_label = 'Geometric mean expression of replicates (log2-scale)'
    #     y_range = (3.9, 16.1)
    
    expr_fig = bokeh.plotting.figure(
                                     background_fill_color=background,
                                     # x_axis_label='Phase or condition',
                                     y_axis_label=y_axis_label,
                                     x_range=x_heatmap_profile, 
                                     y_range=y_range,
                                     tooltips=[
                                         ('TTHERM_ID', '@TTHERM_ID'),
                                         ('module', '@module'),
                                               ],
                                     sizing_mode=plot_sizing_mode,
                                     output_backend="webgl",
                                     tools=["ypan", "xpan", "box_zoom", "ywheel_zoom", "xwheel_zoom", "undo", "reset", "save"]
                                    )

    expr_fig.multi_line('expr_xs', 
                        'expr_ys', 
                        source=expr_source, 
                        alpha='alpha', 
                        line_width=3, 
                        line_join='round',
                        line_color="color",
                        line_dash = 'line_dash',
                       )

    expr_fig.xaxis.major_label_orientation = np.pi/4
    expr_fig.xaxis.major_label_text_font_size = '12pt'
    expr_fig.yaxis.major_label_text_font_size = '12pt'
    expr_fig.yaxis.axis_label_text_font_size = '12pt'
    expr_fig.xgrid.grid_line_color='black'
    expr_fig.xgrid.grid_line_alpha=0.2
    expr_fig.ygrid.grid_line_color='black'
    expr_fig.ygrid.grid_line_alpha=0.2



    all_columns = [
        'module',
        'common_name',
        'TGD2021_description',
        'Description',
        'InterPro_description',
        'Preferred_name',
        'peptide',
        'max_annot_lvl',
        'COG_category',
        'EC',
        'GOs',
        'PFAMs',
        'KEGG_ko',
        'InterPro',
        'KEGG_Pathway',
        'KEGG_Module',
        'KEGG_Reaction',
        'KEGG_rclass',
        'BRITE',
        'KEGG_TC',
        'CAZy',
        'BiGG_Reaction'
    ]

    if 'TTHERM_IDs' in list(embedding_df.columns):
        all_columns.insert(0, 'TTHERM_IDs')

    # For data table
    s2 = bokeh.plotting.ColumnDataSource(data=dict(ID=[]))

    columns_with_widths = {
        "module": 100,
        "common_name": 100,
        "peptide": 100,
        "TGD2021_description": 300,
        "Description": 300,
        'InterPro_description': 300,
        "Preferred_name": 150,
        "max_annot_lvl": 100,
        "COG_category": 100,
        "EC": 100,
        "GOs": 300,
        "PFAMs": 300,
        "KEGG_ko": 300,
        "InterPro": 300,
        "KEGG_Pathway": 300,
        "KEGG_Module": 300,
        "KEGG_Reaction": 300,
        "KEGG_rclass": 300,
        "BRITE": 300
    }

    table_columns = [c for c in all_columns if c in columns_with_widths]

    columns = [TableColumn(field="ID",  title="TTHERM_ID", formatter=HTMLTemplateFormatter(template='<a href="http://tet.ciliate.org/index.php/feature/details/feature_details.php?feature_name=<%= ID %>"target="_blank"><%= ID %></a>'), width=150),
            #    TableColumn(field="YF_ID", title="YF_ID", width=300),

            #    TableColumn(field='KEGG_TC', title='KEGG_TC', width=300),
            #    TableColumn(field='CAZy', title='CAZy', width=300),
            #    TableColumn(field='BiGG_Reaction', title='BiGG_Reaction', width=300),
#                    TableColumn(field="x",  title="x", width=300),
#                    TableColumn(field="y",  title="y")
              ] + [
                  TableColumn(field=c, title=c, width=columns_with_widths[c]) for c in table_columns
              ]
    
    if 'TTHERM_IDs' in list(embedding_df.columns):
        columns.insert(1, TableColumn(field="TTHERM_IDs",  title="TTHERM_IDs", width=600))

    table = DataTable(source=s2, 
                      columns=columns, 
                      editable=True,
                      selectable=True,
                      sortable=True,
                    #   index_width=10,
                      fit_columns=False,
                      sizing_mode=table_sizing_mode
                     )
    

    callback = CustomJS(args=dict(s_expr=expr_source, s2=s2), code="""
    var inds = cb_obj.indices;
    // console.log(inds);
    var d_expr = s_expr.data;
    var d2 = s2.data;

    d_expr['line_dash'] = Array(d_expr['line_dash'].length).fill('solid');
    
    d_expr['alpha'] = Array(d_expr['alpha'].length).fill(Math.min(1, Math.max(7/(d2['ID'].length), 0.05)));
    
    
    if (inds.length > 0) {
        d_expr['line_dash'] = Array(d_expr['line_dash'].length).fill('solid');
        d_expr['line_dash'][inds[0]+1] = 'solid';
                        
        d_expr['alpha'] = Array(d_expr['alpha'].length).fill(Math.min(1/3, Math.max(7/(d2['ID'].length)/3, 0.05)));
        d_expr['alpha'][inds[0]+1] = 1;
    }
                        
    s_expr.change.emit();    
    """)

    s2.selected.js_on_change('indices', callback)
    
    enrich_df = enrich_df
    colors = [color_key[int(m) % len(color_key)] for m in enrich_df['module'].values]
    enrich_df['color'] = colors
    enrich_df['alpha'] = enrich_df.shape[0] * [0.3]
    enrich_df['size'] = enrich_df.shape[0] * [7]
    enrich_df['line_color'] = enrich_df.shape[0] * ['black']

    enrich_df['info'] = enrich_df['info'].astype('str')

    enrich_df2 = enrich_df.copy()
    
    
    enrich_cds = bokeh.models.ColumnDataSource(enrich_df)

    # For the hovertools to not get overwhelmed by the length of annotation text.
    enrich_df2['info'] = enrich_df2['info'].apply(textwrap.shorten, width=50)
    enrich_cds2 = bokeh.models.ColumnDataSource(enrich_df2)
    
    enrich_p = plot_enrichment(enrich_cds2, plot_sizing_mode=plot_sizing_mode)


    if interactive_text_search:
        text_input = TextInput(value="", placeholder=f'Comma-separated descriptive terms: module(s), ID(s), names, or descriptions', sizing_mode=search_sizing_mode)
        text_input2 = TextInput(value="", placeholder=f'Comma-separated functional terms: PFAM names or InterPro/GO/KEGG/EC codes', sizing_mode=search_sizing_mode)

        if interactive_text_search_columns is None:
            interactive_text_search_columns = []
            if hover_data is not None:
                interactive_text_search_columns.extend(hover_data.columns)
            if labels is not None:
                interactive_text_search_columns.append("label")

        if len(interactive_text_search_columns) == 0:
            warn(
                "interactive_text_search_columns set to True, but no hover_data or labels provided."
                "Please provide hover_data or labels to use interactive text search."
            )

        else:
            
            # Descriptive search
            callback = CustomJS(
                args=dict(
                    s1=data_source,
                    s2=s2,
                    table=table,
                    matching_alpha=interactive_text_search_alpha_contrast,
                    non_matching_alpha=1 - interactive_text_search_alpha_contrast,
                    search_columns=interactive_text_search_columns,
                    default_radius=radius,
                    default_alpha=alpha,

                    s_expr=expr_source,

                    s_hm=hm_cds,
                    cols=x,

                    s_enrich=enrich_cds,
                    s_enrich2=enrich_cds2,
                    s_avg=avg_data_source,
                ),
                code="""
                var d1 = s1.data; // embedding
                var d2 = s2.data; // table
                var d_avg = s_avg.data

                var d_expr = s_expr.data; // expression plot
                var d_hm = s_hm.data; // heatmap
                var d_enrich = s_enrich.data; // enrichment table
                var d_enrich2 = s_enrich2.data // enrichment plot

                var selected_ttherm_id = "";

                var ttids = d_hm['TTHERM_ID'].slice(0, """+str(num_genes)+""");
                const num_cols = cols.length;

                var text_search = cb_obj.value;
                var search_terms = text_search.split(',');

                d2['module'] = []
                d2['ID'] = []


                // JS INITIALIZE

                // EMBEDDING
                // Start by making everything tiny and pale
                d1['alpha'] = Array(d1['ID'].length).fill(0.01)
                d1['line_alpha'] = Array(d1['ID'].length).fill(0.01)
                d1['radius'] = Array(d1['ID'].length).fill(0.0001)

                // TABLE
                d2['ID'] = []
                // d2['YF_ID'] = []

                \n"""+'\n'.join([f"d2['{tc}'] = []" for tc in table_columns])+"""\n
                
                // d2['KEGG_TC'] = []
                // d2['CAZy'] = []
                // d2['BiGG_Reaction'] = []

                // EXPRESSION
                d_expr['TTHERM_ID'] = ['blah']
                d_expr['module'] = ['blah']
                d_expr['ID'] = [['blah']]
                d_expr['expr_xs'] = [['Ll']]
                d_expr['expr_ys'] = [[0]]
                d_expr['alpha'] = [0]
                d_expr['color'] = ['black']
                d_expr['line_dash'] = ['solid']

                // HEATMAP
                d_hm['fill_alpha'] = []
                d_hm['line_alpha'] = []
                
                d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)
                d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)
                
                s_avg.selected.indices = []

                d_avg['alpha'] = Array(d_avg['alpha'].length).fill(default_radius)
                d_avg['radius'] = Array(d_avg['radius'].length).fill(default_radius)
                d_avg['line_color'] = Array(d_avg['line_color'].length).fill("black")

                // JS SEARCH

                var search_columns_dict = {}
                for (var col in search_columns){
                    search_columns_dict[col] = search_columns[col]
                }

                // ENRICHMENT
                s_enrich.selected.indices = []
                s_enrich2.selected.indices = []

                d_enrich2['alpha'] = Array(d_enrich2['alpha'].length).fill(0.3)
                d_enrich2['size'] = Array(d_enrich2['size'].length).fill(7)
                d_enrich2['line_color'] = Array(d_enrich2['line_color'].length).fill("black")


                s1.selected.indices = []

                // Run search
                if (text_search.length > 0){
                    
                    // HEATMAP deselect all
                    d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)
                    d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)

                    // Loop over columns and values
                    // If there is no match for any column for a given row, change the alpha value
                    var string_match = false;
                    for (var i = 0; i < d1.x.length; i++) {
                        string_match = false
                        for (var j in search_columns_dict) {
                            if (search_terms.some(t => String(d1[search_columns_dict[j]][i]).includes(t.trim()))) {
                                string_match = true
                            }
                        }
                        if (string_match){
                            // d1['alpha'][i] = matching_alpha
                            // d1['radius'][i] = 1
                            // d2['YF_ID'].push(d1['YF_ID'][i])

                            // d3['xs'].push(ref_expr['xs'][i])
                            // d3['ys'].push(ref_expr['ys'][i])

                            // So that these points are actually considered selected
                            s1.selected.indices.push(i)

                            // TABLE
                            d2['ID'].push(d1['ID'][i])
                            // d2['YF_ID'].push(d1['YF_ID'][i])

                            \n"""+'\n'.join([f"d2['{tc}'].push(d1['{tc}'][i])" for tc in table_columns])+"""\n

                            // d2['KEGG_TC'].push(d1['KEGG_TC'][i])
                            // d2['CAZy'].push(d1['CAZy'][i])
                            // d2['BiGG_Reaction'].push(d1['BiGG_Reaction'][i])
                            
                            // EMBEDDING
                            // d1['alpha'][i] = 1
                            // d1['line_alpha'][i] = 1
                            // d1['radius'][i] = 100

                            // EXPRESSION
                            d_expr['TTHERM_ID'].push(d1['ID'][i])
                            d_expr['module'].push(d1['module'][i])
                            d_expr['ID'].push(Array(18).fill(d1['ID'][i]))
                            d_expr['expr_xs'].push(d1['expr_xs'][i])
                            d_expr['expr_ys'].push(d1['expr_ys'][i])
                            d_expr['color'].push(d1['color'][i])
                            d_expr['line_dash'].push('solid')
                            // console.log(d_expr)
                            // console.log(i)
                            // console.log(
                            //     d_expr['ID'].length, 
                            //     d_expr['expr_xs'].length, 
                            //     d_expr['expr_ys'].length
                            // )

                            // HEATMAP
                            // selected_ttherm_id = d1['ID'][i];
                            // var match = (element) => element == selected_ttherm_id;
                            var gene_index = i;

                            for (var k = 0; k < num_cols; k++) {
                                d_hm['fill_alpha'][gene_index] = 0.7
                                d_hm['line_alpha'][gene_index] = 0.7

                                gene_index = gene_index + ttids.length
                            }

                        }else{
                            // d1['alpha'][i] = non_matching_alpha
                            // d1['radius'][i] = 0.01
                        }
                    }
                }

                d_expr['alpha'].push.apply(d_expr['alpha'],
                    Array(d2['ID'].length).fill(Math.min(1, Math.max(7/(d2['ID'].length), 0.05)))
                );

                var avg_mods = d_avg['label'].slice(0);
                var selected_mods = d2['module'].slice(0);
                
                var avg_nmod_str = ""
                var avg_nmod = -1
                
                for (let mod of selected_mods){
                    let avg_nmod_str = mod.slice(1);
                    let avg_nmod = +avg_nmod_str;
                    // console.log(avg_nmod);
                    avg_mods.forEach((item, index) => {
                        if (item === avg_nmod) {
                            // console.log("IN");
                            // console.log(index);
                            s_avg.selected.indices.push(index);
                        }
                    });
                }

                if (selected_mods.length > 0 && s_avg.selected.indices.length == 0){
                    d_avg['alpha'] = Array(d_avg['alpha'].length).fill(0.05)
                    // d_avg['radius'] = Array(d_avg['radius'].length).fill(default_radius/20)
                    d_avg['line_color'] = Array(d_avg['line_color'].length).fill(null)
                }


                var enrich_mods = d_enrich['module'].slice(0);

                // console.log("selected_mods")
                // console.log(selected_mods)
                
                var enrich_mod_idx = -1
                var nmod_str = ""
                var nmod = -1

                // console.log(s_enrich.selected.indices);
                
                for (let mod of selected_mods){
                    let nmod_str = mod.slice(1);
                    let nmod = +nmod_str;
                    // console.log(nmod);
                    enrich_mods.forEach((item, index) => {
                        if (item === nmod) {
                            // console.log("IN");
                            // console.log(index);
                            s_enrich.selected.indices.push(index);
                            s_enrich2.selected.indices.push(index);
                        }
                    });
                }
                
                // console.log(s_enrich.selected.indices.length);
                // console.log(s_enrich.selected.indices);

                if (selected_mods.length > 0 && s_enrich2.selected.indices.length == 0){
                    // console.log("NONE");
                    d_enrich2['alpha'] = Array2(d_enrich['alpha'].length).fill(0.05)
                    // d_enrich2['size'] = Array(d_enrich2['size'].length).fill(1)
                    d_enrich2['line_color'] = Array(d_enrich2['line_color'].length).fill(null)
                }

                // console.log(s_enrich.selected.indices);

                s1.change.emit();
                s2.change.emit();
                table.change.emit();


                s_expr.change.emit();
                    
                s_hm.change.emit();

                s_enrich.change.emit();
                s_enrich2.change.emit();

                s_avg.change.emit();

                console.log("RAN search");
                // console.log(s1.selected.indices.length);
                // console.log(s1.selected.indices);

            """,
            )

            # text_input.js_on_change("value", callback)
            text_input.js_on_event(events.ValueSubmit, callback)

            # Functional term search
            callback2 = CustomJS(
                args=dict(
                    s1=data_source,
                    s2=s2,
                    table=table,
                    matching_alpha=interactive_text_search_alpha_contrast,
                    non_matching_alpha=1 - interactive_text_search_alpha_contrast,
                    search_columns=interactive_text_search_columns2,
                    default_radius=radius,
                    default_alpha=alpha,

                    s_expr=expr_source,

                    s_hm=hm_cds,
                    cols=x,
                    s_enrich=enrich_cds,
                    s_enrich2=enrich_cds2,
                    s_avg=avg_data_source,
                ),
                code="""
                var d1 = s1.data; // embedding
                var d2 = s2.data; // table
                var d_avg = s_avg.data

                var d_expr = s_expr.data; // expression plot
                var d_hm = s_hm.data; // heatmap
                var d_enrich = s_enrich.data; // enrichment table
                var d_enrich2 = s_enrich2.data; // enrichment plot

                var selected_ttherm_id = "";

                var ttids = d_hm['TTHERM_ID'].slice(0, """+str(num_genes)+""");
                const num_cols = cols.length;

                var text_search = cb_obj.value;
                var search_terms = text_search.split(',');

                d2['module'] = []
                d2['ID'] = []


                // JS INITIALIZE

                // EMBEDDING
                // Start by making everything tiny and pale
                d1['alpha'] = Array(d1['ID'].length).fill(0.01)
                d1['line_alpha'] = Array(d1['ID'].length).fill(0.01)
                d1['radius'] = Array(d1['ID'].length).fill(0.0001)

                // TABLE
                d2['ID'] = []
                // d2['YF_ID'] = []

                \n"""+'\n'.join([f"d2['{tc}'] = []" for tc in table_columns])+"""\n
                
                // d2['KEGG_TC'] = []
                // d2['CAZy'] = []
                // d2['BiGG_Reaction'] = []

                // EXPRESSION
                d_expr['TTHERM_ID'] = ['blah']
                d_expr['module'] = ['blah']
                d_expr['ID'] = [['blah']]
                d_expr['expr_xs'] = [['Ll']]
                d_expr['expr_ys'] = [[0]]
                d_expr['alpha'] = [0]
                d_expr['color'] = ['black']
                d_expr['line_dash'] = ['solid']

                // HEATMAP
                d_hm['fill_alpha'] = []
                d_hm['line_alpha'] = []
                
                d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)
                d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)
                
                s_avg.selected.indices = []

                d_avg['alpha'] = Array(d_avg['alpha'].length).fill(default_radius)
                d_avg['radius'] = Array(d_avg['radius'].length).fill(default_radius)
                d_avg['line_color'] = Array(d_avg['line_color'].length).fill("black")

                // JS SEARCH

                var search_columns_dict = {}
                for (var col in search_columns){
                    search_columns_dict[col] = search_columns[col]
                }

                // ENRICHMENT  
                s_enrich.selected.indices = []
                s_enrich2.selected.indices = []

                d_enrich2['alpha'] = Array(d_enrich2['alpha'].length).fill(0.3)
                d_enrich2['size'] = Array(d_enrich2['size'].length).fill(7)
                d_enrich2['line_color'] = Array(d_enrich2['line_color'].length).fill("black")


                s1.selected.indices = []

                // Run search
                if (text_search.length > 0){
                    
                    // HEATMAP deselect all
                    d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)
                    d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)

                    // Loop over columns and values
                    // If there is no match for any column for a given row, change the alpha value
                    var string_match = false;
                    for (var i = 0; i < d1.x.length; i++) {
                        string_match = false
                        for (var j in search_columns_dict) {
                            if (search_terms.some(t => String(d1[search_columns_dict[j]][i]).includes(t.trim()))) {
                                string_match = true
                            }
                        }
                        if (string_match){
                            // d1['alpha'][i] = matching_alpha
                            // d1['radius'][i] = 1
                            // d2['YF_ID'].push(d1['YF_ID'][i])

                            // d3['xs'].push(ref_expr['xs'][i])
                            // d3['ys'].push(ref_expr['ys'][i])

                            // So that these points are actually considered selected
                            s1.selected.indices.push(i)

                            // TABLE
                            d2['ID'].push(d1['ID'][i])
                            // d2['YF_ID'].push(d1['YF_ID'][i])

                            \n"""+'\n'.join([f"d2['{tc}'].push(d1['{tc}'][i])" for tc in table_columns])+"""\n

                            // d2['KEGG_TC'].push(d1['KEGG_TC'][i])
                            // d2['CAZy'].push(d1['CAZy'][i])
                            // d2['BiGG_Reaction'].push(d1['BiGG_Reaction'][i])
                            
                            // EMBEDDING
                            // d1['alpha'][i] = 1
                            // d1['line_alpha'][i] = 1
                            // d1['radius'][i] = 100

                            // EXPRESSION
                            d_expr['TTHERM_ID'].push(d1['ID'][i])
                            d_expr['module'].push(d1['module'][i])
                            d_expr['ID'].push(Array(18).fill(d1['ID'][i]))
                            d_expr['expr_xs'].push(d1['expr_xs'][i])
                            d_expr['expr_ys'].push(d1['expr_ys'][i])
                            d_expr['color'].push(d1['color'][i])
                            d_expr['line_dash'].push('solid')
                            // console.log(d_expr)
                            // console.log(i)
                            // console.log(
                            //     d_expr['ID'].length, 
                            //     d_expr['expr_xs'].length, 
                            //     d_expr['expr_ys'].length
                            // )

                            // HEATMAP
                            // selected_ttherm_id = d1['ID'][i];
                            // var match = (element) => element == selected_ttherm_id;
                            var gene_index = i;

                            for (var k = 0; k < num_cols; k++) {
                                d_hm['fill_alpha'][gene_index] = 0.7
                                d_hm['line_alpha'][gene_index] = 0.7

                                gene_index = gene_index + ttids.length
                            }

                        }else{
                            // d1['alpha'][i] = non_matching_alpha
                            // d1['radius'][i] = 0.01
                        }
                    }
                }

                d_expr['alpha'].push.apply(d_expr['alpha'],
                    Array(d2['ID'].length).fill(Math.min(1, Math.max(7/(d2['ID'].length), 0.05)))
                );

                var avg_mods = d_avg['label'].slice(0);
                var selected_mods = d2['module'].slice(0);
                
                var avg_nmod_str = ""
                var avg_nmod = -1
                
                for (let mod of selected_mods){
                    let avg_nmod_str = mod.slice(1);
                    let avg_nmod = +avg_nmod_str;
                    // console.log(avg_nmod);
                    avg_mods.forEach((item, index) => {
                        if (item === avg_nmod) {
                            // console.log("IN");
                            // console.log(index);
                            s_avg.selected.indices.push(index);
                        }
                    });
                }

                if (selected_mods.length > 0 && s_avg.selected.indices.length == 0){
                    d_avg['alpha'] = Array(d_avg['alpha'].length).fill(0.05)
                    // d_avg['radius'] = Array(d_avg['radius'].length).fill(default_radius/20)
                    d_avg['line_color'] = Array(d_avg['line_color'].length).fill(null)
                }


                var enrich_mods = d_enrich['module'].slice(0);

                // console.log("selected_mods")
                // console.log(selected_mods)
                
                var enrich_mod_idx = -1
                var nmod_str = ""
                var nmod = -1

                // console.log(s_enrich.selected.indices);
                
                for (let mod of selected_mods){
                    let nmod_str = mod.slice(1);
                    let nmod = +nmod_str;
                    // console.log(nmod);
                    enrich_mods.forEach((item, index) => {
                        if (item === nmod) {
                            // console.log("IN");
                            // console.log(index);
                            s_enrich.selected.indices.push(index);
                            s_enrich2.selected.indices.push(index);
                        }
                    });
                }
                
                // console.log(s_enrich.selected.indices.length);
                // console.log(s_enrich.selected.indices);

                if (selected_mods.length > 0 && s_enrich2.selected.indices.length == 0){
                    // console.log("NONE");
                    d_enrich2['alpha'] = Array(d_enrich2['alpha'].length).fill(0.05)
                    // d_enrich2['size'] = Array(d_enrich2['size'].length).fill(1)
                    d_enrich2['line_color'] = Array(d_enrich2['line_color'].length).fill(null)
                }

                // console.log(s_enrich.selected.indices);

                s1.change.emit();
                s2.change.emit();
                table.change.emit();


                s_expr.change.emit();
                    
                s_hm.change.emit();

                s_enrich.change.emit();
                s_enrich2.change.emit();

                s_avg.change.emit();

                console.log("RAN search");
                // console.log(s1.selected.indices.length);
                // console.log(s1.selected.indices);

            """,
            )

            # text_input.js_on_change("value", callback)
            text_input2.js_on_event(events.ValueSubmit, callback2)

    module_list = list(hover_data['module'].values)
    sorted_module_list = sorted(module_list)
    max_label_num_str = (sorted_module_list[len(module_list) - 1]).replace('m', '')
    max_label_num_len = len(max_label_num_str)

    spinner = Spinner(title="", low=0, high=int(max_label_num_str), step=1, value=0, sizing_mode=search_sizing_mode)

    interactive_text_search_columns_spinner = []
    if labels is not None:
        interactive_text_search_columns_spinner.append('module')

    if len(interactive_text_search_columns_spinner) == 0:
        warn(
            "interactive_text_search_columns_spinner set to True, but no hover_data or labels provided."
            "Please provide hover_data or labels to use interactive text search."
        )

    else:
        pass


    # Lifted from https://stackoverflow.com/questions/31824124/is-there-a-way-to-save-bokeh-data-table-content   + "\t" + 
    download_button1 = Button(label='💾 Data Table', button_type="success")
    download_button1.js_on_click(
        CustomJS(
            args=dict(source_data=data_source),
            code="""
            var inds = source_data.selected.indices;
            var data = source_data.data;
            """
            +
            'var out = "TTHERM_ID\t'+'\t'.join([tc for tc in table_columns])+'\\n";'
            +
            """
            for (var i = 0; i < inds.length; i++) {
            """
                f"""\tout += data['ID'][inds[i]] + "\t" + """+ ' + "\t" + '.join([f"data['{tc}'][inds[i]]" for tc in table_columns]) +""" + "\\n";"""
            """
            }
            var file = new Blob([out], {type: 'text/plain'});
            var elem = window.document.createElement('a');
            elem.href = window.URL.createObjectURL(file);
            elem.download = 'selected-annotation-data.tsv';
            document.body.appendChild(elem);
            elem.click();
            document.body.removeChild(elem);
            """))  
    
    download_button2 = Button(label='💾 Enrichment', button_type="success")
    download_button2.js_on_click(
        CustomJS(
            args=dict(source_data=enrich_cds),
            code="""
            var inds = source_data.selected.indices;
            var data = source_data.data;

            var out = "module\tterm\tinfo\tfold_change\tbonferroni\\n";
            for (var i = 0; i < inds.length; i++) {
                out += data['module'][inds[i]] + "\t" + data['term'][inds[i]] + "\t" + data['info'][inds[i]] + "\t" + data['fold_change'][inds[i]] + "\t" + data['bonferroni'][inds[i]] + "\\n";
            }
            var file = new Blob([out], {type: 'text/plain'});
            var elem = window.document.createElement('a');
            elem.href = window.URL.createObjectURL(file);
            elem.download = 'enrichment-data.tsv';
            document.body.appendChild(elem);
            elem.click();
            document.body.removeChild(elem);
            """))  
    
    
    
    
    if interactive_text_search:


        # FIXME: IMPLEMENT
        # text_search2 = TextInput(value="ENRICHMENT TERM SEARCH (UNDER CONSTRUCTION)", sizing_mode='stretch_width')

        # # FIXME: IMPLEMENT
        # data_module_stats = {'Module': [1, 2, 3, 4],
        #                     'Stats': ['A', 'B', 'C', 'D']}
        # source_module_stats = ColumnDataSource(data=data_module_stats)
        # columns_module_stats = [TableColumn(field="Module", title="Module"),
        #                         TableColumn(field="Stats", title="Stats")]
        # module_stats_table = DataTable(source=source_module_stats, columns=columns_module_stats)


        # spinner_text_html = """
        # <div style="display: flex; justify-content: center; align-items: center; height: 100%;">
        #     <p style="text-align: center;">MODULE #:</p>
        # </div>
        # """
        # spinner_custom_text = Div(text=spinner_text_html)

        cola_sizing_mode = 'stretch_width'
        colb_sizing_mode = 'stretch_width'

        plot_tabs = Tabs(tabs=[TabPanel(child=plot_avg, title='UMAP of clusters'), TabPanel(child=plot, title='UMAP of genes')], sizing_mode=plot_sizing_mode)

        col1a = column(hm, height=800)
        col1a.sizing_mode = cola_sizing_mode
        col2a = column(enrich_p, height=800)
        col2a.sizing_mode = cola_sizing_mode
        # col3a = column(plot, expr_fig, height=800)
        col3a = column(
            plot_tabs,
            expr_fig, height=800) 
        col3a.sizing_mode = cola_sizing_mode

        # col1b = column(module_stats_table, height=800, max_width=450)
        # col1b.sizing_mode = colb_sizing_mode
        col2b = column(table, height=800, max_width=100000)
        col2b.sizing_mode = colb_sizing_mode

        rows_sizing_mode = 'stretch_width'

        row_search = row(text_input, text_input2, download_button1, download_button2, sizing_mode='stretch_width')
        rowa = row(row(col1a, col2a, sizing_mode=rows_sizing_mode), col3a)
        rowa.sizing_mode = rows_sizing_mode
        rowb = row(col2b)
        rowb.sizing_mode = rows_sizing_mode

        plot = column(row_search, rowa, rowb)
        plot.sizing_mode = 'stretch_width'

    else:
        plot = column(row(column(plot, expr_fig), hm, enrich_p), row(download_button1, download_button2), table)


    a_s_callback = CustomJS(args=dict(
                s1=data_source,
                s2=s2,
                table=table,
                matching_alpha=interactive_text_search_alpha_contrast,
                non_matching_alpha=1 - interactive_text_search_alpha_contrast,
                default_radius=radius,
                default_alpha=alpha,

                s_expr=expr_source,

                s_hm=hm_cds,
                cols=x,

                s_enrich=enrich_cds,
                s_enrich2=enrich_cds2,

                s_avg=avg_data_source,

                plot_tabs=plot_tabs,
                ), code="""
            // console.log(plot_tabs.active);
            if (plot_tabs.active == 0){
            var avg_idxs = cb_obj.indices;
            var d1 = s1.data; // embedding
            var d2 = s2.data; // table
            var d_avg = s_avg.data

            var d_expr = s_expr.data; // expression plot
            var d_hm = s_hm.data; // heatmap
            var d_enrich = s_enrich.data; // enrichment table
            var d_enrich2 = s_enrich2.data; // enrichment plot

            var selected_ttherm_id = "";

            var ttids = d_hm['TTHERM_ID'].slice(0, """+str(num_genes)+""");
            const num_cols = cols.length;

            d2['module'] = []
            d2['ID'] = []


            // JS INITIALIZE

            // EMBEDDING
            // Start by making everything tiny and pale
            // d1['alpha'] = Array(d1['ID'].length).fill(0.01)
            // d1['line_alpha'] = Array(d1['ID'].length).fill(0.01)
            // d1['radius'] = Array(d1['ID'].length).fill(0.0001)

            // TABLE
            d2['ID'] = []
            // d2['YF_ID'] = []
            \n"""+'\n'.join([f"d2['{tc}'] = []" for tc in table_columns])+"""\n
            // d2['KEGG_TC'] = []
            // d2['CAZy'] = []
            // d2['BiGG_Reaction'] = []

            // EXPRESSION
            d_expr['TTHERM_ID'] = ['blah']
            d_expr['module'] = ['blah']
            d_expr['ID'] = [['blah']]
            d_expr['expr_xs'] = [['Ll']]
            d_expr['expr_ys'] = [[0]]
            d_expr['alpha'] = [0]
            d_expr['color'] = ['black']
            d_expr['line_dash'] = ['solid']
            

            // HEATMAP
            d_hm['fill_alpha'] = []
            d_hm['line_alpha'] = []
            
            d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)
            d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)


            // ENRICHMENT
            s_enrich.selected.indices = []
            s_enrich2.selected.indices = []

            d_enrich2['alpha'] = Array(d_enrich2['alpha'].length).fill(0.3)
            d_enrich2['size'] = Array(d_enrich2['size'].length).fill(7)
            d_enrich2['line_color'] = Array(d_enrich2['line_color'].length).fill("black")

            // HEATMAP deselect all
            d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)
            d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)


            s1.selected.indices = []

            var avg_selected_mods = (avg_idxs.map(index => d_avg['label'][index]))

            // JS
            for (var i = 0; i < d1.x.length; i++) {
            
                var mod_at_i = +((d1['module'][i]).slice(1))

                if (avg_selected_mods.includes(mod_at_i)) { 
                
                    // console.log(i);
                    // console.log(d1['module'][i]);
                    // console.log(d1['ID'][i]);
                    // d1['alpha'][i] = matching_alpha
                    // d1['radius'][i] = 1
                    // d2['YF_ID'].push(d1['YF_ID'][i])

                    // d3['xs'].push(ref_expr['xs'][i])
                    // d3['ys'].push(ref_expr['ys'][i])

                    s1.selected.indices.push(i)

                    // TABLE
                    d2['ID'].push(d1['ID'][i])
                    // d2['YF_ID'].push(d1['YF_ID'][i])
                    \n"""+'\n'.join([f"d2['{tc}'].push(d1['{tc}'][i])" for tc in table_columns])+"""\n
                    // d2['KEGG_TC'].push(d1['KEGG_TC'][i])
                    // d2['CAZy'].push(d1['CAZy'][i])
                    // d2['BiGG_Reaction'].push(d1['BiGG_Reaction'][i])
                    
                    // EMBEDDING
                    // d1['alpha'][i] = 1
                    // d1['line_alpha'][i] = 1
                    // d1['radius'][i] = 100

                    // EXPRESSION
                    d_expr['TTHERM_ID'].push(d1['ID'][i])
                    d_expr['module'].push(d1['module'][i])
                    d_expr['ID'].push(Array(18).fill(d1['ID'][i]))
                    d_expr['expr_xs'].push(d1['expr_xs'][i])
                    d_expr['expr_ys'].push(d1['expr_ys'][i])
                    d_expr['color'].push(d1['color'][i])
                    d_expr['line_dash'].push('solid')
                    // console.log(d_expr)
                    // console.log(i)
                    // console.log(
                    //     d_expr['ID'].length, 
                    //     d_expr['expr_xs'].length, 
                    //     d_expr['expr_ys'].length
                    // )

                    // HEATMAP
                    // selected_ttherm_id = d1['ID'][i]
                    // var match = (element) => element == selected_ttherm_id
                    var gene_index = i;

                    for (var k = 0; k < num_cols; k++) {
                        d_hm['fill_alpha'][gene_index] = 0.7
                        d_hm['line_alpha'][gene_index] = 0.7

                        gene_index = gene_index + ttids.length
                    }

                }else{
                    // d1['alpha'][i] = non_matching_alpha
                    // d1['radius'][i] = 0.01
                }
            }

            d_expr['alpha'].push.apply(d_expr['alpha'],
                Array(d2['ID'].length).fill(Math.min(1, Math.max(7/(d2['ID'].length), 0.05)))
            );
            
            var enrich_mods = d_enrich['module'].slice(0);
            var selected_mods = d2['module'].slice(0);
            
            var nmod_str = ""
            var nmod = -1
            
            for (let mod of selected_mods){
                let nmod_str = mod.slice(1);
                let nmod = +nmod_str;
                // console.log(nmod);
                enrich_mods.forEach((item, index) => {
                    if (item === nmod) {
                        s_enrich.selected.indices.push(index);
                        s_enrich2.selected.indices.push(index);
                    }
                });
            }

            if (selected_mods.length > 0 && s_enrich2.selected.indices.length == 0){
                d_enrich2['alpha'] = Array(d_enrich2['alpha'].length).fill(0.05)
                // d_enrich2['size'] = Array(d_enrich2['size'].length).fill(1)
                d_enrich2['line_color'] = Array(d_enrich2['line_color'].length).fill(null)
            }


            s1.change.emit();
            s2.change.emit();
            table.change.emit();


            s_expr.change.emit();
                
            s_hm.change.emit();

            s_enrich.change.emit();
            s_enrich2.change.emit();

            // s_avg.change.emit();

            console.log("RAN AVG selection");
            // console.log(idxs.length);
            // console.log(s1.selected.indices);


            // console.log(avg_selected_mods);
            }
            """
        )        

    avg_data_source.selected.js_on_change('indices', a_s_callback)





    s_callback = CustomJS(args=dict(
                s1=data_source,
                s2=s2,
                table=table,
                matching_alpha=interactive_text_search_alpha_contrast,
                non_matching_alpha=1 - interactive_text_search_alpha_contrast,
                default_radius=radius,
                default_alpha=alpha,

                s_expr=expr_source,

                s_hm=hm_cds,
                cols=x,
                s_enrich=enrich_cds,
                s_enrich2=enrich_cds2,

                s_avg=avg_data_source,
                plot_tabs=plot_tabs,
                ), code="""
            if (plot_tabs.active == 1){
            var idxs = cb_obj.indices;
            var d1 = s1.data; // embedding
            var d2 = s2.data; // table
            var d_avg = s_avg.data

            var d_expr = s_expr.data; // expression plot
            var d_hm = s_hm.data; // heatmap
            var d_enrich = s_enrich.data; // enrichment table
            var d_enrich2 = s_enrich2.data; // enrichment plot

            var selected_ttherm_id = "";

            var ttids = d_hm['TTHERM_ID'].slice(0, """+str(num_genes)+""");
            const num_cols = cols.length;

            d2['module'] = []
            d2['ID'] = []


            // JS INITIALIZE

            // EMBEDDING
            // Start by making everything tiny and pale
            // d1['alpha'] = Array(d1['ID'].length).fill(0.01)
            // d1['line_alpha'] = Array(d1['ID'].length).fill(0.01)
            // d1['radius'] = Array(d1['ID'].length).fill(0.0001)

            // TABLE
            d2['ID'] = []
            // d2['YF_ID'] = []
            \n"""+'\n'.join([f"d2['{tc}'] = []" for tc in table_columns])+"""\n
            // d2['KEGG_TC'] = []
            // d2['CAZy'] = []
            // d2['BiGG_Reaction'] = []

            // EXPRESSION
            d_expr['TTHERM_ID'] = ['blah']
            d_expr['module'] = ['blah']
            d_expr['ID'] = [['blah']]
            d_expr['expr_xs'] = [['Ll']]
            d_expr['expr_ys'] = [[0]]
            d_expr['alpha'] = [0]
            d_expr['color'] = ['black']
            d_expr['line_dash'] = ['solid']
            

            // HEATMAP
            d_hm['fill_alpha'] = []
            d_hm['line_alpha'] = []
            
            d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)
            d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)


            // ENRICHMENT
            s_enrich.selected.indices = []
            s_enrich2.selected.indices = []

            d_enrich2['alpha'] = Array(d_enrich2['alpha'].length).fill(0.3)
            d_enrich2['size'] = Array(d_enrich2['size'].length).fill(7)
            d_enrich2['line_color'] = Array(d_enrich2['line_color'].length).fill("black")

            // HEATMAP deselect all
            d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)
            d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)

            
            s_avg.selected.indices = []

            d_avg['alpha'] = Array(d_avg['alpha'].length).fill(default_radius)
            d_avg['radius'] = Array(d_avg['radius'].length).fill(default_radius)
            d_avg['line_color'] = Array(d_avg['line_color'].length).fill("black")

            // JS 
            for (var i = 0; i < d1.x.length; i++) {
                if (idxs.includes(i)) {
                    // console.log(i);
                    // console.log(d1['module'][i]);
                    // console.log(d1['ID'][i]);
                    // d1['alpha'][i] = matching_alpha
                    // d1['radius'][i] = 1
                    // d2['YF_ID'].push(d1['YF_ID'][i])

                    // d3['xs'].push(ref_expr['xs'][i])
                    // d3['ys'].push(ref_expr['ys'][i])

                    // TABLE
                    d2['ID'].push(d1['ID'][i])
                    // d2['YF_ID'].push(d1['YF_ID'][i])
                    \n"""+'\n'.join([f"d2['{tc}'].push(d1['{tc}'][i])" for tc in table_columns])+"""\n
                    // d2['KEGG_TC'].push(d1['KEGG_TC'][i])
                    // d2['CAZy'].push(d1['CAZy'][i])
                    // d2['BiGG_Reaction'].push(d1['BiGG_Reaction'][i])
                    
                    // EMBEDDING
                    // d1['alpha'][i] = 1
                    // d1['line_alpha'][i] = 1
                    // d1['radius'][i] = 100

                    // EXPRESSION
                    d_expr['TTHERM_ID'].push(d1['ID'][i])
                    d_expr['module'].push(d1['module'][i])
                    d_expr['ID'].push(Array(18).fill(d1['ID'][i]))
                    d_expr['expr_xs'].push(d1['expr_xs'][i])
                    d_expr['expr_ys'].push(d1['expr_ys'][i])
                    d_expr['color'].push(d1['color'][i])
                    d_expr['line_dash'].push('solid')
                    // console.log(d_expr)
                    // console.log(i)
                    // console.log(
                    //     d_expr['ID'].length, 
                    //     d_expr['expr_xs'].length, 
                    //     d_expr['expr_ys'].length
                    // )

                    // HEATMAP
                    // selected_ttherm_id = d1['ID'][i]
                    // var match = (element) => element == selected_ttherm_id
                    var gene_index = i;

                    for (var k = 0; k < num_cols; k++) {
                        d_hm['fill_alpha'][gene_index] = 0.7
                        d_hm['line_alpha'][gene_index] = 0.7

                        gene_index = gene_index + ttids.length
                    }

                }else{
                    // d1['alpha'][i] = non_matching_alpha
                    // d1['radius'][i] = 0.01
                }
            }

            d_expr['alpha'].push.apply(d_expr['alpha'],
                Array(idxs.length).fill(Math.min(1, Math.max(7/(idxs.length), 0.05)))
            );


            var avg_mods = d_avg['label'].slice(0);
            var selected_mods = d2['module'].slice(0);
            
            var avg_nmod_str = ""
            var avg_nmod = -1
            
            for (let mod of selected_mods){
                let avg_nmod_str = mod.slice(1);
                let avg_nmod = +avg_nmod_str;
                // console.log(avg_nmod);
                avg_mods.forEach((item, index) => {
                    if (item === avg_nmod) {
                        // console.log("IN");
                        // console.log(index);
                        s_avg.selected.indices.push(index);
                    }
                });
            }

            if (selected_mods.length > 0 && s_avg.selected.indices.length == 0){
                d_avg['alpha'] = Array(d_avg['alpha'].length).fill(0.05)
                d_avg['radius'] = Array(d_avg['radius'].length).fill(default_radius/20)
                d_avg['line_color'] = Array(d_avg['line_color'].length).fill(null)
            }
            
            var enrich_mods = d_enrich['module'].slice(0);
            
            var nmod_str = ""
            var nmod = -1
            
            for (let mod of selected_mods){
                let nmod_str = mod.slice(1);
                let nmod = +nmod_str;
                // console.log(nmod);
                enrich_mods.forEach((item, index) => {
                    if (item === nmod) {
                        s_enrich.selected.indices.push(index);
                        s_enrich2.selected.indices.push(index);
                    }
                });
            }

            if (selected_mods.length > 0 && s_enrich2.selected.indices.length == 0){
                d_enrich2['alpha'] = Array(d_enrich2['alpha'].length).fill(0.05)
                // d_enrich2['size'] = Array(d_enrich2['size'].length).fill(1)
                d_enrich2['line_color'] = Array(d_enrich2['line_color'].length).fill(null)
            }


            s1.change.emit();
            s2.change.emit();
            table.change.emit();


            s_expr.change.emit();
                
            s_hm.change.emit();

            s_enrich.change.emit();
            s_enrich2.change.emit();

            s_avg.change.emit();

            console.log("RAN selection");
            // console.log(idxs.length);
            // console.log(s1.selected.indices);
            
            }
            """
        )

    data_source.selected.js_on_change('indices', s_callback)

    return plot

def get_centroid(module_df):
    
    # get rid of ttherm_ids
    data_cols = [c for c in module_df.columns if ('TTHERM' not in c) and ('label' not in c)]
    data = module_df[data_cols]
    
    centroid = data.apply(np.mean, axis=0).values
    
    return centroid

def get_module_centroid_df(expr_df, cluster_label_df):

    # print(expr_df.columns)
    # print(cluster_label_df.columns)
    
    merge = expr_df.merge(cluster_label_df, on='TTHERM_ID')
    
    grouped = merge.groupby('label')
    
    centroid_rows = []
    
    for label, grp_df in grouped:
        # print(grp_df.head())
        centroid = get_centroid(grp_df)
        centroid_rows.append(centroid)
    
    data_cols = [c for c in merge.columns if ('TTHERM' not in c) and ('label' not in c)]
    
    centroid_df = pd.DataFrame(centroid_rows)
    centroid_df.columns = data_cols
    centroid_df.index.rename('module', inplace=True)
        
    return centroid_df

def get_all_module_centroids(expr_df, cluster_label_df):
    
    merge = expr_df.merge(cluster_label_df, on='TTHERM_ID')
    
    grouped = merge.groupby('label')
    
    module_centroid_list = []
    
    for label, grp_df in grouped:
        
        centroid = get_centroid(grp_df)
        module_centroid_list.append( (label, centroid) )
        
    return module_centroid_list

def arrange_modules(expr_df, cluster_label_df, phases):
    
    if phases == 'full':
        
        x = ['Ll', 
             'Lm', 
             'Lh', 
             'S0', 
             'S3', 
             'S6', 
             'S9', 
             # 'S12',
             'S15', 
             'S24', 
             'C0', 
             # 'C2', 
             'C4', 
             'C6', 
             'C8', 
             'C10', 
             'C12', 
             'C14', 
             'C16', 
             'C18']        
        
    elif phases == 'veg':
        
        x = ['Ll', 
             'Lm', 
             'Lh', 
             'S0', 
             'S3', 
             'S6', 
             'S9', 
             # 'S12', 
             'S15', 
             'S24']
                
    elif phases == 'sex':
        
        x = ['C0', 
             # 'C2', 
             'C4', 
             'C6', 
             'C8', 
             'C10',
             'C12',
             'C14', 
             'C16', 
             'C18']
        
    elif phases == 'rna_seq':
        x = ['000min', '030min', '060min', '090min', '120min', '150min',
       '180min', '210min', '240min']
                
    elif phases == 'rna_seq_single_cycle':
        x = ['000min', '030min', '060min', '090min', '120min', '150min',
       '180min']
                
    else:
        raise(ValueError(f'\"{phases}\" is an invalid choice for parameter \"phases.\"'))
    
    cols = ['TTHERM_ID'] + [c for c in expr_df.columns if c.split('_')[0] in x]
        
    module_centroid_df = get_module_centroid_df(expr_df[cols], cluster_label_df)

    # print(module_centroid_df.shape)
    # print(module_centroid_df.columns)

    # print(np.where(np.isfinite(module_centroid_df) == True))

    # non_finite_indices = np.where(~np.isfinite(module_centroid_df))

    # print(non_finite_indices)

    # print(np.all(np.isfinite(module_centroid_df)))

    # nan_values = module_centroid_df.isna().sum()
    # print("NaN values in each column:\n", nan_values)

    # inf_values = np.isinf(module_centroid_df).sum()
    # print("Infinite values in each column:\n", inf_values)

    
    linkage = scipy.cluster.hierarchy.linkage(module_centroid_df, method='average', metric='correlation', optimal_ordering=True)
    r_cophcorre, ar_copdist = scipy.cluster.hierarchy.cophenet(linkage, scipy.spatial.distance.pdist(module_centroid_df, metric='correlation'))
    
    # print(f'The Copheretic correlation is: {r_cophcorre}')
    
    d_dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True)
    cat_sorted = list(module_centroid_df.iloc[d_dendro['leaves'],:].index)
    
    sorter_index = dict(zip(cat_sorted, range(len(cat_sorted))))
    
    reassigned_df = cluster_label_df.copy(deep=True)
    
    
    
    reassigned_df['label'] = reassigned_df['label'].map(sorter_index)
    # print(len(reassigned_df))
    
    arranged_dfs = []
    
    for cat in cat_sorted:
        
        mini_df = reassigned_df.loc[reassigned_df['label'] == cat]
        # gene_count += len(mini_df)

        arranged_dfs.append(mini_df)
        
#     gene_count = 0
    
#     for mdf in arranged_dfs:
#         gene_count += len(mdf)
    
#     print(gene_count)
        
    arranged_df = pd.concat(arranged_dfs)
    
    
    return arranged_df


def plot_embedding(expression_df, enrich_df, embedding_df, annotation_df, label_df, phases, palette, n_components=2, n_neighbors=15, title=None, random_state=42, radius=0.01, expr_min=0, expr_max=1, yf_to_ttherm_map_df=None, avg_df=None, avg_radius=None):
    
    """
    Function to plot the UMAP of expression data.
    
    z : Bool
        Whether the normalization is by z-score
    """

    if expression_df.shape[0] != embedding_df.shape[0]:
        raise(ValueError(f'The number of rows (genes) in the input dataframes are not equal.\nexpression_df: {expression_df.shape[0]} rows\nembedding_df: {expression_df.shape[0]} rows'))

    # if phases in ['full', 'veg', 'sex']:
    #     mean_expression_df = get_arith_mean_expression(expression_df)

    # elif phases == 'rna_seq':
    #     mean_expression_df = ari_mean_df_of_duplicates(expression_df)

    mean_expression_df = expression_df

    num_genes = mean_expression_df.shape[0]
    
    embedding_df['TTHERM_ID'] = mean_expression_df['TTHERM_ID'].values
    
    merge_unsorted = mean_expression_df.merge(embedding_df, on='TTHERM_ID')

    merge_all = label_df.merge(merge_unsorted, on='TTHERM_ID')

    merge_all_sorted = merge_all

    labels = merge_all_sorted['label'].values

    merge = merge_all_sorted.loc[: , merge_all_sorted.columns]
    
    # take part of annotation df that shared TTHERM_IDs with expression df
    relevant_annot = annotation_df.iloc[np.in1d(annotation_df['TTHERM_ID'].values, merge['TTHERM_ID'].values)]

    if relevant_annot.shape[0] < merge.shape[0]:
           warnings.warn(f'The number of rows (genes) in the relevant annotation dataframe and the merged input dataframe are not equal. If you are using a new genome model, please rerun eggnog to ensure all annotations are populated.\nrelevant_annot: {relevant_annot.shape[0]} rows\nmerge: {merge.shape[0]} rows')

    merge = merge.merge(relevant_annot, on='TTHERM_ID', how='left')

    ttherm_ids = merge['TTHERM_ID'].values
    
    if phases == 'full':
        
        x = ['Ll', 
             'Lm', 
             'Lh', 
             'S0', 
             'S3', 
             'S6', 
             'S9', 
             # 'S12',
             'S15', 
             'S24', 
             'C0', 
             # 'C2', 
             'C4', 
             'C6', 
             'C8', 
             'C10', 
             'C12', 
             'C14', 
             'C16', 
             'C18']
        
        
    elif phases == 'veg':
        
        x = ['Ll', 
             'Lm', 
             'Lh', 
             'S0', 
             'S3', 
             'S6', 
             'S9', 
             # 'S12', 
             'S15', 
             'S24']
        
    elif phases == 'sex':
        
        x = ['C0', 
             # 'C2', 
             'C4', 
             'C6', 
             'C8', 
             'C10',
             'C12',
             'C14', 
             'C16', 
             'C18']
        
    elif phases == 'rna_seq':
        x = ['000min', '030min', '060min', '090min', '120min', '150min',
       '180min', '210min', '240min']
        
    elif phases == 'rna_seq_single_cycle':
        x = ['000min', '030min', '060min', '090min', '120min', '150min', '180min']
        
    else:
        raise(ValueError(f'\"{phases}\" is an invalid choice for parameter \"phases.\"'))
    
    if yf_to_ttherm_map_df is not None:
        merge = merge.merge(yf_to_ttherm_map_df, on='TTHERM_ID', how='left')

    
    # print(len(merge['TTHERM_ID'].values))
    # print(len([f'm{int(l):04d}' for l in labels]))
    # print(merge['TTHERM_ID'].values)
    # print([f'm{int(l):04d}' for l in labels])

    if phases == 'rna_seq':
        rna_seq_phase_dict = {
            '000min': '(000min) G1',
            '030min': '(030min) G1',
            '060min': '(060min) S',
            '090min': '(090min) S',
            '120min': '(120min) G2',
            '150min': '(150min) M',
            '180min': '(180min) M',
            '210min': '(210min) G1',
            '240min': '(240min) S',
        }

        merge.rename(columns=rna_seq_phase_dict, inplace=True)
        # print(merge.columns)

        x = [rna_seq_phase_dict[t] for t in x]
        
    xs = [x for ttid in ttherm_ids]


    ys = [merge.loc[merge['TTHERM_ID'] == ttid, x].values[0] for ttid in ttherm_ids]

    merge['expr_xs'] = xs
    merge['expr_ys'] = ys

    len(str(int(max(labels))))
    
#     pdb.set_trace()
    hover_data = pd.DataFrame({
                               # 'index':np.arange(len(data)),
                               'ID':merge['TTHERM_ID'].values,
                               'module':[f'm{str(int(l)).zfill(len(str(int(max(labels)))))}' for l in labels]})
    
    # print(x)
    # print([rna_seq_phase_dict[t] for t in x])
            
    p = interactive(merge,
                    enrich_df,
                    num_genes,
                    x,
                    # mean_expression_df,
                    title=title,
                    hover_data=hover_data, 
                    labels=labels, 
                    color_key=palette, 
#                     color_key_cmap='Paired',
                    background='white', 
                    radius=radius,
                    alpha=0.3,
#                     width=600, 
#                     height=500,
                    interactive_text_search=True,
                    expr_min=expr_min,
                    expr_max=expr_max,
                    avg_df=avg_df,
                    avg_radius=avg_radius,
                   )
    
    #p.children[1].title = title
    
    return p


# ((range of x values)^2 * (range of y values)^2)^(0.5) * (const) = (optimal radius value)
def compute_2d_embedding_point_radius(embedding_df, const=339.30587926495537):
    return ((((max(embedding_df['x'].values) - min(embedding_df['x'].values))**2) + ((max(embedding_df['y'].values) - min(embedding_df['y'].values))**2))**(0.5)) / const

def generate_umap(expression_df, enrich_df, annotation_df, label_df, phase, palette, title, n_neighbors=5, n_components=2, random_state=42, expr_min=0, expr_max=1, embedding_metric='euclidean', yf_to_ttherm_map_df=None, avg_df=None):
       
    data = expression_df[list(expression_df.columns)[1:]].values
    
    umap_mapper = umap.UMAP(random_state=random_state, n_components=n_components, n_neighbors=n_neighbors, metric=embedding_metric).fit(data)
    embedding = _get_umap_embedding(umap_mapper)
    
    umap_df = pd.DataFrame(np.array(embedding), columns=('x', 'y'))

    radius = compute_2d_embedding_point_radius(umap_df)

    avg_data = avg_df[list(avg_df.columns)[1:]].values

    avg_umap_mapper = umap.UMAP(random_state=random_state, n_components=n_components, n_neighbors=n_neighbors, metric=embedding_metric).fit(avg_data)
    avg_embedding = _get_umap_embedding(avg_umap_mapper)

    avg_umap_df = pd.DataFrame(np.array(avg_embedding), columns=('x', 'y'))

    avg_umap_df['label'] = avg_df['label'].values

    avg_umap_df['num_genes'] = [np.count_nonzero(label_df['label'].values == m) for m in avg_umap_df['label'].values]

    avg_radius = compute_2d_embedding_point_radius(avg_umap_df)
    
    p = plot_embedding(expression_df, enrich_df, umap_df, annotation_df, label_df, phase, palette, title=title, n_neighbors=n_neighbors, radius=radius, expr_min=expr_min, expr_max=expr_max, yf_to_ttherm_map_df=yf_to_ttherm_map_df, avg_df=avg_umap_df, avg_radius=avg_radius)

    return p


def generate_and_save_umap(outfile_name, expression_df, enrich_df, annotation_df, label_df, phase, palette, title, n_neighbors=5, n_components=2, random_state=42, expr_min=0, expr_max=1, embedding_metric='euclidean', yf_to_ttherm_map_df=None, avg_df=None):
    
    p = generate_umap(expression_df, enrich_df, annotation_df, label_df, phase, palette, title, n_neighbors=n_neighbors, n_components=n_components, random_state=random_state, expr_min=expr_min, expr_max=expr_max, embedding_metric=embedding_metric, yf_to_ttherm_map_df=yf_to_ttherm_map_df, avg_df=avg_df)
    
    bokeh.plotting.output_file(filename=outfile_name, title=title, mode='inline')
    bokeh.plotting.save(p)
    print(outfile_name)
    return p    

def generate_and_save_umap_tabbed(outfile_name: str, expression_dfs: list, tab_labels: list, enrich_dfs: list, annotation_df: pd.DataFrame, label_dfs: list, phase, palettes, title, n_neighbors=5, n_components=2, random_state=42, expr_mins=[], expr_maxs=[], embedding_metric='euclidean', yf_to_ttherm_map_df=None, avg_dfs=None):
        if avg_dfs is None:
            avg_dfs = [None for _ in range(len(expression_dfs))]

        num_elements_list = [len(_) for _ in [expression_dfs, tab_labels, enrich_dfs, label_dfs, avg_dfs, expr_mins, expr_maxs]]

        if len(np.unique(num_elements_list)) != 1:
            raise ValueError('The following parameters must all have the same length: expression_dfs, tab_labels, enrich_dfs, label_dfs, and avg_dfs.')

        tabs = []

        for idx in range(num_elements_list[0]):
            expression_df = expression_dfs[idx]
            enrich_df = enrich_dfs[idx]
            label_df = label_dfs[idx]
            avg_df = avg_dfs[idx]
            expr_min = expr_mins[idx]
            expr_max = expr_maxs[idx]
            palette = palettes[idx]

            tab_label = tab_labels[idx]

            p = generate_umap(expression_df, enrich_df, annotation_df, label_df, phase, palette, title, n_neighbors=n_neighbors, n_components=n_components, random_state=random_state, expr_min=expr_min, expr_max=expr_max, embedding_metric=embedding_metric, yf_to_ttherm_map_df=yf_to_ttherm_map_df, avg_df=avg_df)

            tabs.append(TabPanel(child=p, title=tab_label))

        tabbed_plot = Tabs(tabs=tabs, sizing_mode='stretch_both')

        bokeh.plotting.output_file(filename=outfile_name, title=title, mode='inline')

        bokeh.plotting.save(tabbed_plot)

        print(outfile_name)

        return tabbed_plot

def generate_umap_tabbed(expression_dfs: list, tab_labels: list, enrich_dfs: list, annotation_df: pd.DataFrame, label_dfs: list, phase, palettes, title, n_neighbors=5, n_components=2, random_state=42, expr_mins=[], expr_maxs=[], embedding_metric='euclidean', yf_to_ttherm_map_df=None, avg_dfs=None):
        if avg_dfs is None:
            avg_dfs = [None for _ in range(len(expression_dfs))]

        num_elements_list = [len(_) for _ in [expression_dfs, tab_labels, enrich_dfs, label_dfs, avg_dfs, expr_mins, expr_maxs]]

        if len(np.unique(num_elements_list)) != 1:
            raise ValueError('The following parameters must all have the same length: expression_dfs, tab_labels, enrich_dfs, label_dfs, and avg_dfs.')

        tabs = []

        for idx in range(num_elements_list[0]):
            expression_df = expression_dfs[idx]
            enrich_df = enrich_dfs[idx]
            label_df = label_dfs[idx]
            avg_df = avg_dfs[idx]
            expr_min = expr_mins[idx]
            expr_max = expr_maxs[idx]
            palette = palettes[idx]

            tab_label = tab_labels[idx]

            p = generate_umap(expression_df, enrich_df, annotation_df, label_df, phase, palette, title, n_neighbors=n_neighbors, n_components=n_components, random_state=random_state, expr_min=expr_min, expr_max=expr_max, embedding_metric=embedding_metric, yf_to_ttherm_map_df=yf_to_ttherm_map_df, avg_df=avg_df)

            tabs.append(TabPanel(child=p, title=tab_label))

        tabbed_plot = Tabs(tabs=tabs, sizing_mode='stretch_both')

        return tabbed_plot

def generate_and_save_mds(outfile_name, expression_df, annotation_df, label_df, phase, palette, title, n_neighbors=5, n_components=2, random_state=42, expr_min=0, expr_max=1, embedding_metric='euclidean', yf_to_ttherm_map_df=None):
    
    data = expression_df[list(expression_df.columns)[1:]].values

    data_dists = compute_pairwise_distance_matrix(data, embedding_metric)

    mds_mapper = MDS(n_components=2, normalized_stress='auto', dissimilarity='precomputed', n_jobs=-1)
    embedding = mds_mapper.fit_transform(data_dists)
    
    umap_df = pd.DataFrame(np.array(embedding), columns=('x', 'y'))

    radius = compute_2d_embedding_point_radius(umap_df)
    
    bokeh.plotting.output_file(filename=outfile_name, title=title, mode='inline')
    p = plot_embedding(expression_df, umap_df, annotation_df, label_df, phase, palette, title=title, n_neighbors=n_neighbors, radius=radius, expr_min=expr_min, expr_max=expr_max, yf_to_ttherm_map_df=yf_to_ttherm_map_df)
    bokeh.plotting.save(p)
    print(outfile_name)
    return p    

def generate_and_save_pca(outfile_name, expression_df, annotation_df, label_df, phase, palette, title, n_neighbors=5, n_components=2, random_state=42, expr_min=0, expr_max=1):
    
    data = expression_df[list(expression_df.columns)[1:]].values

    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(data)    

    pca_df = pd.DataFrame({
        'x': principal_components[:, 0],
        'y': principal_components[:, 1]
    })

    radius = compute_2d_embedding_point_radius(pca_df)
    
    bokeh.plotting.output_file(filename=outfile_name, title=title, mode='inline')
    p = plot_embedding(expression_df, pca_df, annotation_df, label_df, phase, palette, title=title, n_neighbors=n_neighbors, radius=radius, expr_min=expr_min, expr_max=expr_max)
    bokeh.plotting.save(p)
    print(outfile_name)
    return p

def generate_and_save_tsne(outfile_name, expression_df, annotation_df, label_df, phase, palette, title, n_neighbors=5, n_components=2, random_state=42, expr_min=0, expr_max=1):
    
    data = expression_df[list(expression_df.columns)[1:]].values

    # tsne = TSNE(n_components=n_components, perplexity=30.0, n_iter=1000)
    tsne = TSNE(n_components=n_components, perplexity=50.0, n_iter=5000, n_jobs=-1)
    tsne_components = tsne.fit_transform(data)

    tsne_df = pd.DataFrame({
        'x': tsne_components[:, 0],
        'y': tsne_components[:, 1]
    })

    radius = compute_2d_embedding_point_radius(tsne_df)
    
    bokeh.plotting.output_file(filename=outfile_name, title=title, mode='inline')
    p = plot_embedding(expression_df, tsne_df, annotation_df, label_df, phase, palette, title=title, n_neighbors=n_neighbors, radius=radius, expr_min=expr_min, expr_max=expr_max)
    bokeh.plotting.save(p)
    print(outfile_name)
    return p


# def generate_server_data(expression_df, annotation_df, label_df, phases, palette, n_neighbors=3, random_state=42, embedding_metric='euclidean'):
    
#     data = expression_df[list(expression_df.columns)[1:]].values
    
#     umap_mapper = umap.UMAP(random_state=random_state, n_components=2, n_neighbors=n_neighbors, metric=embedding_metric).fit(data)
#     embedding = _get_umap_embedding(umap_mapper)
    
#     embedding_df = pd.DataFrame(np.array(embedding), columns=('x', 'y'))
    
#     embedding_df['TTHERM_ID'] = expression_df['TTHERM_ID'].values
    
#     merge_unsorted = expression_df.merge(embedding_df, on='TTHERM_ID')

#     merge_all = label_df.merge(merge_unsorted, on='TTHERM_ID')

#     merge_all_sorted = merge_all

#     labels = merge_all_sorted['label'].values

#     merge = merge_all_sorted.loc[: , merge_all_sorted.columns]

#     relevant_annot = annotation_df.iloc[np.in1d(annotation_df['TTHERM_ID'].values, merge['TTHERM_ID'].values)]
#     merge = merge.merge(relevant_annot, on='TTHERM_ID')
    
#     # if phases in ['full', 'veg', 'sex']:
#     #     mean_expression_df = get_arith_mean_expression(expression_df)

#     # elif phases == 'rna_seq':
#     #     mean_expression_df = ari_mean_df_of_duplicates(expression_df)

#     mean_expression_df = expression_df
    
#     embedding_df['TTHERM_ID'] = mean_expression_df['TTHERM_ID'].values
    
#     merge_unsorted = mean_expression_df.merge(embedding_df, on='TTHERM_ID')

#     merge_all = label_df.merge(merge_unsorted, on='TTHERM_ID')

#     merge_all_sorted = merge_all

#     labels = merge_all_sorted['label'].values

#     merge = merge_all_sorted.loc[: , merge_all_sorted.columns]
    
#     # take part of annotation df that shared TTHERM_IDs with expression df
#     relevant_annot = annotation_df.iloc[np.in1d(annotation_df['TTHERM_ID'].values, merge['TTHERM_ID'].values)]
#     merge = merge.merge(relevant_annot, on='TTHERM_ID')

#     ttherm_ids = merge['TTHERM_ID'].values
    
#     if phases == 'full':
        
#         x = ['Ll', 
#              'Lm', 
#              'Lh', 
#              'S0', 
#              'S3', 
#              'S6', 
#              'S9', 
#              # 'S12',
#              'S15', 
#              'S24', 
#              'C0', 
#              # 'C2', 
#              'C4', 
#              'C6', 
#              'C8', 
#              'C10', 
#              'C12', 
#              'C14', 
#              'C16', 
#              'C18']
        
        
#     elif phases == 'veg':
        
#         x = ['Ll', 
#              'Lm', 
#              'Lh', 
#              'S0', 
#              'S3', 
#              'S6', 
#              'S9', 
#              # 'S12', 
#              'S15', 
#              'S24']
        
#     elif phases == 'sex':
        
#         x = ['C0', 
#              # 'C2', 
#              'C4', 
#              'C6', 
#              'C8', 
#              'C10',
#              'C12',
#              'C14', 
#              'C16', 
#              'C18']
        
#     elif phases == 'rna_seq':
#         x = ['000min', '030min', '060min', '090min', '120min', '150min',
#        '180min', '210min', '240min']
        
#     elif phases == 'rna_seq_single_cycle':
#         x = ['000min', '030min', '060min', '090min', '120min', '150min', '180min']
        
#     else:
#         raise(ValueError(f'\"{phases}\" is an invalid choice for parameter \"phases.\"'))

#     xs = [x for ttid in ttherm_ids]

#     ys = [merge.loc[merge['TTHERM_ID'] == ttid, x].values[0] for ttid in ttherm_ids]

#     merge['expr_xs'] = xs
#     merge['expr_ys'] = ys
    
#     merge['module'] = [f'm{int(l):04d}' for l in labels]

#     pr = compute_2d_embedding_point_radius(embedding_df)

#     merge['radius'] = [pr for _ in range(merge.shape[0])]

#     return merge

def main():
    pass


if __name__ == '__main__':
    main()