import pandas as pd
import numpy as np
import bokeh
import umap
from bokeh.layouts import column
from bokeh.models import Button
from bokeh.palettes import RdYlBu3
from bokeh.plotting import figure, curdoc
from bokeh.plotting import show as show_interactive
from bokeh.plotting import output_file, output_notebook
from bokeh.layouts import column, row
from bokeh.models import Div, ColumnDataSource, CustomJS, TextInput, LassoSelectTool, Select, MultiSelect, ColorBar, Legend, LegendItem, Spinner
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn, Button, HTMLTemplateFormatter
from bokeh.events import SelectionGeometry
from bokeh.transform import linear_cmap, jitter
import pickle

with open('./server_data.pkl', 'rb') as f:
    merge = pickle.load(f)

num_genes = merge.shape[0]

num_clusters = max(merge['label'].unique()) + 1

labels = merge['labels'].values

with open(('colors_2000_1'), 'rb') as file:
    color_palette_raw = pickle.load(file)

palette65 = """
white\ngainsboro\n#FA002E\n#22FC22\n#221CFA\n#FF3DD6\n#FFDA00\n#00FEFB\n#F48684\n#CEB4FE\n#FFFFE5\n#0D933D\n#CC00F8\n#800D5D\n#F10084\n#22267A\n#0DADFF\n#CBFD71\n#9A761C\n#F96C00\n#6399A6\n#FFBCDA\n#8D0DA3\n#F79F26\n#00FFBF\n#A37CFB\n#F68EEB\n#720D0D\n#F163AA\n#7E926A\n#826386\n#B41C32\n#9BEBCE\n#E2DB83\n#56D4FA\n#E6E2FB\n#925D58\n#F7C3A7\n#62E970\n#220DBD\n#5583BB\n#7EA01C\n#CDFDB6\n#FD00FB\n#B30D97\n#F5FF00\n#DD77FD\n#4282FC\n#BBA6A4\n#0D8068\n#AB5F26\n#F7C26E\n#9EFE00\n#9B2EFD\n#C56887\n#FD3D68\n#ABF2FD\n#835FAC\n#FF16B1\n#325371\n#CA16CA\n#D26322\n#AFCFFE\n#91A1FA\nfloralwhite
""".split()

color_palette = palette65

if len(color_palette_raw) >= num_clusters:
    color_palette = color_palette_raw[: num_clusters]

palette = color_palette

x = ['000min', '030min', '060min', '090min', '120min', '150min', '180min', '210min', '240min']

title = 'RNA-SEQ SERVER'

expr_min = -3

expr_max = 3

hover_data = None

radius = None

embedding_df,
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
radius=None, 
#     subset_points=None,
interactive_text_search=False,
interactive_text_search_columns=['TTHERM_ID', 'PFAMs', 'Description', 'TGD2021_description', 'module'],
interactive_text_search_alpha_contrast=0.999,
alpha=None,
expr_min = 0,
expr_max = 1,
plot_sizing_mode='stretch_both',
table_sizing_mode = 'stretch_both',
search_sizing_mode = 'stretch_both'

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
            data['color'] = ['green']*len(labels)
        else:

            new_color_key = {k: color_key[i] for i, k in enumerate(unique_labels)}
            data["color"] = pd.Series(labels).map(new_color_key)

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

plot = bokeh.plotting.figure(
    tooltips=tooltips,
    tools="tap,box_select,pan,wheel_zoom,box_zoom,reset,save",
    background_fill_color=background,
    title=title,
    sizing_mode=plot_sizing_mode,
    output_backend="webgl"
#             x_range=(np.floor(min(points[:,0])), np.ceil(max(points[:,0]))), # Get axes
#             y_range=(np.floor(min(points[:,1])), np.ceil(max(points[:,1])))
)

if point_size is not None:

    plot.circle(
        x="x",
        y="y",
        source=data_source,
        color=colors,
        size=point_size,
        alpha="alpha",
        line_color='black'
    )

elif radius is not None:
    plot.circle(
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


x_heatmap_profile = x

# ['Ll', 
#      'Lm', 
#      'Lh', 
#      'S0', 
#      'S3', 
#      'S6', 
#      'S9', 
#      # 'S12', 
#      'S15', 
#      'S24', 
#      'C0', 
#      # 'C2', 
#      'C4', 
#      'C6', 
#      'C8', 
#      'C10', 
#      'C12', 
#      'C14', 
#      'C16', 
#      'C18']

# if normalized:
hm_min = expr_min
hm_max = expr_max
    
# else:
#     hm_min = 2
#     hm_max = 16

# For companion heatmap plot
ttherm_ids = embedding_df['TTHERM_ID'].values
hm_df = embedding_df[['TTHERM_ID'] + x_heatmap_profile]
hm_df['module'] = hover_data['module'].values
hm_df_tidy = hm_df.melt(id_vars=['TTHERM_ID', 'module'], var_name='phase', value_name='normalized_expression')
hm_cds = bokeh.plotting.ColumnDataSource(hm_df_tidy)
hm_cds.data['fill_alpha'] = [0.7]*len(hm_df_tidy)
hm_cds.data['line_alpha'] = [0.7]*len(hm_df_tidy)
# hm_cds.data['y_axis'] = ttherm_ids


hm = heatmap(hm_cds, bokeh.palettes.Inferno256, hm_min, hm_max, x_heatmap_profile, ttherm_ids, plot_sizing_mode=plot_sizing_mode)



# For companion expression plot

expr_source = bokeh.plotting.ColumnDataSource(dict(
    ID=['blah'], 
    expr_xs=[['Ll']], 
    expr_ys=[[0]],
    alpha=[0],
    color=['black']))

# if normalized:
y_axis_label = 'Mean expression z-score'
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
                    line_color="color"
                    )

expr_fig.xaxis.major_label_orientation = np.pi/4
expr_fig.xaxis.major_label_text_font_size = '12pt'
expr_fig.yaxis.major_label_text_font_size = '12pt'
expr_fig.yaxis.axis_label_text_font_size = '12pt'
expr_fig.xgrid.grid_line_color='whitesmoke'
expr_fig.xgrid.grid_line_alpha=0.2
expr_fig.ygrid.grid_line_color='whitesmoke'
expr_fig.ygrid.grid_line_alpha=0.2

# For data table
s2 = bokeh.plotting.ColumnDataSource(data=dict(ID=[]))

columns = [TableColumn(field="ID",  title="TTHERM_ID", formatter=HTMLTemplateFormatter(template='<a href="http://tet.ciliate.org/index.php/feature/details/feature_details.php?feature_name=<%= ID %>"target="_blank"><%= ID %></a>'), width=150),
        #    TableColumn(field="YF_ID", title="YF_ID", width=300),
            TableColumn(field="module",  title="Module", width=100),
            TableColumn(field='TGD2021_description', title='TGD2021_description', width=300),
            TableColumn(field="Description", title="eggNOG_description", width=300),
            TableColumn(field="Preferred_name", title="eggNOG_preferred_name", width=150),
            TableColumn(field="max_annot_lvl", title="max_annot_lvl", width=100),
            TableColumn(field="COG_category", title="COG_category", width=100),
            TableColumn(field='EC', title='EC', width=100),
            TableColumn(field='GOs', title='GOs', width=300),
            TableColumn(field='PFAMs', title='PFAMs', width=300),
            TableColumn(field='KEGG_ko', title='KEGG_ko', width=300),
            TableColumn(field='KEGG_Pathway', title='KEGG_Pathway', width=300),
            TableColumn(field='KEGG_Module', title='KEGG_Module', width=300),
            TableColumn(field='KEGG_Reaction', title='KEGG_Reaction', width=300),
            TableColumn(field='KEGG_rclass', title='KEGG_rclass', width=300),
            TableColumn(field='BRITE', title='BRITE', width=300),
        #    TableColumn(field='KEGG_TC', title='KEGG_TC', width=300),
        #    TableColumn(field='CAZy', title='CAZy', width=300),
        #    TableColumn(field='BiGG_Reaction', title='BiGG_Reaction', width=300),
#                    TableColumn(field="x",  title="x", width=300),
#                    TableColumn(field="y",  title="y")
            ]
table = DataTable(source=s2, 
                    columns=columns, 
                    editable=True,
                    selectable=True,
                    sortable=True,
                #   index_width=10,
                    fit_columns=False,
                    sizing_mode=table_sizing_mode
                    )

heatmap_callback = CustomJS(
    args=dict(
        s1=data_source,
        s_hm=hm_cds,
        cols=x_heatmap_profile
    ),
    code="""
    var d1 = s1.data;
    var d_hm = s_hm.data;
    
    var inds = s1.selected.indices;
    const num_cols = cols.length;
    
    //d_hm['TTHERM_ID'] = []
    //d_hm['normalized_expression'] = []
    d_hm['fill_alpha'] = []
    d_hm['line_alpha'] = []
    
    var selected_ttherm_ids = [];
    
    var ttids = d_hm['TTHERM_ID'].slice(0, """+str(num_genes)+""");
    
    if (inds.length == 0) {
        d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)
        d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.7)
    }else{
    
        // Start with everything deselected
        d_hm['fill_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)
        d_hm['line_alpha'] = Array(d_hm['TTHERM_ID'].length).fill(0.01)
    
        // Get the selected indices
        for (var i = 0; i < inds.length; i++) {
            selected_ttherm_ids.push(d1['ID'][inds[i]])
        }
        console.log(selected_ttherm_ids);
        
        // iterate over the selected ttherm ids
        for (var j = 0; j < selected_ttherm_ids.length; j++) {
        
            // var selected_gene = selected_ttherm_ids[j];
            // console.log(selected_gene);

            // ad hoc function to find if ttherm ids match
            var match = (element) => element == selected_ttherm_ids[j];
        
            // get index of matching ttherm id in heatmap
            var gene_index = ttids.findIndex(match);
            console.log(gene_index);
            
            // loop over the columns and highlight the selected genes
            for (var k = 0; k < num_cols; k++) {
            
                d_hm['fill_alpha'][gene_index] = 0.7
                d_hm['line_alpha'][gene_index] = 0.7

                gene_index = gene_index + ttids.length
            
            }
        
        }
        
    }
    
    console.log(d_hm);
    
    s_hm.change.emit();
    
    """
)

expression_callback = CustomJS(
    args=dict(
        s1=data_source,
        s_expr=expr_source,
        alpha=alpha,
    ),
    code="""
    var d1 = s1.data;
    var d_expr = s_expr.data;

    var inds = s1.selected.indices;
    // console.log(inds)

    // console.log(d1['ID'].length)

    // d1['alpha'] = Array(d1['ID'].length).fill(0.2)

    // console.log(d_expr['ID'].length, d_expr['expr_xs'].length, d_expr['expr_ys'].length)

    d_expr['ID'] = [['blah']]
    d_expr['expr_xs'] = [['Ll']]
    d_expr['expr_ys'] = [[0]]
    d_expr['alpha'] = [0]
    d_expr['color'] = ['black']
    // s_expr.change.emit();

    // debugger;

    for (var i = 0; i < inds.length; i++) {
        // d_expr['alpha'][inds[i]] = 1/(inds.length * 2)
        // console.log(inds[i], i)
        d_expr['ID'].push(Array(18).fill(d1['ID'][inds[i]]))
        d_expr['expr_xs'].push(d1['expr_xs'][inds[i]])
        d_expr['expr_ys'].push(d1['expr_ys'][inds[i]])
        d_expr['alpha'].push(Math.min(1, Math.max(7/(inds.length), 0.05)))
        d_expr['color'].push(d1['color'][inds[i]])
        // console.log(d_expr)
        // console.log(i)
        // console.log(
        //     d_expr['ID'].length, 
        //     d_expr['expr_xs'].length, 
        //     d_expr['expr_ys'].length
        // )
    }

    // s1.change.emit();
    s_expr.change.emit();
    // console.log(s_expr.data)

    """

)

selection_callback =  CustomJS(args=dict(
                                        s1=data_source, 
                                        s2=s2,
                                        table=table), 
                                            code="""

    var d1 = s1.data;
    var d2 = s2.data;


    var inds = s1.selected.indices;

    // Start by making everything tiny and pale
    d1['alpha'] = Array(d1['ID'].length).fill(0.01)
    d1['line_alpha'] = Array(d1['ID'].length).fill(0.01)
    d1['radius'] = Array(d1['ID'].length).fill(0.0001)

    d2['module'] = []
    d2['ID'] = []
    // d2['YF_ID'] = []
    d2['TGD2021_description'] = []
    d2['Description'] = []
    d2['Preferred_name'] = []
    d2['max_annot_lvl'] = []
    d2['COG_category'] = []
    d2['EC'] = []
    d2['GOs'] = []
    d2['PFAMs'] = []
    d2['KEGG_ko'] = []
    d2['KEGG_Pathway'] = []
    d2['KEGG_Module'] = []
    d2['KEGG_Reaction'] = []
    d2['KEGG_rclass'] = []
    d2['BRITE'] = []
    // d2['KEGG_TC'] = []
    // d2['CAZy'] = []
    // d2['BiGG_Reaction'] = []

    for (var i = 0; i < inds.length; i++) {
        d2['module'].push(d1['module'][inds[i]])
        d2['ID'].push(d1['ID'][inds[i]])
        // d2['YF_ID'].push(d1['YF_ID'][inds[i]])
        d2['TGD2021_description'].push(d1['TGD2021_description'][inds[i]])
        d2['Description'].push(d1['Description'][inds[i]])
        d2['Preferred_name'].push(d1['Preferred_name'][inds[i]])
        d2['max_annot_lvl'].push(d1['max_annot_lvl'][inds[i]])
        d2['COG_category'].push(d1['COG_category'][inds[i]])
        d2['EC'].push(d1['EC'][inds[i]])
        d2['GOs'].push(d1['GOs'][inds[i]])
        d2['PFAMs'].push(d1['PFAMs'][inds[i]])
        d2['KEGG_ko'].push(d1['KEGG_ko'][inds[i]])
        d2['KEGG_Pathway'].push(d1['KEGG_Pathway'][inds[i]])
        d2['KEGG_Module'].push(d1['KEGG_Module'][inds[i]])
        d2['KEGG_Reaction'].push(d1['KEGG_Reaction'][inds[i]])
        d2['KEGG_rclass'].push(d1['KEGG_rclass'][inds[i]])
        d2['BRITE'].push(d1['BRITE'][inds[i]])
        // d2['KEGG_TC'].push(d1['KEGG_TC'][inds[i]])
        // d2['CAZy'].push(d1['CAZy'][inds[i]])
        // d2['BiGG_Reaction'].push(d1['BiGG_Reaction'][inds[i]])

        d1['alpha'][inds[i]] = 1
        d1['line_alpha'][inds[i]] = 1
        d1['radius'][inds[i]] = 20
    }
    s1.change.emit();
    s2.change.emit();
    table.change.emit();
""")

data_source.selected.js_on_change('indices', selection_callback, expression_callback, heatmap_callback)

if interactive_text_search:
    text_input = TextInput(value="Search comma-separated module(s), TTHERM_ID(s), or functional term(s), for example: m0003, TTHERM_01207610, Histone", sizing_mode=search_sizing_mode)

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
        callback = CustomJS(
            args=dict(
                source=data_source,
                s2=s2,
                table=table,
                matching_alpha=interactive_text_search_alpha_contrast,
                non_matching_alpha=1 - interactive_text_search_alpha_contrast,
                search_columns=interactive_text_search_columns,
                default_radius=radius,
                default_alpha=alpha
            ),
            code="""
            var data = source.data;
            var text_search = cb_obj.value;
            var d2 = s2.data;

            // var ref_expr = ref_e_s.data;
            // var d3 = sel_e_s.data;

            var search_terms = text_search.split(',');

            d2['module'] = []
            d2['ID'] = []
            // d2['YF_ID'] = []

            // d3['xs'] = []
            // d3['ys'] = []

            var search_columns_dict = {}
            for (var col in search_columns){
                search_columns_dict[col] = search_columns[col]
            }

            // First, clear the data table and selection
            // data['alpha'] = []
            // data['radius'] = []
            source.selected.indices = []

            // source.change.emit();
            s2.change.emit();
            // sel_e_s.change.emit();
            table.change.emit();

            // Run search
            if (text_search.length > 0){
                // Loop over columns and values
                // If there is no match for any column for a given row, change the alpha value
                var string_match = false;
                for (var i = 0; i < data.x.length; i++) {
                    string_match = false
                    for (var j in search_columns_dict) {
                        if (search_terms.some(t => String(data[search_columns_dict[j]][i]).includes(t.trim()))) {
                            string_match = true
                        }
                    }
                    if (string_match){
                        // data['alpha'][i] = matching_alpha
                        // data['radius'][i] = 1
                        d2['module'].push(data['module'][i])
                        d2['ID'].push(data['ID'][i])
                        // d2['YF_ID'].push(data['YF_ID'][i])

                        // d3['xs'].push(ref_expr['xs'][i])
                        // d3['ys'].push(ref_expr['ys'][i])

                        // So that these points are actually considered selected
                        source.selected.indices.push(i)

                    }else{
                        // data['alpha'][i] = non_matching_alpha
                        // data['radius'][i] = 0.01
                    }
                }
                source.change.emit();
                s2.change.emit();
                // sel_e_s.change.emit();
                table.change.emit();

            } else {

                // Loop over columns and values
                // If there is no match for any column for a given row, change the alpha value
                var string_match = false;
                for (var i = 0; i < data.x.length; i++) {
                    string_match = false
                    for (var j in search_columns_dict) {
                        if (search_terms.some(t => String(data[search_columns_dict[j]][i]).includes(t.trim()))) {
                            string_match = true
                        }
                    }
                    if (string_match){
                        // data['alpha'][i] = default_alpha
                        // data['radius'][i] = default_radius
                        d2['module'].push()
                        d2['ID'].push()
                        // d2['YF_ID'].push()

                        // d3['xs'].push()
                        // d3['ys'].push()

                    }else{
                        // data['alpha'][i] = non_matching_alpha
                        // data['radius'][i] = 0.01
                    }
                }
                source.change.emit();
                s2.change.emit();
                // sel_e_s.change.emit();
                table.change.emit();

            }




        """,
        )
        
        text_input.js_on_change("value", callback, selection_callback, expression_callback, heatmap_callback)

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
    spinner_callback = CustomJS(
        args=dict(
            source=data_source,
            s2=s2,
            table=table,
            matching_alpha=interactive_text_search_alpha_contrast,
            non_matching_alpha=1 - interactive_text_search_alpha_contrast,
            search_columns=interactive_text_search_columns_spinner,
            max_label_num_len=max_label_num_len,
            default_radius=radius,
            default_alpha=alpha
        ),
        code="""
function zfill(str, width) {
const diff = width - str.length;
if (diff > 0) {
    return '0'.repeat(diff) + str;
}
return str;
}

var data = source.data;
var spinner_num = cb_obj.value;
var spinner_str = "m" + zfill(spinner_num.toString(), max_label_num_len);
var d2 = s2.data;

console.log(spinner_str);

d2['module'] = [];
d2['ID'] = [];

var search_columns_dict = {}
for (var col in search_columns){
search_columns_dict[col] = search_columns[col]
}

// First, clear the data table and selection
// data['alpha'] = [];
// data['radius'] = [];
source.selected.indices = [];

s2.change.emit();
table.change.emit();

// Run search
if (spinner_str.length > 0) {
// Loop over rows and check if search term is present in any column
for (var i = 0; i < data.x.length; i++) {
    var string_match = false;
    for (var col in search_columns_dict) {
        if (String(data[search_columns_dict[col]][i]) == (spinner_str)) {
            string_match = true;
            break; // Exit the loop if match is found in any column
        }
    }
    if (string_match) {
        // data['alpha'][i] = matching_alpha;
        // data['radius'][i] = 1;
        d2['module'].push(data['module'][i]);
        d2['ID'].push(data['ID'][i]);

        // So that these points are actually considered selected
        source.selected.indices.push(i);
    } else {
        // data['alpha'][i] = non_matching_alpha;
        // data['radius'][i] = 0.01;
    }
}
} else {
// Loop over rows and set default alpha and radius values
for (var i = 0; i < data.x.length; i++) {
    // data['alpha'][i] = default_alpha;
    // data['radius'][i] = default_radius;
}
}

source.change.emit();
s2.change.emit();
table.change.emit();
    """,
    )
    
    spinner.js_on_change("value", spinner_callback, selection_callback, expression_callback, heatmap_callback)


# Lifted from https://stackoverflow.com/questions/31824124/is-there-a-way-to-save-bokeh-data-table-content
download_button1 = Button(label='💾 Annotations', button_type="success")
download_button1.js_on_click(
    CustomJS(
        args=dict(source_data=data_source),
        code="""
        var inds = source_data.selected.indices;
        var data = source_data.data;
        var out = "TTHERM_ID\tmodule\tTGD2021_description\teggNOG_description\teggNOG_preferred_name\tmax_annot_lvl\tCOG_category\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\\n";
        for (var i = 0; i < inds.length; i++) {
            out += data['ID'][inds[i]] + "\t" + data['module'][inds[i]] + "\t" + data['TGD2021_description'][inds[i]] + "\t" + data['Description'][inds[i]] + "\t" + data['Preferred_name'][inds[i]] + "\t" + data['max_annot_lvl'][inds[i]] + "\t" + data['COG_category'][inds[i]] + "\t" + data['GOs'][inds[i]] + "\t" + data['EC'][inds[i]] + "\t" + data['KEGG_ko'][inds[i]] + "\t" + data['KEGG_Pathway'][inds[i]] + "\t" + data['KEGG_Module'][inds[i]] + "\t" + data['KEGG_Reaction'][inds[i]] + "\t" + data['KEGG_rclass'][inds[i]] + "\t" + data['BRITE'][inds[i]] + "\\n";
        }
        var file = new Blob([out], {type: 'text/plain'});
        var elem = window.document.createElement('a');
        elem.href = window.URL.createObjectURL(file);
        elem.download = 'selected-annotation-data.tsv';
        document.body.appendChild(elem);
        elem.click();
        document.body.removeChild(elem);
        """))        

# NEED TO STOP HARDCODING THIS FILE
enrich_df = pd.read_csv(os.path.join(file_dir, main_dir, 'enrichment/test_nn3_full_enrichment.csv'))
colors = [color_key[int(m) % len(color_key)] for m in enrich_df['module'].values]
enrich_df['color'] = colors

enrich_cds = bokeh.models.ColumnDataSource(enrich_df)
enrich_p = plot_enrichment(enrich_cds, plot_sizing_mode=plot_sizing_mode)

download_button2 = Button(label='💾 Enrichment', button_type="success")
download_button2.js_on_click(
    CustomJS(
        args=dict(source_data=enrich_cds),
        code="""
        // var inds = source_data.selected.indices;
        var data = source_data.data;
        var out = "module\tterm\tinfo\tfold_change\tbonferroni\\n";
        for (var i = 0; i < data['module'].length; i++) {
            out += data['module'][i] + "\t" + data['term'][i] + "\t" + data['info'][i] + "\t" + data['fold_change'][i] + "\t" + data['bonferroni'][i] + "\\n";
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
    text_search2 = TextInput(value="ENRICHMENT TERM SEARCH (UNDER CONSTRUCTION)", sizing_mode='stretch_width')

    # FIXME: IMPLEMENT
    data_module_stats = {'Module': [1, 2, 3, 4],
                        'Stats': ['A', 'B', 'C', 'D']}
    source_module_stats = ColumnDataSource(data=data_module_stats)
    columns_module_stats = [TableColumn(field="Module", title="Module"),
                            TableColumn(field="Stats", title="Stats")]
    module_stats_table = DataTable(source=source_module_stats, columns=columns_module_stats)


    spinner_text_html = """
    <div style="display: flex; justify-content: center; align-items: center; height: 100%;">
        <p style="text-align: center;">MODULE #:</p>
    </div>
    """
    spinner_custom_text = Div(text=spinner_text_html)

    cola_sizing_mode = 'stretch_width'
    colb_sizing_mode = 'stretch_width'

    col1a = column(hm, height=800)
    col1a.sizing_mode = cola_sizing_mode
    col2a = column(enrich_p, height=800)
    col2a.sizing_mode = cola_sizing_mode
    col3a = column(plot, expr_fig, height=800)
    col3a.sizing_mode = cola_sizing_mode

    col1b = column(module_stats_table, height=800, max_width=450)
    col1b.sizing_mode = colb_sizing_mode
    col2b = column(table, height=800, max_width=100000)
    col2b.sizing_mode = colb_sizing_mode

    rows_sizing_mode = 'stretch_width'

    row_search = row(text_input, download_button1, download_button2, sizing_mode='stretch_width')
    rowa = row(row(col1a, col2a, sizing_mode=rows_sizing_mode), col3a)
    rowa.sizing_mode = rows_sizing_mode
    rowb = row(col2b)
    rowb.sizing_mode = rows_sizing_mode

    plot = column(row_search, rowa, rowb)
    plot.sizing_mode = 'stretch_width'

else:
    plot = column(row(column(plot, expr_fig), hm, enrich_p), row(download_button1, download_button2), table)