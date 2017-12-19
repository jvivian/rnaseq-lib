import holoviews as hv

gene_curves_opts = {
    'Curve': {'plot': dict(height=120, width=600, tools=['hover'], invert_xaxis=True, yrotation=45, yaxis='left'),
              'style': dict(line_width=1.5)},
    'Curve.Percentage_of_Normal_Samples': {'plot': dict(xaxis=None, invert_yaxis=True),
                                           'style': dict(color='Blue')},
    'Curve.Gene_Expression': {'plot': dict(xaxis=None),
                              'style': dict(color='Green')},
    'Curve.Log2_Fold_Change': {'plot': dict(height=150),
                               'style': dict(color='Purple')},
    'Scatter': {'style': dict(color='red', size=3)}}

gene_kde_opts = {
    'Distribution': {'plot': dict(filled=False),
                    'style': dict(line_color=hv.Cycle())}
}

gene_distribution_opts = {
    'BoxWhisker': {'plot': dict(width=900, xrotation=45)}
}

gene_de_opts = {
    'Scatter': {'plot': dict(color_index='Tissue', legend_position='left', width=700, height=500, tools=['hover']),
                'style': dict(cmap='tab20', size=10, alpha=0.5)}
}
