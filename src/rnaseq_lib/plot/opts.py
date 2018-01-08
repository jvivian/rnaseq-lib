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

gene_kde_opts = {'Overlay': {'plot': dict(width=500, legend_position='left')}}

gene_distribution_opts = {'BoxWhisker': {'plot': dict(width=875, xrotation=70)}}

gene_de_opts = {
    'Scatter': {'plot': dict(color_index='Tissue', legend_position='left', width=700, height=500, tools=['hover']),
                'style': dict(cmap='tab20', size=10, alpha=0.5)}}

sample_count_opts = {
    'Bars': {'plot': dict(width=875, xrotation=70, tools=['hover'], show_legend=False)}
}

l2fc_by_perc_samples_opts = {
    'Curve': {'plot': dict(tools=['hover'])},
    'Overlay': {'plot': dict(legend_position='left', width=500)},
    'Spikes': {'plot': dict(spike_length=100),
               'style': dict(line_alpha=0.4, line_width=5)}
}
