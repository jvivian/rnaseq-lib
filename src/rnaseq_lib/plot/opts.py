gene_curves = {
    'Curve': {'plot': dict(height=120, width=600, tools=['hover'], invert_xaxis=True, yrotation=45, yaxis='left'),
              'style': dict(line_width=1.5)},
    'Curve.Percentage_of_Normal_Samples': {'plot': dict(xaxis=None, invert_yaxis=True),
                                           'style': dict(color='Blue')},
    'Curve.Gene_Expression': {'plot': dict(xaxis=None),
                              'style': dict(color='Green')},
    'Curve.Log2_Fold_Change': {'plot': dict(height=150),
                               'style': dict(color='Purple')},
    'Scatter': {'style': dict(color='red', size=3)}}

gene_kde = {}
