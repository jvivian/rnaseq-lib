"""
This code is copied from Philippjfr's notebook:
https://anaconda.org/philippjfr/sankey/notebook
"""

from functools import cmp_to_key

import holoviews as hv
import numpy as np
import param
from holoviews import Operation

hv.extension('bokeh')

from collections import defaultdict

from holoviews.core.util import basestring, max_range
from holoviews.element.graphs import Graph, Nodes, EdgePaths, Dataset, redim_graph


class Sankey(Graph):
    group = param.String(default='Sankey', constant=True)

    def __init__(self, data, kdims=None, vdims=None, compute=True, **params):
        if isinstance(data, tuple):
            data = data + (None,) * (3 - len(data))
            edges, nodes, edgepaths = data
        else:
            edges, nodes, edgepaths = data, None, None
        if nodes is not None:
            if not isinstance(nodes, Dataset):
                if nodes.ndims == 3:
                    nodes = Nodes(nodes)
                else:
                    nodes = Dataset(nodes)
                    nodes = nodes.clone(kdims=nodes.kdims[0],
                                        vdims=nodes.kdims[1:])
        node_info = nodes
        super(Graph, self).__init__(edges, kdims=kdims, vdims=vdims, **params)
        if compute:
            self._nodes = nodes
            chord = layout_sankey(self)
            self._nodes = chord.nodes
            self._edgepaths = chord.edgepaths
            self._sankey = chord._sankey
        else:
            if not isinstance(nodes, Nodes):
                raise TypeError("Expected Nodes object in data, found %s."
                                % type(nodes))
            self._nodes = nodes
            if not isinstance(edgepaths, EdgePaths):
                raise TypeError("Expected EdgePaths object in data, found %s."
                                % type(edgepaths))
            self._edgepaths = edgepaths
            self._sankey = None
        self._validate()
        self.redim = redim_graph(self, mode='dataset')


from holoviews.plotting.bokeh import GraphPlot
from bokeh.models import Patches


class SankeyPlot(GraphPlot):
    label_index = param.ClassSelector(default=None, class_=(basestring, int),
                                      allow_None=True, doc="""
        Index of the dimension from which the node labels will be drawn""")

    filled = True

    _style_groups = dict(GraphPlot._style_groups, quad='nodes', text='label')

    _draw_order = ['patches', 'multi_line', 'text', 'quad']

    style_opts = GraphPlot.style_opts + ['edge_fill_alpha', 'nodes_line_color', 'label_text_font_size']

    def _init_glyphs(self, plot, element, ranges, source):
        ret = super(SankeyPlot, self)._init_glyphs(plot, element, ranges, source)
        renderer = plot.renderers.pop(plot.renderers.index(self.handles['glyph_renderer']))
        plot.renderers = [renderer] + plot.renderers
        return ret

    def get_data(self, element, ranges, style):
        data, mapping, style = super(SankeyPlot, self).get_data(element, ranges, style)
        quad_mapping = {'left': 'x0', 'right': 'x1', 'bottom': 'y0', 'top': 'y1'}
        quad_data = data['scatter_1']
        quad_data.update({'x0': [], 'x1': [], 'y0': [], 'y1': []})
        for node in element._sankey['nodes']:
            quad_data['x0'].append(node['x0'])
            quad_data['y0'].append(node['y0'])
            quad_data['x1'].append(node['x1'])
            quad_data['y1'].append(node['y1'])
        data['quad_1'] = quad_data
        quad_mapping['fill_color'] = mapping['scatter_1']['node_fill_color']
        mapping['quad_1'] = quad_mapping
        style['nodes_line_color'] = 'black'

        lidx = element.nodes.get_dimension(self.label_index)
        if lidx is None:
            if self.label_index is not None:
                dims = element.nodes.dimensions()[2:]
                self.warning("label_index supplied to Chord not found, "
                             "expected one of %s, got %s." %
                             (dims, self.label_index))
            return data, mapping, style
        if element.vdims:
            edges = Dataset(element)[element[element.vdims[0].name] > 0]
            nodes = list(np.unique([edges.dimension_values(i) for i in range(2)]))
            nodes = element.nodes.select(**{element.nodes.kdims[2].name: nodes})
        else:
            nodes = element
        labels = [lidx.pprint_value(v) for v in nodes.dimension_values(lidx)]
        ys = nodes.dimension_values(1)
        nodes = element._sankey['nodes']
        offset = (nodes[0]['x1'] - nodes[0]['x0']) / 4.
        xs = np.array([node['x1'] for node in nodes])
        data['text_1'] = dict(x=xs + offset, y=ys, text=[str(l) for l in labels])
        mapping['text_1'] = dict(text='text', x='x', y='y', text_baseline='middle', text_align='left')
        return data, mapping, style

    def get_extents(self, element, ranges):
        """
        A Chord plot is always drawn on a unit circle.
        """
        xdim, ydim = element.nodes.kdims[:2]
        xpad = .05 if self.label_index is None else 0.25
        x0, x1 = ranges[xdim.name]
        y0, y1 = ranges[ydim.name]
        xdiff = (x1 - x0)
        ydiff = (y1 - y0)
        x0, x1 = max_range([xdim.range, (x0 - (0.05 * xdiff), x1 + xpad * xdiff)])
        y0, y1 = max_range([ydim.range, (y0 - (0.05 * ydiff), y1 + (0.05 * ydiff))])
        return (x0, y0, x1, y1)

    def _postprocess_hover(self, renderer, source):
        if self.inspection_policy == 'edges':
            if not isinstance(renderer.glyph, Patches):
                return
        else:
            if isinstance(renderer.glyph, Patches):
                return
        super(SankeyPlot, self)._postprocess_hover(renderer, source)


hv.Store.register({Sankey: SankeyPlot}, 'bokeh')

options = hv.Store.options('bokeh')

options.Sankey = hv.Options('plot', xaxis=None, yaxis=None, inspection_policy='edges',
                            selection_policy='nodes', width=1000, height=600, show_frame=False)
options.Sankey = hv.Options('style', node_line_alpha=0, node_nonselection_alpha=0.2, node_size=10,
                            edge_nonselection_alpha=0.2, edge_line_alpha=0, edge_fill_alpha=0.8,
                            label_text_font_size='8pt')


def weightedSource(link):
    return nodeCenter(link['source']) * link['value']


def weightedTarget(link):
    return nodeCenter(link['target']) * link['value']


def nodeCenter(node):
    return (node['y0'] + node['y1']) / 2


def ascendingBreadth(a, b):
    return int(a['y0'] - b['y0'])


def ascendingSourceBreadth(a, b):
    return ascendingBreadth(a['source'], b['source']) | a['index'] - b['index']


def ascendingTargetBreadth(a, b):
    return ascendingBreadth(a['target'], b['target']) | a['index'] - b['index']


def quadratic_bezier(start, end, c0=(0, 0), c1=(0, 0), steps=25):
    """
    Compute quadratic bezier spline given start and end coordinate and
    two control points.
    """
    steps = np.linspace(0, 1, steps)
    sx, sy = start
    ex, ey = end
    cx0, cy0 = c0
    cx1, cy1 = c1
    xs = ((1 - steps) ** 3 * sx + 3 * ((1 - steps) ** 2) * steps * cx0 +
          3 * (1 - steps) * steps ** 2 * cx1 + steps ** 3 * ex)
    ys = ((1 - steps) ** 3 * sy + 3 * ((1 - steps) ** 2) * steps * cy0 +
          3 * (1 - steps) * steps ** 2 * cy1 + steps ** 3 * ey)
    return np.column_stack([xs, ys])


class layout_sankey(Operation):
    """
    Computes a Sankey diagram from a Graph element.

    Adapted from d3-sankey under BSD-3 license.
    """

    bounds = param.NumericTuple(default=(0, 0, 1000, 500))

    node_width = param.Number(default=15)

    node_padding = param.Integer(default=10)

    iterations = param.Integer(32)

    def _process(self, element, key=None):
        graph = {'nodes': [], 'links': []}
        self.computeNodeLinks(element, graph)
        self.computeNodeValues(graph)
        self.computeNodeDepths(graph)
        self.computeNodeBreadths(graph)
        self.computeLinkBreadths(graph)

        paths = []
        for link in graph['links']:
            source, target = link['source'], link['target']
            x0, y0 = source['x1'], link['y0']
            x1, y1 = target['x0'], link['y1']
            start = np.array([(x0, link['width'] + y0),
                              (x0, y0)])
            src = (x0, y0)
            ctr1 = ((x0 + x1) / 2., y0)
            ctr2 = ((x0 + x1) / 2., y1)
            tgt = (x1, y1)
            bottom = quadratic_bezier(src, tgt, ctr1, ctr2)
            mid = np.array([(x1, y1),
                            (x1, y1 + link['width'])])

            xmid = (x0 + x1) / 2.
            y0 = y0 + link['width']
            y1 = y1 + link['width']
            src = (x1, y1)
            ctr1 = (xmid, y1)
            ctr2 = (xmid, y0)
            tgt = (x0, y0)
            top = quadratic_bezier(src, tgt, ctr1, ctr2)
            spline = np.concatenate([start, bottom, mid, top])
            paths.append(spline)

        node_data = []
        for node in graph['nodes']:
            node_data.append((np.mean([node['x0'], node['x1']]),
                              np.mean([node['y0'], node['y1']]),
                              node['index']) + tuple(node['values']))
        nodes = Nodes(node_data, vdims=element.nodes.vdims)
        edges = EdgePaths(paths)
        sankey = Sankey((element.data, nodes, edges), compute=False)
        sankey._sankey = graph
        return sankey

    def computeNodeLinks(self, element, graph):
        """
        Populate the sourceLinks and targetLinks for each node.
        Also, if the source and target are not objects, assume they are indices.
        """
        index = element.nodes.kdims[-1]
        node_map = {}
        values = element.nodes.array(element.nodes.vdims)
        for node, vals in zip(element.nodes.dimension_values(index), values):
            node = {'index': node, 'sourceLinks': [], 'targetLinks': [], 'values': vals}
            graph['nodes'].append(node)
            node_map[node['index']] = node

        links = [element.dimension_values(d) for d in element.dimensions()[:3]]
        for i, (src, tgt, value) in enumerate(zip(*links)):
            source, target = node_map[src], node_map[tgt]
            link = dict(index=i, source=source, target=target, value=value)
            graph['links'].append(link)
            source['sourceLinks'].append(link)
            target['targetLinks'].append(link)

    def computeNodeValues(self, graph):
        """
        Compute the value (size) of each node by summing the associated links.
        """
        for node in graph['nodes']:
            source_val = np.sum([l['value'] for l in node['sourceLinks']])
            target_val = np.sum([l['value'] for l in node['targetLinks']])
            node['value'] = max([source_val, target_val])

    def computeNodeDepths(self, graph):
        """
        Iteratively assign the depth (x-position) for each node.
        Nodes are assigned the maximum depth of incoming neighbors plus one;
        nodes with no incoming links are assigned depth zero, while
        nodes with no outgoing links are assigned the maximum depth.
        """
        nodes = graph['nodes']
        depth = 0
        while nodes:
            next_nodes = []
            for node in nodes:
                node['depth'] = depth
                for link in node['sourceLinks']:
                    if link['target'] not in next_nodes:
                        next_nodes.append(link['target'])
            nodes = next_nodes
            depth += 1

        nodes = graph['nodes']
        depth = 0
        while nodes:
            next_nodes = []
            for node in nodes:
                node['height'] = depth
                for link in node['targetLinks']:
                    if link['source'] not in next_nodes:
                        next_nodes.append(link['source'])
            nodes = next_nodes
            depth += 1

        x0, _, x1, _ = self.p.bounds
        dx = self.p.node_width
        kx = (x1 - x0 - dx) / (depth - 1)
        for node in graph['nodes']:
            d = node['depth'] if node['sourceLinks'] else depth - 1
            node['x0'] = x0 + max([0, min([depth - 1, np.floor(d)]) * kx])
            node['x1'] = node['x0'] + dx

    def computeNodeBreadths(self, graph):
        node_map = hv.OrderedDict()
        for n in graph['nodes']:
            if n['x0'] not in node_map:
                node_map[n['x0']] = []
            node_map[n['x0']].append(n)

        _, y0, _, y1 = self.p.bounds
        py = self.p.node_padding

        def initializeNodeBreadth():
            kys = []
            for nodes in node_map.values():
                nsum = np.sum([node['value'] for node in nodes])
                ky = (y1 - y0 - (len(nodes) - 1) * py) / nsum
                kys.append(ky)
            ky = np.min(kys)

            for nodes in node_map.values():
                for i, node in enumerate(nodes):
                    node['y0'] = i
                    node['y1'] = i + node['value'] * ky

            for link in graph['links']:
                link['width'] = link['value'] * ky

        def relaxLeftToRight(alpha):
            for nodes in node_map.values():
                for node in nodes:
                    if not node['targetLinks']:
                        continue
                    weighted = sum([weightedSource(l) for l in node['targetLinks']])
                    tsum = sum([l['value'] for l in node['targetLinks']])
                    center = nodeCenter(node)
                    dy = (weighted / tsum - center) * alpha
                    node['y0'] += dy
                    node['y1'] += dy

        def relaxRightToLeft(alpha):
            for nodes in list(node_map.values())[::-1]:
                for node in nodes:
                    if not node['sourceLinks']:
                        continue
                    weighted = sum([weightedTarget(l) for l in node['sourceLinks']])
                    tsum = sum([l['value'] for l in node['sourceLinks']])
                    center = nodeCenter(node)
                    dy = (weighted / tsum - center) * alpha
                    node['y0'] += dy
                    node['y1'] += dy

        def resolveCollisions():
            for nodes in node_map.values():
                y = y0
                n = len(nodes)
                nodes.sort(key=cmp_to_key(ascendingBreadth))
                for node in nodes:
                    dy = y - node['y0']
                    if dy > 0:
                        node['y0'] += dy
                        node['y1'] += dy
                    y = node['y1'] + py

                dy = y - py - y1
                if dy > 0:
                    node['y0'] -= dy
                    node['y1'] -= dy
                    y = node['y0']
                    for node in nodes[:-1][::-1]:
                        dy = node['y1'] + py - y;
                        if dy > 0:
                            node['y0'] -= dy
                            node['y1'] -= dy
                        y = node['y0']

        initializeNodeBreadth()
        resolveCollisions()
        alpha = 1
        for _ in range(self.p.iterations):
            alpha = alpha * 0.99
            relaxRightToLeft(alpha)
            resolveCollisions()
            relaxLeftToRight(alpha)
            resolveCollisions()

    def computeLinkBreadths(self, graph):
        for node in graph['nodes']:
            node['sourceLinks'].sort(key=cmp_to_key(ascendingTargetBreadth))
            node['targetLinks'].sort(key=cmp_to_key(ascendingSourceBreadth))

        for node in graph['nodes']:
            y0 = y1 = node['y0']
            for link in node['sourceLinks']:
                link['y0'] = y0
                y0 += link['width']
            for link in node['targetLinks']:
                link['y1'] = y1
                y1 += link['width']