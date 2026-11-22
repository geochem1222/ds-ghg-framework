"""Review candidates: publication-scale styling, with scientific content protected."""
from __future__ import annotations

import ast
import csv
import hashlib
import json
from pathlib import Path
import shutil

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, PathPatch, Polygon, Rectangle, RegularPolygon
from matplotlib.path import Path as ShapePath
import numpy as np

HERE = Path(__file__).resolve().parent
CURRENT_SOURCE = HERE.parent / 'editable_figure_sources'
SOURCE = HERE / 'source_baseline' if (HERE / 'source_baseline').exists() else CURRENT_SOURCE
PACKAGE = HERE.parent.parent
OUT = HERE / 'candidates'
INK = '#272727'
MUTED = '#666B70'
RULE = '#CCD1D5'
BLUE = '#0F4D92'
WATER_EDGE = '#3775BA'
GREEN = '#327B69'
PURPLE = '#79627F'
GOLD = '#B28116'
GAS = '#B36A92'
WATER = '#EFF5FA'
SEDIMENT = '#E8E7E2'
OXIC = '#D9EBE1'
REDUCED = '#E8E1E9'
DEPTH_COLORS = ['#EDF5FA', '#D2E3EF', '#A3C8DF', '#649BBF', '#286083']

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7.2, 'text.color': INK, 'axes.labelcolor': INK,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': .65, 'legend.frameon': False,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'savefig.facecolor': 'white', 'savefig.transparent': False,
})


class Sheet:
    def __init__(self, number, height):
        self.number, self.width, self.height = number, 183, height
        self.fig = plt.figure(figsize=(183 / 25.4, height / 25.4))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set(xlim=(0, 183), ylim=(0, height), aspect='equal')
        self.ax.axis('off')
        self.labels = []

    def text(self, x, y, s, size=7.2, color=INK, ha='left', bold=False, **kwargs):
        obj = self.ax.text(x, y, s, fontsize=size, color=color, ha=ha, va='center',
                           weight='bold' if bold else 'normal', linespacing=1.2, **kwargs)
        self.labels.append(obj)
        return obj

    def arrow(self, a, b, color=BLUE, lw=1., style='-|>', curve=0, dashed=False, scale=7):
        patch = FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=scale,
                                connectionstyle=f'arc3,rad={curve}', linewidth=lw,
                                linestyle=(0, (3, 2)) if dashed else '-', color=color,
                                shrinkA=1, shrinkB=1)
        self.ax.add_patch(patch)
        return patch

    def line(self, xs, ys, color=RULE, lw=.65, **kwargs):
        return self.ax.plot(xs, ys, color=color, lw=lw, **kwargs)

    def rect(self, x, y, w, h, color, edge='none', lw=.65, **kwargs):
        patch = Rectangle((x, y), w, h, facecolor=color, edgecolor=edge, lw=lw, **kwargs)
        self.ax.add_patch(patch)
        return patch

    def panel(self, letter, x, y, title):
        self.text(x, y, letter, 9, bold=True)
        self.text(x + 6, y, title, 8.2, bold=True)

    def finish(self, stem, semantic_checks):
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        outside, collisions = [], []
        for i, label in enumerate(self.labels):
            box = label.get_window_extent(renderer)
            if (box.x0 < 0 or box.y0 < 0 or box.x1 > self.fig.bbox.width
                    or box.y1 > self.fig.bbox.height):
                outside.append(label.get_text())
            for other in self.labels[i+1:]:
                if box.overlaps(other.get_window_extent(renderer)):
                    collisions.append([label.get_text(), other.get_text()])
        stats = {'width_mm': self.width, 'height_mm': self.height,
                 'minimum_label_pt': min(t.get_fontsize() for t in self.labels),
                 'labels': len(self.labels), 'off_canvas_labels': outside,
                 'text_collisions': collisions, 'empirical_data': False,
                 'model_predictions': False, **semantic_checks}
        (OUT / f'{stem}_qa.json').write_text(json.dumps(stats, indent=2), encoding='utf-8')
        assert not outside, outside
        assert not collisions, collisions
        for ext in ('svg', 'pdf'):
            self.fig.savefig(OUT / f'{stem}.{ext}')
        for ext in ('png', 'tiff'):
            kwargs = {'pil_kwargs': {'compression': 'tiff_lzw'}} if ext == 'tiff' else {}
            self.fig.savefig(OUT / f'{stem}.{ext}', dpi=600, **kwargs)
        self.fig.savefig(OUT / f'{stem}_preview.png', dpi=220)
        plt.close(self.fig)
        assert '<text' in (OUT / f'{stem}.svg').read_text(encoding='utf-8')
        return stats


def depth_value(x, y, mode):
    if mode == 'shallow':
        return 0.16 + 0.20 * (0.55 + 0.25 * np.sin(2.8 * x) + 0.20 * np.cos(3.4 * y))
    radius = np.sqrt(((x - 0.52) / 0.72) ** 2 + ((y - 0.48) / 0.82) ** 2)
    return 0.18 + 0.78 * np.clip(1 - radius, 0, 1)


def source_depth_function():
    tree = ast.parse((SOURCE / 'figure1_depth_state_ghg.py').read_text(encoding='utf-8'))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'depth_value')
    namespace = {'np': np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), '<protected depth function>', 'exec'), namespace)
    return namespace['depth_value']


def gas_marker(sheet, x, y, color, shape):
    if shape == 'circle':
        patch = Circle((x, y), .85, fc=color, ec='white', lw=.4)
    elif shape == 'square':
        patch = Rectangle((x-.8, y-.8), 1.6, 1.6, fc=color, ec='white', lw=.4)
    else:
        patch = RegularPolygon((x, y), numVertices=3, radius=1.1, fc=color, ec='white', lw=.4)
    sheet.ax.add_patch(patch)


def figure1():
    s = Sheet(1, 112)
    s.panel('a', 3, 107, 'Depth distribution')
    s.text(6, 100, 'Same area, different bathymetry', 7, MUTED)
    old_depth = source_depth_function()
    maps = []
    for x, name, mode in [(6, 'Shallow-dominated', 'shallow'), (40, 'Deeper basin', 'deep')]:
        values = []
        for i in range(8):
            for j in range(6):
                v = float(np.clip(depth_value((i+.5)/8, (j+.5)/6, mode), 0, .999))
                assert v == float(np.clip(old_depth((i+.5)/8, (j+.5)/6, mode), 0, .999))
                values.append(v)
                s.rect(x+i*29/8, 77+j*16/6, 29/8, 16/6,
                       DEPTH_COLORS[min(int(v*5), 4)], edge='white', lw=.12)
        s.rect(x, 77, 29, 16, 'none', WATER_EDGE, .7)
        s.line([x+1, x+28], [85, 85], color=MUTED, lw=.6, ls=(0, (3, 2)))
        s.text(x+14.5, 73, name, 7, ha='center', bold=True)
        counts = np.histogram(values, bins=np.linspace(0, 1, 6))[0] / 48
        assert np.isclose(counts.sum(), 1)
        s.line([x+2, x+27], [54, 54], color=MUTED, lw=.55)
        for i, count in enumerate(counts):
            s.rect(x+2+i*5, 54, 4.7, 12*count, DEPTH_COLORS[i])
        maps.append({'mode': mode, 'values': values, 'area_fraction': counts.tolist()})
        xx = np.linspace(0, 1, 120)
        bed = 37 - 12*depth_value(xx, .5, mode)
        assert np.array_equal(bed, 37 - 12*old_depth(xx, .5, mode))
        s.ax.fill_between(x+29*xx, 22, bed, color='#D5D0C7', lw=0)
        s.ax.fill_between(x+29*xx, bed, 37, color=WATER, lw=0)
        s.line(x+29*xx, bed, color='#8F8B83', lw=.75)
        s.line([x, x+29], [37, 37], color=WATER_EDGE, lw=.8)
        xb = x+29*xx[len(xx)//2]
        yb = bed[len(xx)//2]
        s.line([xb, xb], [yb, 37], color=WATER_EDGE, lw=.7)
        for y in (yb, 37):
            s.line([xb-.7, xb+.7], [y, y], color=WATER_EDGE, lw=.7)
    s.text(37.5, 48.5, r'Area fraction $p(d)$', 7, MUTED, ha='center')
    s.text(37.5, 42, 'Matched bed transects', 7, MUTED, ha='center')
    s.text(20.5, 17.5, 'Shorter path', 7, MUTED, ha='center')
    s.text(54.5, 17.5, 'Longer path', 7, MUTED, ha='center')
    s.text(37.5, 8, 'Equal area, different storage capacity\nand transfer distance', 7.2, ha='center', bold=True)
    s.line([74, 74], [6, 108], lw=.5)

    s.panel('b', 79, 107, 'Vertical coupling')
    s.text(81, 100, 'Same local depth, different exchange', 7, MUTED)
    for x, color, word in [(82, GOLD, 'Energy'), (111, BLUE, 'Oxidants'), (147, GAS, 'Dissolved gas')]:
        s.line([x, x+4], [93.5, 93.5], color=color, lw=1.45)
        s.text(x+5.3, 93.5, word, 7, MUTED)
    for x, state, name in [(81, 'mixed', 'Mixed'), (116, 'stratified', 'Stratified'),
                            (151, 'ventilating', 'Ventilating')]:
        w, bed, surface, interface = 25, 31, 78, 54
        s.rect(x, bed, w, surface-bed, WATER)
        if state != 'mixed':
            s.rect(x, bed, w, interface-bed, '#AAC6DB' if state == 'stratified' else '#CDDEEB')
            s.line([x, x+w], [interface, interface], color=WATER_EDGE, lw=.7, ls=(0, (3, 2)))
        s.rect(x, bed, w, surface-bed, 'none', WATER_EDGE, .7)
        s.rect(x, 26, w, bed-26, '#D5D0C7')
        s.line([x, x+w], [bed, bed], color='#8F8B83', lw=.75)
        xx = np.linspace(x, x+w, 80)
        s.line(xx, surface+.25*np.sin(np.linspace(0, 4*np.pi, 80)), color=WATER_EDGE, lw=.8)
        s.text(x+w/2, 85, name, 7.4, ha='center', bold=True)
        for dx, color, shape in [(7, '#626B70', 'circle'), (12.5, GAS, 'square'), (18, GREEN, 'triangle')]:
            gas_marker(s, x+dx, bed+2, color, shape)
        if state == 'mixed':
            s.arrow((x+6, 75), (x+6, 36), GOLD, 1.15)
            s.arrow((x+12, 75), (x+12, 36), BLUE, 1.15)
            s.arrow((x+17, 36), (x+17, 81), GAS, 1.15)
            note = 'Continuous\ncoupling'
        elif state == 'stratified':
            s.arrow((x+6, 75), (x+6, 57), GOLD, 1.15)
            s.arrow((x+12, 75), (x+12, 57), BLUE, 1.15)
            s.line([x+17, x+17], [37, 51], color=GAS, lw=1.15)
            for y in (38, 44, 50):
                s.rect(x+16.3, y-.45, 1.4, .9, GAS)
            s.line([x+17, x+17], [52, 75], color=GAS, lw=.7, alpha=.65, ls=(0, (2, 2)))
            note = 'Dissolved-gas\nretention'
        else:
            s.ax.add_patch(Polygon([(x+3, 75), (x+9, 75), (x+14, 36), (x+9, 36)],
                                   fc='white', ec='none', alpha=.7))
            s.arrow((x+8, 75), (x+12, 36), BLUE, 1.35)
            s.arrow((x+17, 36), (x+17, 81), GAS, 1.35)
            s.arrow((x+4, 75), (x+6, 57), GOLD, .9)
            note = 'Reoxygenation\nand pulse release'
        s.arrow((x+22, 36), (x+22, 81), MUTED, .55, scale=5)
        for y in (39, 58, 71):
            s.ax.add_patch(Circle((x+22, y), .65, fc='white', ec=MUTED, lw=.65))
        s.text(x+w/2, 19.5, note, 7, ha='center')
    for x, color, shape, word in [(84, '#626B70', 'circle', 'CO$_2$'),
                                  (106, GAS, 'square', 'CH$_4$'),
                                  (128, GREEN, 'triangle', 'N$_2$O')]:
        gas_marker(s, x, 8, color, shape)
        s.text(x+2, 8, word, 7, MUTED)
    s.ax.add_patch(Circle((151, 8), .65, fc='white', ec=MUTED, lw=.65))
    s.text(154, 8, 'CH$_4$ bubbles', 7, MUTED)
    (OUT/'figure1_schematic_depths.json').write_text(json.dumps(maps, indent=2), encoding='utf-8')
    return s.finish('figure1_refined', {
        'bathymetric_values_identical': True, 'histograms_normalized': True,
        'transects_identical': True, 'bubble_path_in_all_states': True,
        'bubble_radius_mm': .65, 'gas_identity_encoded_by_color_and_shape': True,
        'depth_thresholds_assigned': False})


def figure2():
    s = Sheet(2, 174)
    s.panel('a', 3, 168, 'Light and oxidant access')
    for x, surface, oxic, title, subtitle, illuminated in [
        (10, 141, 7, 'Illuminated, renewed oxygen', 'Shallower-water example', True),
        (109, 149, 2, 'Shaded, restricted renewal', 'Deeper-water example', False),
    ]:
        bed, bottom, w = 122, 100, 64
        s.text(x+w/2, 158, title, 7.6, ha='center', bold=True)
        s.text(x+w/2, 153, subtitle, 7, MUTED, ha='center')
        s.rect(x, bed, w, surface-bed, WATER)
        s.rect(x, bottom, w, bed-bottom, SEDIMENT)
        s.rect(x, bed-oxic, w, oxic, OXIC)
        s.rect(x, bed-oxic-6, w, 6, REDUCED)
        s.line([x, x+w], [surface, surface], color=WATER_EDGE, lw=.9)
        s.line([x, x+w], [bed, bed], color=MUTED, lw=.65)
        s.line([x, x+w], [bed-oxic, bed-oxic], color=GREEN, lw=.6, ls=(0, (2, 2)))
        s.line([x, x+w], [bed-oxic-6, bed-oxic-6], color=PURPLE, lw=.6, ls=(0, (2, 2)))
        beam_end = max(bed, surface-22)
        s.ax.add_patch(Polygon([(x+7, surface), (x+17, surface),
                               (x+21, beam_end), (x+3, beam_end)], fc='#E5CB79', ec='none', alpha=.35))
        s.arrow((x+12, surface-1), (x+12, beam_end+.8), GOLD, 1.05)
        s.text(x+12, surface+2.5, 'Light', 7, GOLD, ha='center')
        if illuminated:
            for dx in (7, 10, 13, 16, 19):
                s.ax.add_patch(Ellipse((x+dx, bed+.5), 1.6, .55, fc=GREEN, ec='white', lw=.3))
            s.arrow((x+12, bed+.8), (x+4, bed+5.5), GREEN, .85)
            s.text(x+3, bed+8.5, 'O$_2$', 7.2, GREEN, ha='center')
            s.arrow((x+20, bed+.8), (x+33, bed+6.5), GOLD, .85)
            s.text(x+36, bed+12, 'Photochemical\nCH$_4$', 7, GOLD, ha='center')
            s.arrow((x+53, surface-2), (x+53, bed+2), GREEN, 1.05)
            s.text(x+57, bed+10, 'O$_2$', 7.2, GREEN)
            s.text(x+32, bed-oxic/2, 'Aerobic oxidation', 7.1, GREEN, ha='center')
        else:
            s.line([x+2, x+w-2], [surface-14, surface-14], color=MUTED, lw=.7, ls=(0, (4, 2)))
            s.arrow((x+53, surface-2), (x+53, surface-11.5), GREEN, .85)
            s.text(x+57, surface-7, 'O$_2$', 7.2, GREEN)
            s.text(x+32, surface-19, 'Restricted renewal', 7.1, MUTED, ha='center')
            s.text(x+31, bed+5, 'Oxic interface', 7, GREEN, ha='center')
            s.arrow((x+48, bed+4), (x+55, bed-oxic/2), GREEN, .65, scale=5)
        s.text(x+32, bed-oxic-3, 'AOM where acceptors persist', 7, PURPLE, ha='center')
        s.text(x+32, bottom+(bed-oxic-6-bottom)/2, 'CH$_4$ production', 7.2, ha='center')
    s.arrow((78, 141), (105, 141), MUTED, .8, style='<->')
    s.text(91.5, 146.5, 'Water-level\nchange', 7.1, MUTED, ha='center')
    s.text(91.5, 123, 'Same\nbed', 7.1, MUTED, ha='center')
    s.text(91.5, 95.5, 'Conditional states; sediment profiles not to scale', 7, MUTED, ha='center')
    s.line([4, 179], [91, 91], lw=.55)

    s.panel('b', 3, 86, 'Methane production, consumption and transport')
    s.rect(17, 31, 150, 41, WATER)
    s.rect(17, 5, 150, 26, SEDIMENT)
    s.rect(17, 28, 150, 3, OXIC)
    s.line([17, 167], [72, 72], color=WATER_EDGE, lw=.8)
    s.line([17, 167], [31, 31], color=MUTED, lw=.65)
    for x in (17, 167):
        s.line([x, x], [5, 72], color=MUTED, lw=.55, ls=(0, (2, 3)))
    s.text(18, 76, 'Atmosphere', 7, MUTED)
    s.text(18, 68, 'Water', 7, MUTED)
    s.text(18, 8, 'Sediment', 7, MUTED)
    for x, y, w, h, word in [(91, 50, 37, 11, 'Dissolved CH$_4$\ninventory'),
                              (91, 15, 34, 10, 'Pore-water CH$_4$')]:
        s.ax.add_patch(Ellipse((x, y), w, h, fc='#DCE8F3', ec=BLUE, lw=.8))
        s.text(x, y, word, 7.4, BLUE, ha='center')
    s.arrow((91, 56), (91, 81), BLUE, 1., style='<->')
    s.text(96, 77, 'Air-water\nexchange', 7.1, BLUE)
    s.text(46, 64, 'Water-column\nproduction', 7.1, ha='center')
    s.text(46, 56.8, 'Where supported', 7, MUTED, ha='center')
    s.arrow((57, 62), (77, 55), BLUE, .9)
    s.arrow((107, 54), (123, 62), GREEN, 1.)
    s.text(137, 63, 'Oxidation\nto DIC', 7.2, GREEN, ha='center')
    s.arrow((2, 49), (70, 49), BLUE, .9)
    s.text(21, 53.5, 'Imported CH$_4$', 7, BLUE)
    s.text(47, 42, 'Interfacial\nphotochemical CH$_4$', 7, GOLD, ha='center')
    s.arrow((62, 42), (75, 46), GOLD, .9)
    s.arrow((2, 36), (42, 36), PURPLE, .9)
    s.text(2, 39.7, 'Resources', 7, PURPLE)
    s.text(21, 33, 'Organic matter, nutrients, acceptors', 7, PURPLE)
    s.arrow((109, 47), (179, 47), MUTED, .85)
    s.text(137, 42, 'Lateral export', 7, MUTED, ha='center')
    s.arrow((91, 20.5), (91, 43.5), BLUE, 1., style='<->')
    s.text(96, 24.5, 'Sediment-water\nexchange', 7, BLUE)
    s.arrow((77, 19), (57, 30), GREEN, .95, curve=.25)
    s.text(48, 22, 'Aerobic oxidation\nto DIC', 7.1, GREEN, ha='center')
    s.arrow((73, 15), (62, 15), PURPLE, .85)
    s.text(49, 15, 'AOM to DIC', 7.1, PURPLE, ha='center')
    s.arrow((48, 11), (72, 14), BLUE, .9)
    s.text(47, 6.5, 'CH$_4$ production', 7.1, ha='center')
    s.arrow((109, 15), (144, 15), BLUE, .9, style='<->')
    s.text(125, 9.5, 'Gas-phase formation\nand dissolution', 7, BLUE, ha='center')
    # Sediment-constrained cavities and a release conduit differ from water-column bubbles.
    cavities = [
        [(150, 18), (150.6, 20.3), (152.7, 21), (155, 20.1),
         (156, 18.9), (155.1, 17.3), (153.3, 17.7), (152, 17)],
        [(146.8, 14.3), (147.6, 15.6), (149.4, 15.9), (150.2, 14.5),
         (149.5, 13.5), (147.6, 13.4)],
        [(154.2, 14.3), (155.1, 15.8), (157.1, 15.4), (158, 14.3),
         (156.6, 13.2), (155.1, 13.4)],
    ]
    for vertices in cavities:
        points = np.asarray(vertices)
        curve, codes = [points[0]], [ShapePath.MOVETO]
        for i in range(len(points)):
            p0, p1, p2, p3 = [points[j % len(points)] for j in (i-1, i, i+1, i+2)]
            curve.extend([p1+(p2-p0)/6, p2-(p3-p1)/6, p2])
            codes.extend([ShapePath.CURVE4] * 3)
        curve.append(points[0])
        codes.append(ShapePath.CLOSEPOLY)
        s.ax.add_patch(PathPatch(ShapePath(curve, codes), fc='white', ec=BLUE, lw=.65))
    s.text(155, 7.5, 'Sediment gas\nstorage', 7, BLUE, ha='center')
    s.line([153, 153.6, 152.5, 154, 154], [21, 24, 26, 28, 31], color=BLUE, lw=1.)
    s.arrow((154, 31), (159, 81), BLUE, 1.)
    for y in (34, 51, 65):
        s.ax.add_patch(Circle((153+(y-21)*6/60, y), .65, fc='white', ec=BLUE, lw=.65))
    s.text(162, 78, 'Ebullition', 7.1, BLUE)
    s.arrow((156.5, 56), (113, 54), BLUE, .85, curve=-.15, dashed=True)
    s.text(131, 54.5, 'Dissolution', 7, BLUE, ha='center')
    return s.finish('figure2_refined', {
        'conditional_light_oxygen_states': True, 'reaction_zones_not_to_scale': True,
        'bubble_radius_mm': .65, 'air_water_exchange_bidirectional': True,
        'bubble_radius_applies_to_water_column_only': True,
        'sediment_gas_cavities_irregular': True,
        'pore_water_gas_exchange_bidirectional': True,
        'sediment_cavities_and_release_conduit_schematic': True,
        'sediment_oxidation_and_AOM_retained': True,
        'water_column_sources_sinks_and_external_fluxes_retained': True,
        'interfacial_photochemical_source_retained': True,
        'resource_loading_distinct_from_imported_methane': True})


def conceptual_series():
    t = np.linspace(0, 1, 2001)
    dt = t[1]-t[0]

    def pulse(start, end, height=1):
        v = np.zeros_like(t)
        m = (t >= start) & (t <= end)
        v[m] = height*np.sin(np.pi*(t[m]-start)/(end-start))**2
        return v

    def integrate(v):
        return np.r_[0, np.cumsum((v[1:]+v[:-1])*.5*dt)]

    temp_shallow = pulse(.15, .43)
    temp_deep = pulse(.20, .62, .50)
    flux_ref = pulse(.17, .48)
    flux_thermal = pulse(.22, .64, .48)
    flux_redox = pulse(.22, .83, 1.22)
    flux_delayed = pulse(.50, .92)
    flux_delayed *= np.trapezoid(flux_ref, t)/np.trapezoid(flux_delayed, t)
    depth = np.interp(t, [0, .10, .28, .47, .67, 1], [0, 0, 1, 1, 0, 0])
    return locals()


def check_series(data):
    source = (SOURCE/'figure3_depth_state_framework.py').read_text(encoding='utf-8')
    start = source.index('t=np.linspace(0,1,2001)')
    stop = source.index('\ndef plot_row(', start)
    namespace = {'np': np}
    exec(compile(source[start:stop], '<protected conceptual arrays>', 'exec'), namespace)
    names = ['t', 'temp_shallow', 'temp_deep', 'flux_ref', 'flux_thermal',
             'flux_redox', 'flux_delayed', 'depth']
    for name in names:
        assert np.array_equal(data[name], namespace[name]), name
    return names


def figure3():
    s = Sheet(3, 154)
    s.panel('a', 3, 149, 'Depth history changes the response conditions')
    s.text(5, 139, 'External forcing $U(t)$', 7.5, BLUE)
    s.text(5, 133.5, 'Climate, light and resources', 7)
    s.text(5, 123, r'Depth history $d(t)$, $\dot d(t)$', 7.5, BLUE)
    s.text(5, 116, 'Atmospheric pressure', 7, MUTED)
    s.arrow((48, 135), (60, 135), BLUE)
    s.arrow((48, 124), (60, 124), BLUE)
    s.rect(61, 115, 66, 26, '#F1F5F8')
    s.text(94, 137, 'Temperature, light, oxidants', 7.1, ha='center')
    s.text(94, 131.5, 'Resources, gas stores, mixing, legacy', 7, ha='center')
    for x, word, term, color in [(72, 'Physics', r'$P_{\mathrm{phys}}$', BLUE),
                                  (94, 'Reactions', r'$R_{\mathrm{rxn}}$', GREEN),
                                  (116, 'Transport', r'$J_{\mathrm{trans}}$', MUTED)]:
        s.text(x, 121.5, word+'\n'+term, 7, color, ha='center')
    s.arrow((128, 133), (139, 133), BLUE)
    s.text(159, 136, 'Local flux $f_g$ (3)', 7.2, ha='center')
    s.arrow((158, 132), (158, 126), BLUE, .7)
    s.text(159, 122.5, 'Area-integrated $F_g$ (1)', 7.2, ha='center')
    s.text(159, 116.5, 'CO$_2$  CH$_4$  N$_2$O', 7.5, ha='center', bold=True)
    s.text(94, 111.7, 'State budget (2)', 7, MUTED, ha='center')
    s.line([48, 54, 54, 135, 135], [123, 123, 107, 107, 130], color=MUTED, lw=.65)
    s.arrow((48, 116), (54, 116), MUTED, .65)
    s.arrow((135, 130), (140, 130), MUTED, .65)
    s.text(91, 107, 'Pressure-sensitive release', 7, MUTED, ha='center',
           bbox={'fc': 'white', 'ec': 'none', 'pad': 1})
    s.line([4, 179], [102.5, 102.5], lw=.55)
    s.line([7, 15], [98, 98], color=MUTED, lw=1.2, ls=(0, (3, 2)))
    s.text(17, 98, 'Shallower / stable depth', 7, MUTED)
    s.line([104, 112], [98, 98], color=BLUE, lw=1.45)
    s.text(114, 98, 'Deeper / depth cycle', 7, BLUE)
    data = conceptual_series()
    identical = check_series(data)
    t = data['t']
    integrate = data['integrate']
    total_scale = integrate(data['flux_redox'])[-1]
    all_fluxes = [data['flux_thermal'], data['flux_redox'], data['flux_delayed']]
    xs = [8, 68, 128]
    titles = ['Thermal buffering', 'Oxidant limitation', 'Retention and release']
    subtitles = ['Oxidation fraction maintained', 'Oxidation loss offsets buffering', 'Excess net CH$_4$ supply matched']

    def row(x, y, h, v1, v2, scale, shade=False, fill=False):
        w = 48
        if shade:
            s.rect(x+w*.13, y, w*.19, h, '#EBEEF0')
        s.line([x, x+w], [y, y], color=RULE, lw=.55)
        if fill:
            s.ax.fill_between(x+w*t, y, y+h*v2/scale, color=BLUE, alpha=.065, lw=0)
        s.line(x+w*t, y+h*v1/scale, color=MUTED, lw=1.1, ls=(0, (3, 2)))
        s.line(x+w*t, y+h*v2/scale, color=BLUE, lw=1.45)

    for x in (62, 122):
        s.line([x, x], [8, 93], lw=.45)
    for j, (x, title, subtitle, flux) in enumerate(zip(xs, titles, subtitles, all_fluxes)):
        s.panel('bcd'[j], x-4, 91.5, title)
        s.text(x, 85.7, subtitle, 7, MUTED)
        s.text(x, 79.8, 'Bed-temperature anomaly' if j < 2 else 'Water-depth change', 7.2)
        v1, v2 = ((data['temp_shallow'], data['temp_deep']) if j < 2
                  else (np.zeros_like(t), data['depth']))
        row(x, 66, 10, v1, v2, 1, shade=True)
        s.text(x, 59.3, 'CH$_4$ emission anomaly', 7.2)
        row(x, 37, 17, data['flux_ref'], flux, 1.25, shade=True, fill=True)
        s.text(x, 30.5, 'Cumulative excess CH$_4$', 7.2)
        row(x, 10, 16, integrate(data['flux_ref']), integrate(flux), total_scale)
        s.arrow((x, 7), (x+48, 7), MUTED, .6)
        s.text(x+24, 3, 'Time', 7, MUTED, ha='center')
    inventory = integrate(data['flux_ref']-data['flux_delayed'])
    assert inventory.min() > -1e-10
    assert abs(inventory[-1]) < 1e-10
    columns = ['t', 'temp_shallow', 'temp_deep', 'depth', 'flux_ref', 'flux_thermal', 'flux_redox', 'flux_delayed']
    with (OUT/'figure3_conceptual_curves.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(zip(*(data[key] for key in columns)))
    return s.finish('figure3_refined', {
        'arrays_identical_to_current': identical, 'arbitrary_relative_scales': True,
        'shading': 'common warming interval; pale pulse fills are not uncertainty intervals',
        'cumulative_scale_shared': True, 'excess_inventory_min': float(inventory.min()),
        'final_excess_inventory': float(inventory[-1]),
        'equal_cumulative_release_in_d': bool(abs(inventory[-1]) < 1e-10),
        'oxidation_conditions_preserved': True})


def compare(number, old, new):
    fig, axs = plt.subplots(1, 2, figsize=(14.4, 14.4*max(old.shape[0]/old.shape[1], new.shape[0]/new.shape[1])/2 + .45))
    for ax, pixels, title in zip(axs, (old, new), ('Current', 'Refined candidate')):
        ax.imshow(pixels)
        ax.axis('off')
        ax.set_title(title, fontsize=11, loc='left', weight='bold')
    fig.subplots_adjust(left=.01, right=.99, bottom=.01, top=.94, wspace=.025)
    fig.savefig(HERE/f'comparison_figure{number}.png', dpi=170)
    plt.close(fig)
    pixels = new[:, :, :3]
    grey = np.dot(pixels, [.2126, .7152, .0722])
    h, w = grey.shape
    fig = plt.figure(figsize=(183/25.4, 183*h/w/25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(grey, cmap='gray', vmin=0, vmax=1)
    ax.axis('off')
    fig.savefig(OUT/f'figure{number}_grayscale.png', dpi=220)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    snapshot = HERE/'baseline_hashes.json'
    protected = list(CURRENT_SOURCE.glob('*')) + list((PACKAGE/'01_upload_now').glob('*'))
    originals_available = CURRENT_SOURCE.is_dir() and (PACKAGE/'01_upload_now').is_dir()
    frozen = HERE/'source_baseline'
    frozen.mkdir(exist_ok=True)
    for name in ('figure1_depth_state_ghg.py', 'figure2_reaction_environments.py',
                 'figure3_depth_state_framework.py', 'Figure_Captions_Current.md',
                 'figure1_depth_state_ghg_framework_preview.png',
                 'figure2_depth_aware_validation_preview.png',
                 'figure3_depth_state_framework_preview.png'):
        if not (frozen/name).exists():
            shutil.copy2(SOURCE/name, frozen/name)
    if not snapshot.exists():
        baseline = {str(p.relative_to(PACKAGE)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in protected if p.is_file()}
        snapshot.write_text(json.dumps(baseline, indent=2), encoding='utf-8')
    results = [figure1(), figure2(), figure3()]
    for i, stem in [(1, 'figure1_depth_state_ghg_framework'),
                    (2, 'figure2_depth_aware_validation'), (3, 'figure3_depth_state_framework')]:
        compare(i, plt.imread(SOURCE/f'{stem}_preview.png'), plt.imread(OUT/f'figure{i}_refined_preview.png'))
    baseline = json.loads(snapshot.read_text(encoding='utf-8'))
    if originals_available:
        assert all(hashlib.sha256((PACKAGE/p).read_bytes()).hexdigest() == digest for p, digest in baseline.items())
    (HERE/'qa_summary.json').write_text(json.dumps({
        'figures': results, 'current_submission_and_artwork_unchanged': True if originals_available else None,
        'matplotlib': mpl.__version__, 'numpy': np.__version__,
        'source': 'figures4papers general design principles; new manuscript-specific code',
        'manual_visual_review': 'pending'}, indent=2), encoding='utf-8')
    print(json.dumps({'status': 'Candidates generated', 'original_files_unchanged': True if originals_available else None,
                      'minimum_text_pt': 7, 'figures': 3}, indent=2))


if __name__ == '__main__':
    main()
