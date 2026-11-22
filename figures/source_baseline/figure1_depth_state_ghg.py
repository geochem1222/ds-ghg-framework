from __future__ import annotations

from pathlib import Path
import json

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figure1_depth_state_ghg_framework"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.transparent": False,
    }
)

INK = "#202428"
MUTED = "#626B70"
RULE = "#C8CDD0"
WATER = "#DDEFF3"
WATER_MID = "#B5D9E2"
WATER_DEEP = "#78B4C3"
WATER_EDGE = "#3E8FA3"
SEDIMENT = "#C7B8A6"
ORANGE = "#E69F00"
BLUE = "#0072B2"
GREEN = "#009E73"
PURPLE = "#CC79A7"
CO2 = "#59636A"
DEPTH_COLORS = ["#E8F4F7", "#C7E5EC", "#94CAD6", "#5CA9BC", "#287B91"]


def setup(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def panel_label(ax, label):
    ax.text(
        -0.035,
        1.01,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        ha="left",
        va="top",
        color=INK,
    )


def arrow(ax, start, end, color=INK, width=1.2, mutation=8, style="-|>", connection="arc3"):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=mutation,
        linewidth=width,
        color=color,
        shrinkA=0,
        shrinkB=0,
        connectionstyle=connection,
        zorder=6,
    )
    ax.add_patch(patch)
    return patch


def depth_value(x, y, mode):
    if mode == "shallow":
        return 0.16 + 0.20 * (0.55 + 0.25 * np.sin(2.8 * x) + 0.20 * np.cos(3.4 * y))
    radius = np.sqrt(((x - 0.52) / 0.72) ** 2 + ((y - 0.48) / 0.82) ** 2)
    return 0.18 + 0.78 * np.clip(1 - radius, 0, 1)


def draw_depth_map(ax, x0, y0, width, height, mode):
    nx, ny = 8, 6
    values = []
    for i in range(nx):
        for j in range(ny):
            xx = (i + 0.5) / nx
            yy = (j + 0.5) / ny
            value = float(np.clip(depth_value(xx, yy, mode), 0, 0.999))
            values.append(value)
            color = DEPTH_COLORS[min(int(value * len(DEPTH_COLORS)), len(DEPTH_COLORS) - 1)]
            ax.add_patch(
                Rectangle(
                    (x0 + i * width / nx, y0 + j * height / ny),
                    width / nx,
                    height / ny,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.18,
                )
            )
    ax.add_patch(Rectangle((x0, y0), width, height, facecolor="none", edgecolor=WATER_EDGE, lw=0.9))
    ax.plot(
        [x0 + 0.04 * width, x0 + 0.96 * width],
        [y0 + 0.50 * height, y0 + 0.50 * height],
        color=INK,
        lw=0.65,
        ls=(0, (3, 2)),
        alpha=0.75,
    )
    return np.asarray(values)


def draw_histogram(ax, x0, y0, width, height, values):
    counts, _ = np.histogram(values, bins=np.linspace(0, 1, 6))
    counts = counts / len(values)
    gap = 0.006
    bar_width = (width - gap * 4) / 5
    ax.plot([x0, x0 + width], [y0, y0], color=MUTED, lw=0.55)
    for i, count in enumerate(counts):
        ax.add_patch(
            Rectangle(
                (x0 + i * (bar_width + gap), y0),
                bar_width,
                height * count,
                facecolor=DEPTH_COLORS[i],
                edgecolor="white",
                linewidth=0.2,
            )
        )


def draw_cross_section(ax, x0, y0, width, height, mode):
    x = np.linspace(0, 1, 120)
    xs = x0 + width * x
    surface = y0 + height
    if mode == "shallow":
        note = "shorter transfer path"
    else:
        note = "longer path; greater capacity"
    bed = surface - height * depth_value(x, 0.5, mode)

    ax.fill_between(xs, y0 - 0.025, bed, color=SEDIMENT, linewidth=0)
    ax.fill_between(xs, bed, surface, color=WATER, linewidth=0)
    ax.plot(xs, bed, color="#7C6956", lw=0.75)
    ax.plot([x0, x0 + width], [surface, surface], color=WATER_EDGE, lw=0.9)
    middle = len(x) // 2
    dimension_x = xs[middle]
    dimension_bottom = bed[middle] + 0.003
    dimension_top = surface - 0.003
    cap = 0.018 * width
    ax.plot([dimension_x, dimension_x], [dimension_bottom, dimension_top], color=WATER_EDGE, lw=0.8)
    ax.plot([dimension_x - cap, dimension_x + cap], [dimension_bottom, dimension_bottom], color=WATER_EDGE, lw=0.8)
    ax.plot([dimension_x - cap, dimension_x + cap], [dimension_top, dimension_top], color=WATER_EDGE, lw=0.8)
    ax.text(x0 + width / 2, y0 - 0.045, note, ha="center", va="top", fontsize=5.7, color=MUTED)


def draw_geometry(ax):
    setup(ax)
    panel_label(ax, "a")
    ax.text(0.50, 0.975, "Geometry sets the depth distribution", ha="center", va="top", fontsize=8, fontweight="bold")
    ax.text(0.50, 0.90, "same inundated area; different bathymetry", ha="center", va="center", fontsize=6.2, color=MUTED)

    positions = [(0.06, "Shallow-dominated", "shallow"), (0.54, "Deeper basin", "deep")]
    for x0, title, mode in positions:
        values = draw_depth_map(ax, x0, 0.64, 0.39, 0.17, mode)
        ax.text(x0 + 0.195, 0.605, title, ha="center", va="center", fontsize=6.2, fontweight="bold")
        draw_histogram(ax, x0 + 0.025, 0.46, 0.34, 0.09, values)
        ax.text(x0 + 0.195, 0.425, "$p(d)$", ha="center", va="top", fontsize=6.2, color=INK)

    ax.text(0.50, 0.365, "matched bed transects", ha="center", va="center", fontsize=5.8, color=MUTED)
    draw_cross_section(ax, 0.08, 0.22, 0.35, 0.105, "shallow")
    draw_cross_section(ax, 0.56, 0.22, 0.35, 0.105, "deep")
    ax.text(
        0.50,
        0.075,
        "Equal area can conceal different volume and transfer distance",
        ha="center",
        va="center",
        fontsize=6.1,
        fontweight="bold",
        color=INK,
    )


def gas_markers(ax, x, y, spacing=0.027, radius=0.009):
    for offset, color in zip((-spacing, 0, spacing), (CO2, PURPLE, GREEN)):
        ax.add_patch(Circle((x + offset, y), radius, facecolor=color, edgecolor="white", lw=0.45, zorder=8))


def draw_surface(ax, x0, width, y):
    xs = np.linspace(x0, x0 + width, 80)
    ys = y + 0.0035 * np.sin(np.linspace(0, 4 * np.pi, xs.size))
    ax.plot(xs, ys, color=WATER_EDGE, lw=1.0, zorder=7)


def draw_state(ax, x0, title, state):
    width = 0.22
    bed = 0.24
    surface = 0.76
    interface = 0.49

    if state == "mixed":
        ax.add_patch(Rectangle((x0, bed), width, surface - bed, facecolor=WATER, edgecolor="none"))
    else:
        ax.add_patch(Rectangle((x0, interface), width, surface - interface, facecolor=WATER, edgecolor="none"))
        ax.add_patch(Rectangle((x0, bed), width, interface - bed, facecolor=WATER_DEEP if state == "stratified" else WATER_MID, edgecolor="none"))
        ax.plot([x0, x0 + width], [interface, interface], color=WATER_EDGE, lw=0.9, ls=(0, (3, 2)))

    ax.add_patch(Rectangle((x0, bed), width, surface - bed, facecolor="none", edgecolor=WATER_EDGE, lw=0.9))
    draw_surface(ax, x0, width, surface)
    ax.add_patch(Rectangle((x0, 0.18), width, bed - 0.18, facecolor=SEDIMENT, edgecolor="none"))
    ax.plot([x0, x0 + width], [bed, bed], color="#7C6956", lw=0.9)
    gas_markers(ax, x0 + width / 2, bed + 0.022)

    ax.text(x0 + width / 2, 0.825, title, ha="center", va="center", fontsize=6.6, fontweight="bold", color=INK)

    if state == "mixed":
        arrow(ax, (x0 + 0.055, surface - 0.035), (x0 + 0.055, bed + 0.06), ORANGE, 1.3, 8)
        arrow(ax, (x0 + 0.105, surface - 0.035), (x0 + 0.105, bed + 0.055), BLUE, 1.3, 8)
        arrow(ax, (x0 + 0.15, bed + 0.06), (x0 + 0.15, surface + 0.035), PURPLE, 1.25, 8)
        arrow(ax, (x0 + 0.035, 0.34), (x0 + 0.035, 0.66), WATER_EDGE, 0.8, 6)
        note = "continuous coupling"
    elif state == "stratified":
        arrow(ax, (x0 + 0.055, surface - 0.035), (x0 + 0.055, interface + 0.035), ORANGE, 1.3, 8)
        arrow(ax, (x0 + 0.105, surface - 0.035), (x0 + 0.105, interface + 0.04), BLUE, 1.3, 8)
        ax.plot([x0 + 0.15, x0 + 0.15], [bed + 0.07, interface - 0.035], color=PURPLE, lw=1.25)
        for y in (0.31, 0.37, 0.43):
            ax.add_patch(Rectangle((x0 + 0.14, y - 0.006), 0.02, 0.012, facecolor=PURPLE, edgecolor="none", zorder=8))
        ax.plot([x0 + 0.15, x0 + 0.15], [interface - 0.03, surface - 0.03], color=PURPLE, lw=0.8, ls=(0, (2, 2)), alpha=0.6)
        note = "dissolved-gas retention"
    else:
        ax.add_patch(
            Polygon(
                [(x0 + 0.03, surface - 0.025), (x0 + 0.10, surface - 0.025), (x0 + 0.125, bed + 0.055), (x0 + 0.075, bed + 0.055)],
                closed=True,
                facecolor="#EDF7F9",
                edgecolor="none",
                zorder=2,
            )
        )
        arrow(ax, (x0 + 0.07, surface - 0.04), (x0 + 0.10, bed + 0.07), BLUE, 1.55, 9)
        arrow(ax, (x0 + 0.15, bed + 0.065), (x0 + 0.15, surface + 0.035), PURPLE, 1.65, 9)
        arrow(ax, (x0 + 0.035, surface - 0.035), (x0 + 0.055, interface + 0.02), ORANGE, 1.0, 7)
        note = "reoxygenation + pulse release"

    # Bubble transport can cross a stratified water column independently of mixing.
    arrow(ax, (x0 + 0.195, bed + 0.06), (x0 + 0.195, surface + 0.035), MUTED, 0.55, 5)
    for y in (0.33, 0.52, 0.68):
        ax.add_patch(Circle((x0 + 0.195, y), 0.0075, facecolor="white", edgecolor=MUTED, lw=0.65, zorder=9))
    ax.text(x0 + width / 2, 0.12, note, ha="center", va="center", fontsize=5.8, color=INK)


def draw_states(ax):
    setup(ax)
    panel_label(ax, "b")
    ax.text(0.50, 0.975, "Hydrodynamic state gates vertical transfer", ha="center", va="top", fontsize=8, fontweight="bold")
    ax.text(0.50, 0.91, "the same local depth can connect, isolate or reconnect sediment and atmosphere", ha="center", va="center", fontsize=6.0, color=MUTED)

    legend = [(0.25, ORANGE, "energy"), (0.43, BLUE, "oxidants"), (0.63, PURPLE, "dissolved gas")]
    for x, color, label in legend:
        ax.plot([x - 0.035, x + 0.015], [0.865, 0.865], color=color, lw=2.0, solid_capstyle="round")
        ax.text(x + 0.028, 0.865, label, ha="left", va="center", fontsize=5.8, color=MUTED)

    draw_state(ax, 0.055, "Mixed", "mixed")
    draw_state(ax, 0.39, "Stratified", "stratified")
    draw_state(ax, 0.725, "Ventilating", "ventilating")

    gas_legend = [(0.18, CO2, "CO$_2$"), (0.33, PURPLE, "CH$_4$"), (0.48, GREEN, "N$_2$O")]
    for x, color, label in gas_legend:
        ax.add_patch(Circle((x, 0.045), 0.008, facecolor=color, edgecolor="none"))
        ax.text(x + 0.018, 0.045, label, ha="left", va="center", fontsize=5.8, color=MUTED)
    ax.add_patch(Circle((0.655, 0.045), 0.0075, facecolor="white", edgecolor=MUTED, lw=0.65))
    ax.text(0.673, 0.045, "CH$_4$ bubbles", ha="left", va="center", fontsize=5.8, color=MUTED)


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(183 / 25.4, 100 / 25.4), facecolor="white")
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[0.38, 0.62],
        left=0.045,
        right=0.985,
        bottom=0.07,
        top=0.96,
        wspace=0.13,
    )
    draw_geometry(fig.add_subplot(gs[0, 0]))
    draw_states(fig.add_subplot(gs[0, 1]))

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    labels = [t for ax in fig.axes for t in ax.texts]
    outside = []
    overlaps = []
    for i, label in enumerate(labels):
        box = label.get_window_extent(renderer)
        if box.x0 < 0 or box.y0 < 0 or box.x1 > fig.bbox.width or box.y1 > fig.bbox.height:
            outside.append(label.get_text())
        for other in labels[i+1:]:
            if box.overlaps(other.get_window_extent(renderer)):
                overlaps.append([label.get_text(), other.get_text()])
    assert not outside, outside
    assert not overlaps, overlaps
    (ROOT / 'figure1_layout_checks.json').write_text(json.dumps({
        'off_canvas_labels': outside, 'text_collisions': overlaps,
        'histogram_normalization': 'area fraction; shared vertical scale',
        'transects': 'sampled from the same schematic depth field as each map',
        'bubble_transport': 'separate pathway across all three hydrodynamic states',
        'empirical_data': False,
    }, indent=2), encoding='utf-8')
    fig.savefig(OUT.with_suffix(".svg"), format="svg")
    fig.savefig(OUT.with_suffix(".pdf"), format="pdf")
    fig.savefig(OUT.with_suffix(".png"), format="png", dpi=600)
    fig.savefig(str(OUT) + "_preview.png", format="png", dpi=180)
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
