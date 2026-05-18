#!/usr/bin/env python3
#
# This file is part of "Polynomials with Small Mahler Measure" (PSMM).
# Copyright (C) 2020-2026, Advanpix LLC.
# License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
#
# Author: Pavel Holoborodko <pavel@advanpix.com>
#
"""
plot_root_heatmap.py — heatmap of polynomial roots from the precomputed
roots/ store.

Default subset: the 100 polynomials in AllKnownAdvanpix with smallest M.
Use --subset maxu to restrict to the (a=2, d=3, sign=-1) Max-U family
(the new finds emitted by scan_pn_convergence; per-degree files in
roots/deg-NNNN.txt are required).

Three coordinate systems for the heatmap:
  polar     2D histogram in (arg(r), log|r|).
            Spreads unit-circle clusters along arg, reveals radial
            departure from |r|=1 vertically.
  annulus   2D histogram on a cartesian zoom |r| in [1-w, 1+w].
            Shows the same data in geometric coordinates; useful for
            judging real-world departures from the circle.
  disk      2D histogram on the full unit disk, with optional radial
            magnification of the near-circle band.

Outputs a PNG at --output (default images/root-heatmap-<subset>-<coord>.png).
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm


def gaussian_filter(H: np.ndarray, sigma: float, truncate: float = 4.0) -> np.ndarray:
    """Pure-numpy separable 2D Gaussian blur (avoid the scipy dependency).
    Reflective edge handling via np.convolve(mode='same'); good enough for
    histograms that fade to zero near the borders, which is our case."""
    if sigma <= 0:
        return H
    radius = int(truncate * sigma + 0.5)
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-(x ** 2) / (2.0 * sigma ** 2))
    kernel /= kernel.sum()
    H = np.apply_along_axis(lambda r: np.convolve(r, kernel, mode="same"),
                            axis=1, arr=H)
    H = np.apply_along_axis(lambda c: np.convolve(c, kernel, mode="same"),
                            axis=0, arr=H)
    return H

REPO = Path(__file__).resolve().parent.parent


# Parse a single AllKnownAdvanpix data row into (N, M_str, half_tuple, raw).
def parse_db_row(line: str):
    line = line.rstrip("\n")
    if not line or line.startswith("#"):
        return None
    toks = line.split()
    if len(toks) < 10:
        return None
    try:
        N = int(toks[0])
        M_str = toks[1]
        half = tuple(int(t) for t in toks[9:])
    except ValueError:
        return None
    if len(half) != N // 2 + 1:
        return None
    return (N, M_str, half, line)


def load_db_entries(db_path: Path):
    out = []
    with db_path.open() as f:
        for line in f:
            row = parse_db_row(line)
            if row is not None:
                out.append(row)
    return out


def select_subset(entries, subset: str, limit: int, family_files: list[Path]):
    """Return a sub-list of entries matching the chosen subset."""
    if subset == "top":
        # Sort by M lexicographically — all DB rows have "1." prefix so
        # string comparison agrees with numeric comparison.
        return sorted(entries, key=lambda e: e[1])[:limit]
    if subset == "maxu":
        family_keys = set()
        for fp in family_files:
            if not fp.exists():
                print(f"warn: {fp} missing", file=sys.stderr)
                continue
            with fp.open() as f:
                for line in f:
                    row = parse_db_row(line)
                    if row is not None:
                        family_keys.add((row[0], row[2]))
        matched = [e for e in entries if (e[0], e[2]) in family_keys]
        return matched
    raise SystemExit(f"unknown subset: {subset}")


def load_roots_for_subset(subset_entries, roots_dir: Path):
    """Read each polynomial's roots from roots/deg-NNNN.txt. Returns a flat
    list of complex numbers."""
    # Group by N to scan each per-degree file once.
    by_n = defaultdict(set)
    for (N, _M, half, _raw) in subset_entries:
        by_n[N].add(half)
    all_roots = []
    n_found = 0
    n_missing_file = 0
    n_missing_block = 0
    for N, halves in sorted(by_n.items()):
        fp = roots_dir / f"deg-{N:04d}.txt"
        if not fp.exists():
            n_missing_file += len(halves)
            continue
        remaining = set(halves)
        with fp.open() as f:
            in_block = False
            block_roots = []
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    if in_block and block_roots:
                        all_roots.extend(block_roots)
                        n_found += 1
                    header = line[1:].split()
                    if len(header) >= 10:
                        try:
                            h_half = tuple(int(t) for t in header[9:])
                        except ValueError:
                            h_half = None
                        in_block = h_half in remaining
                        if in_block:
                            remaining.discard(h_half)
                            block_roots = []
                        else:
                            block_roots = []
                    else:
                        in_block = False
                        block_roots = []
                elif not line:
                    if in_block and block_roots:
                        all_roots.extend(block_roots)
                        n_found += 1
                    in_block = False
                    block_roots = []
                elif in_block:
                    toks = line.split()
                    if len(toks) == 2:
                        try:
                            re = float(toks[0])
                            im = float(toks[1])
                            block_roots.append(complex(re, im))
                        except ValueError:
                            pass
            if in_block and block_roots:
                all_roots.extend(block_roots)
                n_found += 1
        n_missing_block += len(remaining)
    print(f"loaded roots from {n_found} polynomials "
          f"({n_missing_file} missing per-degree files, "
          f"{n_missing_block} blocks not found in their files)",
          file=sys.stderr)
    return np.array(all_roots, dtype=np.complex128)


def filter_near_unit(roots: np.ndarray, eps: float):
    """Drop roots with ||r|-1| < eps."""
    if eps <= 0:
        return roots
    return roots[np.abs(np.abs(roots) - 1.0) >= eps]


def load_roots_from_block_file(path: Path) -> np.ndarray:
    """Read a roots-block file (one block per polynomial: '# header' then
    'real imag' rows then blank line). Returns a flat numpy array of
    complex roots."""
    roots = []
    with path.open() as f:
        for line in f:
            line = line.rstrip("\n").strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) == 2:
                try:
                    re = float(toks[0])
                    im = float(toks[1])
                    roots.append(complex(re, im))
                except ValueError:
                    pass
    return np.array(roots, dtype=np.complex128)


def _draw_family_panel(ax, roots: np.ndarray, label: str, color: str,
                       extent: float, point_size: float, alpha: float,
                       limit_M: float | None = None,
                       formula: str | None = None):
    """One side-by-side panel: dashed unit-circle backdrop (drawn first
    so it sits underneath), then scatter the off-circle roots on top
    (zorder explicitly set so scatter wins over the grid and circle)."""
    # Faint background grid first.
    ax.grid(True, alpha=0.2, linewidth=0.4, zorder=0)
    # Faint dashed unit circle, under everything else.
    theta = np.linspace(0, 2 * np.pi, 720)
    ax.plot(np.cos(theta), np.sin(theta),
            color="black", linestyle="--", linewidth=0.7, alpha=0.25,
            zorder=1)
    # Roots on top.
    ax.scatter(roots.real, roots.imag,
               s=point_size, c=color, alpha=alpha, linewidths=0,
               zorder=3)
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\Re(r)$", fontsize=14)
    ax.set_ylabel(r"$\Im(r)$", fontsize=14)
    ax.tick_params(labelsize=12)
    # Two-line title: family-name + formula on top, count + L below.
    top = f"{label}:  {formula}" if formula else label
    bottom = f"{len(roots):,} off-circle roots"
    if limit_M is not None:
        bottom += rf",  $L = {limit_M:.4f}$"
    ax.set_title(f"{top}\n{bottom}", fontsize=15, pad=12)


def plot_zoom_panel(roots: np.ndarray, label: str, color: str,
                    formula: str | None, limit_M: float | None,
                    args, out_path: Path):
    """Two side-by-side panels of ONE family: full unit-disk on the left,
    a zoomed-in window on the right (showing the off-circle root pocket
    in detail). Used for the Max-U family's tight cluster at Re ~ -1.

    Layout: figure-level suptitle has the family identity, formula, count,
    and L. Each panel just gets a short caption ("Full disk" / "Zoom: ...").
    """
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(17, 9),
                                            dpi=args.dpi)
    full_extent = args.disk_extent
    theta = np.linspace(0, 2 * np.pi, 720)

    # Left panel: full unit-disk
    ax_left.grid(True, alpha=0.2, linewidth=0.4, zorder=0)
    ax_left.plot(np.cos(theta), np.sin(theta),
                 color="black", linestyle="--", linewidth=0.7, alpha=0.25,
                 zorder=1)
    ax_left.scatter(roots.real, roots.imag,
                    s=args.point_size, c=color, alpha=args.alpha,
                    linewidths=0, zorder=3)
    ax_left.set_xlim(-full_extent, full_extent)
    ax_left.set_ylim(-full_extent, full_extent)
    ax_left.set_aspect("equal")
    ax_left.set_xlabel(r"$\Re(r)$", fontsize=14)
    ax_left.set_ylabel(r"$\Im(r)$", fontsize=14)
    ax_left.tick_params(labelsize=12)
    ax_left.set_title("Full unit disk", fontsize=15, pad=10)
    if args.zoom_rect:
        rect_x = [args.zoom_xmin, args.zoom_xmax, args.zoom_xmax,
                  args.zoom_xmin, args.zoom_xmin]
        rect_y = [args.zoom_ymin, args.zoom_ymin, args.zoom_ymax,
                  args.zoom_ymax, args.zoom_ymin]
        ax_left.plot(rect_x, rect_y, color="black", linewidth=0.9,
                     alpha=0.6, zorder=2)

    # Right panel: zoom into the off-circle region
    ax_right.grid(True, alpha=0.2, linewidth=0.4, zorder=0)
    ax_right.plot(np.cos(theta), np.sin(theta),
                  color="black", linestyle="--", linewidth=0.7, alpha=0.25,
                  zorder=1)
    ax_right.scatter(roots.real, roots.imag,
                     s=args.zoom_point_size, c=color, alpha=args.alpha,
                     linewidths=0, zorder=3)
    ax_right.set_xlim(args.zoom_xmin, args.zoom_xmax)
    ax_right.set_ylim(args.zoom_ymin, args.zoom_ymax)
    ax_right.set_aspect("equal")
    ax_right.set_xlabel(r"$\Re(r)$", fontsize=14)
    ax_right.set_ylabel(r"$\Im(r)$", fontsize=14)
    ax_right.tick_params(labelsize=12)
    zoom_title = (rf"Zoom:  $\Re(r) \in [{args.zoom_xmin}, {args.zoom_xmax}]$,  "
                  rf"$\Im(r) \in [{args.zoom_ymin}, {args.zoom_ymax}]$")
    ax_right.set_title(zoom_title, fontsize=15, pad=10)

    # Figure-level title carries the family identity.
    if args.title:
        suptitle = args.title
    else:
        top = f"{label}:  {formula}" if formula else label
        bottom = f"{len(roots):,} off-circle roots"
        if limit_M is not None:
            bottom += rf",  $L = {limit_M:.4f}$"
        suptitle = f"{top}\n{bottom}"
    fig.suptitle(suptitle, fontsize=16, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(out_path)
    plt.close()


def plot_two_family(roots_a: np.ndarray, label_a: str, color_a: str,
                    roots_b: np.ndarray, label_b: str, color_b: str,
                    args, out_path: Path):
    """Side-by-side textbook panels: one family per axis, each with its
    own dashed unit-circle backdrop and Boyd-Lawton limit annotation.
    Assumes the caller has already applied --filter-unit if requested."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(17, 9), dpi=args.dpi)
    extent = args.disk_extent
    _draw_family_panel(ax_a, roots_a, label_a, color_a, extent,
                       args.point_size, args.alpha,
                       limit_M=args.family_a_limit,
                       formula=args.family_a_formula or None)
    _draw_family_panel(ax_b, roots_b, label_b, color_b, extent,
                       args.point_size, args.alpha,
                       limit_M=args.family_b_limit,
                       formula=args.family_b_formula or None)
    if args.title:
        fig.suptitle(args.title, fontsize=17, y=0.99)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def plot_polar(roots: np.ndarray, args, out_path: Path):
    arg = np.angle(roots)
    radius_log = np.log10(np.abs(roots))
    # Symmetric y-limit so the unit circle (log = 0) is centered.
    y_abs = max(abs(radius_log.min()), abs(radius_log.max()), 1e-2)
    fig, ax = plt.subplots(figsize=(12, 6), dpi=160)
    H, xedges, yedges, img = ax.hist2d(
        arg, radius_log,
        bins=(args.bins_arg, args.bins_radius),
        range=[[-np.pi, np.pi], [-y_abs, y_abs]],
        norm=LogNorm() if args.log_counts else None,
        cmap=args.cmap,
    )
    ax.axhline(0.0, color="white", linewidth=0.6, alpha=0.7)
    ax.set_xlabel(r"$\arg(r)$")
    ax.set_ylabel(r"$\log_{10}|r|$")
    ax.set_title(args.title or f"Root density — {args.subset} "
                                f"({len(roots):,} roots, {args.coord} coords)")
    plt.colorbar(img, ax=ax, label="root count per bin")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def plot_annulus(roots: np.ndarray, args, out_path: Path):
    width = args.annulus_width
    re = roots.real
    im = roots.imag
    keep = np.abs(np.abs(roots) - 1.0) <= width
    re = re[keep]
    im = im[keep]
    fig, ax = plt.subplots(figsize=(8, 8), dpi=160)
    H, xedges, yedges, img = ax.hist2d(
        re, im,
        bins=args.bins,
        range=[[-1.0 - width, 1.0 + width], [-1.0 - width, 1.0 + width]],
        norm=LogNorm() if args.log_counts else None,
        cmap=args.cmap,
    )
    theta = np.linspace(0, 2 * np.pi, 720)
    ax.plot(np.cos(theta), np.sin(theta),
            color="white", linewidth=0.6, alpha=0.7)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\Re(r)$")
    ax.set_ylabel(r"$\Im(r)$")
    ax.set_title(args.title or f"Root density — {args.subset} "
                                f"({len(roots):,} roots, annulus width {width:.3f})")
    plt.colorbar(img, ax=ax, label="root count per bin")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def get_aesthetic_cmap(name: str):
    """Custom dark-red-to-cream colormap used when --aesthetic is on.
    The reference image Pavel pointed to runs roughly black -> deep blood red
    -> red -> orange -> peach -> cream -> white."""
    if name == "redhot":
        return LinearSegmentedColormap.from_list("redhot", [
            (0.00, "#000000"),
            (0.18, "#1a0000"),
            (0.32, "#5a0a05"),
            (0.50, "#a01810"),
            (0.65, "#e04020"),
            (0.78, "#f08850"),
            (0.88, "#ffcc99"),
            (1.00, "#ffffff"),
        ])
    return name


def plot_disk(roots: np.ndarray, args, out_path: Path):
    # Optional radial magnification: rho_view = sign(rho-1)|rho-1|^alpha + 1
    # where rho = |r| and alpha < 1 stretches near-unit-circle values away
    # from the circle for visibility. With --aesthetic we keep alpha=1
    # (the reference plot doesn't magnify; the colormap+smoothing reveals
    # the structure).
    rho = np.abs(roots)
    theta = np.angle(roots)
    if args.disk_alpha < 1.0:
        rho_view = np.sign(rho - 1.0) * np.abs(rho - 1.0) ** args.disk_alpha + 1.0
    else:
        rho_view = rho
    re = rho_view * np.cos(theta)
    im = rho_view * np.sin(theta)

    extent = args.disk_extent
    bins = args.bins
    H, xedges, yedges = np.histogram2d(
        re, im, bins=bins,
        range=[[-extent, extent], [-extent, extent]],
    )
    if args.smooth_sigma > 0:
        H = gaussian_filter(H, sigma=args.smooth_sigma)

    cmap = get_aesthetic_cmap(args.cmap)
    norm = LogNorm(vmin=max(H.min(), args.vmin_floor)) if args.log_counts else None

    if args.aesthetic:
        # Black background, no axes/labels/colorbar/title -- pure data image.
        fig = plt.figure(figsize=(10, 10), dpi=args.dpi, facecolor="black")
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_facecolor("black")
        ax.imshow(
            H.T, origin="lower",
            extent=[-extent, extent, -extent, extent],
            cmap=cmap, norm=norm,
            interpolation=args.interp,
            aspect="equal",
        )
        ax.set_xlim(-extent, extent)
        ax.set_ylim(-extent, extent)
        ax.axis("off")
        plt.savefig(out_path, facecolor="black",
                    bbox_inches="tight", pad_inches=0)
        plt.close()
        return

    # Scientific style (legacy default for diagnostic runs)
    fig, ax = plt.subplots(figsize=(9, 9), dpi=160)
    img = ax.imshow(
        H.T, origin="lower",
        extent=[-extent, extent, -extent, extent],
        cmap=cmap, norm=norm,
        interpolation=args.interp,
        aspect="equal",
    )
    t = np.linspace(0, 2 * np.pi, 720)
    ax.plot(np.cos(t), np.sin(t),
            color="white", linewidth=0.6, alpha=0.6)
    ax.set_aspect("equal")
    ax.set_xlabel("Re")
    ax.set_ylabel("Im")
    title = (args.title or f"Root density — {args.subset} "
             f"({len(roots):,} roots, disk")
    if args.disk_alpha < 1.0:
        title += f", radial^{args.disk_alpha}"
    title += ")"
    ax.set_title(title)
    plt.colorbar(img, ax=ax, label="root count per bin")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", choices=["top", "maxu"], default="top",
                    help="top = 100 smallest-M polynomials; "
                         "maxu = (a=2,d=3,sign=-1) Max-U family new finds")
    ap.add_argument("--limit", type=int, default=100,
                    help="for --subset top, how many to pick (default 100)")
    ap.add_argument("--coord", choices=["polar", "annulus", "disk"],
                    default="polar")
    ap.add_argument("--filter-unit", type=float, default=0.0,
                    help="drop roots with ||r|-1| < this (default 0 = keep all)")
    ap.add_argument("--annulus-width", type=float, default=0.05,
                    help="half-width around |r|=1 for --coord annulus")
    ap.add_argument("--disk-alpha", type=float, default=0.5,
                    help="radial magnification exponent for --coord disk "
                         "(1.0 = no magnification; smaller = more spread)")
    ap.add_argument("--bins", type=int, default=400,
                    help="bins per axis for cartesian coords "
                         "(use 1500-2500 for aesthetic disk plots)")
    ap.add_argument("--bins-arg", type=int, default=720,
                    help="arg-axis bins for polar coords")
    ap.add_argument("--bins-radius", type=int, default=300,
                    help="radius-axis bins for polar coords")
    ap.add_argument("--cmap", default="magma",
                    help="matplotlib colormap name OR 'redhot' for the "
                         "built-in black->dark-red->cream->white aesthetic")
    ap.add_argument("--aesthetic", action="store_true",
                    help="black background, no axes/title/colorbar, "
                         "tight crop -- for README / paper figures")
    ap.add_argument("--disk-extent", type=float, default=1.3,
                    help="half-width of the disk view (default 1.3)")
    ap.add_argument("--smooth-sigma", type=float, default=0.0,
                    help="Gaussian blur sigma (in bins) applied to the "
                         "2D histogram before colormapping; 1-3 gives a "
                         "smooth glow, 0 keeps hard bin edges (default 0)")
    ap.add_argument("--interp", default="bilinear",
                    help="imshow interpolation (default bilinear; 'nearest' "
                         "for hard pixels)")
    ap.add_argument("--vmin-floor", type=float, default=1.0,
                    help="LogNorm vmin floor so empty bins don't blow up "
                         "the dynamic range (default 1.0)")
    ap.add_argument("--dpi", type=int, default=200,
                    help="output DPI (default 200; bump to 300+ for print)")
    ap.add_argument("--log-counts", action="store_true", default=True,
                    help="use log scale for the colorbar (default on)")
    ap.add_argument("--linear-counts", dest="log_counts", action="store_false",
                    help="use linear scale for the colorbar")
    ap.add_argument("--title", default=None)
    ap.add_argument("--output", default=None,
                    help="output PNG (default: images/root-heatmap-<subset>-<coord>.png)")
    # Two-family overlay mode (textbook style, scatter, dashed unit circle).
    ap.add_argument("--two-family", action="store_true",
                    help="overlay two families' off-circle roots in two "
                         "colors. Use --family-a-* and --family-b-* below.")
    ap.add_argument("--family-a-subset", choices=["top", "maxu"],
                    default="maxu",
                    help="family A is one of the DB-resident subsets "
                         "(default maxu = (2,3,-1) Max-U family)")
    ap.add_argument("--family-a-label", default="(2,3,-1) Max-U family")
    ap.add_argument("--family-a-color", default="#cc1010")
    ap.add_argument("--family-a-limit", type=float, default=1.2554340377272518,
                    help="Boyd-Lawton limit L for family A, annotated under "
                         "the panel title")
    ap.add_argument("--family-a-formula", default=
                    r"$P_N(x) = (x+1)(x^{N-1}+1) - x^{(N-2)/2}\,\Phi_3(x)$",
                    help="LaTeX formula for family A (panel subtitle)")
    ap.add_argument("--family-b-roots", default="",
                    help="path to a roots-block file for family B "
                         "(produced by scan_pn_convergence --roots-output)")
    ap.add_argument("--family-b-label", default="(2,12,-1) family")
    ap.add_argument("--family-b-color", default="#1060d8")
    ap.add_argument("--family-b-limit", type=float, default=1.3091,
                    help="Boyd-Lawton limit L for family B, annotated under "
                         "the panel title")
    ap.add_argument("--family-b-formula", default=
                    r"$P_N(x) = (x+1)(x^{N-1}+1) - x^{(N-4)/2}\,\Phi_{12}(x)$",
                    help="LaTeX formula for family B (panel subtitle)")
    ap.add_argument("--point-size", type=float, default=3.5)
    ap.add_argument("--alpha", type=float, default=0.55)
    # --zoom-panel: one family in two panels (full disk + zoom). For
    # Max-U the default window crops to the Re ~ -1 crescent.
    ap.add_argument("--zoom-panel", action="store_true",
                    help="render one family in two side-by-side panels: "
                         "full unit-disk on the left, zoomed-in window "
                         "(--zoom-xmin/xmax/ymin/ymax) on the right.")
    ap.add_argument("--zoom-xmin", type=float, default=-1.10)
    ap.add_argument("--zoom-xmax", type=float, default=-0.90)
    ap.add_argument("--zoom-ymin", type=float, default=-0.45)
    ap.add_argument("--zoom-ymax", type=float, default=+0.45)
    ap.add_argument("--zoom-point-size", type=float, default=1.5)
    ap.add_argument("--zoom-rect", action="store_true",
                    help="draw the zoom window outline on the left panel "
                         "(default off; the right-panel title carries the "
                         "range, which is usually enough)")
    args = ap.parse_args()

    entries = load_db_entries(REPO / "AllKnownAdvanpix")
    print(f"DB: {len(entries)} entries", file=sys.stderr)

    family_files = [
        REPO / "doc" / "new_finds_d_only.txt",
        REPO / "doc" / "new_finds_d_only_pn_extended.txt",
        # part-2 entries don't have roots in roots/ yet — skip
    ]

    # ---- Single-family zoom-panel branch ---------------------------
    if args.zoom_panel:
        sub = select_subset(entries, args.subset, args.limit, family_files)
        roots = load_roots_for_subset(sub, REPO / "roots")
        if args.filter_unit > 0:
            before = len(roots)
            roots = filter_near_unit(roots, args.filter_unit)
            print(f"unit-circle filter eps={args.filter_unit}: "
                  f"{before} -> {len(roots)} roots", file=sys.stderr)
        if len(roots) == 0:
            raise SystemExit("no roots after filtering")
        out_path = (Path(args.output) if args.output
                    else REPO / "images" /
                         f"root-heatmap-{args.subset}-zoom.png")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Use family-A label/color/limit/formula by default for the
        # Max-U case; CLI override possible via those flags.
        plot_zoom_panel(roots, args.family_a_label, args.family_a_color,
                        args.family_a_formula or None,
                        args.family_a_limit, args, out_path)
        print(f"wrote {out_path}", file=sys.stderr)
        return

    # ---- Two-family overlay branch ---------------------------------
    if args.two_family:
        if not args.family_b_roots:
            raise SystemExit("--two-family requires --family-b-roots PATH")

        # Family A: DB-resident subset.
        a_subset = select_subset(entries, args.family_a_subset,
                                 args.limit, family_files)
        a_roots = load_roots_for_subset(a_subset, REPO / "roots")
        # Family B: read from roots-block file.
        b_path = Path(args.family_b_roots)
        if not b_path.exists():
            raise SystemExit(f"family-b roots file missing: {b_path}")
        b_roots = load_roots_from_block_file(b_path)

        if args.filter_unit > 0:
            a_before, b_before = len(a_roots), len(b_roots)
            a_roots = filter_near_unit(a_roots, args.filter_unit)
            b_roots = filter_near_unit(b_roots, args.filter_unit)
            print(f"unit-circle filter eps={args.filter_unit}: "
                  f"A {a_before}->{len(a_roots)}, "
                  f"B {b_before}->{len(b_roots)}",
                  file=sys.stderr)
        if len(a_roots) == 0 and len(b_roots) == 0:
            raise SystemExit("no roots in either family after filtering")

        out_path = (Path(args.output) if args.output
                    else REPO / "work" / "root-heatmap-two-family.png")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plot_two_family(a_roots, args.family_a_label, args.family_a_color,
                        b_roots, args.family_b_label, args.family_b_color,
                        args, out_path)
        print(f"wrote {out_path}", file=sys.stderr)
        return
    # ---- End two-family branch -------------------------------------

    subset_entries = select_subset(entries, args.subset, args.limit, family_files)
    print(f"selected subset {args.subset}: {len(subset_entries)} entries",
          file=sys.stderr)
    if not subset_entries:
        raise SystemExit("empty subset")

    roots = load_roots_for_subset(subset_entries, REPO / "roots")
    if args.filter_unit > 0:
        before = len(roots)
        roots = filter_near_unit(roots, args.filter_unit)
        print(f"unit-circle filter eps={args.filter_unit}: "
              f"{before} -> {len(roots)} roots", file=sys.stderr)
    if len(roots) == 0:
        raise SystemExit("no roots after filtering")

    out_path = (Path(args.output) if args.output
                else REPO / "work" /
                     f"root-heatmap-{args.subset}-{args.coord}.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.coord == "polar":
        plot_polar(roots, args, out_path)
    elif args.coord == "annulus":
        plot_annulus(roots, args, out_path)
    elif args.coord == "disk":
        plot_disk(roots, args, out_path)

    print(f"wrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
