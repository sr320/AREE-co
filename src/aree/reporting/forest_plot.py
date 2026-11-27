"""Pure-Python SVG forest plots for evidence cards, Quarto synthesis pages, and the Streamlit app."""

import html
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple


_MIN_ROW_HEIGHT = 24
_HEADER_HEIGHT = 44
_AXIS_HEIGHT = 40
_LEFT_COL_WIDTH = 230
_PLOT_WIDTH = 290
_RIGHT_COL_WIDTH = 220
_PADDING_X = 16


def _is_finite(value) -> bool:
    if value is None or value == "":
        return False
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _truncate(text: str, max_len: int = 32) -> str:
    text = str(text or "")
    if len(text) <= max_len:
        return text
    return text[: max_len - 1] + "…"


def _nice_ticks(x_min: float, x_max: float) -> Tuple[float, float, List[float]]:
    if x_min >= x_max:
        span = abs(x_min) * 0.5 or 1.0
        x_min, x_max = x_min - span, x_max + span
    pad = (x_max - x_min) * 0.08
    lo = min(x_min - pad, 0.0)
    hi = max(x_max + pad, 0.0)
    if hi - lo < 1e-6:
        lo, hi = -1.0, 1.0

    raw_step = (hi - lo) / 4.0
    mag = 10.0 ** math.floor(math.log10(raw_step))
    norm = raw_step / mag
    step_unit = 1.0 if norm <= 1.5 else 2.0 if norm <= 3.5 else 5.0 if norm <= 7.5 else 10.0
    step = step_unit * mag
    start = math.floor(lo / step) * step
    end = math.ceil(hi / step) * step
    ticks = []
    cur = start
    while cur <= end + step * 1e-6:
        val = 0.0 if abs(cur) < step * 1e-9 else round(cur, 6)
        ticks.append(val)
        cur += step
    return start, end, ticks


def _row_bounds(row: Dict[str, object]) -> Optional[Tuple[float, Optional[float], Optional[float], Optional[float]]]:
    if not _is_finite(row.get("effect")):
        return None
    eff = float(row["effect"])
    ci_lo = float(row["ci_low"]) if _is_finite(row.get("ci_low")) else None
    ci_hi = float(row["ci_high"]) if _is_finite(row.get("ci_high")) else None
    se = float(row["se"]) if _is_finite(row.get("se")) and float(row["se"]) > 0 else None
    if (ci_lo is None or ci_hi is None) and se is not None:
        ci_lo = eff - 1.96 * se
        ci_hi = eff + 1.96 * se
    return eff, ci_lo, ci_hi, se


def _format_stat(eff: float, ci_lo: Optional[float], ci_hi: Optional[float], subtext: str = "") -> str:
    if ci_lo is not None and ci_hi is not None:
        base = "{:+.2f} [{:+.2f}, {:+.2f}]".format(eff, ci_lo, ci_hi)
    else:
        base = "{:+.2f} (SE missing)".format(eff)
    return "{} {}".format(base, subtext).strip() if subtext else base


def _effect_color(eff: float, has_ci: bool) -> str:
    if not has_ci:
        return "#6b7280"
    if eff > 0:
        return "#b91c1c"
    if eff < 0:
        return "#1d4ed8"
    return "#4b5563"


def render_forest_svg(
    rows: List[Dict[str, object]],
    pooled: Optional[Dict[str, object]] = None,
    title: str = "Effect-size forest plot",
    xlabel: str = "Effect size (log2 fold change)",
    width: int = 760,
) -> str:
    """Render a deterministic standalone SVG forest plot for a feature across studies or strata."""
    parsed = []
    values = [0.0]
    for row in rows:
        bounds = _row_bounds(row)
        if bounds is None:
            continue
        eff, ci_lo, ci_hi, se = bounds
        values.append(eff)
        if ci_lo is not None and ci_hi is not None:
            values.extend([ci_lo, ci_hi])
        parsed.append((str(row.get("label", "")), eff, ci_lo, ci_hi, se, str(row.get("subtext", ""))))

    pooled_bounds = _row_bounds(pooled) if pooled else None
    if pooled_bounds is not None:
        values.append(pooled_bounds[0])
        if pooled_bounds[1] is not None and pooled_bounds[2] is not None:
            values.extend([pooled_bounds[1], pooled_bounds[2]])

    if not parsed and pooled_bounds is None:
        return (
            '<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="64" role="img" '
            'aria-label="{t}">'
            '<rect width="{w}" height="64" fill="#ffffff" stroke="#e5e7eb"/>'
            '<text x="16" y="37" font-family="system-ui, -apple-system, sans-serif" font-size="12" fill="#6b7280">'
            'No numeric effect sizes available for plotting.</text></svg>'
        ).format(w=width, t=html.escape(title))

    x_min, x_max, ticks = _nice_ticks(min(values), max(values))
    plot_x0 = _PADDING_X + _LEFT_COL_WIDTH
    plot_x1 = max(plot_x0 + 180, width - _PADDING_X - _RIGHT_COL_WIDTH)
    right_x = plot_x1 + 14

    def scale_x(val: float) -> float:
        return plot_x0 + ((val - x_min) / (x_max - x_min)) * (plot_x1 - plot_x0)

    n_rows = len(parsed) + (1 if pooled_bounds is not None else 0)
    height = _HEADER_HEIGHT + max(1, n_rows) * _MIN_ROW_HEIGHT + _AXIS_HEIGHT + (8 if pooled_bounds else 0)

    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        'role="img" aria-label="{title}">'.format(w=width, h=height, title=html.escape(title)),
        "<title>{}</title>".format(html.escape(title)),
        '<style>'
        'text { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }'
        '.hdr { font-size: 13px; font-weight: 600; fill: #111827; }'
        '.colhdr { font-size: 11px; font-weight: 600; fill: #4b5563; }'
        '.lbl { font-size: 11px; fill: #1f2937; }'
        '.stat { font-size: 11px; fill: #374151; font-variant-numeric: tabular-nums; }'
        '.tick { font-size: 10px; fill: #6b7280; }'
        '</style>',
        '<rect width="{w}" height="{h}" rx="6" fill="#ffffff" stroke="#e5e7eb"/>'.format(w=width, h=height),
        '<text x="{x}" y="20" class="hdr">{t}</text>'.format(x=_PADDING_X, t=html.escape(_truncate(title, 76))),
        '<text x="{x}" y="36" class="colhdr">Contrast / Stratum</text>'.format(x=_PADDING_X),
        '<text x="{x}" y="36" class="colhdr" text-anchor="middle">{lbl}</text>'.format(
            x=(plot_x0 + plot_x1) / 2.0, lbl=html.escape(xlabel)
        ),
        '<text x="{x}" y="36" class="colhdr">Effect [95% CI]</text>'.format(x=right_x),
        '<line x1="{p}" y1="42" x2="{r}" y2="42" stroke="#d1d5db" stroke-width="1"/>'.format(
            p=_PADDING_X, r=width - _PADDING_X
        ),
    ]

    zero_x = scale_x(0.0)
    axis_y = height - _AXIS_HEIGHT
    lines.append(
        '<line x1="{zx:.2f}" y1="44" x2="{zx:.2f}" y2="{ay}" stroke="#9ca3af" '
        'stroke-width="1" stroke-dasharray="3,3"/>'.format(zx=zero_x, ay=axis_y)
    )

    se_vals = [item[4] for item in parsed if item[4] is not None and item[4] > 0]
    max_inv_var = max((1.0 / (s * s) for s in se_vals), default=1.0)

    y = _HEADER_HEIGHT + (_MIN_ROW_HEIGHT / 2.0)
    for idx, (label, eff, ci_lo, ci_hi, se, subtext) in enumerate(parsed):
        if idx % 2 == 1:
            lines.append(
                '<rect x="{x}" y="{ry:.1f}" width="{rw}" height="{rh}" fill="#f9fafb"/>'.format(
                    x=_PADDING_X, ry=y - _MIN_ROW_HEIGHT / 2.0, rw=width - 2 * _PADDING_X, rh=_MIN_ROW_HEIGHT
                )
            )
        has_ci = ci_lo is not None and ci_hi is not None
        color = _effect_color(eff, has_ci)
        cx = scale_x(eff)
        lines.append(
            '<text x="{x}" y="{ty:.1f}" class="lbl">{t}</text>'.format(
                x=_PADDING_X + 2, ty=y + 4, t=html.escape(_truncate(label, 31))
            )
        )
        if has_ci:
            x1 = scale_x(ci_lo)
            x2 = scale_x(ci_hi)
            lines.append(
                '<line x1="{x1:.2f}" y1="{y:.1f}" x2="{x2:.2f}" y2="{y:.1f}" '
                'stroke="{c}" stroke-width="1.6"/>'.format(x1=x1, x2=x2, y=y, c=color)
            )
            inv_var = (1.0 / (se * se)) if se else max_inv_var * 0.4
            half = 3.2 + 3.0 * math.sqrt(min(1.0, inv_var / max_inv_var))
            lines.append(
                '<rect x="{rx:.2f}" y="{ry:.2f}" width="{sz:.2f}" height="{sz:.2f}" '
                'fill="{c}" rx="1"/>'.format(rx=cx - half, ry=y - half, sz=half * 2.0, c=color)
            )
        else:
            lines.append(
                '<circle cx="{cx:.2f}" cy="{y:.1f}" r="3.5" fill="#ffffff" '
                'stroke="#6b7280" stroke-width="1.5"/>'.format(cx=cx, y=y)
            )
        stat_str = _format_stat(eff, ci_lo, ci_hi, subtext=subtext)
        lines.append(
            '<text x="{x}" y="{ty:.1f}" class="stat">{t}</text>'.format(
                x=right_x, ty=y + 4, t=html.escape(_truncate(stat_str, 34))
            )
        )
        y += _MIN_ROW_HEIGHT

    if pooled_bounds is not None and pooled is not None:
        p_eff, p_lo, p_hi, _ = pooled_bounds
        div_y = y - _MIN_ROW_HEIGHT / 2.0 + 2
        lines.append(
            '<line x1="{x1}" y1="{dy:.1f}" x2="{x2}" y2="{dy:.1f}" stroke="#d1d5db" stroke-width="1"/>'.format(
                x1=_PADDING_X, x2=width - _PADDING_X, dy=div_y
            )
        )
        py = y + 6
        p_label = str(pooled.get("label", "Pooled RE"))
        p_sub = str(pooled.get("subtext", ""))
        px = scale_x(p_eff)
        pl = scale_x(p_lo if p_lo is not None else p_eff)
        pr = scale_x(p_hi if p_hi is not None else p_eff)
        p_color = "#7c2d12" if p_eff > 0 else "#1e3a8a" if p_eff < 0 else "#374151"
        lines.append(
            '<text x="{x}" y="{ty:.1f}" class="lbl" font-weight="600">{t}</text>'.format(
                x=_PADDING_X + 2, ty=py + 4, t=html.escape(_truncate(p_label, 31))
            )
        )
        pts = "{l:.2f},{y:.1f} {c:.2f},{t:.1f} {r:.2f},{y:.1f} {c:.2f},{b:.1f}".format(
            l=pl, r=pr, c=px, y=py, t=py - 6, b=py + 6
        )
        lines.append('<polygon points="{pts}" fill="{c}" opacity="0.9"/>'.format(pts=pts, c=p_color))
        lines.append(
            '<text x="{x}" y="{ty:.1f}" class="stat" font-weight="600">{t}</text>'.format(
                x=right_x, ty=py + 4, t=html.escape(_truncate(_format_stat(p_eff, p_lo, p_hi, p_sub), 34))
            )
        )

    lines.append(
        '<line x1="{x0}" y1="{ay}" x2="{x1}" y2="{ay}" stroke="#9ca3af" stroke-width="1"/>'.format(
            x0=plot_x0, x1=plot_x1, ay=axis_y
        )
    )
    for tv in ticks:
        tx = scale_x(tv)
        if tx < plot_x0 - 4 or tx > plot_x1 + 4:
            continue
        lines.append(
            '<line x1="{tx:.2f}" y1="{ay}" x2="{tx:.2f}" y2="{ay2}" stroke="#9ca3af" stroke-width="1"/>'.format(
                tx=tx, ay=axis_y, ay2=axis_y + 4
            )
        )
        lbl = "0" if tv == 0 else "{:.2g}".format(tv)
        lines.append(
            '<text x="{tx:.2f}" y="{ty}" class="tick" text-anchor="middle">{lbl}</text>'.format(
                tx=tx, ty=axis_y + 16, lbl=html.escape(lbl)
            )
        )
    lines.append(
        '<text x="{cx:.1f}" y="{ty}" class="tick" text-anchor="middle">{lbl}</text>'.format(
            cx=(plot_x0 + plot_x1) / 2.0, ty=axis_y + 30, lbl=html.escape(xlabel)
        )
    )
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def write_forest_plot_svg(
    output_path: Path,
    rows: List[Dict[str, object]],
    pooled: Optional[Dict[str, object]] = None,
    title: str = "Effect-size forest plot",
    xlabel: str = "Effect size (log2 fold change)",
    width: int = 760,
) -> Path:
    """Write a pure-Python SVG forest plot to ``output_path`` and return the path."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    svg = render_forest_svg(rows=rows, pooled=pooled, title=title, xlabel=xlabel, width=width)
    output_path.write_text(svg, encoding="utf-8")
    return output_path
