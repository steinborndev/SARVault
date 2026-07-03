"""Header branding: the SARVault logo rendered as a flush-left inline <img>.

The logo is an SVG data-URI. Streamlit places the <img> element at the page's
content margin, but if the SVG has empty space baked into the left of its viewBox
(a common artefact of re-exporting a logo), the *visible* icon ends up inset from
that margin even though the element isn't — so it no longer lines up with the text
below it. ``flush_left_svg`` trims that whitespace so the leftmost drawn element sits
at the SVG's left edge, making the icon flush regardless of how the asset was exported.
"""

import base64
import re
from pathlib import Path

_LOGO = Path(__file__).resolve().parents[1] / "assets" / "logo.svg"


def flush_left_svg(svg: str) -> str:
    """Return the SVG with its leftmost geometry pinned to the left edge.

    Moves the viewBox's min-x to the leftmost drawn coordinate (across rect/text ``x``,
    circle ``cx − r`` and polygon points) and drops the root ``width``/``height`` so the
    ``<img height=…>`` plus the trimmed viewBox drive the aspect ratio cleanly. Handles
    both padding styles: content shifted right within a ``0 …`` viewBox, and a negative
    viewBox min-x. Returns the input unchanged (bar the width/height strip) when it can't
    find geometry, so an unparseable asset degrades gracefully rather than breaking.
    """
    m = re.search(r'viewBox="([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)"', svg)
    if not m:
        return svg
    vx, vy, vw, vh = map(float, m.groups())

    xs: list[float] = [float(v) for v in re.findall(r'\bx="([\d.\-]+)"', svg)]
    for cx, r in re.findall(r'cx="([\d.\-]+)"[^>]*?\br="([\d.\-]+)"', svg):
        xs.append(float(cx) - float(r))
    for pts in re.findall(r'points="([^"]+)"', svg):
        xs += [float(a) for a, _ in re.findall(r"(-?[\d.]+)[, ](-?[\d.]+)", pts)]
    if not xs:
        return svg

    # Drop root width/height so a trimmed viewBox isn't fighting a stale intrinsic size.
    def _strip_wh(mm: re.Match) -> str:
        tag = re.sub(r'\s+width="[^"]*"', "", mm.group(0), count=1)
        return re.sub(r'\s+height="[^"]*"', "", tag, count=1)

    svg = re.sub(r"<svg\b[^>]*>", _strip_wh, svg, count=1)

    minx = min(xs)
    if minx > vx + 0.5:
        new_vx, new_vw = minx, vw - (minx - vx)
        svg = re.sub(
            r'viewBox="[^"]+"', f'viewBox="{new_vx:g} {vy:g} {new_vw:g} {vh:g}"', svg, count=1
        )
    return svg


def logo_html(height: int = 84) -> str:
    """SARVault logo as a flush-left, block-level inline <img> (SVG data-URI)."""
    svg = flush_left_svg(_LOGO.read_text())
    b64 = base64.b64encode(svg.encode()).decode()
    return (
        f'<img src="data:image/svg+xml;base64,{b64}" height="{height}" '
        'style="display:block;margin:0">'
    )
