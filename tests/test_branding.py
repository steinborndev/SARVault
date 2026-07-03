"""Tests for the header logo rendering (dashboard/branding.py).

The key behaviour: the logo must sit flush-left regardless of whether the source SVG
has empty space baked into the left of its viewBox, so the icon lines up with the text
below it. flush_left_svg pins the leftmost drawn element to the SVG's left edge.
"""

import re

from dashboard import branding


def _viewbox(svg: str) -> str:
    return re.search(r'viewBox="([^"]+)"', svg).group(1)


def test_flush_logo_viewbox_unchanged_but_wh_stripped():
    # A logo whose content already starts at x=0 keeps its viewBox; the root
    # width/height are stripped so a later <img height=…> sizes it cleanly.
    svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 244 64" width="244" height="64">'\
          '<rect x="0" y="8" width="48" height="48"/></svg>'
    out = branding.flush_left_svg(svg)
    assert _viewbox(out) == "0 0 244 64"
    assert "width=" not in out.split(">", 1)[0] and "height=" not in out.split(">", 1)[0]


def test_content_shifted_right_is_trimmed():
    # Content shifted right within a 0-based viewBox leaves left whitespace; trim it.
    svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 244 64" width="244" height="64">'\
          '<rect x="21" y="8" width="48" height="48"/>'\
          '<circle cx="45" cy="32" r="4.5"/>'\
          '<polygon points="35,17 48,24 48,40 35,47"/></svg>'
    out = branding.flush_left_svg(svg)
    # leftmost element is the rect at x=21 -> new min-x 21, width 244-21=223
    assert _viewbox(out) == "21 0 223 64"


def test_negative_viewbox_padding_is_trimmed():
    # Padding expressed as a negative viewBox min-x (content at x=0) is also removed.
    svg = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="-20 0 264 64" width="264" height="64">'\
          '<rect x="0" y="8" width="48" height="48"/></svg>'
    assert _viewbox(branding.flush_left_svg(svg)) == "0 0 244 64"


def test_no_viewbox_returned_unchanged():
    svg = '<svg xmlns="http://www.w3.org/2000/svg"><rect x="5" y="0" width="1" height="1"/></svg>'
    assert branding.flush_left_svg(svg) == svg


def test_real_logo_is_flush_and_renders():
    # The shipped asset is already flush; it must survive round-trip and render to an <img>.
    html = branding.logo_html(height=84)
    assert html.startswith('<img src="data:image/svg+xml;base64,')
    assert 'height="84"' in html and "display:block" in html
