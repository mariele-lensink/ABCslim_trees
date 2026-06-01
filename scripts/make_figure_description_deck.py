#!/usr/bin/env python3
"""Build a PowerPoint deck from figure_png_inventory_descriptions.txt.

This intentionally uses only the Python standard library so the deck can be
created on systems without python-pptx installed.
"""

from __future__ import annotations

import html
import re
import struct
import zipfile
from pathlib import Path
from textwrap import wrap


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "figure_png_inventory_descriptions.txt"
OUT = ROOT / "figure_png_inventory_descriptions_fixed.pptx"

EMU_PER_INCH = 914400
SLIDE_W = 13.333333
SLIDE_H = 7.5


def emu(inches: float) -> int:
    return int(round(inches * EMU_PER_INCH))


def xml_text(text: str) -> str:
    return html.escape(text, quote=False)


def parse_report() -> list[tuple[str, str]]:
    text = REPORT.read_text(encoding="utf-8")
    pattern = re.compile(
        r"^Path: (?P<path>.+?)\nDescription: (?P<desc>.*?)(?=\n\nPath: |\n\nFigures in |\Z)",
        re.MULTILINE | re.DOTALL,
    )
    entries: list[tuple[str, str]] = []
    for match in pattern.finditer(text):
        path = match.group("path").strip()
        desc = " ".join(match.group("desc").split())
        entries.append((path, desc))
    return entries


def png_size(path: Path) -> tuple[int, int]:
    with path.open("rb") as fh:
        sig = fh.read(24)
    if sig[:8] != b"\x89PNG\r\n\x1a\n" or sig[12:16] != b"IHDR":
        raise ValueError(f"Not a PNG file: {path}")
    return struct.unpack(">II", sig[16:24])


def fit_box(img_w: int, img_h: int, box_w: float, box_h: float) -> tuple[float, float]:
    ratio = min(box_w / img_w, box_h / img_h)
    return img_w * ratio, img_h * ratio


def paragraph_xml(text: str, size: int = 1200, bold: bool = False) -> str:
    safe = xml_text(text)
    b = "<a:b/>" if bold else ""
    return (
        "<a:p><a:r><a:rPr lang=\"en-US\" sz=\""
        f"{size}\">{b}</a:rPr><a:t>{safe}</a:t></a:r></a:p>"
    )


def textbox_xml(shape_id: int, name: str, x: float, y: float, w: float, h: float, paragraphs: list[str]) -> str:
    body = "".join(paragraphs) or paragraph_xml("")
    return f"""
      <p:sp>
        <p:nvSpPr>
          <p:cNvPr id="{shape_id}" name="{xml_text(name)}"/>
          <p:cNvSpPr txBox="1"/>
          <p:nvPr/>
        </p:nvSpPr>
        <p:spPr>
          <a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(w)}" cy="{emu(h)}"/></a:xfrm>
          <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
          <a:noFill/>
          <a:ln><a:noFill/></a:ln>
        </p:spPr>
        <p:txBody>
          <a:bodyPr wrap="square" anchor="t"/>
          <a:lstStyle/>
          {body}
        </p:txBody>
      </p:sp>
    """


def picture_xml(shape_id: int, rel_id: str, name: str, x: float, y: float, w: float, h: float) -> str:
    return f"""
      <p:pic>
        <p:nvPicPr>
          <p:cNvPr id="{shape_id}" name="{xml_text(name)}"/>
          <p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr>
          <p:nvPr/>
        </p:nvPicPr>
        <p:blipFill>
          <a:blip r:embed="{rel_id}"/>
          <a:stretch><a:fillRect/></a:stretch>
        </p:blipFill>
        <p:spPr>
          <a:xfrm><a:off x="{emu(x)}" y="{emu(y)}"/><a:ext cx="{emu(w)}" cy="{emu(h)}"/></a:xfrm>
          <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
        </p:spPr>
      </p:pic>
    """


def slide_xml(content: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
       xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
       xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
  <p:cSld>
    <p:spTree>
      <p:nvGrpSpPr>
        <p:cNvPr id="1" name=""/>
        <p:cNvGrpSpPr/>
        <p:nvPr/>
      </p:nvGrpSpPr>
      <p:grpSpPr>
        <a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm>
      </p:grpSpPr>
      {content}
    </p:spTree>
  </p:cSld>
  <p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr>
</p:sld>
"""


def rels_xml(rels: list[tuple[str, str, str]]) -> str:
    items = "\n".join(
        f'  <Relationship Id="{rid}" Type="{typ}" Target="{target}"/>'
        for rid, typ, target in rels
    )
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
{items}
</Relationships>
"""


def make_content_types(slide_count: int) -> str:
    overrides = [
        '<Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/>',
        '<Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/>',
        '<Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/>',
        '<Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/>',
        '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>',
        '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>',
    ]
    for idx in range(1, slide_count + 1):
        overrides.append(
            f'<Override PartName="/ppt/slides/slide{idx}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/>'
        )
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Default Extension="png" ContentType="image/png"/>
  {chr(10).join(overrides)}
</Types>
"""


THEME_XML = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<a:theme xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" name="Office Theme">
  <a:themeElements>
    <a:clrScheme name="Office"><a:dk1><a:srgbClr val="000000"/></a:dk1><a:lt1><a:srgbClr val="FFFFFF"/></a:lt1><a:dk2><a:srgbClr val="1F1F1F"/></a:dk2><a:lt2><a:srgbClr val="F2F2F2"/></a:lt2><a:accent1><a:srgbClr val="4472C4"/></a:accent1><a:accent2><a:srgbClr val="ED7D31"/></a:accent2><a:accent3><a:srgbClr val="A5A5A5"/></a:accent3><a:accent4><a:srgbClr val="FFC000"/></a:accent4><a:accent5><a:srgbClr val="5B9BD5"/></a:accent5><a:accent6><a:srgbClr val="70AD47"/></a:accent6><a:hlink><a:srgbClr val="0563C1"/></a:hlink><a:folHlink><a:srgbClr val="954F72"/></a:folHlink></a:clrScheme>
    <a:fontScheme name="Office"><a:majorFont><a:latin typeface="Aptos Display"/><a:ea typeface=""/><a:cs typeface=""/></a:majorFont><a:minorFont><a:latin typeface="Aptos"/><a:ea typeface=""/><a:cs typeface=""/></a:minorFont></a:fontScheme>
    <a:fmtScheme name="Office"><a:fillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:fillStyleLst><a:lnStyleLst><a:ln w="6350"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln><a:ln w="12700"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln><a:ln w="19050"><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:ln></a:lnStyleLst><a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst><a:bgFillStyleLst><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill><a:solidFill><a:schemeClr val="phClr"/></a:solidFill></a:bgFillStyleLst></a:fmtScheme>
  </a:themeElements>
</a:theme>
"""


SLIDE_LAYOUT_XML = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldLayout xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" type="blank" preserve="1">
  <p:cSld name="Blank"><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr></p:spTree></p:cSld>
  <p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr>
</p:sldLayout>
"""


SLIDE_MASTER_XML = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:sldMaster xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
  <p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr></p:spTree></p:cSld>
  <p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" folHlink="folHlink"/>
  <p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/></p:sldLayoutIdLst>
  <p:txStyles><p:titleStyle/><p:bodyStyle/><p:otherStyle/></p:txStyles>
</p:sldMaster>
"""


def build() -> None:
    entries = parse_report()
    if not entries:
        raise SystemExit(f"No figure entries found in {REPORT}")

    slide_count = len(entries) + 1
    with zipfile.ZipFile(OUT, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", make_content_types(slide_count))
        zf.writestr(
            "_rels/.rels",
            rels_xml(
                [
                    ("rId1", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument", "ppt/presentation.xml"),
                    ("rId2", "http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties", "docProps/core.xml"),
                    ("rId3", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties", "docProps/app.xml"),
                ]
            ),
        )
        zf.writestr(
            "docProps/core.xml",
            """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"
  xmlns:dc="http://purl.org/dc/elements/1.1/"
  xmlns:dcterms="http://purl.org/dc/terms/"
  xmlns:dcmitype="http://purl.org/dc/dcmitype/"
  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <dc:title>ABCslim Trees Figure Inventory</dc:title>
  <dc:creator>Codex</dc:creator>
</cp:coreProperties>
""",
        )
        zf.writestr(
            "docProps/app.xml",
            f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
  xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">
  <Application>Codex OpenXML generator</Application>
  <Slides>{slide_count}</Slides>
</Properties>
""",
        )

        slide_rels = []
        for idx in range(1, slide_count + 1):
            slide_rels.append(
                (
                    f"rId{idx + 1}",
                    "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide",
                    f"slides/slide{idx}.xml",
                )
            )
        slide_rels.insert(
            0,
            (
                "rId1",
                "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster",
                "slideMasters/slideMaster1.xml",
            ),
        )
        zf.writestr("ppt/_rels/presentation.xml.rels", rels_xml(slide_rels))
        sld_ids = "\n".join(
            f'    <p:sldId id="{255 + idx}" r:id="rId{idx + 1}"/>'
            for idx in range(1, slide_count + 1)
        )
        zf.writestr(
            "ppt/presentation.xml",
            f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<p:presentation xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
  xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
  xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">
  <p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst>
  <p:sldIdLst>
{sld_ids}
  </p:sldIdLst>
  <p:sldSz cx="{emu(SLIDE_W)}" cy="{emu(SLIDE_H)}" type="wide"/>
  <p:notesSz cx="{emu(10)}" cy="{emu(7.5)}"/>
  <p:defaultTextStyle>
    <a:defPPr><a:defRPr lang="en-US"/></a:defPPr>
  </p:defaultTextStyle>
</p:presentation>
""",
        )
        zf.writestr("ppt/theme/theme1.xml", THEME_XML)
        zf.writestr("ppt/slideLayouts/slideLayout1.xml", SLIDE_LAYOUT_XML)
        zf.writestr("ppt/slideMasters/slideMaster1.xml", SLIDE_MASTER_XML)
        zf.writestr(
            "ppt/slideMasters/_rels/slideMaster1.xml.rels",
            rels_xml(
                [
                    ("rId1", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout", "../slideLayouts/slideLayout1.xml"),
                    ("rId2", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/theme", "../theme/theme1.xml"),
                ]
            ),
        )
        zf.writestr(
            "ppt/slideLayouts/_rels/slideLayout1.xml.rels",
            rels_xml(
                [
                    ("rId1", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideMaster", "../slideMasters/slideMaster1.xml"),
                ]
            ),
        )

        cover = textbox_xml(
            2,
            "Title",
            0.7,
            1.2,
            12.0,
            1.0,
            [paragraph_xml("ABCslim Trees PNG Figure Inventory", 2800, True)],
        ) + textbox_xml(
            3,
            "Subtitle",
            0.75,
            2.4,
            11.8,
            3.4,
            [
                paragraph_xml(f"{len(entries)} PNG figures from figures/ and output/ABC_results/", 1700),
                paragraph_xml("Each slide contains the figure path, the PNG, and its provenance/description text.", 1500),
                paragraph_xml("Descriptions were generated from figure_png_inventory_descriptions.txt.", 1300),
            ],
        )
        zf.writestr("ppt/slides/slide1.xml", slide_xml(cover))
        zf.writestr(
            "ppt/slides/_rels/slide1.xml.rels",
            rels_xml(
                [
                    ("rIdLayout", "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout", "../slideLayouts/slideLayout1.xml"),
                ]
            ),
        )

        for entry_idx, (path_str, desc) in enumerate(entries, start=2):
            img_path = ROOT / path_str
            if not img_path.exists():
                raise FileNotFoundError(img_path)
            media_name = f"image{entry_idx - 1}.png"
            zf.write(img_path, f"ppt/media/{media_name}")

            img_w, img_h = png_size(img_path)
            box_x, box_y, box_w, box_h = 0.55, 0.95, 8.0, 5.95
            draw_w, draw_h = fit_box(img_w, img_h, box_w, box_h)
            img_x = box_x + (box_w - draw_w) / 2
            img_y = box_y + (box_h - draw_h) / 2

            title = path_str
            title_lines = wrap(title, width=82)
            desc_lines = wrap(desc, width=62)
            desc_text = "\n".join(desc_lines[:24])
            if len(desc_lines) > 24:
                desc_text += "\n..."

            content = (
                textbox_xml(2, "Path", 0.45, 0.15, 12.45, 0.55, [paragraph_xml("\n".join(title_lines), 1050, True)])
                + picture_xml(3, "rIdImg", Path(path_str).name, img_x, img_y, draw_w, draw_h)
                + textbox_xml(4, "Description", 8.75, 0.95, 4.1, 5.95, [paragraph_xml(desc_text, 820)])
                + textbox_xml(5, "Slide number", 11.75, 7.05, 1.1, 0.25, [paragraph_xml(f"{entry_idx - 1}/{len(entries)}", 750)])
            )
            zf.writestr(f"ppt/slides/slide{entry_idx}.xml", slide_xml(content))
            zf.writestr(
                f"ppt/slides/_rels/slide{entry_idx}.xml.rels",
                rels_xml(
                    [
                        (
                            "rIdImg",
                            "http://schemas.openxmlformats.org/officeDocument/2006/relationships/image",
                            f"../media/{media_name}",
                        ),
                        (
                            "rIdLayout",
                            "http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout",
                            "../slideLayouts/slideLayout1.xml",
                        ),
                    ]
                ),
            )

    print(OUT)


if __name__ == "__main__":
    build()
