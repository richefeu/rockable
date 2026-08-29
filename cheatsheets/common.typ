// ============================================================================
//  Shared styling for the Rockable cheatsheets
//  ---------------------------------------------------------------------------
//  Every sheet starts with:
//      #import "common.typ": *
//      #show: rockable-sheet.with(title: [...])
//  and is meant to fit on a *single* landscape A4 page.
// ============================================================================

#import "@preview/boxed-sheet:0.1.0": *

#let rockable-homepage = "https://github.com/richefeu/rockable"
#let rockable-author = "Vincent Richefeu"

// The page setup shared by all the sheets. Individual sheets may override
// `font-size`, `num-columns`, ... when their content needs it.
#let rockable-sheet(
  title: [],
  write-title: false,
  font-size: 6.5pt,
  line-skip: 5.5pt,
  num-columns: 4,
  column-gutter: 9pt,
  body,
) = {
  show: cheatsheet.with(
    title: title,
    homepage: rockable-homepage,
    authors: rockable-author,
    write-title: write-title,
    title-align: left,
    title-number: true,
    title-delta: 2pt,
    scaling-size: true,
    font-size: font-size,
    line-skip: line-skip,
    x-margin: 10pt,
    y-margin: 30pt,
    num-columns: num-columns,
    column-gutter: column-gutter,
    numbered-units: false,
  )
  body
}

// ---------------------------------------------------------------------------
// A keyword, an option value, or a command line, typeset in a grey box.
#let command(body, fill: luma(90%)) = {
  set text(black, font: "Courier New", weight: "semibold")
  box(fill: fill, outset: 2pt, radius: 2pt, [#body])
}

// A file name (input files, output files, companion files).
#let file(name) = {
  set text(rgb("#1c4f7c"), font: "Courier New", weight: "semibold")
  box(fill: rgb("#e3eefa"), outset: 2pt, radius: 2pt, [#name])
}

// ---------------------------------------------------------------------------
// Marker for the features that are only available when Rockable has been
// compiled with the corresponding CMake option (see the "CLI & Files" sheet).
#let flag(name) = {
  box(
    fill: rgb("#ffe9b8"),
    stroke: 0.4pt + rgb("#c8912a"),
    outset: 1.6pt,
    radius: 2pt,
    text(fill: rgb("#7a5606"), size: 0.82em, weight: "bold", font: "Courier New", raw(name)),
  )
}

#let periodic-flag = flag("ROCKABLE_ENABLE_PERIODIC")
#let soft-flag = flag("ROCKABLE_ENABLE_SOFT_PARTICLES")
#let boundary-flag = flag("ROCKABLE_ENABLE_BOUNDARY")
#let ftcorr-flag = flag("ROCKABLE_USE_FT_CORR")

// ---------------------------------------------------------------------------
// A short note that does not deserve a full paragraph.
#let note(body) = text(style: "italic", size: 0.95em)[#body]

// A compact two-column table used for keyword references.
#let kwtable(..rows) = {
  set text(size: 0.97em)
  table(
    columns: (auto, 1fr),
    stroke: none,
    inset: (x: 1.5pt, y: 1.6pt),
    align: (left + top, left + top),
    ..rows
  )
}
