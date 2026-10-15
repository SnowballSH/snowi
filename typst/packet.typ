// Packet layout. `packet(manifest, variant, language)` renders one PDF body.
#import "notes.typ": variant as variant-state, lang as lang-state

#let mono = "Inconsolata"
#let serif = "New Computer Modern"
#let text-size = 7.5pt
#let muted = luma(90)
#let keep-together = 45 // lines
#let margin = (x: 0.35in, top: 0.42in, bottom: 0.28in)
#let gutter = 0.16in

/// Plain text with `code` spans and $typst math$ spans.
#let rich(source) = {
  let pattern = regex("`[^`]+`|\$[^$]+\$")
  let pieces = ()
  let pos = 0
  for m in source.matches(pattern) {
    pieces.push(source.slice(pos, m.start))
    if m.text.starts-with("`") {
      pieces.push(raw(m.text.slice(1, -1)))
    } else {
      pieces.push(eval(m.text, mode: "markup"))
    }
    pos = m.end
  }
  pieces.push(source.slice(pos))
  pieces.join()
}

#let raw-lang(language) = if language == "cpp" { "cpp" } else { "py" }

#let label-run(name, body) = [#text(weight: "bold")[#name:] #body]

#let usage-block(source, language) = block(
  width: 100%,
  inset: (left: 3pt, y: 1pt),
  stroke: (left: 0.6pt + luma(150)),
  spacing: 0.35em,
  raw(source, lang: raw-lang(language), block: true),
)

// Records where an entry starts (for tools/build.py --report). The
// placed context takes no space; it sits at the current flow position.
#let mark-position(entry) = place(context {
  let pos = here().position()
  [#metadata((id: entry.id, page: pos.page, x: pos.x.pt(), y: pos.y.pt())) <snowi-pos>]
})

// Label of an entry's heading, unique across the four bundled packets.
#let entry-label(entry, variant, language) = label("snowi:" + variant + ":" + language + ":" + entry.id)

#let item(entry, language, variant, by-id) = {
  mark-position(entry)
  let snippet = entry.code.at(language)
  let fields = snippet.fields
  let is-full = variant == "full"
  let ext = snippet.file.split(".").last()
  block(sticky: true, spacing: 0.4em, below: 0.3em)[#heading(level: 2, entry.title)#entry-label(entry, variant, language)]
  let info = ()
  info.push(rich(fields.Description))
  if "Precondition" in fields { info.push(label-run("Requires", rich(fields.Precondition))) }
  if "Time" in fields { info.push(label-run("Time", rich(fields.Time))) }
  if snippet.uses.len() > 0 {
    let names = snippet.uses.map(u => by-id.at(u, default: u))
    info.push(label-run("Uses", names.join(", ") + "."))
  }
  if "Tweak" in fields { info.push(label-run("Tweak", rich(fields.Tweak))) }
  if is-full and "Notes" in fields { info.push(label-run("Notes", rich(fields.Notes))) }
  block(sticky: true, spacing: 0.35em, par(justify: true, info.join(" ")))
  block(sticky: true, usage-block(fields.Usage, language))
  if is-full and "Example" in fields {
    block(sticky: true, usage-block(fields.Example, language))
  }
  block(spacing: 0.35em, width: 100%)[
    #set text(size: 0.92em, fill: muted)
    #snippet.file #h(1fr) #raw(snippet.hash)
    #if is-full and "Source" in fields [ · #rich(fields.Source)]
  ]
  // Short snippets stay in one column; long ones may break.
  block(breakable: snippet.lines > keep-together, raw(snippet.code, lang: raw-lang(language), block: true))
}

#let note(entry, language, variant) = {
  mark-position(entry)
  block(sticky: true, spacing: 0.4em, below: 0.3em)[#heading(level: 2, entry.title)#entry-label(entry, variant, language)]
  include entry.path
}

/// One run-in paragraph per section: "Section 3: Item 3, Item 4, ...".
#let compact-toc(sections, variant, language) = context {
  set par(justify: false, leading: 0.4em, spacing: 0.45em)
  let page-of(entry) = counter(page).at(entry-label(entry, variant, language)).first()
  for (title, entries) in sections {
    let items = entries.map(e => link(entry-label(e, variant, language))[#e.title~#text(fill: muted, str(page-of(e)))])
    par[*#title*~#text(weight: "bold", str(page-of(entries.first())))#h(0.4em)#items.join([, ])]
  }
}

#let selected(entry, variant, language) = {
  let in-variant = variant == "full" or entry.icpc
  let has-lang = entry.kind == "note" or language in entry.code
  in-variant and has-lang
}

#let packet(manifest, variant, language) = {
  let config = manifest.config
  let layout-config = config.layout
  let lang-name = manifest.langs.at(language)
  let title = config.at(variant).title + " (" + lang-name + ")"
  let by-id = (:)
  for section in manifest.sections {
    for entry in section.entries { by-id.insert(entry.id, entry.title) }
  }

  set document(title: title)
  set page(
    paper: layout-config.paper,
    flipped: true,
    columns: layout-config.columns,
    margin: margin,
    header-ascent: 35%,
    header: context block(
      width: 100%,
      inset: (bottom: 2.5pt),
      stroke: (bottom: 0.4pt),
      {
        set text(size: 8pt, font: serif)
        grid(
          columns: (1fr, auto, 1fr),
          align(left)[*#config.team.university* · #config.team.team],
          text(fill: muted, title),
          align(right, text(weight: "bold", counter(page).display())),
        )
      },
    ),
  )
  set columns(gutter: gutter)
  set text(font: serif, size: text-size, lang: "en", hyphenate: true)
  set par(justify: true, leading: 0.38em, spacing: 0.5em, first-line-indent: 0pt)
  set list(indent: 0.4em, body-indent: 0.35em, spacing: 0.35em)
  set enum(indent: 0.4em, body-indent: 0.35em, spacing: 0.35em)
  set table(inset: (x: 2pt, y: 1.5pt), stroke: 0.4pt + luma(140))
  set math.equation(numbering: none)
  show math.equation: set text(size: text-size)
  // Tall inline math (roots, fractions) grows its own line instead of
  // touching the line above.
  show math.equation.where(block: false): set text(top-edge: "bounds", bottom-edge: "bounds")
  set raw(theme: "theme.tmTheme")
  show raw: set text(font: mono, weight: 500)
  show raw.where(block: true): set par(leading: 0.28em, justify: false)
  show raw.where(block: true): set block(spacing: 0.35em, width: 100%)

  show heading.where(level: 1): it => block(
    width: 100%,
    above: 0.8em,
    below: 0.45em,
    sticky: true,
    {
      set text(size: 10pt, weight: "bold")
      it.body
      v(-0.55em)
      line(length: 100%, stroke: 0.8pt)
    },
  )
  show heading.where(level: 2): it => block(
    width: 100%,
    above: 0.75em,
    below: 0.3em,
    sticky: true,
    text(size: 8.5pt, weight: "bold", it.body),
  )
  show heading.where(level: 3): set text(size: text-size, weight: "bold")
  show heading.where(level: 3): set block(above: 0.55em, below: 0.3em, sticky: true)

  variant-state.update(variant)
  lang-state.update(language)

  // Size the code font so `line-limit` characters fill one column.
  context {
    let probe = measure(text(font: mono, size: 10pt, weight: 500, "0" * layout-config.line-limit))
    let column = (page.width - 2 * margin.x - (layout-config.columns - 1) * gutter) / layout-config.columns
    let code-size = 10pt * (column / probe.width) * 0.995
    show raw: set text(size: code-size)
    show raw.where(block: false): set text(size: text-size * 0.98)

    block(below: 0.6em)[
      #text(size: 13pt, weight: "bold", title) \
      #text(fill: muted)[Generated from source. Each item shows its file and the
      first 6 hex digits of its typing hash (see Contest setup). Section and
      item page numbers follow.]
    ]
    let sections = manifest
      .sections
      .map(section => (section.title, section.entries.filter(e => selected(e, variant, language))))
      .filter(((title, entries)) => entries.len() > 0)
    compact-toc(sections, variant, language)

    for (section-title, entries) in sections {
      heading(level: 1, section-title)
      for entry in entries {
        if entry.kind == "note" { note(entry, language, variant) } else { item(entry, language, variant, by-id) }
      }
    }
  }

  // Page geometry for tools/build.py --report.
  place(context [#metadata((
    pages: counter(page).final().first(),
    columns: layout-config.columns,
    width: page.width.pt(),
    height: page.height.pt(),
    margin-x: margin.x.pt(),
    top: margin.top.pt(),
    bottom: margin.bottom.pt(),
    gutter: gutter.pt(),
  )) <snowi-pages>])

  if variant == "icpc" {
    let limit = config.icpc.page-limit
    context {
      let pages = counter(page).final().first()
      if pages > limit {
        panic("ICPC packet (" + lang-name + ") has " + str(pages) + " pages; the limit is " + str(limit))
      }
    }
  }
}
