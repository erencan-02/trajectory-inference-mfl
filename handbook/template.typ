// Styling and environments for the handbook.

#let accent = rgb("#1d4f7a")
#let accent-soft = rgb("#e9f0f7")
#let connect-color = rgb("#b5541c")
#let connect-soft = rgb("#fbf0e8")

#let argmin = math.op("arg min", limits: true)
#let argmax = math.op("arg max", limits: true)
#let OT = math.op("OT")
#let Fit = math.op("Fit")

// ---------------------------------------------------------------- environments
#let env-counter = counter("env")

#let env(kind, color, fill, title: none, numbered: true, body) = {
  if numbered { env-counter.step() }
  block(
    width: 100%, inset: (x: 10pt, y: 8pt), radius: 2pt, fill: fill,
    stroke: (left: 2.2pt + color), breakable: true,
    {
      text(fill: color, weight: "bold")[
        #kind#if numbered [ #context env-counter.display()]
      ]
      if title != none [ #text(fill: color)[(#title)]]
      [. ]
      body
    },
  )
}

#let definition(title: none, body) = env("Definition", accent, accent-soft, title: title, body)
#let theorem(title: none, body) = env("Theorem", accent, accent-soft, title: title, body)
#let proposition(title: none, body) = env("Proposition", accent, accent-soft, title: title, body)
#let lemma(title: none, body) = env("Lemma", accent, accent-soft, title: title, body)
#let example(title: none, body) = env("Example", luma(90), luma(245), title: title, body)
#let remark(title: none, body) = env("Remark", luma(90), luma(245), title: title, numbered: false, body)
#let connection(title: none, body) = env("Connection", connect-color, connect-soft, title: title, body)

#let proof(body) = block(width: 100%, inset: (x: 10pt, y: 2pt))[
  _Proof sketch._ #body #h(1fr) $square$
]

// ---------------------------------------------------------------- document setup
#let handbook(title: "", subtitle: "", author: "", date: "", abstract: [], body) = {
  set document(title: title, author: author)
  set text(font: "New Computer Modern", size: 10.5pt, lang: "en")
  set par(justify: true, leading: 0.62em, spacing: 0.95em)
  set page(
    paper: "a4", margin: (x: 2.4cm, top: 2.6cm, bottom: 2.4cm),
    numbering: "1",
    header: context {
      if counter(page).get().first() > 1 {
        set text(size: 8.5pt, fill: luma(110))
        smallcaps(title)
        h(1fr)
        author
        v(-4pt)
        line(length: 100%, stroke: 0.4pt + luma(190))
      }
    },
  )
  set heading(numbering: "1.1")
  show heading.where(level: 1): it => {
    v(1.2em)
    set text(size: 16pt, fill: accent, weight: "bold")
    if it.numbering != none {
      block[#counter(heading).display(it.numbering)#h(0.6em)#it.body]
    } else { block(it.body) }
    v(0.5em)
  }
  show heading.where(level: 2): it => {
    v(0.4em)
    set text(size: 12pt, fill: accent)
    it
    v(0.15em)
  }
  set math.equation(numbering: "(1)", supplement: [Eq.])
  show math.equation: set block(spacing: 0.9em)
  show figure.caption: set text(size: 9pt)
  show figure.caption: set par(justify: true)
  show figure.caption: it => align(left, it)
  set figure(gap: 0.6em)
  show link: set text(fill: accent)
  show ref: set text(fill: accent)
  set list(indent: 0.8em)
  set enum(indent: 0.8em)

  // ---- title page
  v(3.2cm)
  align(center)[
    #text(size: 23pt, weight: "bold", fill: accent)[#title]
    #v(0.5em)
    #text(size: 13pt, style: "italic")[#subtitle]
    #v(2.2em)
    #text(size: 12pt)[#author]
    #v(0.3em)
    #text(size: 10pt, fill: luma(90))[#date]
  ]
  v(2.2cm)
  block(inset: (x: 1.2cm))[
    #set par(justify: true)
    #align(center, text(weight: "bold", fill: accent)[Abstract])
    #v(0.2em)
    #abstract
  ]
  v(1fr)
  pagebreak()
  {
    show outline.entry.where(level: 1): set text(weight: "bold")
    outline(title: text(fill: accent)[Contents], depth: 2, indent: auto)
  }
  pagebreak()
  body
}
