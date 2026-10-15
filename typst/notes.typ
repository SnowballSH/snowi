// Helpers for notes (content/**/*.typ) and the packet template.
// The packet sets these states at the start of each document.

#let variant = state("snowi-variant", "full")
#let lang = state("snowi-lang", "cpp")

/// Shown only in the full packet.
#let full(body) = context if variant.get() == "full" { body }
/// Shown only in the ICPC packet.
#let icpc(body) = context if variant.get() == "icpc" { body }
/// Shown only in the C++ packets.
#let cpp(body) = context if lang.get() == "cpp" { body }
/// Shown only in the PyPy packets.
#let py(body) = context if lang.get() == "py" { body }

/// A named fact: bold run-in name, then the statement.
#let fact(name, body) = block(spacing: 0.5em, breakable: true)[*#name.* #body]

/// A short code block inside a note. lang: "cpp", "py", "sh" or "".
#let code(language, source) = raw(source, lang: language, block: true)

/// Code in the packet's language: #lang-code(cpp: "...", py: "...").
#let lang-code(cpp: none, py: none) = context {
  let source = if lang.get() == "cpp" { cpp } else { py }
  if source != none {
    raw(source, lang: if lang.get() == "cpp" { "cpp" } else { "py" }, block: true)
  }
}
