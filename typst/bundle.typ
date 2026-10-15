// All four packets in one bundle compile (typst --features bundle).
#import "packet.typ": packet

#let manifest = json(sys.inputs.at("manifest", default: "/build/manifest.json"))

#for variant in ("icpc", "full") {
  for language in ("cpp", "py") {
    document("snowi-" + variant + "-" + language + ".pdf", packet(manifest, variant, language))
  }
}
