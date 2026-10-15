// One packet, chosen with --input variant=icpc|full --input lang=cpp|py.
#import "packet.typ": packet

#packet(
  json("/build/manifest.json"),
  sys.inputs.at("variant", default: "full"),
  sys.inputs.at("lang", default: "cpp"),
)
