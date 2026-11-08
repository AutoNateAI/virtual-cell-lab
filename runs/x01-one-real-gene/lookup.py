#!/usr/bin/env python3
"""Run x01: look one real gene up across public records, and check the records against each other.

Everything comes from public APIs and is saved, unedited, under snapshots/. Nothing here uses a model. Standard library
only, so anyone can re-run it:

    python3 -I runs/x01-one-real-gene/lookup.py

The chain of checks mirrors S01E03 and S01E04: find the gene, read its letters, translate them three at a time, compare
with the protein another database reports, then remove the three letters that delta F508 removes and see what changes.
"""
import json
import re
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
SNAPSHOTS = HERE / "snapshots"
SYMBOL = "CFTR"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "tool=research-intelligence-lab"

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = {a + b + c: AMINO[i] for i, (a, b, c) in enumerate((a, b, c) for a in BASES for b in BASES for c in BASES)}


def fetch(name, url, accept="*/*"):
    """Download a record once per run and keep exactly the bytes the server sent."""
    request = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "research-intelligence-lab/x01"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                body = response.read()
            break
        except OSError:
            if attempt == 3:
                raise
            time.sleep(2 * (attempt + 1))
    (SNAPSHOTS / name).write_bytes(body)
    time.sleep(0.4)  # NCBI asks for at most three requests a second without a key
    return body.decode()


def translate(cds):
    return "".join(CODONS[cds[i : i + 3]] for i in range(0, len(cds) - len(cds) % 3, 3))


def genbank_parts(text):
    """The coding-sequence range, the record's own translation and the nucleotide sequence of a GenBank record."""
    start, end = map(int, re.search(r"^ {5}CDS {13}(\d+)\.\.(\d+)", text, re.M).groups())
    protein = re.sub(r"\s+", "", re.search(r'/translation="([^"]+)"', text).group(1))
    origin = text[text.index("ORIGIN") :]
    letters = re.sub(r"[^a-z]", "", origin.split("\n", 1)[1]).upper()
    return start, end, protein, letters


def uniprot_features(record, wanted):
    return [f for f in record["features"] if f["type"] == wanted]


def main():
    SNAPSHOTS.mkdir(exist_ok=True)
    out = {"symbol": SYMBOL}

    hgnc = json.loads(fetch("hgnc.json", f"https://rest.genenames.org/fetch/symbol/{SYMBOL}", "application/json"))
    doc = hgnc["response"]["docs"][0]
    out["hgnc"] = {k: doc.get(k) for k in ("hgnc_id", "symbol", "name", "location", "entrez_id", "ensembl_gene_id", "uniprot_ids", "mane_select")}
    entrez, ensembl_id, accession = doc["entrez_id"], doc["ensembl_gene_id"], doc["uniprot_ids"][0]
    refseq = next(m for m in doc["mane_select"] if m.startswith("NM_"))

    gene = json.loads(fetch("ncbi_gene.json", f"{EUTILS}/esummary.fcgi?db=gene&id={entrez}&retmode=json&{TOOL}"))["result"][str(entrez)]
    info = gene["genomicinfo"][0]
    out["ncbi_gene"] = {"id": entrez, "name": gene["name"], "description": gene["description"], "chromosome": gene["chromosome"], "map_location": gene["maplocation"],
                        "start": min(info["chrstart"], info["chrstop"]), "stop": max(info["chrstart"], info["chrstop"]), "accession": info["chraccver"], "exon_count": info["exoncount"]}
    out["ncbi_gene"]["span_letters"] = out["ncbi_gene"]["stop"] - out["ncbi_gene"]["start"] + 1

    ensembl = json.loads(fetch("ensembl.json", f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"))
    out["ensembl"] = {k: ensembl.get(k) for k in ("id", "display_name", "seq_region_name", "start", "end", "strand", "biotype", "assembly_name")}

    uniprot = json.loads(fetch("uniprot.json", f"https://rest.uniprot.org/uniprotkb/{accession}.json"))
    protein = uniprot["sequence"]["value"]
    out["uniprot"] = {"accession": accession, "name": uniprot["proteinDescription"]["recommendedName"]["fullName"]["value"], "length": uniprot["sequence"]["length"], "reviewed": uniprot["entryType"]}
    out["uniprot"]["function"] = [t["value"] for c in uniprot["comments"] if c["commentType"] == "FUNCTION" for t in c["texts"]]
    out["uniprot"]["subcellular_locations"] = [loc["location"]["value"] for c in uniprot["comments"] if c["commentType"] == "SUBCELLULAR LOCATION" for loc in c.get("subcellularLocations", [])]
    out["uniprot"]["domains"] = [{"description": f["description"], "start": f["location"]["start"]["value"], "end": f["location"]["end"]["value"]} for f in uniprot_features(uniprot, "Domain")]
    out["uniprot"]["regions_at_508"] = [{"type": f["type"], "description": f["description"], "start": f["location"]["start"]["value"], "end": f["location"]["end"]["value"]}
                                        for f in uniprot["features"] if f["location"]["start"]["value"] <= 508 <= f["location"]["end"]["value"] and f["type"] != "Natural variant"]
    out["uniprot"]["natural_variants_at_508"] = [{"description": f["description"], "original": f.get("alternativeSequence", {}).get("originalSequence"), "alternatives": f.get("alternativeSequence", {}).get("alternativeSequences"),
                                                  "start": f["location"]["start"]["value"], "end": f["location"]["end"]["value"], "id": f.get("featureId")}
                                                 for f in uniprot_features(uniprot, "Natural variant") if f["location"]["start"]["value"] <= 508 <= f["location"]["end"]["value"]]
    out["uniprot"]["natural_variant_count"] = len(uniprot_features(uniprot, "Natural variant"))

    record = fetch("refseq.gb", f"{EUTILS}/efetch.fcgi?db=nuccore&id={refseq}&rettype=gb&retmode=text&{TOOL}")
    start, end, refseq_protein, letters = genbank_parts(record)
    cds = letters[start - 1 : end]
    translated = translate(cds)
    assert translated.endswith("*"), "the coding sequence should end on a stop codon"
    translated = translated[:-1]
    out["refseq"] = {"transcript": refseq, "transcript_letters": len(letters), "cds_start": start, "cds_end": end, "cds_letters": len(cds), "codons": len(cds) // 3}
    out["checks"] = {"our_translation_matches_refseq_translation": translated == refseq_protein, "our_translation_matches_uniprot": translated == protein, "protein_length": len(translated), "position_508": translated[507]}

    # Protein position p sits at coding letters 3p-2 to 3p, counted from 1.
    codons = {str(p): cds[3 * p - 3 : 3 * p] for p in range(505, 512)}
    edited = cds[:1520] + cds[1523:]
    edited_protein = translate(edited).rstrip("*")
    out["delta_f508"] = {"removed_letters_1521_to_1523": cds[1520:1523], "codons_before": codons, "amino_acids_before": {p: CODONS[c] for p, c in codons.items()},
                         "codons_after": {str(p): edited[3 * p - 3 : 3 * p] for p in range(505, 511)}, "edited_protein_length": len(edited_protein),
                         "same_as_normal_minus_position_508": edited_protein == translated[:507] + translated[508:], "cds_letters_after": len(edited)}

    found = json.loads(fetch("clinvar_search.json", f"{EUTILS}/esearch.fcgi?db=clinvar&term=CFTR[gene]+AND+%22c.1521_1523del%22&retmode=json&{TOOL}"))
    ids = found["esearchresult"]["idlist"]
    out["clinvar_search_ids"] = ids
    summaries = json.loads(fetch("clinvar.json", f"{EUTILS}/esummary.fcgi?db=clinvar&id={','.join(ids)}&retmode=json&{TOOL}"))["result"]
    matches = [summaries[i] for i in ids if "c.1521_1523del" in summaries[i]["title"] and "Phe508del" in summaries[i]["title"]]
    out["clinvar"] = [{"id": m["uid"], "title": m["title"], "accession": m["accession"], "germline_classification": m["germline_classification"]["description"], "review_status": m["germline_classification"]["review_status"], "last_evaluated": m["germline_classification"]["last_evaluated"]} for m in matches]

    (HERE / "results.json").write_text(json.dumps(out, indent=2) + "\n")
    (HERE / "results.md").write_text(report(out))
    print(report(out))


def report(o):
    """The human-readable result, generated from results.json so the two cannot disagree."""
    d, c, u, n, e = o["delta_f508"], o["checks"], o["uniprot"], o["ncbi_gene"], o["ensembl"]
    clin = o["clinvar"][0] if o["clinvar"] else None
    at_508 = next((r["description"] for r in u["regions_at_508"] if r["type"] == "Domain"), "none")
    same_gene = o["hgnc"]["symbol"] == n["name"] == e["display_name"] and n["chromosome"] == e["seq_region_name"] == "7"
    rows = [
        ("HGNC, NCBI Gene, Ensembl and UniProt point at one gene on chromosome 7", same_gene),
        ("Our translation of the RefSeq letters equals UniProt's sequence", c["our_translation_matches_uniprot"]),
        ("Position 508 is phenylalanine, inside a named domain", c["position_508"] == "F" and at_508 != "none"),
        ("Removing the three letters changes only that one amino acid", d["same_as_normal_minus_position_508"]),
        ("ClinVar classifies the variant as pathogenic", bool(clin and clin["germline_classification"] == "Pathogenic")),
    ]
    lines = ["# Run x01: one real gene, CFTR", "", "Generated by `lookup.py` from the public records saved in `snapshots/`.", "", "| Prediction, registered before the lookup | Held? |", "|---|---|"]
    lines += [f"| {t} | {'yes' if ok else 'NO'} |" for t, ok in rows]
    lines += ["", "## The records", "", f"- HGNC {o['hgnc']['hgnc_id']}, {o['hgnc']['name']}, at {o['hgnc']['location']}",
              f"- NCBI Gene {n['id']}: chromosome {n['chromosome']}, {n['start']:,} to {n['stop']:,} on {n['accession']}, {n['span_letters']:,} letters, {n['exon_count']} exons",
              f"- Ensembl {e['id']}: chromosome {e['seq_region_name']}, {e['start']:,} to {e['end']:,} on {e['assembly_name']}, {e['end'] - e['start'] + 1:,} letters",
              f"- RefSeq {o['refseq']['transcript']}: {o['refseq']['transcript_letters']:,} letters, of which {o['refseq']['cds_letters']:,} code for protein ({o['refseq']['cds_start']} to {o['refseq']['cds_end']})",
              f"- UniProt {u['accession']}: {u['length']:,} amino acids; {u['natural_variant_count']} natural variants annotated; position 508 is in {at_508}",
              "", "## Delta F508, letter by letter", "", "| Codon | Before | Amino acid | After |", "|---|---|---|---|"]
    for p in map(str, range(505, 511)):
        lines.append(f"| {p} | {d['codons_before'][p]} | {d['amino_acids_before'][p]} | {d['codons_after'][p]} |")
    lines += ["", f"Coding letters 1521 to 1523 are `{d['removed_letters_1521_to_1523']}`. Removing them leaves {d['cds_letters_after']:,} coding letters and a {d['edited_protein_length']:,} amino acid protein.", ""]
    if clin:
        lines += [f"ClinVar {clin['accession']} ({clin['title']}): {clin['germline_classification']}, review status: {clin['review_status']}, last evaluated {clin['last_evaluated']}.", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
