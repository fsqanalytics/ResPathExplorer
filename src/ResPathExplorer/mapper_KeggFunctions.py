import re
from Bio.KEGG import REST
from typing import Optional

def search_gene_id_kegg(gene_name: str, org_code: Optional[str] = None) -> Optional[str]:
    """Return the first matching KEGG ID, optionally restricted to an organism."""
    if not isinstance(gene_name, str) or not gene_name.strip():
        raise ValueError("gene_name must be a non-empty string.")
    if org_code is not None and (not isinstance(org_code, str) or not org_code.strip()):
        raise ValueError("org_code must be a non-empty string if provided.")

    result = REST.kegg_find("genes", gene_name.strip()).read()
    for line in result.splitlines():
        identifier = line.split("\t", 1)[0].strip()
        if ":" not in identifier:
            continue
        if org_code is None or identifier.split(":", 1)[0] == org_code.strip():
            return identifier
    return None


def get_gene_name_by_kegg_id(kegg_id: str) -> Optional[str]:
    """Retrieve the gene symbol from a KEGG gene ID."""
    if not isinstance(kegg_id, str) or ":" not in kegg_id:
        raise ValueError(f"'{kegg_id}' is not a valid KEGG gene ID (expected format: 'eco:b0002').")

    try:
        result = REST.kegg_get(kegg_id).read()
    except Exception as e:
        raise ValueError(f"Error accessing KEGG API: {e}")

    match = re.search(r"^SYMBOL\s+(.+)", result, re.MULTILINE)
    return match.group(1).strip() if match else None
