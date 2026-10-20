import re
from Bio.KEGG import REST
from typing import Optional

def search_gene_id_kegg(gene_name: str, org_code: Optional[str] = None) -> Optional[str]:
    """Convert a gene name to a KEGG gene ID."""
    org = org_code
    result = REST.kegg_find("genes", gene_name).read()
    match = re.findall(rf'\b{org}:\w+\b', result)
    return match[0] if match else None


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