from bioservices import KEGG
from typing import Dict, Optional
from .validate_color_code import validate_color_code


def get_url_pathway(
    target_path: str,
    gene_color_dict: Optional[Dict[str, str]] = None
) -> str:
    """
    Generate a KEGG pathway visualization URL with specific genes highlighted using custom colors.

    Args:
        target_path (str): KEGG pathway identifier (e.g., "hsa04110").
        gene_color_dict (Optional[Dict[str, str]]): Dictionary where keys are gene identifiers and values are
                                                    color strings in the format "background,border".

    Returns:
        str: URL to the KEGG pathway visualization with colored genes.

    Raises:
        ValueError: If the target_path is invalid or if gene_color_dict contains invalid entries.
    """
    if not target_path or not isinstance(target_path, str):
        raise ValueError("target_path must be a non-empty string.")

    if gene_color_dict is not None:
        if not isinstance(gene_color_dict, dict):
            raise ValueError("gene_color_dict must be a dictionary if provided.")

        for gene_id, color in gene_color_dict.items():
            if not isinstance(gene_id, str) or not isinstance(color, str):
                raise ValueError(f"Each gene ID and color must be a string. Found: {gene_id} -> {color}")
            if "," not in color:
                raise ValueError(f"Color format must be 'background,border'. Found: {color}")
            bg, border = color.split(",", 1)
            validate_color_code(bg.strip())
            validate_color_code(border.strip())

    s = KEGG()
    url = s.show_pathway(target_path, scale=None, keggid=gene_color_dict, show=True)
    return url
