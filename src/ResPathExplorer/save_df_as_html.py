import pandas as pd

def save_df_as_html(df: pd.DataFrame, filename: str, title: str = "Table") -> None:
    """
    Save a pandas DataFrame as a styled HTML file.

    Args:
        df (pd.DataFrame): DataFrame to export.
        filename (str): Output filename (.html).
        title (str): Title for the HTML document and visible heading.

    Returns:
        None
    """
    try:
        html_table = df.to_html(index=False, border=1)

        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title>{title}</title>
            <style>
                body {{
                    font-family: Arial, sans-serif;
                }}
                table {{
                    width: 100%;
                    border-collapse: collapse;
                    margin-top: 20px;
                }}
                th, td {{
                    border: 1px solid black;
                    padding: 8px;
                    text-align: center;
                }}
                th {{
                    background-color: #f2f2f2;
                }}
                h1 {{
                    text-align: center;
                    margin-top: 20px;
                }}
            </style>
        </head>
        <body>
            <h1>{title}</h1>
            {html_table}
        </body>
        </html>
        """

        with open(filename, "w", encoding="utf-8") as f:
            f.write(html_content)

        print(f"HTML file saved: {filename}")
    except Exception as e:
        print(f"Error saving HTML: {e}")


