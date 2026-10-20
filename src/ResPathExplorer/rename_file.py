import os

def rename_file(output_directory: str, old_name: str, new_name: str) -> None:
    """
    Rename a file in the given directory.

    Args:
        output_directory (str): Path to the directory where the file exists.
        old_name (str): Current file name.
        new_name (str): New file name.

    Returns:
        None
    """
    old_path = os.path.join(output_directory, old_name)
    new_path = os.path.join(output_directory, new_name)

    if os.path.exists(old_path):
        os.rename(old_path, new_path)
        print(f"✅ File renamed to: {new_name}")
    else:
        raise FileNotFoundError(f"File '{old_name}' not found in {output_directory}.")