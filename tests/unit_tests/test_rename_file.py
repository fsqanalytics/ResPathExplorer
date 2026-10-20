import os
import tempfile
from unittest import TestCase
from src.ResPathExplorer.rename_file import rename_file


class TestRenameFile(TestCase):
    def setUp(self):

        self.test_dir = tempfile.TemporaryDirectory()
        self.dir_path = self.test_dir.name

        self.old_name = "old_file.txt"
        self.new_name = "new_file.txt"

        old_file_path = os.path.join(self.dir_path, self.old_name)
        with open(old_file_path, "w") as f:
            f.write("Test file.")

    def tearDown(self):

        self.test_dir.cleanup()

    def test_successful_rename(self):
        rename_file(self.dir_path, self.old_name, self.new_name)

        self.assertFalse(os.path.exists(os.path.join(self.dir_path, self.old_name)))
        self.assertTrue(os.path.exists(os.path.join(self.dir_path, self.new_name)))

    def test_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            rename_file(self.dir_path, "non-existent_file.txt", "new.txt")
