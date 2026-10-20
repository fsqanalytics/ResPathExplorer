from unittest import TestCase
from src.ResPathExplorer.validate_color_code import validate_color_code


class TestValidateColorCode(TestCase):
    def test_valid_named_colors(self):
        valid_colors = ['red', 'blue', 'DarkGreen', 'lightblue']
        for color in valid_colors:
            try:
                validate_color_code(color)
            except ValueError:
                self.fail(f"validate_color_code raised ValueError unexpectedly for color: {color}")

    def test_valid_hex_colors(self):
        valid_hex = ['#FF0000', '#00ff00', '#abc', '#ABC']
        for color in valid_hex:
            try:
                validate_color_code(color)
            except ValueError:
                self.fail(f"validate_color_code raised ValueError unexpectedly for color: {color}")

    def test_invalid_colors(self):
        invalid_colors = ['notacolor', '#GGGGGG', '123456', 'bluegreen', '#12345G', '', None]
        for color in invalid_colors:
            with self.assertRaises(ValueError, msg=f"Expected ValueError for invalid color: {color}"):
                validate_color_code(color)

