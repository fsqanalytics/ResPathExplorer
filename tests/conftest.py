"""Keep unit tests independent of the network, working directory, and display."""

import socket

import matplotlib
import pytest

matplotlib.use("Agg")


@pytest.fixture(autouse=True)
def isolated_test_environment(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def reject_network(*args, **kwargs):
        raise AssertionError("Unit tests must mock external network requests.")

    monkeypatch.setattr(socket.socket, "connect", reject_network)
    yield

    import matplotlib.pyplot as plt
    plt.close("all")
