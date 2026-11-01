from sp_utils import read_ccd_gain_file


def write_gain_file(tmp_path, contents):
    gain_file = tmp_path / "gain.txt"
    gain_file.write_text(contents, encoding="utf-8")
    return gain_file


def test_gain_zero_returns_without_checking_file(monkeypatch):
    def fail_if_file_checked(*args):
        raise AssertionError("File should not be checked when gain is zero")

    monkeypatch.setattr("sp_utils.os.path.isfile", fail_if_file_checked)

    assert read_ccd_gain_file("missing.txt", 0) == (None, None)


def test_reads_requested_gain_and_explicit_zero_gain(tmp_path):
    gain_file = write_gain_file(
        tmp_path,
        "# CCD gain table\n"
        "GAIN SYSTEM_GAIN ORIGINAL\n"
        "0 7.3440 1\n"
        "20 5.9488 1\n",
    )

    assert read_ccd_gain_file(gain_file, 20) == (5.9488, 7.344)


def test_interpolates_missing_gain_within_table_range(tmp_path):
    gain_file = write_gain_file(
        tmp_path,
        "GAIN SYSTEM_GAIN\n"
        "0 8.0\n"
        "10 4.0\n"
        "20 2.0\n",
    )

    interpolated_gain, system_gain_0 = read_ccd_gain_file(gain_file, 5)

    assert 4.0 < interpolated_gain < 8.0
    assert system_gain_0 == 8.0


def test_requires_explicit_gain_zero_row(tmp_path, capsys):
    gain_file = write_gain_file(tmp_path, "10 4.0\n20 2.0\n")

    assert read_ccd_gain_file(gain_file, 15) == (None, None)
    assert "Gain 0 must be present" in capsys.readouterr().out


def test_does_not_extrapolate_outside_table_range(tmp_path, capsys):
    gain_file = write_gain_file(tmp_path, "0 8.0\n10 4.0\n20 2.0\n")

    assert read_ccd_gain_file(gain_file, 21) == (None, None)
    assert "outside the available range" in capsys.readouterr().out


def test_missing_file_returns_none_pair(tmp_path):
    gain_file = tmp_path / "missing.txt"

    assert read_ccd_gain_file(gain_file, 10) == (None, None)
