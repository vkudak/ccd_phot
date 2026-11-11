import shutil
import subprocess
import sys
from pathlib import Path

from sp_utils import get_star_horizontal, read_config_sat


def test_star_horizontal_accepts_hour_and_degree_ra():
    common_args = ("+20:00:00", 40, -70, 0, "2024-01-10T18:07:00.0000000")
    hour_ra = get_star_horizontal("10:00:00", *common_args)
    degree_ra = get_star_horizontal(150, *common_args, ra_unit="degrees")

    assert hour_ra == degree_ra
    assert 0 <= hour_ra[0] < 360
    assert -90 <= hour_ra[1] <= 90


def test_resolves_orbit_file_relative_to_config_directory():
    config_path = Path(__file__).parent / "fits_19574" / "config_sat.ini"

    config = read_config_sat(config_path)

    assert config["orbit_file"] == str(config_path.parent / "tle_20231122.txt")
    assert Path(config["orbit_file"]).is_file()


def test_sat_phot_sample_run_removes_generated_result(tmp_path):
    project_dir = Path(__file__).resolve().parents[1]
    sample_dir = Path(__file__).parent / "fits_19574"
    run_dir = tmp_path / "fits_19574"
    shutil.copytree(sample_dir, run_dir)

    try:
        completed = subprocess.run(
            [sys.executable, str(project_dir / "sat_phot.py"), str(run_dir)],
            cwd=project_dir,
            capture_output=True,
            text=True,
            check=False,
        )

        assert completed.returncode == 0, completed.stdout + completed.stderr
        assert "Photometry is DONE." in completed.stdout

        result_files = list(run_dir.glob("result_6275_*.phR"))
        assert len(result_files) == 1
        assert "ops7323_00030.fits" in result_files[0].read_text(encoding="utf-8")
    finally:
        for result_file in run_dir.glob("result_*.phR"):
            result_file.unlink()

    assert not list(run_dir.glob("result_*.phR"))