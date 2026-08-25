#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: launch_job.py PIPELINE_DIR CONFIG_PATH FAILED_MARKER")

    pipeline_dir = Path(sys.argv[1])
    config_path = Path(sys.argv[2])
    failed_marker = Path(sys.argv[3])

    cmd = ["bash", str(pipeline_dir / "run_pipeline.sh"), "--json", str(config_path)]

    result = subprocess.run(cmd, cwd=pipeline_dir, check=False)
    if result.returncode != 0:
        failed_marker.write_text(f"{result.returncode}\n", encoding="utf-8")
    else:
        failed_marker.unlink(missing_ok=True)
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
