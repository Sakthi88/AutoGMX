#!/usr/bin/env python3
import os
import sys


def steps_from_env(override_key, time_key, timestep_key, scale=1.0):
    override = os.environ.get(override_key, "").strip()
    if override:
        return str(int(override))
    return str(int(round(float(os.environ[time_key]) * scale / float(os.environ[timestep_key]))))


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: render_mdp.py template.mdp output.mdp")

    template_path, output_path = sys.argv[1], sys.argv[2]
    with open(template_path, "r", encoding="utf-8") as handle:
        text = handle.read()

    values = {
        "TEMPERATURE_K": os.environ["TEMPERATURE_K"],
        "PRESSURE_BAR": os.environ["PRESSURE_BAR"],
        "POSRES_FC": os.environ["POSRES_FC"],
        "EM_STEPS": os.environ["EM_STEPS"],
        "NVT_STEPS": steps_from_env("NVT_STEPS", "NVT_TIME_PS", "TIMESTEP_PS"),
        "NPT_STEPS": steps_from_env("NPT_STEPS", "NPT_TIME_PS", "TIMESTEP_PS"),
        "MD_STEPS": steps_from_env("MD_STEPS", "MD_TIME_NS", "TIMESTEP_PS", scale=1000.0),
        "TIMESTEP_PS": os.environ["TIMESTEP_PS"],
        "MD_CONTINUATION": os.environ.get("MD_CONTINUATION") or "yes",
        "MD_GEN_VEL": os.environ.get("MD_GEN_VEL") or "no",
        "MD_GEN_SEED": os.environ.get("MD_GEN_SEED") or "-1",
    }

    rendered = text.format(**values)
    with open(output_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(rendered)


if __name__ == "__main__":
    main()
