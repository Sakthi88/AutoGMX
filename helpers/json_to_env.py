#!/usr/bin/env python3
import json
import shlex
import sys


def shell_quote(value):
    return shlex.quote(str(value))


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: json_to_env.py config.json")
    path = sys.argv[1]
    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    for key, value in data.items():
        if isinstance(value, list):
            value = " ".join(str(item) for item in value)
        elif isinstance(value, bool):
            value = "yes" if value else "no"
        print(f"export {key}={shell_quote(value)}")


if __name__ == "__main__":
    main()
