#!/usr/bin/env python3

import subprocess
import re

# Get git diff
print("Running git diff on modules/...")
diff = subprocess.run(
    ["git", "diff", "main", "modules/"], capture_output=True, text=True
).stdout
print("Git diff fetched.\n")

lines = diff.splitlines()

table = []
versions = {}

for i, line in enumerate(lines):
    # Look for conda dependency version change
    if line.startswith("-  - ") and "::" in line and "=" in line:
        match = re.search(r"-  - [^:]+::([^=]+)=([\w\.\-]+)", line)
        if match:
            tool = match.group(1)
            old_version = match.group(2)
            print(f"[{i}] Found old conda version: {tool} {old_version}")
            if tool not in versions.keys():
                versions = {**versions, **{tool: {"old": None, "new": None}}}
            versions[tool]["old"] = old_version
            continue

    if line.startswith("+  - ") and "::" in line and "=" in line:
        match = re.search(r"\+  - [^:]+::([^=]+)=([\w\.\-]+)", line)
        if match:
            tool = match.group(1)
            new_version = match.group(2)
            print(f"[{i}] Found new conda version: {tool} {new_version}")
            if tool not in versions.keys():
                print(f"[{i}] Tool {tool} not found in versions, skipping...")
                continue
            versions[tool]["new"] = new_version
            continue

    if (
        line.startswith("-")
        and (
            "biocontainers/" in line
            or "docker.io/" in line
            or "community.wave.seqera.io" in line
        )
        and "mulled" not in line
    ):
        match = re.search(r"/([^/:]+):([\w\.]+)", line)
        if match:
            tool = match.group(1)
            old_version = match.group(2)
            print(f"[{i}] Found old container version: {tool} {old_version}")
            if tool not in versions.keys():
                versions = {**versions, **{tool: {"old": None, "new": None}}}
            versions[tool]["old"] = old_version
            continue

    if (
        line.startswith("+")
        and (
            "biocontainers/" in line
            or "docker.io/" in line
            or "community.wave.seqera.io" in line
        )
        and "mulled" not in line
    ):
        match = re.search(r"/([^/:]+):([\w\.]+)", line)
        if match:
            tool = match.group(1)
            new_version = match.group(2)
            print(f"[{i}] Found new container version: {tool} {new_version}")
            if tool not in versions.keys():
                print(f"[{i}] Tool {tool} not found in versions, skipping...")
                continue
            versions[tool]["new"] = new_version
            continue

# Print markdown table
print("\n## Tool Updates Detected\n")
print("| Tool | Old Version | New Version |")
print("|------|-------------|-------------|")
for tool, version in versions.items():
    if version["old"] == version["new"]:
        continue
    print(f"| {tool} | {version["old"]} | {version["new"]} |")
