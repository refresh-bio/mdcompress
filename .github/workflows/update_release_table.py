#!/usr/bin/env python3
"""Append a formatted "Release Assets" table to a GitHub release body.

Reads the live asset list for a given tag from the GitHub API, builds a
markdown table of per-platform download links, and PATCHes it into the
release description. Idempotent: if the table is already present it does
nothing, so it is safe to run repeatedly (and safe against the
`release: edited` event it triggers when it writes).

Usage:
    python update_release_table.py <tag>

Environment (all set automatically inside GitHub Actions):
    GITHUB_TOKEN       required; a token with `contents: write` on the repo
    GITHUB_REPOSITORY  "owner/repo" -- selects the repo, so this works
                       unchanged on a fork and on upstream
    GITHUB_API_URL     API base URL (defaults to the public GitHub API)

Uses only the Python standard library, so it needs no `pip install`.
"""

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request

TABLE_HEADER = "### Release Assets"

# Asset filenames produced by build_release.py look like:
#   mdcompress-1.0.0.linux.arm64.tar.gz
ASSET_RE = re.compile(r"\.(linux|mac|windows)\.(x64|arm64)\.tar\.gz$")

OS_NAMES = {"linux": "Linux", "mac": "macOS", "windows": "Windows"}
OS_ORDER = ["linux", "mac", "windows"]      # stable display order
ARCH_ORDER = ["x64", "arm64"]


def _env(name):
    value = os.getenv(name)
    if not value:
        sys.exit(f"ERROR: required environment variable {name} is not set")
    return value


def _api_base():
    return os.getenv("GITHUB_API_URL", "https://api.github.com")


def _request(method, url, data=None):
    body = json.dumps(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body, method=method)
    req.add_header("Authorization", f"token {_env('GITHUB_TOKEN')}")
    req.add_header("Accept", "application/vnd.github+json")
    if body is not None:
        req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")
        sys.exit(f"ERROR: {method} {url} -> {exc.code} {exc.reason}\n{detail}")


def get_release(tag):
    repo = _env("GITHUB_REPOSITORY")
    return _request("GET", f"{_api_base()}/repos/{repo}/releases/tags/{tag}")


def update_release_body(release_id, new_body):
    repo = _env("GITHUB_REPOSITORY")
    _request("PATCH", f"{_api_base()}/repos/{repo}/releases/{release_id}",
             {"body": new_body})


def arch_label(os_key, arch):
    if arch == "x64":
        return "Intel/AMD 64-bit"
    if arch == "arm64":
        return "Apple Silicon (ARM 64-bit)" if os_key == "mac" else "ARM 64-bit"
    return arch


def generate_markdown_table(release):
    rows = []
    for asset in release.get("assets", []):
        name = asset["name"]
        if name.endswith(".md5"):
            continue
        match = ASSET_RE.search(name)
        if not match:
            continue
        os_key, arch = match.group(1), match.group(2)
        rows.append((os_key, arch, asset["browser_download_url"]))

    if not rows:
        return ""

    rows.sort(key=lambda r: (
        OS_ORDER.index(r[0]) if r[0] in OS_ORDER else len(OS_ORDER),
        ARCH_ORDER.index(r[1]) if r[1] in ARCH_ORDER else len(ARCH_ORDER),
    ))

    lines = [
        TABLE_HEADER,
        "| OS | Architecture | Link |",
        "|----|--------------|------|",
    ]
    has_mac = False
    for os_key, arch, url in rows:
        has_mac = has_mac or os_key == "mac"
        lines.append(
            f"| {OS_NAMES.get(os_key, os_key)} | "
            f"{arch_label(os_key, arch)} | [Download]({url}) |"
        )

    table = "\n".join(lines) + "\n"
    if has_mac:
        table += (
            "\n> **macOS note:** the binaries are unsigned. If macOS blocks one, "
            "clear the quarantine attribute:\n"
            "> ```\n"
            "> sudo xattr -dr com.apple.quarantine <path>/mdcompress\n"
            "> ```\n"
        )
    return table


def main(tag):
    release = get_release(tag)
    body = release.get("body") or ""

    if TABLE_HEADER in body:
        print(f"Asset table already present in '{tag}'; nothing to do.")
        return

    table = generate_markdown_table(release)
    if not table:
        print(f"No matching assets found for '{tag}'; nothing to do.")
        return

    print(table)
    new_body = (body + "\n\n" + table) if body.strip() else table
    update_release_body(release["id"], new_body)
    print(f"Asset table added to release '{tag}'.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add a release-assets table to a GitHub release body.")
    parser.add_argument("tag", help="Release tag, e.g. v1.0.0")
    main(parser.parse_args().tag)
