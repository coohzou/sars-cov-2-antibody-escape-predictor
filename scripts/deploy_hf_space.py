"""Upload this project to a Hugging Face Docker Space (no HF↔GitHub link required)."""

from __future__ import annotations

import os
from pathlib import Path

from huggingface_hub import HfApi, create_repo

ROOT = Path(__file__).resolve().parents[1]
IGNORE = [
    ".git/*",
    ".venv/*",
    ".idea/*",
    "uploads/*",
    "__pycache__/*",
    "*.py[cod]",
    "*.log",
    "response.json",
    "upload_test.json",
    "u*.json",
    "paper_*.json",
    "assets/*",
    "data/training/drdb_payload/*",
]


def main() -> None:
    token = os.environ.get("HF_TOKEN")
    if not token:
        raise SystemExit("Set HF_TOKEN to a Hugging Face write token.")

    api = HfApi(token=token)
    whoami = api.whoami()
    username = os.environ.get("HF_USERNAME") or whoami.get("name") or whoami.get("fullname")
    if not username:
        raise SystemExit("Could not determine Hugging Face username from HF_TOKEN.")

    space_id = os.environ.get("HF_SPACE", "sars-cov-2-antibody-escape-predictor")
    repo_id = f"{username}/{space_id}"
    print(f"Deploying Space: {repo_id}")
    create_repo(
        repo_id,
        repo_type="space",
        space_sdk="docker",
        private=False,
        exist_ok=True,
    )

    api.upload_folder(
        folder_path=str(ROOT),
        repo_id=repo_id,
        repo_type="space",
        ignore_patterns=IGNORE,
        commit_message="Deploy Flask app from GitHub Actions",
    )

    space_readme = """---
title: SARS-CoV-2 Antibody Escape Predictor
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
app_port: 7860
pinned: false
---

Public demo for https://github.com/coohzou/sars-cov-2-antibody-escape-predictor
"""
    api.upload_file(
        path_or_fileobj=space_readme.encode("utf-8"),
        path_in_repo="README.md",
        repo_id=repo_id,
        repo_type="space",
        commit_message="Add Space README",
    )
    print(f"Deployed: https://huggingface.co/spaces/{repo_id}")


if __name__ == "__main__":
    main()
