#!/usr/bin/env python3
"""
Traverse all XML files in data/europepmc_articles, extract the <body> element
as plain text, and append/update that text in the corresponding metadata JSON
file, keyed by the PMCID found via an element with attribute pub-id-type="pmcid".

Usage:
  python -m src.update_body_from_xml \
      --articles-dir /Users/mahfouz/Code/biohackathon/biohackathon2025/data/europepmc_articles

Defaults to the repository's data/europepmc_articles directory if run from the
repo root or any subdirectory inside it.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import traceback
import xml.etree.ElementTree as ET
from typing import Optional


def find_repo_root(start: str) -> str:
    cur = os.path.abspath(start)
    while True:
        if os.path.isdir(os.path.join(cur, '.git')) or os.path.exists(os.path.join(cur, 'README.md')):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            return start
        cur = parent


def extract_pmcid(root: ET.Element) -> Optional[str]:
    """Find the PMCID text from any element with attribute pub-id-type="pmcid".

    Handles namespaced and non-namespaced tags by checking attributes only.
    """
    for el in root.iter():
        if el.attrib.get('pub-id-type') == 'pmcid' or el.attrib.get('{*}pub-id-type') == 'pmcid':
            text = (el.text or '').strip()
            if text:
                return text
    # Some JATS use <article-id pub-id-type="pmcid"> or <pub-id> or <article-id>.
    # Try common tag names heuristically as a fallback.
    for tag in ('article-id', 'pub-id', 'article-id', 'article-id'):
        for el in root.iter():
            local = el.tag.split('}', 1)[-1]
            if local == tag and (el.attrib.get('pub-id-type') == 'pmcid' or el.attrib.get('{*}pub-id-type') == 'pmcid'):
                text = (el.text or '').strip()
                if text:
                    return text
    return None


def extract_body_text(root: ET.Element) -> str:
    """Extract plain text from the <body> element (any namespace),
    concatenating all text nodes and normalizing whitespace.
    """
    body_el: Optional[ET.Element] = None
    for el in root.iter():
        local = el.tag.split('}', 1)[-1]
        if local == 'body':
            body_el = el
            break
    if body_el is None:
        return ''
    # itertext yields all text content in document order
    text_chunks = [t for t in body_el.itertext() if t is not None]
    text = '\n'.join(line for line in '\n'.join(text_chunks).splitlines())
    # Normalize whitespace: collapse multiple spaces, trim lines, collapse excessive blank lines
    text = re.sub(r'\s+', ' ', text)
    text = text.strip()
    return text


def load_json(path: str) -> Optional[dict]:
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except json.JSONDecodeError:
        # Attempt a light repair: convert single-quoted body text values to double quotes
        try:
            with open(path, 'r', encoding='utf-8') as f:
                raw = f.read()
            # Replace patterns like: "body": {  "text" : '...'} with double quotes
            def _fix_body_text_quotes(match: re.Match) -> str:
                prefix, inner, suffix = match.group(1), match.group(2), match.group(3)
                # Escape double quotes inside the inner text, keep single quotes as-is
                inner_escaped = inner.replace('"', '\\"')
                return prefix + '"' + inner_escaped + '"' + suffix

            repaired = re.sub(
                r'("body"\s*:\s*\{\s*"text"\s*:\s*)\'(.*?)\'(\s*\})',
                _fix_body_text_quotes,
                raw,
                flags=re.DOTALL,
            )
            return json.loads(repaired)
        except Exception:
            return None
    except Exception:
        return None


def save_json(path: str, data: dict) -> None:
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
        f.write('\n')
    os.replace(tmp, path)


def process_xml_file(xml_path: str, articles_dir: str, dry_run: bool = False) -> bool:
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except Exception:
        print(f"[WARN] Failed to parse XML: {xml_path}", file=sys.stderr)
        traceback.print_exc()
        return False

    pmcid = extract_pmcid(root)
    if not pmcid:
        print(f"[WARN] No PMCID found in XML: {xml_path}", file=sys.stderr)
        return False

    body_text = extract_body_text(root)

    meta_filename = f"{pmcid}_metadata.json"
    meta_path = os.path.join(articles_dir, f"PMC{meta_filename}")
    if not os.path.exists(meta_path):
        print(f"[WARN] Missing metadata JSON for {pmcid}: {meta_path}", file=sys.stderr)
        return False

    meta = load_json(meta_path)
    if meta is None:
        print(f"[WARN] Could not load JSON (possibly invalid) for {pmcid}: {meta_path}", file=sys.stderr)
        return False

    # Update/append the body text
    if 'body' not in meta or not isinstance(meta['body'], dict):
        meta['body'] = {}
    meta['body']['text'] = body_text

    if dry_run:
        print(f"[DRY-RUN] Would update {meta_path} with body.text length={len(body_text)}")
        return True

    try:
        save_json(meta_path, meta)
        print(f"[OK] Updated {meta_path} (pmcid={pmcid}, text_len={len(body_text)})")
        return True
    except Exception:
        print(f"[ERROR] Failed to write JSON: {meta_path}", file=sys.stderr)
        traceback.print_exc()
        return False


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Append/update body.text in metadata JSONs from XML bodies.")
    parser.add_argument('--articles-dir', default=None, help='Path to europepmc_articles directory')
    parser.add_argument('--dry-run', action='store_true', help='Do not write changes, only report actions')
    args = parser.parse_args(argv)

    repo_root = find_repo_root(os.getcwd())
    articles_dir = args.articles_dir or os.path.join(repo_root, 'data', 'europepmc_articles')

    if not os.path.isdir(articles_dir):
        print(f"[ERROR] Articles directory not found: {articles_dir}", file=sys.stderr)
        return 2

    successes = 0
    failures = 0

    for name in os.listdir(articles_dir):
        if not name.lower().endswith('.xml'):
            continue
        xml_path = os.path.join(articles_dir, name)
        ok = process_xml_file(xml_path, articles_dir, dry_run=args.dry_run)
        if ok:
            successes += 1
        else:
            failures += 1

    print(f"Done. Success: {successes}, Failures: {failures}")
    return 0 if failures == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
