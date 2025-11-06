"""
Script to download full text articles from Europe PMC with "blood" in title.
Published in the last 5 years and open access.
"""
import os
from os.path import join

import requests
import json
from pathlib import Path
import time
from datetime import datetime


def search_europepmc(query, page_size=100):
    """
    Search Europe PMC for articles matching the query.

    Args:
        query: Search query string
        page_size: Number of results per page

    Returns:
        Tuple of (list of article metadata, total hit count)
    """
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": query,
        "format": "json",
        "pageSize": page_size,
        "resultType": "core",
    }

    response = requests.get(base_url, params=params)
    response.raise_for_status()

    data = response.json()
    results = data.get("resultList", {}).get("result", [])
    hit_count = data.get("hitCount", 0)

    return results, hit_count


def search_europepmc_all(query, page_size=1000):
    """
    Search Europe PMC and fetch all pages of results.

    Args:
        query: Search query string
        page_size: Number of results per page (max 1000)

    Returns:
        Tuple of (list of all article metadata, total hit count)
    """
    all_results = []
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

    # First request to get total count
    params = {
        "query": query,
        "format": "json",
        "pageSize": page_size,
        "cursorMark": "*",
        "resultType": "core",
    }

    response = requests.get(base_url, params=params)
    response.raise_for_status()
    data = response.json()

    hit_count = data.get("hitCount", 0)
    results = data.get("resultList", {}).get("result", [])
    all_results.extend(results)

    next_cursor = data.get("nextCursorMark")
    cursor_mark = data.get("cursorMark")

    # Fetch remaining pages
    while next_cursor and next_cursor != cursor_mark and len(all_results) < hit_count:
        params["cursorMark"] = next_cursor
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        results = data.get("resultList", {}).get("result", [])
        all_results.extend(results)

        cursor_mark = next_cursor
        next_cursor = data.get("nextCursorMark")

        print(f"  Fetched {len(all_results)}/{hit_count} articles...")
        time.sleep(0.3)  # Be nice to the API

    return all_results, hit_count


def get_full_text(pmcid):
    """
    Download full text for an article given its PMC ID.

    Args:
        pmcid: PubMed Central ID (e.g., "PMC1234567")

    Returns:
        Tuple of (full text content or None, size in bytes)
    """
    # Try to get full text XML
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/" f"{pmcid}/fullTextXML"

    try:
        response = requests.get(url)
        if response.status_code == 200:
            content = response.text
            size = len(content.encode("utf-8"))
            return content, size
        else:
            print(
                f"  Full text not available for {pmcid} "
                f"(Status: {response.status_code})"
            )
            return None, 0
    except Exception as e:
        print(f"  Error fetching full text for {pmcid}: {e}")
        return None, 0


def format_size(size_bytes):
    """Format bytes to human-readable size."""
    for unit in ["B", "KB", "MB", "GB"]:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} TB"


def check_download_done(basepath):
    lockpath = join(basepath, "data/europepmc_articles/download_finished.lock")
    return os.path.exists(lockpath)


def save_download_done(basepath):
    lockpath = join(basepath, "data/europepmc_articles/download_finished.lock")
    with open(lockpath, "w") as file:
        file.write("")


def main():
    # Create output directory
    p = Path(__file__).parent.parent
    output_dir = p / "data/europepmc_articles"
    output_dir.mkdir(exist_ok=True)

    if check_download_done(basepath=p):
        print("Download was previously completed.")
        return

    # Calculate year range (last 5 years)
    current_year = datetime.now().year
    start_year = current_year - 5

    print(f"Searching for open access articles in Europe PMC...")
    print(f"Published between {start_year} and {current_year}")
    print("=" * 60)

    # Query searching only in abstracts for peripheral blood immune cells
    # Quotes needed for multi-word phrases to ensure exact phrase matching
    query = f"""OPEN_ACCESS:Y AND PUB_YEAR:[{start_year} TO {current_year}]
    AND ABSTRACT:(human OR "Homo sapiens")
    AND ABSTRACT:(healthy OR "normal donor")
    AND ABSTRACT:("peripheral blood" OR PBMC OR
         "T cell" OR "B cell" OR "NK cell" OR
         "dendritic cell" OR "immune cell" OR
         monocyte OR neutrophil OR eosinophil OR
         basophil OR lymphocyte OR granulocyte OR leukocyte)
    NOT ABSTRACT:(mouse OR murine OR rat OR "cell line"
         OR cancer OR tumor OR tumour OR carcinoma
         OR "COVID-19" OR "SARS-CoV-2" OR coronavirus
         OR leukemia OR leukaemia OR lymphoma
         OR infection OR infected OR sepsis
         OR autoimmune OR diabetes OR "rheumatoid arthritis"
         OR Crohn OR "inflammatory bowel" OR "multiple sclerosis"
         OR "systemic lupus" OR HIV OR hepatitis
         OR Alzheimer OR Parkinson
         OR transplant OR transplantation
         OR stroke OR myocardial OR cardiovascular)"""

    # First, get count and sample
    _, total_count = search_europepmc(query, page_size=100)

    print(f"\nTotal open access articles found: {total_count}")
    print("Fetching all article metadata...")

    # Fetch all articles
    articles, _ = search_europepmc_all(query, page_size=1000)

    print(f"Retrieved metadata for: {len(articles)} articles")

    # Sample articles to estimate average size
    print("\nEstimating total memory requirements...")
    sample_sizes = []
    sample_count = min(5, len(articles))

    for i, article in enumerate(articles[:sample_count]):
        pmcid = article.get("pmcid")
        if not pmcid:
            continue

        print(f"  Sampling article {i+1}/{sample_count}: {pmcid}...")
        _, size = get_full_text(pmcid)
        if size > 0:
            sample_sizes.append(size)
        time.sleep(0.3)

    if sample_sizes:
        avg_size = sum(sample_sizes) / len(sample_sizes)
        estimated_total = avg_size * total_count
        print(f"\nAverage article size: {format_size(avg_size)}")
        print(
            f"Estimated total memory for all {total_count} articles: "
            f"{format_size(estimated_total)}"
        )
    else:
        print("\nCould not estimate size (no full texts available)")

    # Now download all articles
    print("\n" + "=" * 60)
    print(f"Downloading all {total_count} articles...")
    print("=" * 60)

    downloaded_count = 0
    target_count = total_count  # Download all found articles
    total_downloaded_size = 0

    for article in articles:
        if downloaded_count >= target_count:
            break

        # Check if article has PMC ID (required for full text)
        pmcid = article.get("pmcid")
        if not pmcid:
            continue

        filepath = output_dir / f"{pmcid}.xml"
        if os.path.exists(filepath):
            # print("Article already downloaded. Proceeding to next.")
            downloaded_count += 1
            continue

        title = article.get("title", "Unknown Title")
        pub_year = article.get("pubYear", "Unknown")
        print(
            f"\nProcessing ({downloaded_count + 1}/{target_count}): " f"{title[:60]}..."
        )
        print(f"  PMC ID: {pmcid} | Year: {pub_year}")

        # Get full text
        full_text, size = get_full_text(pmcid)

        if full_text:
            # Save full text XML
            with open(filepath, "w", encoding="utf-8") as f:
                f.write(full_text)

            # Save metadata as JSON
            metadata_filename = output_dir / f"{pmcid}_metadata.json"
            with open(metadata_filename, "w", encoding="utf-8") as f:
                json.dump(article, f, indent=2)

            print(f"  ✓ Saved to {filepath} ({format_size(size)})")
            downloaded_count += 1
            total_downloaded_size += size

            # Be nice to the API
            time.sleep(0.5)

    save_download_done(basepath=p)

    print("\n" + "=" * 60)
    print(f"Download complete!")
    print(f"  Articles downloaded: {downloaded_count}")
    print(f"  Total size: {format_size(total_downloaded_size)}")
    print(f"  Location: '{output_dir}' directory")
    print("=" * 60)


if __name__ == "__main__":
    main()
