#!/usr/bin/env python3
"""
Preliminary filtering of downloaded articles.
Combines automated disease/model filtering with optional interactive review.

This script filters out papers with:
- Disease-related terms in titles
- Animal models and non-human studies
- Aging/cohort studies
- Terms ending in "pathy"
"""

import json
import shutil
from pathlib import Path
import sys


def automatic_filter(source_dir, removed_dir, verbose=True):
    """
    Automatically filter articles based on disease terms, animal models, etc.

    Args:
        source_dir: Path to directory containing articles
        removed_dir: Path to directory for removed articles
        verbose: Whether to print progress

    Returns:
        tuple: (kept_count, removed_count)
    """
    # Create removed directory
    removed_dir.mkdir(exist_ok=True)

    # Disease-related terms to filter out (case-insensitive)
    disease_terms = [
        # General terms
        'Patient', 'Patients', 'Disease', 'Disorder', 'Syndrome',
        'Cancer', 'Tumor', 'Tumour', 'Carcinoma', 'Malignancy', 'Neoplasm',
        'Infection', 'Infected', 'Sepsis', 'Inflammation', 'Inflammatory',
        'Autoimmune', 'Immune-mediated',
        'Injury', 'Injured', 'Wound', 'Trauma', 'Burn',
        'Transplant', 'Transplantation',
        'Therapy', 'Therapeutic', 'Treatment',
        'Pathogen', 'Pathology', 'Pathological',
        'Dysfunction', 'Deficiency', 'Defect',
        'Ulcer', 'Lesion', 'Erosion',

        # Study design and demographics to exclude
        'Cohort', 'Aging', 'Ageing', 'Aged', 'Elderly',

        # Animal models and non-human
        'Mouse', 'Mice', 'Murine', 'Rat', 'Rats',
        'Porcine', 'Pig', 'Pigs', 'Swine',
        'Canine', 'Dog', 'Dogs',
        'Bovine', 'Cattle', 'Cow',
        'Ovine', 'Sheep',
        'Equine', 'Horse',
        'Feline', 'Cat', 'Cats',
        'Primate', 'Monkey', 'Macaque',
        'Rabbit', 'Rabbits',
        'Zebrafish',
        'Drosophila',
        'C. elegans', 'Caenorhabditis',

        # Metabolic/Endocrine
        'Diabetes', 'Diabetic', 'Hyperglycemia', 'Hypoglycemia',
        'Obesity', 'Obese', 'Metabolic Syndrome',
        'Hypertension', 'Hyperlipidemia', 'Dyslipidemia',
        'Thyroid', 'Hyperthyroid', 'Hypothyroid',

        # Cardiovascular
        'Stroke', 'Myocardial', 'Infarction', 'Ischemic', 'Ischaemic',
        'Atherosclerosis', 'Thrombosis', 'Embolism',
        'Hypertensive', 'Cardiovascular Disease',
        'Heart Failure', 'Cardiac',

        # Respiratory
        'Asthma', 'Asthmatic', 'COPD', 'Bronchitis',
        'Pneumonia', 'Tuberculosis', 'Fibrosis',
        'Emphysema', 'Respiratory Disease',

        # Neurological
        'Alzheimer', "Alzheimer's", 'Parkinson', "Parkinson's",
        'Dementia', 'Multiple Sclerosis', 'MS',
        'Epilepsy', 'Seizure', 'Neuropathy',
        'Schizophrenia', 'Bipolar', 'Depression', 'Depressive',
        'Autism', 'ADHD', 'Migraine',

        # Rheumatologic/Autoimmune
        'Arthritis', 'Rheumatoid', 'Osteoarthritis',
        'Lupus', 'SLE', 'Systemic Lupus',
        'Scleroderma', 'Sjogren', 'Psoriasis',

        # Hematologic
        'Leukemia', 'Leukaemia', 'Lymphoma', 'Myeloma',
        'Anemia', 'Anaemia', 'Thalassemia',
        'Hemophilia', 'Bleeding Disorder',

        # Infectious
        'COVID', 'SARS', 'Coronavirus', 'COVID-19',
        'HIV', 'AIDS', 'Hepatitis',
        'Influenza', 'Flu', 'Viral',

        # Gastrointestinal
        'Crohn', "Crohn's", 'Colitis', 'IBD',
        'Inflammatory Bowel', 'Cirrhosis',
        'Gastritis', 'Enteritis',

        # Renal
        'Nephropathy', 'Kidney Disease', 'Renal Failure',
        'Dialysis', 'Uremia',

        # Ophthalmic
        'Glaucoma', 'Retinopathy', 'Macular Degeneration',
        'Cataract', 'Uveitis',

        # Allergic/Immunologic
        'Allergy', 'Allergic', 'Atopic', 'Anaphylaxis',
        'Hypersensitivity',

        # Obstetric
        'Pregnancy Loss', 'Miscarriage', 'Abortion',
        'Preeclampsia', 'Eclampsia', 'Gestational Diabetes',

        # Other
        'Necrosis', 'Apoptosis', 'Cachexia',
        'Endometriosis', 'Polycystic',
        'Narcolepsy', 'Sleep Disorder'
    ]

    # Get all JSON metadata files
    json_files = sorted(source_dir.glob("*_metadata.json"))

    removed_count = 0
    kept_count = 0

    if verbose:
        print(f"\nScanning {len(json_files)} articles for disease-related terms...")
        print(f"{'='*80}\n")

    for json_file in json_files:
        try:
            with open(json_file, 'r') as f:
                metadata = json.load(f)

            pmcid = metadata.get('pmcid', 'Unknown')
            title = metadata.get('title', '').lower()

            # Check if any disease term is in the title or if "pathy" is in title
            contains_disease_term = any(term.lower() in title for term in disease_terms)
            contains_pathy = 'pathy' in title

            if contains_disease_term or contains_pathy:
                # Move to removed folder
                xml_file = source_dir / f"{pmcid}.xml"

                if json_file.exists():
                    shutil.move(str(json_file), str(removed_dir / json_file.name))
                if xml_file.exists():
                    shutil.move(str(xml_file), str(removed_dir / xml_file.name))

                removed_count += 1
                if verbose:
                    print(f"REMOVED: {pmcid}")
            else:
                kept_count += 1

        except Exception as e:
            print(f"Error processing {json_file}: {e}")
            continue

    if verbose:
        print(f"\n{'='*80}")
        print(f"Automatic Filtering Complete!")
        print(f"{'='*80}")
        print(f"Total articles processed: {len(json_files)}")
        print(f"Kept (healthy studies): {kept_count}")
        print(f"Removed (disease-related): {removed_count}")
        print(f"\nKept articles remain in: {source_dir}")
        print(f"Removed articles moved to: {removed_dir}")
        print(f"{'='*80}\n")

    return kept_count, removed_count


def interactive_filter(source_dir, keep_dir, removed_dir):
    """
    Interactively review articles one by one.

    Args:
        source_dir: Path to directory containing articles
        keep_dir: Path to directory for kept articles
        removed_dir: Path to directory for removed articles
    """
    # Create directories
    keep_dir.mkdir(exist_ok=True)
    removed_dir.mkdir(exist_ok=True)

    # Get all JSON metadata files
    json_files = sorted(source_dir.glob("*_metadata.json"))

    if not json_files:
        print("No metadata files found!")
        return

    total = len(json_files)
    print(f"\n{'='*80}")
    print(f"Interactive Article Filter")
    print(f"{'='*80}")
    print(f"\nTotal articles to review: {total}")
    print(f"\nFor each article, choose:")
    print(f"  'k' or 'y' = KEEP (move to filtered folder)")
    print(f"  'r' or 'n' = REMOVE (move to removed folder)")
    print(f"  's' = SKIP (leave in original location)")
    print(f"  'q' = QUIT")
    print(f"\n{'='*80}\n")

    kept_count = 0
    removed_count = 0
    skipped_count = 0

    for idx, json_file in enumerate(json_files, 1):
        # Read metadata
        try:
            with open(json_file, 'r') as f:
                metadata = json.load(f)

            pmcid = metadata.get('pmcid', 'Unknown')
            title = metadata.get('title', 'No title available')
            year = metadata.get('pubYear', 'Unknown')
            authors = metadata.get('authorString', 'Unknown authors')

            # Show article info
            print(f"\n[{idx}/{total}] PMC ID: {pmcid}")
            print(f"Year: {year}")
            print(f"Title: {title}")
            print(f"Authors: {authors[:100]}..." if len(authors) > 100 else f"Authors: {authors}")
            print(f"\n{'-'*80}")

            # Get user decision
            while True:
                choice = input(f"Decision (k=keep, r=remove, s=skip, q=quit): ").lower().strip()

                if choice == 'q':
                    print(f"\n\nExiting early...")
                    print(f"Processed: {idx-1}/{total}")
                    print(f"Kept: {kept_count}, Removed: {removed_count}, Skipped: {skipped_count}")
                    return

                if choice in ['k', 'y', 'r', 'n', 's']:
                    break

                print("Invalid choice. Please enter k, r, s, or q.")

            # Get associated XML file
            xml_file = source_dir / f"{pmcid}.xml"

            # Process decision
            if choice in ['k', 'y']:
                # Keep - move to filtered folder
                if json_file.exists():
                    shutil.move(str(json_file), str(keep_dir / json_file.name))
                if xml_file.exists():
                    shutil.move(str(xml_file), str(keep_dir / xml_file.name))
                print(f"✓ KEPT - moved to filtered folder")
                kept_count += 1

            elif choice in ['r', 'n']:
                # Remove - move to removed folder
                if json_file.exists():
                    shutil.move(str(json_file), str(removed_dir / json_file.name))
                if xml_file.exists():
                    shutil.move(str(xml_file), str(removed_dir / xml_file.name))
                print(f"✗ REMOVED - moved to removed folder")
                removed_count += 1

            else:  # 's'
                # Skip - leave in place
                print(f"○ SKIPPED - left in original location")
                skipped_count += 1

        except Exception as e:
            print(f"Error processing {json_file}: {e}")
            continue

    # Final summary
    print(f"\n{'='*80}")
    print(f"Filtering Complete!")
    print(f"{'='*80}")
    print(f"Total reviewed: {total}")
    print(f"Kept: {kept_count} (in {keep_dir})")
    print(f"Removed: {removed_count} (in {removed_dir})")
    print(f"Skipped: {skipped_count} (in {source_dir})")
    print(f"{'='*80}\n")


def main():
    """Main entry point for the filtering script."""
    # Get project root (parent of src directory where this script is located)
    project_root = Path(__file__).parent.parent

    # Default paths relative to project root
    source_dir = project_root / "data" / "europepmc_articles"
    removed_dir = project_root / "data" / "europepmc_articles_disease_removed"
    keep_dir = project_root / "data" / "europepmc_articles_filtered"

    print(f"\n{'='*80}")
    print(f"Preliminary Article Filtering")
    print(f"{'='*80}\n")
    print(f"This script will filter articles based on:")
    print(f"  1. Disease-related terms in titles")
    print(f"  2. Animal models and non-human studies")
    print(f"  3. Aging/cohort studies")
    print(f"  4. Terms ending in 'pathy'")
    print(f"\nSource directory: {source_dir}")
    print(f"\nWhat would you like to do?")
    print(f"  1. Automatic filtering only")
    print(f"  2. Interactive review only (requires TTY)")
    print(f"  3. Both (automatic first, then interactive)")
    print(f"  q. Quit")

    choice = input(f"\nEnter your choice (1/2/3/q): ").strip().lower()

    if choice == 'q':
        print("Exiting...")
        return

    if choice in ['1', '3']:
        # Run automatic filtering
        kept, removed = automatic_filter(source_dir, removed_dir, verbose=True)

        if choice == '1':
            return

    if choice in ['2', '3']:
        # Run interactive filtering
        if not sys.stdin.isatty():
            print("\nError: Interactive mode requires a TTY (terminal).")
            print("This won't work in non-interactive SSH sessions.")
            print("Use option 1 for automatic filtering only.")
            return

        interactive_filter(source_dir, keep_dir, removed_dir)

    if choice not in ['1', '2', '3']:
        print("Invalid choice. Exiting...")


if __name__ == "__main__":
    main()
