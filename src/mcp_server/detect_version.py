#!/usr/bin/env python3
"""
Detect ChromaDB Version

This script attempts to determine what version of ChromaDB created your database.
"""

import sqlite3
import sys
from pathlib import Path

DB_PATH = "/Users/mahfouz/Code/biohackathon/vectors.chromadb/chroma.sqlite3"  # UPDATE THIS!


def check_sqlite_version():
    """Check SQLite version info."""
    print("=" * 60)
    print("ChromaDB Version Detector")
    print("=" * 60)

    # if DB_PATH == "/Users/mahfouz/Code/biohackathon/vectors.chromadb/chroma.sqlite3":
    #     print("❌ ERROR: Please update DB_PATH in the script!")
    #     return

    db_file = Path(DB_PATH)
    if not db_file.exists():
        print(f"❌ ERROR: Database file not found: {DB_PATH}")
        return

    print(f"\n📂 Database: {DB_PATH}")
    print(f"📊 Size: {db_file.stat().st_size / 1024 / 1024:.2f} MB")

    try:
        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()

        # Check user version
        print("\n🔍 Checking database version...")
        cursor.execute("PRAGMA user_version;")
        user_version = cursor.fetchone()[0]
        print(f"   User Version: {user_version}")

        # List all tables
        print("\n📋 Tables in database:")
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;")
        tables = cursor.fetchall()
        for table in tables:
            table_name = table[0]
            cursor.execute(f"SELECT COUNT(*) FROM {table_name};")
            count = cursor.fetchone()[0]
            print(f"   - {table_name}: {count} rows")

        # Check schema of key tables
        print("\n🔍 Checking table schemas...")

        key_tables = ['collections', 'embeddings', 'segments']
        for table_name in key_tables:
            cursor.execute(f"SELECT sql FROM sqlite_master WHERE type='table' AND name='{table_name}';")
            result = cursor.fetchone()
            if result:
                print(f"\n   {table_name} schema:")
                print(f"   {result[0][:200]}...")

        # Try to detect version based on schema
        print("\n🔬 Version Detection:")

        # Check if collections table has topic column (newer versions)
        cursor.execute("PRAGMA table_info(collections);")
        columns = [row[1] for row in cursor.fetchall()]
        print(f"   Collections columns: {columns}")

        if 'topic' in columns:
            print("   ✅ Likely ChromaDB 0.4.x or newer (has 'topic' column)")
        else:
            print("   📌 Likely ChromaDB 0.3.x or older (no 'topic' column)")

        # Check for segments table
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='segments';")
        if cursor.fetchone():
            print("   ✅ Has 'segments' table (0.4.x feature)")
        else:
            print("   📌 No 'segments' table (older version)")

        conn.close()

        print("\n" + "=" * 60)
        print("Recommendations:")
        print("=" * 60)

        if 'topic' in columns:
            print("Your database appears to be from ChromaDB 0.4.x+")
            print("Current ChromaDB should work. The error might be from:")
            print("  1. Corrupted database")
            print("  2. Mid-version (0.4.0-0.4.15 had breaking changes)")
            print("\nTry:")
            print("  pip install --upgrade chromadb")
        else:
            print("Your database appears to be from ChromaDB 0.3.x or older")
            print("\nMigration needed! Try:")
            print("  1. pip install chromadb==0.3.26")
            print("  2. Run the migration script to export data")
            print("  3. pip install --upgrade chromadb")
            print("  4. Import data to new database")

    except Exception as e:
        print(f"\n❌ Error reading database: {e}")
        print("\nThis could mean:")
        print("  - Database is corrupted")
        print("  - Database is from a very old/new version")
        print("  - Database is encrypted")


if __name__ == "__main__":
    check_sqlite_version()