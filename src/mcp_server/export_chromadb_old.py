#!/usr/bin/env python3
"""
Export ChromaDB 0.3.x Database

Run this with: pip install chromadb==0.3.26 pydantic==1.10.13
"""

import chromadb
import json
from pathlib import Path

# CONFIGURATION - UPDATE THESE!
OLD_DB_PATH = "/Users/mahfouz/Code/biohackathon/vectors.chromadb"  # UPDATE THIS
EXPORT_FILE = "chromadb_export.json"

def export_database():
    """Export data from ChromaDB 0.3.x to JSON."""
    print("="*60)
    print("ChromaDB 0.3.x Export Tool")
    print("="*60)
    
    # if OLD_DB_PATH == "/path/to/your/chroma/data":
    #     print("❌ ERROR: Please update OLD_DB_PATH in the script!")
    #     return False
    
    print(f"\n📂 Connecting to: {OLD_DB_PATH}")
    
    try:
        # Connect to old database
        client = chromadb.PersistentClient(path=OLD_DB_PATH)
        print("✅ Successfully connected!")
        
    except Exception as e:
        print(f"❌ Error connecting: {e}")
        return False
    
    # Get all collections
    collections = client.list_collections()
    print(f"\n📚 Found {len(collections)} collections:")
    for col in collections:
        print(f"   - {col.name}")
    
    export_data = {
        "version": "0.3.26",
        "export_date": str(Path(EXPORT_FILE).stat().st_mtime) if Path(EXPORT_FILE).exists() else "new",
        "collections": {}
    }
    
    # Export each collection
    for collection in collections:
        print(f"\n📖 Exporting: {collection.name}")
        
        try:
            # Get all data
            results = collection.get(
                include=["documents", "metadatas", "embeddings"]
            )
            
            num_docs = len(results["ids"])
            print(f"   Documents: {num_docs}")
            
            export_data["collections"][collection.name] = {
                "name": collection.name,
                "metadata": collection.metadata if hasattr(collection, 'metadata') else {},
                "data": {
                    "ids": results["ids"],
                    "documents": results.get("documents", []),
                    "metadatas": results.get("metadatas", []),
                    "embeddings": results.get("embeddings", []),
                }
            }
            
            print(f"   ✅ Exported successfully")
            
        except Exception as e:
            print(f"   ❌ Error: {e}")
            continue
    
    # Save to JSON
    print(f"\n💾 Saving to {EXPORT_FILE}...")
    try:
        with open(EXPORT_FILE, 'w') as f:
            json.dump(export_data, f, indent=2)
        
        file_size = Path(EXPORT_FILE).stat().st_size / 1024 / 1024
        print(f"✅ Export complete!")
        print(f"   File: {EXPORT_FILE}")
        print(f"   Size: {file_size:.2f} MB")
        print(f"   Collections: {len(export_data['collections'])}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error saving file: {e}")
        return False


if __name__ == "__main__":
    import sys
    
    # Check ChromaDB version
    import chromadb
    version = chromadb.__version__
    print(f"ChromaDB version: {version}")
    
    if not version.startswith("0.3"):
        print("\n⚠️  WARNING: This script is designed for ChromaDB 0.3.x")
        print("   Install with: pip install chromadb==0.3.26 pydantic==1.10.13")
        sys.exit(1)
    
    success = export_database()
    
    if success:
        print("\n" + "="*60)
        print("Next Steps:")
        print("="*60)
        print("1. Deactivate this virtual environment")
        print("2. Install latest ChromaDB: pip install chromadb")
        print("3. Run: python import_chromadb.py")
        print("="*60)
    else:
        print("\n❌ Export failed!")
        sys.exit(1)
