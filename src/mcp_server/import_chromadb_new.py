#!/usr/bin/env python3
"""
Import to New ChromaDB

Run this with: pip install chromadb (latest version)
"""

import chromadb
import json
from pathlib import Path

# CONFIGURATION - UPDATE THESE!
NEW_DB_PATH = "/path/to/new/chroma/data"  # UPDATE THIS
IMPORT_FILE = "chromadb_export.json"

def import_database():
    """Import data into new ChromaDB."""
    print("="*60)
    print("ChromaDB Import Tool")
    print("="*60)
    
    if NEW_DB_PATH == "/path/to/new/chroma/data":
        print("❌ ERROR: Please update NEW_DB_PATH in the script!")
        return False
    
    # Check if export file exists
    if not Path(IMPORT_FILE).exists():
        print(f"❌ ERROR: Export file not found: {IMPORT_FILE}")
        print("   Run export_chromadb_old.py first!")
        return False
    
    # Load export data
    print(f"\n📖 Loading: {IMPORT_FILE}")
    try:
        with open(IMPORT_FILE, 'r') as f:
            export_data = json.load(f)
        
        print(f"✅ Loaded successfully")
        print(f"   Exported from version: {export_data.get('version', 'unknown')}")
        print(f"   Collections: {len(export_data['collections'])}")
        
    except Exception as e:
        print(f"❌ Error loading file: {e}")
        return False
    
    # Create new database
    print(f"\n🔄 Creating new database at: {NEW_DB_PATH}")
    try:
        client = chromadb.PersistentClient(path=NEW_DB_PATH)
        print("✅ Database created")
        
    except Exception as e:
        print(f"❌ Error creating database: {e}")
        return False
    
    # Import each collection
    for collection_name, collection_data in export_data["collections"].items():
        print(f"\n📥 Importing: {collection_name}")
        
        try:
            # Create collection
            collection = client.get_or_create_collection(
                name=collection_name,
                metadata=collection_data.get("metadata", {})
            )
            
            data = collection_data["data"]
            num_docs = len(data["ids"])
            print(f"   Documents to import: {num_docs}")
            
            if num_docs == 0:
                print("   ⚠️  No documents to import")
                continue
            
            # Import in batches
            batch_size = 1000
            for i in range(0, num_docs, batch_size):
                end_idx = min(i + batch_size, num_docs)
                
                batch_docs = data["documents"][i:end_idx] if data.get("documents") else None
                batch_metas = data["metadatas"][i:end_idx] if data.get("metadatas") else None
                batch_embeds = data["embeddings"][i:end_idx] if data.get("embeddings") else None
                
                collection.add(
                    ids=data["ids"][i:end_idx],
                    documents=batch_docs,
                    metadatas=batch_metas,
                    embeddings=batch_embeds
                )
                
                print(f"   ✅ Imported {end_idx}/{num_docs} documents")
            
            # Verify
            final_count = collection.count()
            print(f"   ✅ Collection complete! Final count: {final_count}")
            
        except Exception as e:
            print(f"   ❌ Error: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    print("\n" + "="*60)
    print("🎉 Import Complete!")
    print("="*60)
    print(f"New database location: {NEW_DB_PATH}")
    print("\nUpdate your chromadb_mcp_server.py with:")
    print(f'chroma_client = chromadb.PersistentClient(path="{NEW_DB_PATH}")')
    print("="*60)
    
    return True


if __name__ == "__main__":
    import sys
    
    # Check ChromaDB version
    import chromadb
    version = chromadb.__version__
    print(f"ChromaDB version: {version}")
    
    success = import_database()
    
    if not success:
        print("\n❌ Import failed!")
        sys.exit(1)
