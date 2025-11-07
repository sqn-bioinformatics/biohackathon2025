# ChromaDB Migration Guide
## From 0.3.x to Latest Version

Your old ChromaDB database is incompatible with the MCP server due to version differences. 
This guide will help you migrate your data safely.

## The Problem

- **ChromaDB 0.3.26** requires Pydantic v1
- **MCP library** requires Pydantic v2
- They cannot coexist in the same environment

## The Solution

Use two separate virtual environments:
1. **Export env** (old ChromaDB) - to export your data
2. **MCP env** (new ChromaDB + MCP) - to run the server

---

## Step-by-Step Migration

### Step 1: Export Your Old Data

```bash
# Create a temporary virtual environment for export
python3 -m venv export_env
source export_env/bin/activate

# Install old ChromaDB
pip install chromadb==0.3.26 pydantic==1.10.13

# Edit export_chromadb_old.py and update OLD_DB_PATH
# Change this line:
#   OLD_DB_PATH = "/path/to/your/chroma/data"
# To your actual path, for example:
#   OLD_DB_PATH = "/Users/mahfouz/my_chroma_data"

# Run the export
python export_chromadb_old.py

# This creates: chromadb_export.json
# Deactivate the environment
deactivate
```

### Step 2: Import to New Database

```bash
# Use your main environment (or create a new one)
# Make sure you have latest ChromaDB
pip install --upgrade chromadb

# Edit import_chromadb_new.py and update NEW_DB_PATH
# This should be a NEW directory (not your old one!)
#   NEW_DB_PATH = "/Users/mahfouz/chroma_data_new"

# Run the import
python import_chromadb_new.py

# This creates a new database at NEW_DB_PATH
```

### Step 3: Update Your MCP Server

Edit `chromadb_mcp_server.py` line 26-28:

```python
# Change from:
chroma_client = chromadb.Client(Settings(
    anonymized_telemetry=False,
))

# To (use your NEW_DB_PATH):
chroma_client = chromadb.PersistentClient(
    path="/Users/mahfouz/chroma_data_new"
)
```

### Step 4: Test the MCP Server

```bash
# Make sure you have MCP and latest ChromaDB
pip install mcp chromadb

# Test the server
python chromadb_mcp_server.py
```

---

## Quick Reference Commands

### Check Your Current ChromaDB Version
```bash
python -c "import chromadb; print(chromadb.__version__)"
```

### Find Your Old Database Path
Look for a directory containing `chroma.sqlite3`:
```bash
find ~ -name "chroma.sqlite3" 2>/dev/null
```

### Verify Export File
```bash
ls -lh chromadb_export.json
# Should see a file with reasonable size

# Check content
python -c "import json; data=json.load(open('chromadb_export.json')); print(f'Collections: {len(data[\"collections\"])}')"
```

---

## Troubleshooting

### "Cannot import name 'TypeAdapter' from 'pydantic'"
- You're mixing ChromaDB 0.3.x with MCP (which needs Pydantic v2)
- Solution: Use separate virtual environments as shown above

### "spawn python ENOENT" 
- Use `python3` instead of `python` in Claude Desktop config

### Export script fails
- Make sure you're in the export_env with old ChromaDB
- Verify OLD_DB_PATH is correct
- Check that chroma.sqlite3 exists in that path

### Import script fails
- Make sure NEW_DB_PATH is a NEW directory (don't overwrite old one)
- Ensure chromadb_export.json exists in current directory
- Check you have latest ChromaDB installed

### MCP Server won't start after migration
- Verify NEW_DB_PATH in chromadb_mcp_server.py is correct
- Make sure Claude Desktop config uses `python3`
- Restart Claude Desktop completely after config changes

---

## File Summary

- **export_chromadb_old.py** - Run with ChromaDB 0.3.26 to export data
- **import_chromadb_new.py** - Run with latest ChromaDB to import data
- **chromadb_export.json** - Intermediate JSON file with your data
- **chromadb_mcp_server.py** - The MCP server (use with new database)

---

## Alternative: Start Fresh

If migration seems too complex, you can also:

1. Start with a new empty ChromaDB
2. Re-add your documents through the MCP server
3. Keep your old database as backup

To do this:
```python
# In chromadb_mcp_server.py, use a new path:
chroma_client = chromadb.PersistentClient(
    path="/Users/mahfouz/chroma_data_fresh"
)
```

Then use Claude to add documents to your new database!

---

## Need Help?

1. Check what version created your database: `python detect_chromadb_version.py`
2. Make sure paths are correct (use absolute paths, not relative)
3. Check the export JSON file was created successfully
4. Verify new database was created at NEW_DB_PATH

Good luck with the migration! 🚀
