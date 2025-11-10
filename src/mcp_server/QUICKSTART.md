# Quick Start Guide

## Installation & Setup (5 minutes)

### Step 1: Install Dependencies

```bash
pip install -r requirements.txt
```

### Step 2: Load Sample Data (Optional)

```bash
python load_sample_data.py
```

This creates three collections with sample data:
- **programming_languages**: Information about 8 programming languages
- **ai_concepts**: Explanations of 10 AI/ML concepts  
- **quotes**: 8 famous tech quotes

### Step 3: Configure Claude Desktop

1. Find your Claude Desktop config file:
   - **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
   - **Windows**: `%APPDATA%/Claude/claude_desktop_config.json`

2. Edit the file and add the MCP server configuration:

```json
{
  "mcpServers": {
    "chromadb": {
      "command": "python",
      "args": [
        "/full/path/to/chromadb_mcp_server.py"
      ]
    }
  }
}
```

**Important**: Replace `/full/path/to/` with the actual absolute path to where you saved the file.

3. **Restart Claude Desktop completely** (Quit and reopen)

### Step 4: Test It Out

Open Claude Desktop and try these queries:

```
What collections are available in ChromaDB?
```

```
Search the programming_languages collection for languages that are compiled
```

```
Add a document to the quotes collection: "Talk is cheap. Show me the code." by Linus Torvalds
```

```
Show me the first 5 documents in the ai_concepts collection
```

## Common Issues

### "chromadb not found"
- Make sure you installed dependencies: `pip install chromadb mcp`

### Server doesn't show up in Claude
- Verify the path in `claude_desktop_config.json` is correct and absolute
- Restart Claude Desktop **completely** (quit and reopen, don't just close the window)
- Check Claude Desktop logs for errors

### Permission denied
- Make sure the script is executable: `chmod +x chromadb_mcp_server.py`
- On Windows, you may need to use `python.exe` instead of `python` in the config

## Using with Your Own Data

### Persistent Storage

By default, data is stored in memory and lost when the server restarts. To persist data, modify line 16 in `chromadb_mcp_server.py`:

```python
# Change this:
chroma_client = chromadb.Client(Settings(
    anonymized_telemetry=False,
))

# To this:
chroma_client = chromadb.PersistentClient(
    path="./chroma_data"  # Data will be saved in this directory
)
```

### Connecting to an Existing ChromaDB Server

If you're running ChromaDB in client/server mode:

```python
chroma_client = chromadb.HttpClient(
    host="localhost",
    port=8000
)
```

## Next Steps

- Read the full [README.md](README.md) for all available tools and features
- Explore ChromaDB documentation: https://docs.trychroma.com/
- Check out MCP documentation: https://modelcontextprotocol.io/

## Example Workflows

### Building a Knowledge Base

```
Create a collection called "company_docs"

Add these documents to company_docs:
1. "Our Q4 revenue was $2.5M, up 40% YoY"
2. "The new product launch is scheduled for March 15th"
3. "Engineering team has 12 members across 4 time zones"

Search company_docs for information about revenue
```

### Semantic Search

```
Create a collection called "support_tickets"

Add several support ticket descriptions to the collection

Search for tickets similar to "User cannot login to dashboard"
```

Enjoy your ChromaDB MCP server! 🚀
