#!/usr/bin/env python3
"""
ChromaDB MCP Server

An MCP server that exposes ChromaDB vector database functionality.
Allows querying, adding, and managing documents in ChromaDB collections.
"""

import asyncio
import json
from typing import Any, Optional
import chromadb
from chromadb.config import Settings
from mcp.server import Server
from mcp.types import (
    Resource,
    Tool,
    TextContent,
    ImageContent,
    EmbeddedResource,
)
import mcp.server.stdio

# Initialize ChromaDB client
# You can customize this with your own settings
# chroma_client = chromadb.Client(Settings(
#     anonymized_telemetry=False,
# ))

chroma_client = chromadb.PersistentClient(
    path="/Users/mahfouz/Code/biohackathon/vectors.chromadb"
)

# Create MCP server
app = Server("chromadb-mcp")


@app.list_tools()
async def list_tools() -> list[Tool]:
    """List available tools for ChromaDB operations."""
    return [
        Tool(
            name="query_collection",
            description="Search for similar documents in a ChromaDB collection using semantic similarity",
            inputSchema={
                "type": "object",
                "properties": {
                    "collection_name": {
                        "type": "string",
                        "description": "Name of the collection to query"
                    },
                    "query_texts": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Text queries to search for"
                    },
                    "n_results": {
                        "type": "integer",
                        "description": "Number of results to return per query",
                        "default": 5
                    },
                    "where": {
                        "type": "object",
                        "description": "Optional metadata filter",
                        "default": None
                    },
                    "include": {
                        "type": "array",
                        "items": {
                            "type": "string",
                            "enum": ["documents", "metadatas", "distances", "embeddings"]
                        },
                        "description": "What to include in results",
                        "default": ["documents", "metadatas", "distances"]
                    }
                },
                "required": ["collection_name", "query_texts"]
            }
        ),
        Tool(
            name="add_documents",
            description="Add documents to a ChromaDB collection",
            inputSchema={
                "type": "object",
                "properties": {
                    "collection_name": {
                        "type": "string",
                        "description": "Name of the collection to add documents to"
                    },
                    "documents": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Documents to add"
                    },
                    "metadatas": {
                        "type": "array",
                        "items": {"type": "object"},
                        "description": "Optional metadata for each document"
                    },
                    "ids": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Optional IDs for documents (auto-generated if not provided)"
                    }
                },
                "required": ["collection_name", "documents"]
            }
        ),
        Tool(
            name="create_collection",
            description="Create a new ChromaDB collection",
            inputSchema={
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Name of the collection to create"
                    },
                    "metadata": {
                        "type": "object",
                        "description": "Optional metadata for the collection"
                    }
                },
                "required": ["name"]
            }
        ),
        Tool(
            name="delete_collection",
            description="Delete a ChromaDB collection",
            inputSchema={
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Name of the collection to delete"
                    }
                },
                "required": ["name"]
            }
        ),
        Tool(
            name="list_collections",
            description="List all ChromaDB collections",
            inputSchema={
                "type": "object",
                "properties": {}
            }
        ),
        Tool(
            name="get_collection_info",
            description="Get information about a specific collection",
            inputSchema={
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Name of the collection"
                    }
                },
                "required": ["name"]
            }
        ),
        Tool(
            name="delete_documents",
            description="Delete documents from a collection by IDs",
            inputSchema={
                "type": "object",
                "properties": {
                    "collection_name": {
                        "type": "string",
                        "description": "Name of the collection"
                    },
                    "ids": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "IDs of documents to delete"
                    }
                },
                "required": ["collection_name", "ids"]
            }
        ),
        Tool(
            name="peek_collection",
            description="Peek at the first few documents in a collection",
            inputSchema={
                "type": "object",
                "properties": {
                    "collection_name": {
                        "type": "string",
                        "description": "Name of the collection"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Number of documents to return",
                        "default": 10
                    }
                },
                "required": ["collection_name"]
            }
        )
    ]


@app.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    """Handle tool calls for ChromaDB operations."""

    try:
        if name == "query_collection":
            collection = chroma_client.get_collection(arguments["collection_name"])

            results = collection.query(
                query_texts=arguments["query_texts"],
                n_results=arguments.get("n_results", 5),
                where=arguments.get("where"),
                include=arguments.get("include", ["documents", "metadatas", "distances"])
            )

            return [TextContent(
                type="text",
                text=json.dumps(results, indent=2)
            )]

        elif name == "add_documents":
            collection = chroma_client.get_or_create_collection(arguments["collection_name"])

            # Generate IDs if not provided
            ids = arguments.get("ids")
            if not ids:
                import uuid
                ids = [str(uuid.uuid4()) for _ in arguments["documents"]]

            collection.add(
                documents=arguments["documents"],
                metadatas=arguments.get("metadatas"),
                ids=ids
            )

            return [TextContent(
                type="text",
                text=f"Successfully added {len(arguments['documents'])} documents to collection '{arguments['collection_name']}'"
            )]

        elif name == "create_collection":
            collection = chroma_client.create_collection(
                name=arguments["name"],
                metadata=arguments.get("metadata")
            )

            return [TextContent(
                type="text",
                text=f"Successfully created collection '{arguments['name']}'"
            )]

        elif name == "delete_collection":
            chroma_client.delete_collection(arguments["name"])

            return [TextContent(
                type="text",
                text=f"Successfully deleted collection '{arguments['name']}'"
            )]

        elif name == "list_collections":
            collections = chroma_client.list_collections()

            collection_info = [
                {
                    "name": col.name,
                    "metadata": col.metadata,
                    "count": col.count()
                }
                for col in collections
            ]

            return [TextContent(
                type="text",
                text=json.dumps(collection_info, indent=2)
            )]

        elif name == "get_collection_info":
            collection = chroma_client.get_collection(arguments["name"])

            info = {
                "name": collection.name,
                "metadata": collection.metadata,
                "count": collection.count()
            }

            return [TextContent(
                type="text",
                text=json.dumps(info, indent=2)
            )]

        elif name == "delete_documents":
            collection = chroma_client.get_collection(arguments["collection_name"])

            collection.delete(ids=arguments["ids"])

            return [TextContent(
                type="text",
                text=f"Successfully deleted {len(arguments['ids'])} documents from collection '{arguments['collection_name']}'"
            )]

        elif name == "peek_collection":
            collection = chroma_client.get_collection(arguments["collection_name"])

            # Don't include embeddings to avoid serialization issues
            results = collection.peek(limit=arguments.get("limit", 10))

            # Remove embeddings from results if present
            if "embeddings" in results:
                del results["embeddings"]

            # Format nicely
            formatted = {
                "collection": arguments["collection_name"],
                "count": len(results.get("ids", [])),
                "documents": []
            }

            if results.get("ids"):
                for i in range(len(results["ids"])):
                    doc = {"id": results["ids"][i]}

                    if results.get("documents") and i < len(results["documents"]):
                        doc["document"] = results["documents"][i]

                    if results.get("metadatas") and i < len(results["metadatas"]):
                        doc["metadata"] = results["metadatas"][i]

                    formatted["documents"].append(doc)

            return [TextContent(
                type="text",
                text=json.dumps(formatted, indent=2)
            )]

        elif name == "get_documents":
            collection = chroma_client.get_collection(arguments["collection_name"])

            # Build get parameters
            get_params = {}

            if "ids" in arguments and arguments["ids"]:
                get_params["ids"] = arguments["ids"]

            if "where" in arguments and arguments["where"]:
                get_params["where"] = arguments["where"]

            if "limit" in arguments and arguments["limit"]:
                get_params["limit"] = arguments["limit"]

            # Don't include embeddings - they cause serialization issues
            get_params["include"] = ["documents", "metadatas"]

            results = collection.get(**get_params)

            # Format results nicely
            formatted_results = {
                "collection": arguments["collection_name"],
                "count": len(results["ids"]) if results.get("ids") else 0,
                "documents": []
            }

            if results.get("ids"):
                for i in range(len(results["ids"])):
                    doc = {"id": results["ids"][i]}

                    if results.get("documents") and i < len(results["documents"]):
                        doc["document"] = results["documents"][i]

                    if results.get("metadatas") and i < len(results["metadatas"]):
                        doc["metadata"] = results["metadatas"][i]

                    formatted_results["documents"].append(doc)

            return [TextContent(
                type="text",
                text=json.dumps(formatted_results, indent=2)
            )]

        else:
            return [TextContent(
                type="text",
                text=f"Unknown tool: {name}"
            )]

    except Exception as e:
        return [TextContent(
            type="text",
            text=f"Error: {str(e)}"
        )]


@app.list_resources()
async def list_resources() -> list[Resource]:
    """List available ChromaDB resources."""
    resources = []

    # Add a resource for listing all collections
    resources.append(
        Resource(
            uri="chroma://collections",
            name="All Collections",
            mimeType="application/json",
            description="List of all ChromaDB collections"
        )
    )

    # Add a resource for each collection
    try:
        collections = chroma_client.list_collections()
        for col in collections:
            resources.append(
                Resource(
                    uri=f"chroma://collection/{col.name}",
                    name=f"Collection: {col.name}",
                    mimeType="application/json",
                    description=f"Documents in collection '{col.name}' (count: {col.count()})"
                )
            )
    except Exception as e:
        pass

    return resources


@app.read_resource()
async def read_resource(uri: str) -> str:
    """Read a ChromaDB resource."""

    if uri == "chroma://collections":
        collections = chroma_client.list_collections()
        collection_info = [
            {
                "name": col.name,
                "metadata": col.metadata,
                "count": col.count()
            }
            for col in collections
        ]
        return json.dumps(collection_info, indent=2)

    elif uri.startswith("chroma://collection/"):
        collection_name = uri.replace("chroma://collection/", "")
        collection = chroma_client.get_collection(collection_name)

        # Get all documents in the collection
        results = collection.get(include=["documents", "metadatas", "embeddings"])

        return json.dumps({
            "name": collection.name,
            "metadata": collection.metadata,
            "count": collection.count(),
            "documents": results
        }, indent=2)

    else:
        raise ValueError(f"Unknown resource: {uri}")


async def main():
    """Run the MCP server."""
    async with mcp.server.stdio.stdio_server() as (read_stream, write_stream):
        await app.run(
            read_stream,
            write_stream,
            app.create_initialization_options()
        )


if __name__ == "__main__":
    asyncio.run(main())




# #!/usr/bin/env python3
# """
# ChromaDB MCP Server
#
# An MCP server that exposes ChromaDB vector database functionality.
# Allows querying, adding, and managing documents in ChromaDB collections.
# """
#
# import asyncio
# import json
# from typing import Any, Optional
# import chromadb
# from chromadb.config import Settings
# from mcp.server import Server
# from mcp.types import (
#     Resource,
#     Tool,
#     TextContent,
#     ImageContent,
#     EmbeddedResource,
# )
# import mcp.server.stdio
#
# # Initialize ChromaDB client
# # You can customize this with your own settings
# # chroma_client = chromadb.Client(Settings(
# #     anonymized_telemetry=False,
# # ))
#
# chroma_client = chromadb.PersistentClient(
#     path="/Users/mahfouz/Code/biohackathon/vectors.chromadb"
# )
#
# # Create MCP server
# app = Server("chromadb-mcp")
#
#
# @app.list_tools()
# async def list_tools() -> list[Tool]:
#     """List available tools for ChromaDB operations."""
#     return [
#         Tool(
#             name="query_collection",
#             description="Search for similar documents in a ChromaDB collection using semantic similarity",
#             inputSchema={
#                 "type": "object",
#                 "properties": {
#                     "collection_name": {
#                         "type": "string",
#                         "description": "Name of the collection to query"
#                     },
#                     "query_texts": {
#                         "type": "array",
#                         "items": {"type": "string"},
#                         "description": "Text queries to search for"
#                     },
#                     "n_results": {
#                         "type": "integer",
#                         "description": "Number of results to return per query",
#                         "default": 5
#                     },
#                     "where": {
#                         "type": "object",
#                         "description": "Optional metadata filter",
#                         "default": None
#                     },
#                     "include": {
#                         "type": "array",
#                         "items": {
#                             "type": "string",
#                             "enum": ["documents", "metadatas", "distances", "embeddings"]
#                         },
#                         "description": "What to include in results",
#                         "default": ["documents", "metadatas", "distances"]
#                     }
#                 },
#                 "required": ["collection_name", "query_texts"]
#             }
#         ),
#         Tool(
#             name="add_documents",
#             description="Add documents to a ChromaDB collection",
#             inputSchema={
#                 "type": "object",
#                 "properties": {
#                     "collection_name": {
#                         "type": "string",
#                         "description": "Name of the collection to add documents to"
#                     },
#                     "documents": {
#                         "type": "array",
#                         "items": {"type": "string"},
#                         "description": "Documents to add"
#                     },
#                     "metadatas": {
#                         "type": "array",
#                         "items": {"type": "object"},
#                         "description": "Optional metadata for each document"
#                     },
#                     "ids": {
#                         "type": "array",
#                         "items": {"type": "string"},
#                         "description": "Optional IDs for documents (auto-generated if not provided)"
#                     }
#                 },
#                 "required": ["collection_name", "documents"]
#             }
#         ),
#         Tool(
#             name="create_collection",
#             description="Create a new ChromaDB collection",
#             inputSchema={
#                 "type": "object",
#                 "properties": {
#                     "name": {
#                         "type": "string",
#                         "description": "Name of the collection to create"
#                     },
#                     "metadata": {
#                         "type": "object",
#                         "description": "Optional metadata for the collection"
#                     }
#                 },
#                 "required": ["name"]
#             }
#         ),
#         Tool(
#             name="delete_collection",
#             description="Delete a ChromaDB collection",
#             inputSchema={
#                 "type": "object",
#                 "properties": {
#                     "name": {
#                         "type": "string",
#                         "description": "Name of the collection to delete"
#                     }
#                 },
#                 "required": ["name"]
#             }
#         ),
#         Tool(
#             name="list_collections",
#             description="List all ChromaDB collections",
#             inputSchema={
#                 "type": "object",
#                 "properties": {}
#             }
#         ),
#         Tool(
#             name="get_collection_info",
#             description="Get information about a specific collection",
#             inputSchema={
#                 "type": "object",
#                 "properties": {
#                     "name": {
#                         "type": "string",
#                         "description": "Name of the collection"
#                     }
#                 },
#                 "required": ["name"]
#             }
#         ),
#         Tool(
#             name="delete_documents",
#             description="Delete documents from a collection by IDs",
#             inputSchema={
#                 "type": "object",
#                 "properties": {
#                     "collection_name": {
#                         "type": "string",
#                         "description": "Name of the collection"
#                     },
#                     "ids": {
#                         "type": "array",
#                         "items": {"type": "string"},
#                         "description": "IDs of documents to delete"
#                     }
#                 },
#                 "required": ["collection_name", "ids"]
#             }
#         ),
#         Tool(
#             name="peek_collection",
#             description="Peek at the first few documents in a collection",
#             inputSchema={
#                 "type": "object",
#                 "properties": {
#                     "collection_name": {
#                         "type": "string",
#                         "description": "Name of the collection"
#                     },
#                     "limit": {
#                         "type": "integer",
#                         "description": "Number of documents to return",
#                         "default": 10
#                     }
#                 },
#                 "required": ["collection_name"]
#             }
#         )
#     ]
#
#
# @app.call_tool()
# async def call_tool(name: str, arguments: Any) -> list[TextContent]:
#     """Handle tool calls for ChromaDB operations."""
#
#     try:
#         if name == "query_collection":
#             collection = chroma_client.get_collection(arguments["collection_name"])
#
#             results = collection.query(
#                 query_texts=arguments["query_texts"],
#                 n_results=arguments.get("n_results", 5),
#                 where=arguments.get("where"),
#                 include=arguments.get("include", ["documents", "metadatas", "distances"])
#             )
#
#             return [TextContent(
#                 type="text",
#                 text=json.dumps(results, indent=2)
#             )]
#
#         elif name == "add_documents":
#             collection = chroma_client.get_or_create_collection(arguments["collection_name"])
#
#             # Generate IDs if not provided
#             ids = arguments.get("ids")
#             if not ids:
#                 import uuid
#                 ids = [str(uuid.uuid4()) for _ in arguments["documents"]]
#
#             collection.add(
#                 documents=arguments["documents"],
#                 metadatas=arguments.get("metadatas"),
#                 ids=ids
#             )
#
#             return [TextContent(
#                 type="text",
#                 text=f"Successfully added {len(arguments['documents'])} documents to collection '{arguments['collection_name']}'"
#             )]
#
#         elif name == "create_collection":
#             collection = chroma_client.create_collection(
#                 name=arguments["name"],
#                 metadata=arguments.get("metadata")
#             )
#
#             return [TextContent(
#                 type="text",
#                 text=f"Successfully created collection '{arguments['name']}'"
#             )]
#
#         elif name == "delete_collection":
#             chroma_client.delete_collection(arguments["name"])
#
#             return [TextContent(
#                 type="text",
#                 text=f"Successfully deleted collection '{arguments['name']}'"
#             )]
#
#         elif name == "list_collections":
#             collections = chroma_client.list_collections()
#
#             collection_info = [
#                 {
#                     "name": col.name,
#                     "metadata": col.metadata,
#                     "count": col.count()
#                 }
#                 for col in collections
#             ]
#
#             return [TextContent(
#                 type="text",
#                 text=json.dumps(collection_info, indent=2)
#             )]
#
#         elif name == "get_collection_info":
#             collection = chroma_client.get_collection(arguments["name"])
#
#             info = {
#                 "name": collection.name,
#                 "metadata": collection.metadata,
#                 "count": collection.count()
#             }
#
#             return [TextContent(
#                 type="text",
#                 text=json.dumps(info, indent=2)
#             )]
#
#         elif name == "delete_documents":
#             collection = chroma_client.get_collection(arguments["collection_name"])
#
#             collection.delete(ids=arguments["ids"])
#
#             return [TextContent(
#                 type="text",
#                 text=f"Successfully deleted {len(arguments['ids'])} documents from collection '{arguments['collection_name']}'"
#             )]
#
#         elif name == "peek_collection":
#             collection = chroma_client.get_collection(arguments["collection_name"])
#
#             results = collection.peek(limit=arguments.get("limit", 10))
#
#             return [TextContent(
#                 type="text",
#                 text=json.dumps(results, indent=2)
#             )]
#
#         else:
#             return [TextContent(
#                 type="text",
#                 text=f"Unknown tool: {name}"
#             )]
#
#     except Exception as e:
#         return [TextContent(
#             type="text",
#             text=f"Error: {str(e)}"
#         )]
#
#
# @app.list_resources()
# async def list_resources() -> list[Resource]:
#     """List available ChromaDB resources."""
#     resources = []
#
#     # Add a resource for listing all collections
#     resources.append(
#         Resource(
#             uri="chroma://collections",
#             name="All Collections",
#             mimeType="application/json",
#             description="List of all ChromaDB collections"
#         )
#     )
#
#     # Add a resource for each collection
#     try:
#         collections = chroma_client.list_collections()
#         for col in collections:
#             resources.append(
#                 Resource(
#                     uri=f"chroma://collection/{col.name}",
#                     name=f"Collection: {col.name}",
#                     mimeType="application/json",
#                     description=f"Documents in collection '{col.name}' (count: {col.count()})"
#                 )
#             )
#     except Exception as e:
#         pass
#
#     return resources
#
#
# @app.read_resource()
# async def read_resource(uri: str) -> str:
#     """Read a ChromaDB resource."""
#
#     if uri == "chroma://collections":
#         collections = chroma_client.list_collections()
#         collection_info = [
#             {
#                 "name": col.name,
#                 "metadata": col.metadata,
#                 "count": col.count()
#             }
#             for col in collections
#         ]
#         return json.dumps(collection_info, indent=2)
#
#     elif uri.startswith("chroma://collection/"):
#         collection_name = uri.replace("chroma://collection/", "")
#         collection = chroma_client.get_collection(collection_name)
#
#         # Get all documents in the collection
#         results = collection.get(include=["documents", "metadatas", "embeddings"])
#
#         return json.dumps({
#             "name": collection.name,
#             "metadata": collection.metadata,
#             "count": collection.count(),
#             "documents": results
#         }, indent=2)
#
#     else:
#         raise ValueError(f"Unknown resource: {uri}")
#
#
# async def main():
#     """Run the MCP server."""
#     async with mcp.server.stdio.stdio_server() as (read_stream, write_stream):
#         await app.run(
#             read_stream,
#             write_stream,
#             app.create_initialization_options()
#         )
#
#
# if __name__ == "__main__":
#     asyncio.run(main())
