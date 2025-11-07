#!/usr/bin/env python3
"""
Example script to populate ChromaDB with sample data.
Run this before using the MCP server to have some test data available.
"""

import chromadb
from chromadb.config import Settings

# Initialize ChromaDB
client = chromadb.Client(Settings(
    anonymized_telemetry=False,
))

# Create a collection for programming languages
print("Creating 'programming_languages' collection...")
prog_lang_collection = client.get_or_create_collection("programming_languages")

prog_lang_collection.add(
    documents=[
        "Python is a high-level, interpreted programming language known for its simplicity and readability.",
        "JavaScript is the programming language of the web, running in browsers and on servers via Node.js.",
        "Rust is a systems programming language focused on safety, speed, and concurrency.",
        "Go is a statically typed, compiled language designed for simplicity and efficiency.",
        "TypeScript is a typed superset of JavaScript that compiles to plain JavaScript.",
        "Java is a class-based, object-oriented programming language designed to have few implementation dependencies.",
        "C++ is a general-purpose programming language created as an extension of the C programming language.",
        "Ruby is a dynamic, open source programming language with a focus on simplicity and productivity.",
    ],
    metadatas=[
        {"category": "interpreted", "year": 1991, "creator": "Guido van Rossum"},
        {"category": "interpreted", "year": 1995, "creator": "Brendan Eich"},
        {"category": "compiled", "year": 2010, "creator": "Mozilla Research"},
        {"category": "compiled", "year": 2009, "creator": "Google"},
        {"category": "compiled", "year": 2012, "creator": "Microsoft"},
        {"category": "compiled", "year": 1995, "creator": "Sun Microsystems"},
        {"category": "compiled", "year": 1985, "creator": "Bjarne Stroustrup"},
        {"category": "interpreted", "year": 1995, "creator": "Yukihiro Matsumoto"},
    ],
    ids=[f"lang_{i}" for i in range(8)]
)

print(f"Added {prog_lang_collection.count()} documents to 'programming_languages' collection")

# Create a collection for AI/ML concepts
print("\nCreating 'ai_concepts' collection...")
ai_collection = client.get_or_create_collection("ai_concepts")

ai_collection.add(
    documents=[
        "Machine Learning is a subset of AI that enables systems to learn from data without explicit programming.",
        "Deep Learning uses neural networks with multiple layers to learn hierarchical representations of data.",
        "Natural Language Processing enables computers to understand, interpret, and generate human language.",
        "Computer Vision enables machines to derive meaningful information from digital images and videos.",
        "Reinforcement Learning trains agents to make decisions by rewarding desired behaviors.",
        "Supervised Learning uses labeled data to train models to make predictions.",
        "Unsupervised Learning finds patterns in unlabeled data without predefined categories.",
        "Transfer Learning leverages knowledge from one task to improve performance on another related task.",
        "Neural Networks are computing systems inspired by biological neural networks in animal brains.",
        "Transformers are a type of neural network architecture that has revolutionized NLP tasks.",
    ],
    metadatas=[
        {"difficulty": "beginner", "field": "general"},
        {"difficulty": "intermediate", "field": "neural networks"},
        {"difficulty": "intermediate", "field": "language"},
        {"difficulty": "intermediate", "field": "vision"},
        {"difficulty": "advanced", "field": "decision making"},
        {"difficulty": "beginner", "field": "general"},
        {"difficulty": "intermediate", "field": "general"},
        {"difficulty": "advanced", "field": "general"},
        {"difficulty": "beginner", "field": "neural networks"},
        {"difficulty": "advanced", "field": "language"},
    ],
    ids=[f"ai_{i}" for i in range(10)]
)

print(f"Added {ai_collection.count()} documents to 'ai_concepts' collection")

# Create a collection for famous quotes
print("\nCreating 'quotes' collection...")
quotes_collection = client.get_or_create_collection("quotes")

quotes_collection.add(
    documents=[
        "The only way to do great work is to love what you do.",
        "Innovation distinguishes between a leader and a follower.",
        "Stay hungry, stay foolish.",
        "The best way to predict the future is to invent it.",
        "Any sufficiently advanced technology is indistinguishable from magic.",
        "The computer was born to solve problems that did not exist before.",
        "Software is a great combination between artistry and engineering.",
        "Code is like humor. When you have to explain it, it's bad.",
    ],
    metadatas=[
        {"author": "Steve Jobs", "topic": "work"},
        {"author": "Steve Jobs", "topic": "innovation"},
        {"author": "Steve Jobs", "topic": "life"},
        {"author": "Alan Kay", "topic": "future"},
        {"author": "Arthur C. Clarke", "topic": "technology"},
        {"author": "Bill Gates", "topic": "technology"},
        {"author": "Bill Gates", "topic": "software"},
        {"author": "Cory House", "topic": "programming"},
    ],
    ids=[f"quote_{i}" for i in range(8)]
)

print(f"Added {quotes_collection.count()} documents to 'quotes' collection")

# List all collections
print("\n" + "="*60)
print("Summary of collections created:")
print("="*60)

collections = client.list_collections()
for col in collections:
    print(f"  - {col.name}: {col.count()} documents")

print("\n✅ Sample data loaded successfully!")
print("\nYou can now start the MCP server and query these collections through Claude.")
print("\nExample queries to try:")
print("  - 'What programming languages are in the database?'")
print("  - 'Search for AI concepts related to neural networks'")
print("  - 'Find quotes about innovation'")
print("  - 'Show me compiled programming languages'")
