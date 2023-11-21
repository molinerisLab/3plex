#!/bin/env python
import sqlite3
import sys
# Connect to SQLite database (creates a new file if it doesn't exist)
conn = sqlite3.connect('index_tpx_stability.db')
cursor = conn.cursor()

# Create the table
cursor.execute('''
    CREATE TABLE TPX_Stability (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        tfo_start INTEGER,
        tfo_end INTEGER,
        Duplex_ID TEXT,
        TTS_start INTEGER,
        TTS_end INTEGER,
        Score REAL,
        Error_rate REAL,
        Errors TEXT,
        Motif TEXT,
        Strand TEXT,
        Orientation TEXT,
        Guanine_rate REAL,
        Stability REAL
    )
''')

# Commit the changes and close the connection
conn.commit()
first = True
for line_ in sys.stdin:
    if (first):
        first = False 
    else:
        line = line_.split("\t")
        tfo_start = int(line[1])
        tfo_end = int(line[2])
        Duplex_ID = line[3]
        chr_start = int(Duplex_ID.split(":")[-1].split("-")[0])
        TTS_start = int(line[4])+chr_start
        TTS_end = int(line[5])+chr_start
        Score = float(line[6])
        Error_rate = float(line[7])
        Errors = line[8]
        Motif = line[9]
        Strand = line[10]
        Orientation = line[11]
        Guanine_rate = float(line[12])
        Stability = float(line[13])
        # Save the object to the database
        cursor.execute('''
            INSERT INTO YourTableName
            (tfo_start, tfo_end, Duplex_ID, TTS_start, TTS_end, Score, Error_rate, Errors, Motif, Strand, Orientation, Guanine_rate, Stability)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (tfo_start, tfo_end, Duplex_ID, TTS_start, TTS_end,
            Score, Error_rate, Errors, Motif, Strand, Orientation,
            Guanine_rate, Stability))

conn.commit()
conn.close()
