import json
import os
import csv
from collections import defaultdict
import sys

csvfile = os.path.splitext(sys.argv[1])[0]
jsonfile = csvfile + '.json'

FILES = defaultdict(lambda: defaultdict(list))
with open(csvfile+'.csv') as f:
    reader = csv.reader(f, delimiter=";")
    header = next(reader)
    for row in reader:
        sample_name = row[0].strip()
        for i in range(1, len(header)):
            FILES[sample_name][header[i]] = row[i]

with open(jsonfile, 'w') as f:
    json.dump(FILES, f, indent = 4, sort_keys = True)
