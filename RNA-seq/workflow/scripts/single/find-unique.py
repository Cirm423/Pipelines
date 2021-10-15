import argparse
import sys

parser = argparse.ArgumentParser(description='Get the number of uniquely aligned reads from STAR Log.final.out')

parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),default=sys.stdin,help='Log.final.out file to read')

args = parser.parse_args()

with args.infile as f:
    for line in f:
        if line.strip().startswith("Uniquely mapped reads number"):
            res=line.split("|")[1].strip()
            break

print(res)
