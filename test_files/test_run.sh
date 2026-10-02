#!/bin/bash
python3 ../badger2.py --contrasts-a "4-1,6-4,7-6" --contrasts-b "D-W" --transform none --input test_input.csv --min-avg 10 --metadata test_metadata.csv --factorname-a TimePoint --factorname-b Treatment --pair Individual
