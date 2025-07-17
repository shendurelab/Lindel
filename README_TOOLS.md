# Lindel Prediction Tools

This repository contains tools for predicting CRISPR/Cas9-induced indel patterns using the Lindel model.

## Overview

Lindel is a machine learning model that predicts the probability distribution of indel patterns at CRISPR/Cas9 target sites. This repository provides both a web application for Google App Engine and a standalone batch processing script.

## Tools Provided

### 1. Google App Engine Web Application

A Flask-based web service that provides:

- Single sequence prediction via REST API
- Batch prediction endpoint (up to 100 sequences)
- Full web interface with templates and styling

### 2. Standalone Batch Processing Script (`Lindel_prediction.py`)

A command-line tool for processing sequences locally:


#### Usage

```bash
# Single sequence prediction
python Lindel_prediction.py -s ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAGGCGATCGATCGATCGATCG

# Batch processing from file
python Lindel_prediction.py -f example_sequences.txt -o results.tsv

# Output in different formats
python Lindel_prediction.py -f input.txt -o results.json --format json
python Lindel_prediction.py -f input.txt -o results.csv --format csv

# Control number of top predictions
python Lindel_prediction.py -f input.txt -o results.tsv --top 10
```

#### Input File Format

The input file should contain one sequence per line, optionally with a tab-separated name:

```
TAACGTTATCAACGCCTATATTAAAGCGACCGTCGGTTGAACTGCGTGGATCAATGCGTC	seq_1
TAACGTTATCAACGCCTATATCAGAGCGACCGTTGGTAGAACTGCGTCGATCAATGCGTC	seq_2
```

## Input Requirements

- **Sequence length**: Exactly 60 base pairs
- **Sequence composition**: Only A, T, C, G characters
- **PAM sequence**: Must contain a valid PAM (AGG, TGG, CGG, or GGG) at positions 33-36
- **Cut site**: The model assumes the cut site is at position 30 (between positions 30 and 31)

## Output Formats

### TSV Format (Default)
Tab-separated values with columns:
- Name, Sequence, Frameshift_Ratio, Indel_Type, Size, Position, Frequency, Description, Visual

### JSON Format
Structured JSON with complete prediction data including visual representations.

### CSV Format
Comma-separated values, same structure as TSV.

## Model Details

The Lindel model predicts:
- **Deletions**: Size and position relative to the cut site
- **Insertions**: Size and inserted sequence
- **Frameshift probability**: Overall likelihood of causing a frameshift
- **Frequency**: Relative frequency of each indel pattern

## Example Output

For a single sequence, you might see output like:

```
Sequence: ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAGGCGATCGATCGATCGATCG
Frameshift ratio: 0.6234
Number of predictions: 15

Top predictions:
 1. D1 at position -1 - 12.45%
     ATCGATCGATCGATCGATCGATCGATCGATC | -AGATCGATCGATCGATCGATCG
 2. I1+A - 8.73%
     ATCGATCGATCGATCGATCGATCGATCGATC A AGATCGATCGATCGATCGATCG
 3. D2 at position 0 - 7.21%
     ATCGATCGATCGATCGATCGATCGATCGATC | --ATCGATCGATCGATCGATCG
```

## Error Handling

The tools provide comprehensive error handling for:
- Invalid sequence composition
- Incorrect sequence length
- Missing PAM sequences
- File I/O errors
- Model prediction errors

## Support

For questions about the Lindel model or these tools, please refer to the original paper or contact the authors.

## Citation

If you use these tools in your research, please cite the original Lindel paper: 
https://academic.oup.com/nar/article/47/15/7989/5511473
