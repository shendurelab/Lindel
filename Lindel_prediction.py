#!/usr/bin/env python3
"""
Lindel Batch Prediction Script

This script processes DNA sequences for indel prediction using the Lindel model.
It can handle single sequences or batch processing from files.

Author: Wei Chen
Updated for batch processing and improved usability

Usage:
    # Single sequence prediction
    python lindel_batch.py -s TAACGTTATCAACGCCTATATTAAAGCGACCGTCGGTTGAACTGCGTGGATCAATGCGTC

    # Batch processing from file
    python lindel_batch.py -f input_sequences.txt -o output_results.txt

    # Batch processing with custom options
    python lindel_batch.py -f input.txt -o output.txt --format json --top 20

Input file format:
    Each line should contain one 60bp sequence (can include sequence names after a tab)
    Example:
        TAACGTTATCAACGCCTATATTAAAGCGACCGTCGGTTGAACTGCGTGGATCAATGCGTC    seq_1

Requirements:
    - 60bp sequences (30bp upstream + 30bp downstream of cut site)
    - Only ATCG characters allowed
    - PAM sequence required in positions 33-36
"""

import argparse
import sys
import os
import json
import pickle as pkl
import re
from typing import List, Dict, Tuple, Optional

# Add the current directory to Python path to import Lindel
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

try:
    from Lindel.Predictor import gen_prediction
    import Lindel
except ImportError as e:
    print(f"Error: Could not import Lindel module: {e}")
    print("Make sure the Lindel folder is in the same directory as this script.")
    print(f"Current directory: {current_dir}")
    print(f"Looking for: {os.path.join(current_dir, 'Lindel')}")
    print(f"Exists: {os.path.exists(os.path.join(current_dir, 'Lindel'))}")
    sys.exit(1)

class LindelBatchPredictor:
    def __init__(self):
        """Initialize the predictor by loading model weights and prerequisites."""
        try:
            # Load model weights and prerequisites
            lindel_path = Lindel.__path__[0]
            
            with open(os.path.join(lindel_path, "Model_weights.pkl"), 'rb') as f:
                self.weights = pkl.load(f)
            
            with open(os.path.join(lindel_path, 'model_prereq.pkl'), 'rb') as f:
                self.prerequesites = pkl.load(f)
                
            print("Model loaded successfully.")
            
        except FileNotFoundError as e:
            print(f"Error: Required model files not found: {e}")
            print("Make sure Model_weights.pkl and model_prereq.pkl are in the Lindel folder.")
            sys.exit(1)
        except Exception as e:
            print(f"Error loading model: {e}")
            sys.exit(1)
    
    def validate_sequence(self, sequence: str) -> Tuple[bool, str]:
        """
        Validate input sequence.
        
        Args:
            sequence: DNA sequence to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        sequence = sequence.upper().strip()
        
        if not re.match(r"^[ATCG]+$", sequence):
            return False, "Invalid characters in sequence. Only A, T, C, G allowed."
        
        if len(sequence) < 60:
            return False, f"Sequence too short: {len(sequence)}bp. Need at least 60bp."
        
        if len(sequence) > 60:
            print(f"Warning: Sequence longer than 60bp ({len(sequence)}bp), using first 60bp.")
            sequence = sequence[:60]

        # Check for PAM sequence (any of AGG, TGG, CGG, GGG) at positions 33-36
        pam_region = sequence[33:36]
        valid_pams = ['AGG', 'TGG', 'CGG', 'GGG']
        if pam_region not in valid_pams:
            return False, f"No valid PAM sequence found at position 33-36. Found: {pam_region}"
        
        return True, sequence[:60]
    
    def predict_single(self, sequence: str, top_n: int = 20) -> Dict:
        """
        Predict indels for a single sequence.
        
        Args:
            sequence: 60bp DNA sequence
            top_n: Number of top predictions to return
            
        Returns:
            Dictionary with prediction results
        """
        is_valid, result = self.validate_sequence(sequence)
        if not is_valid:
            return {"error": result, "sequence": sequence}
        
        sequence = result
        
        try:
            y_hat, fs = gen_prediction(sequence, self.weights, self.prerequesites)
            
            rev_index = self.prerequesites[1]
            pred_freq = {rev_index[i]: y_hat[i] for i in range(len(y_hat)) if y_hat[i] != 0}
            pred_sorted = sorted(pred_freq.items(), key=lambda kv: kv[1], reverse=True)
            
            # Limit to top N predictions
            pred_sorted = pred_sorted[:top_n]
            
            predictions = self._format_predictions(sequence, pred_sorted)
            
            return {
                "sequence": sequence,
                "frameshift_ratio": round(fs, 4),
                "num_predictions": len(predictions),
                "predictions": predictions
            }
            
        except Exception as e:
            return {"error": f"Prediction failed: {str(e)}", "sequence": sequence}
    
    def _format_predictions(self, seq: str, pred_sorted: List[Tuple]) -> List[Dict]:
        """Format predictions into a readable format."""
        predictions = []
        ss = 13
        cs = ss + 17
        
        for pt, freq in pred_sorted:
            try:
                # Deletion
                idx1, dl = map(int, pt.split('+'))
                indel_type = f"D{dl}"
                position = idx1 - 30
                
                # Create visual representation
                idx1 += cs
                idx2 = idx1 + dl
                if idx1 < cs:
                    if idx2 >= cs:
                        visual = f"{seq[:idx1]}{'-' * (cs - idx1)} | {'-' * (idx2 - cs)}{seq[idx2:]}"
                    else:
                        visual = f"{seq[:idx1]}{'-' * (idx2 - idx1)}{seq[idx2:cs]} | {seq[cs:]}"
                elif idx1 > cs:
                    visual = f"{seq[:cs]} | {seq[cs:idx1]}{'-' * dl}{seq[idx2:]}"
                else:
                    visual = f"{seq[:idx1]} | {'-' * dl}{seq[idx2:]}"
                    
                predictions.append({
                    "type": "deletion",
                    "size": dl,
                    "position": position,
                    "frequency": round(freq * 100, 3),
                    "description": f"{indel_type} at position {position}",
                    "visual": visual
                })
                
            except ValueError:
                # Insertion
                idx1 = int(pt.split('+')[0])
                if pt != '3':
                    bp = pt.split('+')[1]
                    size = len(bp)
                    description = f"I{size}+{bp}"
                else:
                    bp = 'X'  # label any insertion >= 3bp as X
                    size = "≥3"
                    description = f"I3+{bp}"
                
                visual = f"{seq[:cs]} {bp}{' ' * (2 - len(bp))}{seq[cs:]}"
                
                predictions.append({
                    "type": "insertion",
                    "size": size,
                    "sequence": bp,
                    "frequency": round(freq * 100, 3),
                    "description": description,
                    "visual": visual
                })
        
        return predictions
    
    def process_batch_file(self, input_file: str, output_file: str, output_format: str = 'tsv', top_n: int = 20):
        """
        Process batch sequences from a file.
        
        Args:
            input_file: Path to input file with sequences
            output_file: Path to output file
            output_format: Output format ('tsv', 'json', 'csv')
            top_n: Number of top predictions per sequence
        """
        try:
            with open(input_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"Error: Input file '{input_file}' not found.")
            return
        except Exception as e:
            print(f"Error reading input file: {e}")
            return
        
        results = []
        processed = 0
        errors = 0
        
        print(f"Processing {len(lines)} sequences...")
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue
            
            # Parse line (sequence and optional name)
            parts = line.split('\t')
            sequence = parts[0].strip()
            seq_name = parts[1].strip() if len(parts) > 1 else f"seq_{i+1}"
            
            print(f"Processing {i+1}/{len(lines)}: {seq_name}", end=" ... ")
            
            result = self.predict_single(sequence, top_n)
            result['name'] = seq_name
            result['index'] = i + 1
            
            if 'error' in result:
                print(f"ERROR: {result['error']}")
                errors += 1
            else:
                print(f"OK (FS: {result['frameshift_ratio']}, {result['num_predictions']} predictions)")
                processed += 1
            
            results.append(result)
        
        # Write results
        try:
            if output_format.lower() == 'json':
                self._write_json_results(results, output_file)
            elif output_format.lower() == 'csv':
                self._write_csv_results(results, output_file)
            else:  # Default to TSV
                self._write_tsv_results(results, output_file)
            
            print(f"\nResults written to: {output_file}")
            print(f"Successfully processed: {processed}")
            print(f"Errors: {errors}")
            
        except Exception as e:
            print(f"Error writing output file: {e}")
    
    def _write_json_results(self, results: List[Dict], output_file: str):
        """Write results in JSON format."""
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
    
    def _write_tsv_results(self, results: List[Dict], output_file: str):
        """Write results in TSV format."""
        with open(output_file, 'w') as f:
            # Header
            f.write("Name\tSequence\tFrameshift_Ratio\tIndel_Type\tSize\tPosition\tFrequency\tDescription\tVisual\n")
            
            for result in results:
                name = result.get('name', 'unknown')
                sequence = result.get('sequence', 'N/A')
                
                if 'error' in result:
                    f.write(f"{name}\t{sequence}\tERROR\t\t\t\t\t{result['error']}\t\n")
                    continue
                
                fs_ratio = result.get('frameshift_ratio', 0)
                
                if not result.get('predictions'):
                    f.write(f"{name}\t{sequence}\t{fs_ratio}\tNo predictions\t\t\t\t\t\n")
                    continue
                
                for pred in result['predictions']:
                    f.write(f"{name}\t{sequence}\t{fs_ratio}\t{pred['type']}\t"
                           f"{pred['size']}\t{pred.get('position', '')}\t{pred['frequency']}\t"
                           f"{pred['description']}\t{pred['visual']}\n")
    
    def _write_csv_results(self, results: List[Dict], output_file: str):
        """Write results in CSV format."""
        import csv
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Header
            writer.writerow(["Name", "Sequence", "Frameshift_Ratio", "Indel_Type", "Size", "Position", "Frequency", "Description", "Visual"])
            
            for result in results:
                name = result.get('name', 'unknown')
                sequence = result.get('sequence', 'N/A')
                
                if 'error' in result:
                    writer.writerow([name, sequence, "ERROR", "", "", "", "", result['error'], ""])
                    continue
                
                fs_ratio = result.get('frameshift_ratio', 0)
                
                if not result.get('predictions'):
                    writer.writerow([name, sequence, fs_ratio, "No predictions", "", "", "", "", ""])
                    continue
                
                for pred in result['predictions']:
                    writer.writerow([name, sequence, fs_ratio, pred['type'], pred['size'], 
                                   pred.get('position', ''), pred['frequency'], pred['description'], pred['visual']])


def main():
    parser = argparse.ArgumentParser(
        description="Lindel Batch Prediction Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-s', '--sequence', type=str, 
                           help='Single sequence to predict (60bp)')
    input_group.add_argument('-f', '--file', type=str,
                           help='Input file with sequences (one per line)')
    
    # Output options
    parser.add_argument('-o', '--output', type=str,
                       help='Output file (required for batch processing)')
    parser.add_argument('--format', choices=['tsv', 'csv', 'json'], default='tsv',
                       help='Output format (default: tsv)')
    parser.add_argument('--top', type=int, default=20,
                       help='Number of top predictions to include (default: 20)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate arguments
    if args.file and not args.output:
        parser.error("Output file (-o) is required for batch processing")
    
    # Initialize predictor
    predictor = LindelBatchPredictor()
    
    if args.sequence:
        # Single sequence prediction
        print(f"Predicting for sequence: {args.sequence}")
        result = predictor.predict_single(args.sequence, args.top)
        
        if 'error' in result:
            print(f"Error: {result['error']}")
            sys.exit(1)
        
        print(f"\nResults:")
        print(f"Sequence: {result['sequence']}")
        print(f"Frameshift ratio: {result['frameshift_ratio']}")
        print(f"Number of predictions: {result['num_predictions']}")
        print("\nTop predictions:")
        
        for i, pred in enumerate(result['predictions'][:10], 1):
            print(f"{i:2d}. {pred['description']} - {pred['frequency']:.2f}%")
            print(f"     {pred['visual']}")
        
        # Optionally save to file
        if args.output:
            if args.format == 'json':
                with open(args.output, 'w') as f:
                    json.dump(result, f, indent=2)
            else:
                predictor._write_tsv_results([{**result, 'name': 'input_sequence', 'index': 1}], args.output)
            print(f"\nResults saved to: {args.output}")
    
    else:
        # Batch processing
        print(f"Processing batch file: {args.file}")
        predictor.process_batch_file(args.file, args.output, args.format, args.top)


if __name__ == "__main__":
    main()
