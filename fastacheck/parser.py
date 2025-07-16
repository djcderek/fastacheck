import tracemalloc
import re
from pathlib import Path

class FastaParser:
    def __init__(self, file_path: str, streaming:bool = False):
        self.file_path = Path(file_path)
        self.streaming = streaming

        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        self.stats = {
            'total_sequences': 0,
            'total_length': 0,
            'parsing_errors': [],
            'file_size': self.file_path.stat().st_size
        }

    def _open_file(self):
        return open(self.file_path, "r")
    
    def _is_valid_header(self, line: str):
        """Check if a line is a valid FASTA header."""
        return line.startswith('>') and len(line.strip()) > 1
    
    def _is_valid_sequence_char(self, char: str) -> bool:
        """Check if a character is valid in a biological sequence."""

        valid_chars = set('ACGTUWSMKRYBDHVNX*-')
        return char.upper() in valid_chars
    
    def _validate_sequence(self, sequence: str):
        """
        Validate a sequence string.
        
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not sequence:
            return False, "Empty sequence"
        
        invalid_chars = set()
        for char in sequence:
            if not self._is_valid_sequence_char(char):
                invalid_chars.add(char)
        
        if invalid_chars:
            return False, f"Invalid characters found: {', '.join(sorted(invalid_chars))}"
        
        return True, ""
    
    def _parse_header(self, header_line: str):
        """
        Parse FASTA header into components.
        
        Returns:
            Dictionary with header information
        """
        header = header_line[1:].strip()
        
        parts = re.split(r'\s+', header, 1)
        
        result = {
            'full_header': header,
            'sequence_id': parts[0],
            'description': parts[1] if len(parts) > 1 else ''
        }
        
        if '|' in result['sequence_id']:
            # Format like: >gi|12345|ref|NM_123456.1|
            db_parts = result['sequence_id'].split('|')
            if len(db_parts) >= 4:
                result['database'] = db_parts[0]
                result['gi'] = db_parts[1]
                result['ref_type'] = db_parts[2]
                result['accession'] = db_parts[3]
        
        return result
    
    def validate_format(self):
        """
        Validate FASTA format

        Return dictionary with validation results
        """
        validation_result = {
            "is_valid": True,
            "errors": [],
            "warnings": [],
            "sequence_count": 0
        }

        in_sequence = False
        current_sequence_length = 0
        try:
            with self._open_file() as file:
                for line_num, line in enumerate(file, 1):
                    line = line.rstrip('\n\r')
                    
                    if not line.strip():
                        continue
                    
                    if line.startswith('>'):
                        if not self._is_valid_header(line):
                            validation_result['errors'].append(
                                f"Line {line_num}: Invalid header format"
                            )
                            validation_result['is_valid'] = False
                        
                        if in_sequence and current_sequence_length == 0:
                            validation_result['warnings'].append(
                                f"Line {line_num}: Empty sequence before this header"
                            )
                        
                        validation_result['sequence_count'] += 1
                        in_sequence = True
                        current_sequence_length = 0
                    
                    else:
                        if not in_sequence:
                            validation_result['errors'].append(
                                f"Line {line_num}: Sequence data without header"
                            )
                            validation_result['is_valid'] = False
                        
                        # Validate sequence characters
                        is_valid, error_msg = self._validate_sequence(line)
                        if not is_valid:
                            validation_result['errors'].append(
                                f"Line {line_num}: {error_msg}"
                            )
                            validation_result['is_valid'] = False
                        
                        current_sequence_length += len(line)
                
                # Check if file ended with a sequence
                if in_sequence and current_sequence_length == 0:
                    validation_result['warnings'].append(
                        "File ends with empty sequence"
                    )
        
        except UnicodeDecodeError as e:
            validation_result['errors'].append(f"File encoding error: {e}")
            validation_result['is_valid'] = False
        except Exception as e:
            validation_result['errors'].append(f"Unexpected error during validation: {e}")
            validation_result['is_valid'] = False
        
        return validation_result

    
    def parse_sequences(self):
        """
        Parse FASTA sequences from the file.
        
        Yields:
            Tuple of (header_dict, sequence_string)
        """
        try:
            with self._open_file() as file:
                current_header = None
                current_sequence = []
                
                for line_num, line in enumerate(file, 1):
                    self.line_number = line_num
                    line = line.rstrip('\n\r')
                    
                    # Skip empty lines
                    if not line.strip():
                        continue
                    
                    if line.startswith('>'):
                        # New header - yield previous sequence if exists
                        if current_header is not None:
                            sequence = ''.join(current_sequence)
                            self.stats['total_sequences'] += 1
                            self.stats['total_length'] += len(sequence)
                            yield current_header, sequence
                        
                        # Parse new header
                        current_header = self._parse_header(line)
                        current_sequence = []
                    
                    else:
                        # Sequence line
                        if current_header is None:
                            error_msg = f"Line {line_num}: Sequence data without header"
                            self.stats['parsing_errors'].append(error_msg)
                            raise error_msg
                        
                        # Validate and add sequence
                        is_valid, error_msg = self._validate_sequence(line)
                        if not is_valid:
                            error_msg = f"Line {line_num}: {error_msg}"
                            self.stats['parsing_errors'].append(error_msg)
                            if not self.streaming:  # In streaming mode, continue with warnings
                                raise error_msg
                        
                        current_sequence.append(line.upper())
                
                # Yield last sequence
                if current_header is not None:
                    sequence = ''.join(current_sequence)
                    self.stats['total_sequences'] += 1
                    self.stats['total_length'] += len(sequence)
                    yield current_header, sequence
        
        except Exception as e:
                error_msg = f"Unexpected parsing error: {e}"
                self.stats['parsing_errors'].append(error_msg)
                raise error_msg
        
    def get_sequences_list(self) -> list:
        """
        Get all sequences as a list (non-streaming mode).
        
        Returns:
            List of (header_dict, sequence_string) tuples
        """
        if self.streaming:
            raise ValueError("Cannot use get_sequences_list in streaming mode")
        
        return list(self.parse_sequences())
        
    def get_stats(self):
        """Get parsing statistics."""
        return self.stats.copy()

if __name__ == "__main__":
    tracemalloc.start()
    # Test the class
    parser = FastaParser("/Users/DerekChan1/Project/fastacheck/examples/realistic_genome_assembly.fasta")
    print(f"Parser initialized with: {parser.file_path}")
    validation = parser.validate_format()
    print("Validation results:", validation)
    
    if validation['is_valid']:
        # Parse sequences
        sequences = parser.get_sequences_list()
        
        print(f"\nParsed {len(sequences)} sequences:")
        for i, (header, seq) in enumerate(sequences[:3]):  # Show first 3
            print(f"\nSequence {i+1}:")
            print(f"  ID: {header['sequence_id']}")
            print(f"  Description: {header['description']}")
            print(f"  Length: {len(seq)}")
            print(f"  First 50 chars: {seq[:50]}...")
    
    # Get statistics
    stats = parser.get_stats()
    print(f"\nParsing statistics: {stats}")
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage: {current / 1024 / 1024:.2f} MB")
    print(f"Peak memory usage: {peak / 1024 / 1024:.2f} MB")

    tracemalloc.stop()