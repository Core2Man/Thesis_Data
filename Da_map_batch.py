# Interactive Protein ID Converter
# This script converts protein identifiers between UniProt, NCBI RefSeq, and TRO formats.
# It loads all necessary databases on startup for fast, repeated conversions.

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
import tempfile
import os
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from functools import partial
import time

# --- Color codes for output formatting ---
RED = "\033[31m"
BRIGHT_RED = "\033[91m"
CYAN = "\033[36m"
GREEN = "\033[32m"
RESET = "\033[0m"
BOLD = "\033[1m"

# --- Basic Setup ---
# Configure logging to display informational messages
logging.basicConfig(level=logging.INFO, format='%(message)s')

# --- 1. CONFIGURE YOUR DATABASE PATHS HERE ---
# Edit the paths below to point to your local FASTA files and BLAST databases.
# This is the only section you need to modify.
DATABASE_PATHS = {
    "uniprot_fasta": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/Da_uniprot.fasta",
    "refseq_fasta": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/Da_RefSeq.fasta",
    "tro_fasta": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/Da_TRO.fasta",
    "uniprot_blast_db": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/BLAST/Da_BLAST_UniProt",
    "refseq_blast_db": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/BLAST/Da_BLAST_NCBI",
    "tro_blast_db": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/BLAST/Da_BLAST_TRO",
    "deeplocpro_path": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/DeepLockDA.csv",
    "uniprot_diamond_db": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/DIAMOND/Da_DIAMOND_UniProt",
    "refseq_diamond_db": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/DIAMOND/Da_DIAMOND_NCBI",
    "tro_diamond_db": "/Users/sysoevm/Documents/Working_Folder/Da_map/data/DIAMOND/Da_DIAMOND_TRO"}
# --- End of Configuration ---


# --- Core Helper Functions ---

def get_sanitized_path(prompt_message):
    """
    Gets a file path from the user and removes leading/trailing quotes and whitespace.
    """
    path = input(prompt_message).strip()
    # Remove single and double quotes from the start and end of the path
    return path.strip("'\"")

def detect_identifier_type(protein_id):
    """
    Detects the type of a protein identifier based on its prefix.
    """
    if not isinstance(protein_id, str):
        return "Unknown"
    if protein_id.startswith("WP_") or protein_id.startswith("NP_"):
        return "NCBI RefSeq"
    elif protein_id.startswith("A0A") or (len(protein_id) > 5 and protein_id[0].isalpha() and protein_id[1].isdigit()):
        return "UniProt"
    elif protein_id.startswith("TRO"):
        return "TRO"
    else:
        return "Unknown"

def extract_protein_description(header, id_type):
    """
    Extracts protein description from FASTA header based on database type.
    """
    try:
        if id_type == "UniProt":
            # UniProt format: >tr|A0A550J8C9|A0A550J8C9_9BACT Probable dual-specificity RNA methyltransferase RlmN OS=...
            parts = header.split(" ", 1)  # Split on first space
            if len(parts) > 1:
                desc_part = parts[1]
                # Extract description before OS= (organism)
                if " OS=" in desc_part:
                    return desc_part.split(" OS=")[0].strip()
                else:
                    return desc_part.strip()
        elif id_type == "NCBI RefSeq":
            # RefSeq format varies, but typically: >WP_012345678.1 protein description [organism]
            parts = header.split(" ", 1)
            if len(parts) > 1:
                desc_part = parts[1]
                # Extract description before [organism] bracket
                if "[" in desc_part:
                    return desc_part.split("[")[0].strip()
                else:
                    return desc_part.strip()
        elif id_type == "TRO":
            # TRO format: similar to RefSeq
            parts = header.split(" ", 1)
            if len(parts) > 1:
                desc_part = parts[1]
                if "[" in desc_part:
                    return desc_part.split("[")[0].strip()
                else:
                    return desc_part.strip()
    except Exception as e:
        logging.warning(f"Could not parse description from header: {header[:50]}...")
    
    return "No description available"

def load_protein_ids(fasta_file, id_type):
    """
    Loads protein IDs, sequences, and descriptions from a FASTA file into dictionaries.
    Returns tuple: (sequences_dict, descriptions_dict)
    """
    protein_dict = {}
    description_dict = {}
    try:
        logging.info(f"Loading {id_type} database from {fasta_file}...")
        start_time = time.time()
        
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_id = None
            if id_type == "UniProt":
                try:
                    protein_id = record.id.split("|")[1]
                except IndexError:
                    logging.warning(f"Could not parse UniProt ID: {record.id}. Using full ID.")
                    protein_id = record.id
            elif id_type == "NCBI RefSeq":
                protein_id = record.id.split(" ")[0]
            elif id_type == "TRO":
                protein_id = record.id.split(" ")[0]
            
            if protein_id:
                protein_dict[protein_id] = str(record.seq)
                description_dict[protein_id] = extract_protein_description(record.description, id_type)

        load_time = time.time() - start_time
        logging.info(f"{GREEN}Loaded {id_type} protein library: {len(protein_dict)} proteins (took {load_time:.2f}s){RESET}")

    except FileNotFoundError:
        logging.error(f"{RED}FASTA file not found: {fasta_file}{RESET}")
        return None, None
    except Exception as e:
        logging.error(f"{RED}An error occurred while reading {fasta_file}: {e}{RESET}")
        return None, None
    return protein_dict, description_dict

def load_deeplocpro_data(file_path):
    """
    Load DeepLocPro data from a CSV file.
    """
    try:
        logging.info(f"Loading DeepLocPro data from {file_path}...")
        start_time = time.time()
        
        deeplocpro_df = pd.read_csv(file_path)
        
        load_time = time.time() - start_time
        logging.info(f"{GREEN}Loaded DeepLocPro file: {len(deeplocpro_df)} entries found (took {load_time:.2f}s){RESET}")
        return deeplocpro_df
    except FileNotFoundError:
        logging.error(f"{RED}DeepLocPro file not found: {file_path}{RESET}")
        return None
    except Exception as e:
        logging.error(f"{RED}An error occurred while reading DeepLocPro file: {e}{RESET}")
        return None

def find_matching_id(protein_id, source_dict, target_dict_rev):
    """
    Finds a matching ID in the target dictionary by comparing sequences.
    Uses a pre-built reverse dictionary for speed.
    """
    seq = source_dict.get(protein_id)
    if seq:
        return target_dict_rev.get(seq)
    return None

def run_blast(query_seq, db_path, cache):
    """
    Runs a BLASTp search, using a cache to avoid re-running identical searches.
    """
    if query_seq in cache:
        return cache[query_seq]
        
    result = None
    try:
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".fasta") as query_file:
            query_record = SeqRecord(Seq(query_seq), id="query_sequence")
            SeqIO.write(query_record, query_file.name, "fasta")
            query_file_path = query_file.name

        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".xml") as result_file:
            result_file_path = result_file.name

        blastp_cline = ["blastp", "-query", query_file_path, "-db", db_path, "-evalue", "0.001", "-outfmt", "5", "-out", result_file_path]
        
        subprocess.run(blastp_cline, check=True, capture_output=True, text=True)

        with open(result_file_path) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                if blast_record.alignments:
                    result = blast_record.alignments[0].hit_def.split(" ")[0]
                    break # Found the best hit
        
    except FileNotFoundError:
        logging.error(f"{RED}`blastp` command not found. Make sure BLAST+ is installed and in your system's PATH.{RESET}")
    except subprocess.CalledProcessError as e:
        logging.error(f"{RED}BLAST process failed with error: {e.stderr}{RESET}")
    except Exception as e:
        logging.error(f"{RED}An unexpected error occurred during BLAST search: {e}{RESET}")
    finally:
        if 'query_file_path' in locals() and os.path.exists(query_file_path): os.remove(query_file_path)
        if 'result_file_path' in locals() and os.path.exists(result_file_path): os.remove(result_file_path)
    
    cache[query_seq] = result # Cache the result (even if it's None)
    return result

def get_protein_localization(df, protein_id):
    """
    Get protein localization from DeepLocPro data.
    """
    if df is None:
        return None, False
    
    match = df[df["ACC"] == protein_id]
    if not match.empty:
        localization = match.iloc[0]["Localization"]
        return localization, True
    return None, False

def convert_id_comprehensive(source_id, source_dict, source_desc_dict, databases, blast_dbs, blast_cache, deeplocpro_df):
    """
    Converts a single source ID to all possible target IDs and gets localization info.
    Returns a dictionary with all conversions, descriptions, and localization data.
    """
    if not isinstance(source_id, str) or source_id not in source_dict:
        logging.warning(f"ID '{source_id}' not found in the source FASTA file. Skipping.")
        return {
            'original_id': source_id,
            'uniprot_id': source_id,
            'refseq_id': source_id,
            'tro_id': source_id,
            'uniprot_description': None,
            'refseq_description': None,
            'tro_description': None,
            'localization': None,
            'sequence': None
        }

    query_seq = source_dict.get(source_id)
    result = {
        'original_id': source_id,
        'uniprot_id': None,
        'refseq_id': None,
        'tro_id': None,
        'uniprot_description': None,
        'refseq_description': None,
        'tro_description': None,
        'localization': None,
        'sequence': query_seq
    }

    # Set the known ID and description based on source type
    source_type = detect_identifier_type(source_id)
    source_description = source_desc_dict.get(source_id, "No description available")
    
    if source_type == "UniProt":
        result['uniprot_id'] = source_id
        result['uniprot_description'] = source_description
    elif source_type == "NCBI RefSeq":
        result['refseq_id'] = source_id
        result['refseq_description'] = source_description
    elif source_type == "TRO":
        result['tro_id'] = source_id
        result['tro_description'] = source_description

    # Find matches in other databases
    if source_type == "UniProt":
        # Convert to RefSeq
        result['refseq_id'] = find_matching_id(source_id, source_dict, databases['refseq_rev'])
        if result['refseq_id']:
            result['refseq_description'] = databases['refseq_descriptions'].get(result['refseq_id'])
        elif not result['refseq_id']:
            result['refseq_id'] = run_blast(query_seq, blast_dbs['refseq'], blast_cache)
            if result['refseq_id']:
                result['refseq_description'] = databases['refseq_descriptions'].get(result['refseq_id'])
        
        # Convert to TRO
        result['tro_id'] = find_matching_id(source_id, source_dict, databases['tro_rev'])
        if result['tro_id']:
            result['tro_description'] = databases['tro_descriptions'].get(result['tro_id'])
        elif not result['tro_id']:
            result['tro_id'] = run_blast(query_seq, blast_dbs['tro'], blast_cache)
            if result['tro_id']:
                result['tro_description'] = databases['tro_descriptions'].get(result['tro_id'])

    elif source_type == "NCBI RefSeq":
        # Convert to UniProt
        result['uniprot_id'] = find_matching_id(source_id, source_dict, databases['uniprot_rev'])
        if result['uniprot_id']:
            result['uniprot_description'] = databases['uniprot_descriptions'].get(result['uniprot_id'])
        elif not result['uniprot_id']:
            result['uniprot_id'] = run_blast(query_seq, blast_dbs['uniprot'], blast_cache)
            if result['uniprot_id']:
                result['uniprot_description'] = databases['uniprot_descriptions'].get(result['uniprot_id'])
        
        # Convert to TRO
        result['tro_id'] = find_matching_id(source_id, source_dict, databases['tro_rev'])
        if result['tro_id']:
            result['tro_description'] = databases['tro_descriptions'].get(result['tro_id'])
        elif not result['tro_id']:
            result['tro_id'] = run_blast(query_seq, blast_dbs['tro'], blast_cache)
            if result['tro_id']:
                result['tro_description'] = databases['tro_descriptions'].get(result['tro_id'])

    elif source_type == "TRO":
        # Convert to UniProt
        result['uniprot_id'] = find_matching_id(source_id, source_dict, databases['uniprot_rev'])
        if result['uniprot_id']:
            result['uniprot_description'] = databases['uniprot_descriptions'].get(result['uniprot_id'])
        elif not result['uniprot_id']:
            result['uniprot_id'] = run_blast(query_seq, blast_dbs['uniprot'], blast_cache)
            if result['uniprot_id']:
                result['uniprot_description'] = databases['uniprot_descriptions'].get(result['uniprot_id'])
        
        # Convert to RefSeq
        result['refseq_id'] = find_matching_id(source_id, source_dict, databases['refseq_rev'])
        if result['refseq_id']:
            result['refseq_description'] = databases['refseq_descriptions'].get(result['refseq_id'])
        elif not result['refseq_id']:
            result['refseq_id'] = run_blast(query_seq, blast_dbs['refseq'], blast_cache)
            if result['refseq_id']:
                result['refseq_description'] = databases['refseq_descriptions'].get(result['refseq_id'])

    # Get localization info (prefer RefSeq ID for lookup)
    lookup_id = result['refseq_id'] or source_id
    if lookup_id:
        localization, found = get_protein_localization(deeplocpro_df, lookup_id)
        if found:
            result['localization'] = localization

    return result

def convert_id_simple(source_id, source_dict, target_dict_rev, target_blast_db, blast_cache):
    """
    Converts a single source ID to a target ID (original simple function for backward compatibility).
    """
    if not isinstance(source_id, str) or source_id not in source_dict:
        logging.warning(f"ID '{source_id}' not found in the source FASTA file. Skipping.")
        return source_id

    # 1. Try direct sequence match using the pre-built reverse dictionary
    target_id = find_matching_id(source_id, source_dict, target_dict_rev)

    # 2. Fall back to BLAST
    if not target_id:
        logging.info(f"No direct match for {source_id}. Running BLAST...")
        query_seq = source_dict.get(source_id)
        if query_seq:
            target_id = run_blast(query_seq, target_blast_db, blast_cache)

    return target_id if target_id else source_id


def check_diamond_installed():
    """
    Check if DIAMOND is installed and accessible.
    """
    try:
        result = subprocess.run(['diamond', 'version'], capture_output=True, text=True, check=True)
        version_info = result.stdout.strip()
        logging.info(f"{GREEN}DIAMOND found: {version_info}{RESET}")
        return True
    except FileNotFoundError:
        logging.warning(f"{RED}DIAMOND not found. Please install DIAMOND for faster searches.{RESET}")
        logging.info(
            f"{CYAN}To install: conda install -c bioconda diamond OR download from https://github.com/bbuchfink/diamond{RESET}")
        return False
    except subprocess.CalledProcessError:
        logging.warning(f"{RED}DIAMOND found but not working properly.{RESET}")
        return False


def create_diamond_database(fasta_path, diamond_db_path):
    """
    Creates a DIAMOND database from a FASTA file.
    Returns True if successful, False otherwise.
    """
    if not os.path.exists(fasta_path):
        logging.error(f"{RED}FASTA file not found: {fasta_path}{RESET}")
        return False

    # Check if DIAMOND is installed
    if not check_diamond_installed():
        return False

    try:
        # Create directory if it doesn't exist
        db_dir = os.path.dirname(diamond_db_path)
        if db_dir and not os.path.exists(db_dir):
            os.makedirs(db_dir, exist_ok=True)

        # Check if database already exists and is newer than FASTA file
        if os.path.exists(diamond_db_path + ".dmnd"):
            fasta_mtime = os.path.getmtime(fasta_path)
            db_mtime = os.path.getmtime(diamond_db_path + ".dmnd")
            if db_mtime > fasta_mtime:
                logging.info(f"{GREEN}DIAMOND database already up-to-date: {diamond_db_path}.dmnd{RESET}")
                return True

        logging.info(f"Creating DIAMOND database from {fasta_path}...")

        # Run diamond makedb command
        cmd = [
            'diamond', 'makedb',
            '--in', fasta_path,
            '--db', diamond_db_path,
            '--threads', str(os.cpu_count() or 4)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        if os.path.exists(diamond_db_path + ".dmnd"):
            file_size = os.path.getsize(diamond_db_path + ".dmnd") / (1024 * 1024)  # Size in MB
            logging.info(
                f"{GREEN}Successfully created DIAMOND database: {diamond_db_path}.dmnd ({file_size:.1f} MB){RESET}")
            return True
        else:
            logging.error(f"{RED}DIAMOND database creation failed - output file not found{RESET}")
            return False

    except subprocess.CalledProcessError as e:
        logging.error(f"{RED}DIAMOND makedb failed: {e.stderr}{RESET}")
        return False
    except Exception as e:
        logging.error(f"{RED}Error creating DIAMOND database: {e}{RESET}")
        return False


def run_diamond(query_seq, db_path, cache):
    """
    Runs a DIAMOND blastp search, using a cache to avoid re-running identical searches.
    Returns the best hit ID or None.
    """
    if query_seq in cache:
        return cache[query_seq]

    result = None

    # Check if DIAMOND database exists
    if not os.path.exists(db_path + ".dmnd"):
        logging.debug(f"DIAMOND database not found: {db_path}.dmnd - falling back to BLAST")
        return None

    try:
        # Create temporary query file
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".fasta") as query_file:
            query_record = SeqRecord(Seq(query_seq), id="query_sequence")
            SeqIO.write(query_record, query_file.name, "fasta")
            query_file_path = query_file.name

        # Create temporary output file
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".tsv") as result_file:
            result_file_path = result_file.name

        # Run DIAMOND blastp
        cmd = [
            'diamond', 'blastp',
            '--db', db_path,
            '--query', query_file_path,
            '--out', result_file_path,
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
            'send', 'evalue', 'bitscore',
            '--max-target-seqs', '1',
            '--evalue', '0.001',
            '--threads', str(os.cpu_count() or 4),
            '--quiet'
        ]

        subprocess.run(cmd, check=True, capture_output=True, text=True)

        # Parse results
        if os.path.exists(result_file_path) and os.path.getsize(result_file_path) > 0:
            with open(result_file_path, 'r') as f:
                first_line = f.readline().strip()
                if first_line:
                    # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                    fields = first_line.split('\t')
                    if len(fields) >= 2:
                        subject_id = fields[1]
                        # Remove version number if present (e.g., WP_123456.1 -> WP_123456)
                        result = subject_id.split('.')[0] if '.' in subject_id else subject_id

    except subprocess.CalledProcessError as e:
        logging.debug(f"DIAMOND search failed, will fall back to BLAST")
    except FileNotFoundError:
        logging.debug(f"DIAMOND command not found")
    except Exception as e:
        logging.debug(f"Error during DIAMOND search: {e}")
    finally:
        # Clean up temporary files
        if 'query_file_path' in locals() and os.path.exists(query_file_path):
            os.remove(query_file_path)
        if 'result_file_path' in locals() and os.path.exists(result_file_path):
            os.remove(result_file_path)

    cache[query_seq] = result  # Cache the result (even if None)
    return result


def convert_id_simple_diamond(source_id, source_dict, target_dict_rev, target_blast_db, target_diamond_db, blast_cache,
                              use_diamond=True):
    """
    Converts a single source ID to a target ID using DIAMOND or BLAST.
    """
    if not isinstance(source_id, str) or source_id not in source_dict:
        logging.warning(f"ID '{source_id}' not found in the source FASTA file. Skipping.")
        return source_id

    # 1. Try direct sequence match using the pre-built reverse dictionary
    target_id = find_matching_id(source_id, source_dict, target_dict_rev)

    # 2. Fall back to DIAMOND or BLAST
    if not target_id:
        query_seq = source_dict.get(source_id)
        if query_seq:
            if use_diamond and target_diamond_db:
                # Try DIAMOND first
                target_id = run_diamond(query_seq, target_diamond_db, blast_cache)

            # Fall back to BLAST if DIAMOND didn't work or wasn't requested
            if not target_id:
                target_id = run_blast(query_seq, target_blast_db, blast_cache)

    return target_id if target_id else source_id


def convert_id_comprehensive_diamond(source_id, source_dict, source_desc_dict, databases, blast_dbs, diamond_dbs,
                                     blast_cache, deeplocpro_df, use_diamond=True):
    """
    Converts a single source ID to all possible target IDs using DIAMOND or BLAST.
    Returns a dictionary with all conversions, descriptions, and localization data.
    """
    if not isinstance(source_id, str) or source_id not in source_dict:
        logging.warning(f"ID '{source_id}' not found in the source FASTA file. Skipping.")
        return {
            'original_id': source_id,
            'uniprot_id': None,
            'refseq_id': None,
            'tro_id': None,
            'uniprot_description': None,
            'refseq_description': None,
            'tro_description': None,
            'localization': None,
            'sequence': None
        }

    query_seq = source_dict.get(source_id)
    result = {
        'original_id': source_id,
        'uniprot_id': None,
        'refseq_id': None,
        'tro_id': None,
        'uniprot_description': None,
        'refseq_description': None,
        'tro_description': None,
        'localization': None,
        'sequence': query_seq
    }

    # Set the known ID and description based on source type
    source_type = detect_identifier_type(source_id)
    source_description = source_desc_dict.get(source_id, "No description available")

    if source_type == "UniProt":
        result['uniprot_id'] = source_id
        result['uniprot_description'] = source_description
    elif source_type == "NCBI RefSeq":
        result['refseq_id'] = source_id
        result['refseq_description'] = source_description
    elif source_type == "TRO":
        result['tro_id'] = source_id
        result['tro_description'] = source_description

    # Helper function to search with DIAMOND/BLAST fallback
    def search_database(db_name, query_sequence):
        target_id = None
        if use_diamond and diamond_dbs and db_name in diamond_dbs:
            target_id = run_diamond(query_sequence, diamond_dbs[db_name], blast_cache)
        if not target_id and db_name in blast_dbs:
            target_id = run_blast(query_sequence, blast_dbs[db_name], blast_cache)
        return target_id

    # Find matches in other databases (same logic as convert_id_comprehensive but with DIAMOND support)
    if source_type == "UniProt":
        # Convert to RefSeq
        result['refseq_id'] = find_matching_id(source_id, source_dict, databases['refseq_rev'])
        if result['refseq_id']:
            result['refseq_description'] = databases['refseq_descriptions'].get(result['refseq_id'])
        else:
            result['refseq_id'] = search_database('refseq', query_seq)
            if result['refseq_id']:
                result['refseq_description'] = databases['refseq_descriptions'].get(result['refseq_id'])

        # Convert to TRO
        result['tro_id'] = find_matching_id(source_id, source_dict, databases['tro_rev'])
        if result['tro_id']:
            result['tro_description'] = databases['tro_descriptions'].get(result['tro_id'])
        else:
            result['tro_id'] = search_database('tro', query_seq)
            if result['tro_id']:
                result['tro_description'] = databases['tro_descriptions'].get(result['tro_id'])

    elif source_type == "NCBI RefSeq":
        # Convert to UniProt
        result['uniprot_id'] = find_matching_id(source_id, source_dict, databases['uniprot_rev'])
        if result['uniprot_id']:
            result['uniprot_description'] = databases['uniprot_descriptions'].get(result['uniprot_id'])
        else:
            result['uniprot_id'] = search_database('uniprot', query_seq)
            if result['uniprot_id']:
                result['uniprot_description'] = databases['uniprot_descriptions'].get(result['uniprot_id'])

        # Convert to TRO
        result['tro_id'] = find_matching_id(source_id, source_dict, databases['tro_rev'])
        if result['tro_id']:
            result['tro_description'] = databases['tro_descriptions'].get(result['tro_id'])
        else:
            result['tro_id'] = search_database('tro', query_seq)
            if result['tro_id']:
                result['tro_description'] = databases['tro_descriptions'].get(result['tro_id'])

    elif source_type == "TRO":
        # Similar logic for TRO...
        result['uniprot_id'] = find_matching_id(source_id, source_dict, databases['uniprot_rev'])
        if result['uniprot_id']:
            result['uniprot_description'] = databases['uniprot_descriptions'].get(result['uniprot_id'])
        else:
            result['uniprot_id'] = search_database('uniprot', query_seq)
            if result['uniprot_id']:
                result['uniprot_description'] = databases['uniprot_descriptions'].get(result['uniprot_id'])

        result['refseq_id'] = find_matching_id(source_id, source_dict, databases['refseq_rev'])
        if result['refseq_id']:
            result['refseq_description'] = databases['refseq_descriptions'].get(result['refseq_id'])
        else:
            result['refseq_id'] = search_database('refseq', query_seq)
            if result['refseq_id']:
                result['refseq_description'] = databases['refseq_descriptions'].get(result['refseq_id'])

    # Get localization info (prefer RefSeq ID for lookup)
    lookup_id = result['refseq_id'] or source_id
    if lookup_id:
        localization, found = get_protein_localization(deeplocpro_df, lookup_id)
        if found:
            result['localization'] = localization

    return result

# --- Mode-Specific Logic ---

def single_protein_mode_with_diamond(databases):
    """    Handles the logic for converting a single protein ID using pre-loaded databases."""
    print(f"\n{CYAN}--- Single Protein Mode ---{RESET}")
    blast_cache = {}  # Simple cache for this session

    blast_dbs = {
        'uniprot': DATABASE_PATHS["uniprot_blast_db"],
        'refseq': DATABASE_PATHS["refseq_blast_db"],
        'tro': DATABASE_PATHS["tro_blast_db"]
    }
    diamond_dbs = {
        'uniprot': DATABASE_PATHS.get("uniprot_diamond_db"),
        'refseq': DATABASE_PATHS.get("refseq_diamond_db"),
        'tro': DATABASE_PATHS.get("tro_diamond_db")
    }

    # Remove DIAMOND options since they're not implemented
    # If you want to keep the UI for future implementation, you can leave this in
    # but it won't actually use DIAMOND

    # Ask user about search method preference (keeping for future compatibility)
    print(f"\nSearch method:")
    print("1: BLAST (standard method)")
    print("2: Quick mode (sequence matching only, no BLAST)")

    search_choice = input("Enter choice (1-2, default is 1): ").strip()
    use_blast = True  # Default

    if search_choice == "2":
        use_blast = False
        print(f"{GREEN}Using quick mode (sequence matching only).{RESET}")
    else:
        use_blast = True
        print(f"{GREEN}Using BLAST for sequence searches.{RESET}")

    while True:
        user_input = input(f"\nEnter a protein ID (or type 'exit' to quit): ").strip()
        if user_input.lower() == 'exit':
            break

        id_type = detect_identifier_type(user_input)
        print(f"{GREEN}Detected ID type: {id_type}{RESET}")

        if id_type == "Unknown":
            print(f"{RED}Unknown or unsupported protein ID format.{RESET}")
            continue

        # Get the appropriate source dictionary and descriptions
        if id_type == "UniProt":
            source_dict = databases["uniprot"]
            source_desc_dict = databases["uniprot_descriptions"]
        elif id_type == "NCBI RefSeq":
            source_dict = databases["refseq"]
            source_desc_dict = databases["refseq_descriptions"]
        elif id_type == "TRO":
            source_dict = databases["tro"]
            source_desc_dict = databases["tro_descriptions"]
        else:
            print(f"{RED}Unsupported ID type.{RESET}")
            continue

        # Check if the ID exists in the database
        if user_input not in source_dict:
            print(f"{RED}ID '{user_input}' not found in the {id_type} database.{RESET}")
            continue

        # Perform comprehensive conversion using the existing function
        result = convert_id_comprehensive_diamond(
            user_input,
            source_dict,
            source_desc_dict,
            databases,
            blast_dbs,
            diamond_dbs,  # ADD this parameter
            blast_cache,
            databases.get("deeplocpro"),
            use_diamond=True  # ADD this parameter
        )

        # Display results
        print(f"\n{RED}Original input:{RESET} {BOLD}{result['original_id']}{RESET}")

        if result['sequence']:
            # Show first 50 amino acids of sequence (truncate if too long)
            seq_display = result['sequence'][:50]
            if len(result['sequence']) > 50:
                seq_display += "..."
            print(f"Amino acid sequence ({len(result['sequence'])} aa): {seq_display}")

        # Display UniProt result
        if result['uniprot_id'] and result['uniprot_id'] != user_input:
            print(f"\n{GREEN}UniProt reference:{RESET} {result['uniprot_id']}")
            if result['uniprot_description']:
                print(f"  Description: {result['uniprot_description']}")
        elif id_type == "UniProt":
            print(f"\n{GREEN}UniProt reference:{RESET} {result['original_id']} (original)")
            if result['uniprot_description']:
                print(f"  Description: {result['uniprot_description']}")
        else:
            print(f"\n{GREEN}UniProt reference:{RESET} Not found")

        # Display RefSeq result
        if result['refseq_id'] and result['refseq_id'] != user_input:
            print(f"\n{GREEN}RefSeq reference:{RESET} {result['refseq_id']}")
            if result['refseq_description']:
                print(f"  Description: {result['refseq_description']}")
        elif id_type == "NCBI RefSeq":
            print(f"\n{GREEN}RefSeq reference:{RESET} {result['original_id']} (original)")
            if result['refseq_description']:
                print(f"  Description: {result['refseq_description']}")
        else:
            print(f"\n{GREEN}RefSeq reference:{RESET} Not found")

        # Display TRO result
        if result['tro_id'] and result['tro_id'] != user_input:
            print(f"\n{GREEN}TRO reference:{RESET} {result['tro_id']}")
            if result['tro_description']:
                print(f"  Description: {result['tro_description']}")
        elif id_type == "TRO":
            print(f"\n{GREEN}TRO reference:{RESET} {result['original_id']} (original)")
            if result['tro_description']:
                print(f"  Description: {result['tro_description']}")
        else:
            print(f"\n{GREEN}TRO reference:{RESET} Not found")

        # Display localization information
        if result['localization']:
            lookup_id = result['refseq_id'] or user_input
            print(
                f"\n{GREEN}Predicted cellular localization for {RED}{lookup_id}{RESET} is {RED}{result['localization']}{RESET}.")
        else:
            print(f"\n{RED}No localization data available for this protein.{RESET}")

        # Offer to save results
        save_choice = input("\nSave results to file? (y/N): ").strip().lower()
        if save_choice in ['y', 'yes']:
            output_file = get_sanitized_path("Enter output filename (e.g., results.txt): ")
            try:
                with open(output_file, 'w') as f:
                    f.write(f"Protein ID Conversion Results\n")
                    f.write(f"{'=' * 50}\n")
                    f.write(f"Original ID: {result['original_id']}\n")
                    f.write(f"ID Type: {id_type}\n\n")

                    f.write(f"UniProt ID: {result['uniprot_id'] or 'Not found'}\n")
                    if result['uniprot_description']:
                        f.write(f"UniProt Description: {result['uniprot_description']}\n")
                    f.write(f"\n")

                    f.write(f"RefSeq ID: {result['refseq_id'] or 'Not found'}\n")
                    if result['refseq_description']:
                        f.write(f"RefSeq Description: {result['refseq_description']}\n")
                    f.write(f"\n")

                    f.write(f"TRO ID: {result['tro_id'] or 'Not found'}\n")
                    if result['tro_description']:
                        f.write(f"TRO Description: {result['tro_description']}\n")
                    f.write(f"\n")

                    if result['localization']:
                        f.write(f"Predicted Localization: {result['localization']}\n")

                    if result['sequence']:
                        f.write(f"\nSequence ({len(result['sequence'])} aa):\n")
                        # Write sequence in FASTA format (60 characters per line)
                        for i in range(0, len(result['sequence']), 60):
                            f.write(f"{result['sequence'][i:i + 60]}\n")

                print(f"{GREEN}Results saved to {output_file}{RESET}")
            except Exception as e:
                print(f"{RED}Failed to save results: {e}{RESET}")


def bulk_file_mode_with_diamond(databases):
    """    Handles the logic for converting a CSV or Excel file in bulk using pre-loaded databases.
    Now with DIAMOND support for faster searches.    """
    print(f"\n{CYAN}--- Bulk File Mode (CSV/Excel) ---{RESET}")
    input_file = get_sanitized_path("Enter path to the input file (CSV or Excel): ")

    # Determine file type and read accordingly
    file_extension = os.path.splitext(input_file)[1].lower()

    try:
        if file_extension == '.csv':
            df = pd.read_csv(input_file)
            print(f"{GREEN}Successfully loaded CSV file.{RESET}")
        elif file_extension in ['.xlsx', '.xls']:
            # For Excel files, we might need to specify which sheet to use
            try:
                # First, try to read the first sheet
                df = pd.read_excel(input_file)
                print(f"{GREEN}Successfully loaded Excel file (first sheet).{RESET}")
            except Exception as e:
                # If that fails, let user choose the sheet
                try:
                    # Get sheet names
                    excel_file = pd.ExcelFile(input_file)
                    sheet_names = excel_file.sheet_names

                    if len(sheet_names) > 1:
                        print(f"\nAvailable sheets in the Excel file:")
                        for i, sheet in enumerate(sheet_names, 1):
                            print(f"  {i}: {sheet}")

                        while True:
                            try:
                                sheet_choice = input("Enter sheet number or name: ").strip()
                                if sheet_choice.isdigit():
                                    sheet_index = int(sheet_choice) - 1
                                    if 0 <= sheet_index < len(sheet_names):
                                        selected_sheet = sheet_names[sheet_index]
                                        break
                                    else:
                                        print(f"{RED}Invalid sheet number. Please try again.{RESET}")
                                elif sheet_choice in sheet_names:
                                    selected_sheet = sheet_choice
                                    break
                                else:
                                    print(f"{RED}Sheet not found. Please try again.{RESET}")
                            except ValueError:
                                print(f"{RED}Invalid input. Please try again.{RESET}")
                    else:
                        selected_sheet = sheet_names[0]

                    df = pd.read_excel(input_file, sheet_name=selected_sheet)
                    print(f"{GREEN}Successfully loaded Excel sheet '{selected_sheet}'.{RESET}")

                except Exception as inner_e:
                    logging.error(f"{RED}Failed to read Excel file: {inner_e}{RESET}")
                    return
        else:
            logging.error(f"{RED}Unsupported file format. Please use .csv, .xlsx, or .xls files.{RESET}")
            return

        # Check for ProteinID column
        protein_col = 'ProteinID'
        if protein_col not in df.columns:
            # If ProteinID column doesn't exist, show available columns and let user choose
            print(f"\n{RED}Column '{protein_col}' not found.{RESET}")
            print("Available columns:")
            for i, col in enumerate(df.columns, 1):
                print(f"  {i}: {col}")

            while True:
                try:
                    col_choice = input("Enter column number or name containing protein IDs: ").strip()
                    if col_choice.isdigit():
                        col_index = int(col_choice) - 1
                        if 0 <= col_index < len(df.columns):
                            protein_col = df.columns[col_index]
                            break
                        else:
                            print(f"{RED}Invalid column number. Please try again.{RESET}")
                    elif col_choice in df.columns:
                        protein_col = col_choice
                        break
                    else:
                        print(f"{RED}Column not found. Please try again.{RESET}")
                except ValueError:
                    print(f"{RED}Invalid input. Please try again.{RESET}")

            print(f"{GREEN}Using column '{protein_col}' for protein IDs.{RESET}")

    except Exception as e:
        logging.error(f"{RED}Failed to read file: {e}{RESET}")
        return

    # Check if the column has any data
    if df[protein_col].dropna().empty:
        logging.error(f"{RED}No protein IDs found in column '{protein_col}'.{RESET}")
        return

    first_id = df[protein_col].dropna().iloc[0]
    source_type = detect_identifier_type(first_id)
    if source_type == "Unknown":
        print(f"{RED}Could not automatically detect the protein ID format.{RESET}")
        return

    print(f"\n{GREEN}Detected '{source_type}' format in the '{protein_col}' column.{RESET}")

    # Ask user about search method preference
    print(f"\nSearch method:")
    print("1: DIAMOND (faster, recommended)")
    print("2: BLAST (slower, traditional)")
    print("3: Auto (try DIAMOND first, fall back to BLAST)")

    search_choice = input("Enter choice (1-3, default is 3): ").strip()
    use_diamond = True  # Default

    if search_choice == "2":
        use_diamond = False
        print(f"{GREEN}Using BLAST for sequence searches.{RESET}")
    elif search_choice == "1":
        use_diamond = True
        print(f"{GREEN}Using DIAMOND for sequence searches.{RESET}")
    else:
        print(f"{GREEN}Using auto mode (DIAMOND preferred).{RESET}")

    # Ask user what type of conversion they want
    print("\nConversion options:")
    print("1: Simple conversion (convert to one target format)")
    print("2: Comprehensive conversion (convert to all formats + localization)")

    mode_choice = input("Enter choice (1 or 2): ").strip()

    blast_dbs = {
        'uniprot': DATABASE_PATHS["uniprot_blast_db"],
        'refseq': DATABASE_PATHS["refseq_blast_db"],
        'tro': DATABASE_PATHS["tro_blast_db"]
    }
    diamond_dbs = {
        'uniprot': DATABASE_PATHS.get("uniprot_diamond_db"),
        'refseq': DATABASE_PATHS.get("refseq_diamond_db"),
        'tro': DATABASE_PATHS.get("tro_diamond_db")
    }

    if mode_choice == "1":
        # Simple conversion mode
        source_dict, target_dict_rev, target_blast_db = (None, None, None)

        if source_type == "UniProt":
            print("The script will convert IDs to 'NCBI RefSeq'.")
            source_dict = databases["uniprot"]
            target_dict_rev = databases["refseq_rev"]
            target_blast_db = DATABASE_PATHS["refseq_blast_db"]
            # target_diamond_db = DATABASE_PATHS.get("refseq_diamond_db", "")
        elif source_type == "NCBI RefSeq":
            print("The script will convert IDs to 'UniProt'.")
            source_dict = databases["refseq"]
            target_dict_rev = databases["uniprot_rev"]
            target_blast_db = DATABASE_PATHS["uniprot_blast_db"]
            # target_diamond_db = DATABASE_PATHS.get("uniprot_diamond_db", "")
        elif source_type == "TRO":
            print("The script will convert IDs to 'UniProt'.")
            source_dict = databases["tro"]
            target_dict_rev = databases["uniprot_rev"]
            target_blast_db = DATABASE_PATHS["uniprot_blast_db"]
            # target_diamond_db = DATABASE_PATHS.get("uniprot_diamond_db", "")

        output_file = get_sanitized_path("Enter path for the output file (will match input format): ")

        logging.info(f"Starting parallel protein ID conversion...")
        ids_to_convert = df[protein_col].tolist()
        blast_cache = {}  # Cache for this bulk run
        target_diamond_db = diamond_dbs.get('refseq')
        # NOTE: Using regular convert_id_simple since convert_id_simple_diamond doesn't exist
        converter_func = partial(convert_id_simple_diamond,
                                 source_dict=source_dict,
                                 target_dict_rev=target_dict_rev,
                                 target_blast_db=target_blast_db,
                                 target_diamond_db=target_diamond_db,
                                 blast_cache=blast_cache,
                                 use_diamond=use_diamond)

        num_workers = os.cpu_count() or 4
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            results_iterator = executor.map(converter_func, ids_to_convert)
            converted_ids = list(tqdm(results_iterator, total=len(ids_to_convert), desc="Converting IDs"))

        # Update the DataFrame with converted IDs
        df[protein_col] = converted_ids

    elif mode_choice == "2":
        # Comprehensive conversion mode
        output_file = get_sanitized_path("Enter path for the output file (will match input format): ")

        # Get the appropriate source dictionary and descriptions
        if source_type == "UniProt":
            source_dict = databases["uniprot"]
            source_desc_dict = databases["uniprot_descriptions"]
        elif source_type == "NCBI RefSeq":
            source_dict = databases["refseq"]
            source_desc_dict = databases["refseq_descriptions"]
        elif source_type == "TRO":
            source_dict = databases["tro"]
            source_desc_dict = databases["tro_descriptions"]

        logging.info(f"Starting comprehensive parallel protein ID conversion...")
        ids_to_convert = df[protein_col].tolist()
        blast_cache = {}  # Cache for this bulk run
        target_diamond_db = diamond_dbs.get('refseq')

        converter_func = partial(convert_id_comprehensive_diamond,
                                 source_dict=source_dict,
                                 source_desc_dict=source_desc_dict,
                                 databases=databases,
                                 blast_dbs=blast_dbs,
                                 diamond_dbs=diamond_dbs,
                                 blast_cache=blast_cache,
                                 deeplocpro_df=databases.get("deeplocpro"),
                                 use_diamond=use_diamond)

        num_workers = os.cpu_count() or 4
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            results_iterator = executor.map(converter_func, ids_to_convert)
            conversion_results = list(tqdm(results_iterator, total=len(ids_to_convert), desc="Converting IDs"))

        # Create new DataFrame with all results
        results_df = pd.DataFrame(conversion_results)

        # Add original data columns (except the protein ID column which is now expanded)
        for col in df.columns:
            if col != protein_col:
                results_df[col] = df[col].values

        df = results_df

    else:
        print(f"{RED}Invalid choice. Returning to main menu.{RESET}")
        return

    # Save output file (match input format)
    try:
        output_extension = os.path.splitext(output_file)[1].lower()
        input_extension = os.path.splitext(input_file)[1].lower()

        # If no extension specified, use the same as input
        if not output_extension:
            output_file += input_extension
            output_extension = input_extension

        if output_extension == '.csv':
            df.to_csv(output_file, index=False)
            logging.info(f"{GREEN}Conversion complete. Saved to {output_file}{RESET}")
        elif output_extension in ['.xlsx', '.xls']:
            df.to_excel(output_file, index=False)
            logging.info(f"{GREEN}Conversion complete. Saved to {output_file}{RESET}")
        else:
            # Default to CSV if unsupported output format
            logging.warning(f"{RED}Unsupported output format. Saving as CSV.{RESET}")
            output_file = os.path.splitext(output_file)[0] + '.csv'
            df.to_csv(output_file, index=False)
            logging.info(f"{GREEN}Conversion complete. Saved to {output_file}{RESET}")

    except Exception as e:
        logging.error(f"{RED}Failed to save output file: {e}{RESET}")
        # Try to save as CSV as fallback
        try:
            fallback_file = os.path.splitext(output_file)[0] + '_fallback.csv'
            df.to_csv(fallback_file, index=False)
            logging.info(f"{GREEN}Saved as CSV fallback: {fallback_file}{RESET}")
        except:
            logging.error(f"{RED}Could not save file in any format.{RESET}")

# --- Main Execution ---

def main():
    """
    Main function to run the script with DIAMOND support.
    """
    print("=" * 60)
    print(f"{CYAN}           Protein ID Conversion Tool{RESET}")
    print("=" * 60)
    print(f"{GREEN}Program version from 24 August 2025 - Enhanced with DIAMOND{RESET}")

    # --- 2. Load all databases on startup ---
    logging.info(f"{CYAN}Starting program. Loading databases...{RESET}")
    databases = {}
    diamond_tasks = [
        ("UniProt", DATABASE_PATHS["uniprot_fasta"], DATABASE_PATHS["uniprot_diamond_db"]),
        ("RefSeq", DATABASE_PATHS["refseq_fasta"], DATABASE_PATHS["refseq_diamond_db"]),
        ("TRO", DATABASE_PATHS["tro_fasta"], DATABASE_PATHS["tro_diamond_db"])
    ]

    try:
        with ThreadPoolExecutor(max_workers=4) as executor:
            # Submit all database loading tasks
            uniprot_future = executor.submit(load_protein_ids, DATABASE_PATHS["uniprot_fasta"], "UniProt")
            refseq_future = executor.submit(load_protein_ids, DATABASE_PATHS["refseq_fasta"], "NCBI RefSeq")
            tro_future = executor.submit(load_protein_ids, DATABASE_PATHS["tro_fasta"], "TRO")
            deeplocpro_future = executor.submit(load_deeplocpro_data, DATABASE_PATHS["deeplocpro_path"])

            # Get results (sequences and descriptions)
            databases["uniprot"], databases["uniprot_descriptions"] = uniprot_future.result()
            databases["refseq"], databases["refseq_descriptions"] = refseq_future.result()
            databases["tro"], databases["tro_descriptions"] = tro_future.result()
            databases["deeplocpro"] = deeplocpro_future.result()

        # Check if critical databases loaded successfully
        if not all([databases["uniprot"], databases["refseq"], databases["tro"]]):
            raise ValueError("One or more critical databases failed to load. Please check paths and file integrity.")

        # --- 3. Pre-build reverse dictionaries for fast lookups ---
        logging.info(f"{CYAN}Building reverse dictionaries for matching...{RESET}")
        start_time = time.time()

        databases["uniprot_rev"] = {v: k for k, v in databases["uniprot"].items()}
        databases["refseq_rev"] = {v: k for k, v in databases["refseq"].items()}
        databases["tro_rev"] = {v: k for k, v in databases["tro"].items()}

        build_time = time.time() - start_time
        logging.info(f"{GREEN}Reverse dictionaries built (took {build_time:.2f}s){RESET}")

        # --- 4. Create DIAMOND databases if possible ---
        logging.info(f"{CYAN}Checking/creating DIAMOND databases for faster searches...{RESET}")

        diamond_success = {}
        diamond_tasks = [
            ("UniProt", DATABASE_PATHS["uniprot_fasta"], DATABASE_PATHS["uniprot_diamond_db"]),
            ("RefSeq", DATABASE_PATHS["refseq_fasta"], DATABASE_PATHS["refseq_diamond_db"]),
            ("TRO", DATABASE_PATHS["tro_fasta"], DATABASE_PATHS["tro_diamond_db"])
        ]

        for db_name, fasta_path, diamond_path in diamond_tasks:
            success = create_diamond_database(fasta_path, diamond_path)
            diamond_success[db_name] = success

        if any(diamond_success.values()):
            logging.info(f"{GREEN}DIAMOND databases ready. Searches will be much faster!{RESET}")
        else:
            logging.warning(f"{RED}DIAMOND not available. Will use BLAST for sequence searches.{RESET}")
            logging.info(
                f"{CYAN}To install DIAMOND: conda install -c bioconda diamond OR see https://github.com/bbuchfink/diamond{RESET}")

        logging.info(f"{GREEN}All databases loaded and ready.{RESET}")
        logging.info(
            f"{GREEN}DaMap helps matching proteins of Desulfuromonas acetexigens between RefSeq NCBI, TRO NCBI and UniProt databases, provides amino acid sequences and predicted cellular localization from DeepLocPro 1.0{RESET}")

    except Exception as e:
        logging.error(f"{RED}Fatal error during database initialization: {e}{RESET}")
        logging.error(f"{RED}Please check the paths in the DATABASE_PATHS configuration section.{RESET}")
        return  # Exit if databases can't be loaded

    # --- 5. Main application loop ---
    while True:
        print(f"\n{CYAN}Select mode:{RESET}")
        print("  1: Single Protein ID Conversion")
        print("  2: Bulk File Conversion (CSV/Excel)")
        print("  3: Create/Update DIAMOND Databases")
        print("  q: Quit")

        mode = input("Enter choice: ").strip()

        if mode == '1':
            single_protein_mode_with_diamond(databases)
        elif mode == '2':
            bulk_file_mode_with_diamond(databases)
        elif mode == '3':
            # Manual DIAMOND database creation/update
            print(f"\n{CYAN}--- DIAMOND Database Management ---{RESET}")
            print("This will create or update DIAMOND databases from your FASTA files.")
            print("DIAMOND databases enable much faster protein sequence searches.")

            confirm = input("Proceed with DIAMOND database creation/update? (y/N): ").strip().lower()
            if confirm in ['y', 'yes']:
                for db_name, fasta_path, diamond_path in diamond_tasks:
                    print(f"\nProcessing {db_name} database...")
                    success = create_diamond_database(fasta_path, diamond_path)
                    if success:
                        print(f"{GREEN}{db_name} DIAMOND database ready.{RESET}")
                    else:
                        print(f"{RED}Failed to create {db_name} DIAMOND database.{RESET}")
            else:
                print("Database update cancelled.")
        elif mode.lower() == 'q':
            print(f"{GREEN}Exiting program.{RESET}")
            break
        else:
            print(f"{RED}Invalid choice. Please enter 1, 2, 3, or q.{RESET}")


if __name__ == "__main__":
    main()