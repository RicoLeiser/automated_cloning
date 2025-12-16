import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import io
import tempfile
import os
import zipfile
import re

# Hardcoded sequences
BEGIN_SEQ = "CAGCGGCCGCCAAAAAACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTATGCTAAGAAGTTTACTAATACGACTCACTATAGGGATAAT"
END_SEQ = "ATTATCCCTATAGTGAGTCGTATTAGAAAGTTGTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTCTCGAGTA"

# Annotations for BEGIN sequence (1-indexed to 0-indexed)
BEGIN_ANNOTATIONS = [
    ("TT7-antisense", 10, 57),  # base pairs 11 to 57
    ("PT7-sense", 67, 87)       # base pairs 68 to 86
]

# Annotations for END sequence (1-indexed to 0-indexed)
END_ANNOTATIONS = [
    ("PT7-antisense", 5, 25),   # base 1 to 25
    ("TT7-sense", 33, 79)       # base 34 to 80
]


def clean_fasta_content(file_content):
    """Clean FASTA content by removing non-ASCII characters and formatting properly."""
    cleaned_lines = []
    
    for line in file_content.split('\n'):
        if line.strip().startswith(">"):
            # Clean header: remove non-ASCII characters
            clean_header = re.sub(r"[^\x00-\x7F]+", "", line)
            cleaned_lines.append(clean_header.strip())
        else:
            # Clean sequence: keep only letters, convert to uppercase
            clean_seq = re.sub(r"[^A-Za-z]", "", line)
            if clean_seq:  # Only add non-empty lines
                cleaned_lines.append(clean_seq.upper())
    
    return '\n'.join(cleaned_lines)

st.title("🧬 Cloning Manager - Sequence Assembly")
st.write("Assemble sequences with hardcoded begin/end sequences and custom inserts")

# Initialize session state variables
if 'sequence_names' not in st.session_state:
    st.session_state.sequence_names = {}
if 'generated_files' not in st.session_state:
    st.session_state.generated_files = False

# Step 1: Species input
st.header("Step 1: Define Species")
species = st.text_input("Enter species name:", placeholder="e.g., Homo sapiens")

# Step 2: File upload
st.header("Step 2: Upload Sequence File")
uploaded_file = st.file_uploader("Upload FASTA or FNA file", type=["fasta", "fna", "fa"])

if uploaded_file is not None and species:
    try:
        # Try to decode with different encodings
        try:
            file_content = uploaded_file.getvalue().decode("utf-8", errors="ignore")
        except Exception:
            file_content = uploaded_file.getvalue().decode("latin-1", errors="ignore")
        
        # Clean the FASTA content
        cleaned_content = clean_fasta_content(file_content)
        
        # Parse the cleaned file
        sequences = list(SeqIO.parse(io.StringIO(cleaned_content), "fasta"))
        
        if len(sequences) == 0:
            st.error("❌ No valid sequences found in the file. Please check your FASTA format.")
        else:
            st.success(f"✅ Found {len(sequences)} sequence(s) in the file")
            
            # Step 3: Name assignment
            st.header("Step 3: Assign Names to Sequences")
            
            for i, seq_record in enumerate(sequences):
                original_id = seq_record.id
                st.subheader(f"Sequence {i+1}")
                st.write(f"**Original Header:** {seq_record.description}")
                
                # Convert sequence to string safely
                try:
                    # Get the sequence as a string
                    seq_str = str(seq_record.seq)
                    # Clean the sequence - keep only valid nucleotides
                    seq_str = ''.join([c.upper() for c in seq_str if c.upper() in 'ATCGRYSWKMBDHVN'])
                    seq_length = len(seq_str)
                    
                    st.write(f"**Length:** {seq_length} bp")
                    st.write(f"**Preview:** {seq_str[:60]}..." if len(seq_str) > 60 else f"**Sequence:** {seq_str}")
                    
                except Exception as e:
                    st.error(f"Error processing sequence {i+1}: {str(e)}")
                    continue
                
                # Input field for custom name
                custom_name = st.text_input(
                    f"Assign name for sequence {i+1}:",
                    key=f"name_{i}",
                    placeholder=f"e.g., Clone_{i+1}"
                )
                
                if custom_name:
                    st.session_state.sequence_names[i] = {
                        'custom_name': custom_name,
                        'original_id': original_id,
                        'sequence': seq_str,
                        'length': seq_length
                    }
    
    except Exception as e:
        st.error(f"❌ Error reading file: {str(e)}")
        st.info("Please ensure your file is a valid FASTA or FNA format.")
    
    # Step 4: Generate clones
    if len(st.session_state.sequence_names) == len(sequences):
        st.header("Step 4: Generate Clones")
        
        if st.button("🔬 Generate Clone Files"):
            # Store generated files in session state
            st.session_state.generated_files = True
            
            genbank_records = []
            excel_data = []
            gb_files = {}  # Dictionary to store individual GenBank files
            
            for i, seq_info in st.session_state.sequence_names.items():
                custom_name = seq_info['custom_name']
                insert_seq = seq_info['sequence']
                insert_length = seq_info['length']
                original_id = seq_info['original_id']
                
                # Assemble the complete sequence
                full_sequence = BEGIN_SEQ + insert_seq + END_SEQ
                
                # Create SeqRecord with required annotations
                seq_record = SeqRecord(
                    Seq(full_sequence),
                    id=custom_name,
                    description=f"{custom_name} - assembled clone",
                    annotations={
                        "molecule_type": "DNA",
                        "topology": "linear",
                        "organism": species
                    }
                )
                
                # Add features for BEGIN annotations
                for feat_name, start, end in BEGIN_ANNOTATIONS:
                    feature = SeqFeature(
                        FeatureLocation(start, end),
                        type="misc_feature",
                        qualifiers={"label": feat_name}
                    )
                    seq_record.features.append(feature)
                
                # Add feature for INSERT
                insert_start = len(BEGIN_SEQ)
                insert_end = insert_start + insert_length
                insert_feature = SeqFeature(
                    FeatureLocation(insert_start, insert_end),
                    type="CDS",
                    qualifiers={"label": custom_name, "gene": original_id}
                )
                seq_record.features.append(insert_feature)
                
                # Add features for END annotations (offset by begin + insert length)
                end_offset = len(BEGIN_SEQ) + insert_length
                for feat_name, start, end in END_ANNOTATIONS:
                    feature = SeqFeature(
                        FeatureLocation(start + end_offset, end + end_offset),
                        type="misc_feature",
                        qualifiers={"label": feat_name}
                    )
                    seq_record.features.append(feature)
                
                genbank_records.append(seq_record)
                
                # Create individual GenBank file for this clone
                gb_buffer = io.StringIO()
                SeqIO.write([seq_record], gb_buffer, "genbank")
                gb_files[f"{custom_name}.gb"] = gb_buffer.getvalue()
                
                # Prepare Excel data
                # Column 3: name_goi_species_genename_length
                col3_name = f"{custom_name}_GOI_{species}_{original_id}_{insert_length}"
                excel_data.append({
                    'Clone Name': custom_name,
                    'Full Sequence': full_sequence,
                    'Identifier': col3_name
                })
            
            # Create ZIP file containing all GenBank files
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for filename, content in gb_files.items():
                    zip_file.writestr(filename, content)
            zip_content = zip_buffer.getvalue()
            
            # Create Excel file
            df = pd.DataFrame(excel_data)
            excel_buffer = io.BytesIO()
            df.to_excel(excel_buffer, index=False, engine='openpyxl')
            excel_content = excel_buffer.getvalue()
            
            # Store files in session state for persistent download buttons
            st.session_state.zip_content = zip_content
            st.session_state.excel_content = excel_content
            st.session_state.genbank_records = genbank_records
            st.session_state.excel_df = df
            st.session_state.num_files = len(gb_files)
        
        # Display download buttons if files have been generated
        if st.session_state.get('generated_files', False):
            # Display success and download buttons
            st.success("✅ Clone files generated successfully!")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.download_button(
                    label=f"📥 Download GenBank Files (ZIP - {st.session_state.num_files} files)",
                    data=st.session_state.zip_content,
                    file_name="clones_genbank.zip",
                    mime="application/zip",
                    key="download_zip"
                )
            
            with col2:
                st.download_button(
                    label="📥 Download Excel File",
                    data=st.session_state.excel_content,
                    file_name="clones_summary.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    key="download_excel"
                )
            
            # Display preview
            st.header("Preview")
            st.write(f"**Generated {len(st.session_state.genbank_records)} clone(s)**")
            for i, record in enumerate(st.session_state.genbank_records):
                with st.expander(f"Clone: {record.id}"):
                    st.write(f"**Length:** {len(record.seq)} bp")
                    st.write(f"**Number of features:** {len(record.features)}")
                    st.write("**Features:**")
                    for feat in record.features:
                        label = feat.qualifiers.get('label', ['Unknown'])[0]
                        st.write(f"  - {label}: {feat.location.start}-{feat.location.end}")
            
            st.dataframe(st.session_state.excel_df)

else:
    if not species:
        st.info("👆 Please enter a species name to continue")
    if uploaded_file is None:
        st.info("👆 Please upload a FASTA/FNA file to continue")
