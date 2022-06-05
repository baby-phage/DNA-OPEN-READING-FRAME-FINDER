from __future__ import annotations
import streamlit as st
import requests


########################################################################################################################
################################################## BACKEND #############################################################
########################################################################################################################

def fetch_sequence_NCBI(accession_ID: str) -> str | int:
    """
    :param accession_ID: Accession ID of the desired DNA
    :return: DNA sequence in FASTA format
    """

    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={accession_ID}&db=nuccore&report=fasta&extrafeat=null&conwithfeat=on&hide-cdd=on&retmode=html&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000"
    Sequence = requests.get(url)
    if not Sequence.text.startswith(">"):
        return -1
    else:
        return Sequence.text


def dna_seq_validity_checker(dna: str) -> int:
    """
    :param dna: DNA Sequence that needs to be verified
    :return: -1 if DNA sequence is not valid, else +1
    """

    for nt in dna.upper():
        if nt not in "ATGC":
            return -1
    else:
        return 1


def FASTA_Parser(FASTA_seq: str) -> int | str:
    """
    :param FASTA_seq: DNA sequence in FASTA format
    :return: The DNA sequence in one line
    """

    if not FASTA_seq.startswith(">"):
        return -1
    else:
        FASTA_seq = FASTA_seq.split("\n")
        if len(FASTA_seq) < 2:
            return -1
        else:
            dna_seq = "".join(FASTA_seq[1:])

            return dna_seq


def complement(dna: str) -> str:
    """
    :param dna: DNA sequence
    :return: reverse complement of DNA sequence
    """

    dna = dna.upper()
    c_dna = ""
    for nucleotide in dna:
        if nucleotide == "A":
            c_dna += "T"
        elif nucleotide == "T":
            c_dna += "A"
        elif nucleotide == "G":
            c_dna += "C"
        elif nucleotide == "C":
            c_dna += "G"
    return c_dna


def transcribe(dna: str, strand: str) -> str:
    """
    :param dna: DNA sequence
    :param strand : Strand type [ Coding strand / Template strand]
    :return: transcribed mRNA sequence
    """

    dna = dna.upper()

    if strand == "C":
        m_rna = dna.replace("T", "U")

    elif strand == "T":
        m_rna = complement(dna).replace("T", "U")

    return m_rna


def translate(m_rna: str) -> str:
    """
    :param m_rna: mRNA sequence
    :return: translated peptide sequence
    """

    codon_table = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
                   "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                   "UAU": "Y", "UAC": "Y", "UAA": "", "UAG": "",
                   "UGU": "C", "UGC": "C", "UGA": "", "UGG": "W",
                   "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
                   "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                   "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                   "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                   "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
                   "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                   "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                   "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                   "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
                   "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                   "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                   "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }

    rna_codons = []
    for nucleotide in range(0, len(m_rna), 3):
        rna_codons.append(m_rna[nucleotide:nucleotide + 3])

    peptide_seq = ""
    for codon in rna_codons:
        peptide_seq = peptide_seq + codon_table[codon]

    return peptide_seq


def orf_span_finder(dna: str, strand: str) -> list[(int, int)]:
    """
    :param dna: DNA sequence
    :param strand : Strand type [ Coding strand / Template strand]
    :return: A list of tuple(s) containing start codon position and stop codon position of ORF(s)
    """

    if strand == "C":
        start_codon_dna = "ATG"
        stop_codon_dna_lst = ["TAA", "TAG", "TGA"]
    elif strand == "T":
        start_codon_dna = "TAC"
        stop_codon_dna_lst = ["ATT", "ATC", "ACT"]

    start_codon_pos_lst = []
    stop_codon_pos_lst = []

    for codon_pos in range(len(dna) - 2):
        codon_dna = dna[codon_pos:(codon_pos + 3)]

        if codon_dna == start_codon_dna:
            start_codon_pos_lst.append(codon_pos)

        elif codon_dna in stop_codon_dna_lst:
            stop_codon_pos_lst.append(codon_pos)

    orf_span_lst = []
    viability_chk = []
    for start_pos in start_codon_pos_lst:
        for stop_pos in stop_codon_pos_lst:
            if start_pos not in viability_chk:
                if stop_pos > start_pos and (stop_pos - start_pos) % 3 == 0:
                    orf_span_lst.append((start_pos, stop_pos))
                    viability_chk.append(start_pos)

    return orf_span_lst


def orf_seq_generator(dna, span_lst: list[(int, int)]) -> list[str]:
    """
    :param dna: DNA sequence
    :param span_lst: list of tuple(s) containing start codon position and stop codon position of ORF(s)
    :return: A list containing DNA sequence(s) of ORF(s)
    """

    orf_seq = []
    for start, stop in span_lst:
        orf_seq.append(dna[start:stop + 3])

    return orf_seq


def peptide_mw(peptide: str) -> float:
    """
    :param peptide: Amino acid sequence of a peptide
    :return: The peptide's molecular weight
    """

    AA_weights = {'A': 71.04, 'C': 103.01, 'D': 115.03, 'E': 129.04, 'F': 147.07,
                  'G': 57.02, 'H': 137.06, 'I': 113.08, 'K': 128.09, 'L': 113.08,
                  'M': 131.04, 'N': 114.04, 'P': 97.05, 'Q': 128.06, 'R': 156.10,
                  'S': 87.03, 'T': 101.05, 'V': 99.07, 'W': 186.08, 'Y': 163.06}

    peptide_weight = 0
    for AA in peptide:
        peptide_weight += AA_weights[AA]

    return round(peptide_weight, 2)


def preview_seq(seq: str) -> str:
    """
    :param seq: sequence that needs to be previewed in website
    :return: formatted string
    """
    output = []
    n = 75
    for i in range(0, len(seq), n):
        output.append(f"{i + 1}{(5 - len(str(i))) * ' '}{seq[i:i + n]}")
    return "\n".join(output)


def preview_ORF(orf_no, dna_orf, start, stop, strand) -> str:
    """
    :param orf_no: The ORF number
    :param dna_orf: DNA sequence of an ORF
    :param start: ORF start index
    :param stop: ORF stop index
    :param strand: strand type -> Template strand / Coding strand
    :return: formatted string
    """

    if strand == "T":
        dna = (f"Coding   DNA → [{start}] {complement(dna_orf)} [{stop}]\n"
               f"Template DNA → [{start}] {dna_orf} [{stop}]\n")

        rna = transcribe(dna_orf, strand="T")

    elif strand == "C":
        dna = (f"Template DNA → [{start}] {complement(dna_orf)} [{stop}]\n"
               f"Coding   DNA → [{start}] {dna_orf} [{stop}]\n")

        rna = transcribe(dna_orf, strand="C")

    pepide = translate(rna)

    preview = (
        f"| ORF - {orf_no} |\n"
        f"{dna}"
        "\n"
        f"m-RNA        → [{start}] {rna} [{stop}]\n"
        "\n"
        f"PEPTIDE      → [{start}] {pepide} [{stop}]\n"
        f"               Molecular weight = {peptide_mw(pepide)}\n"
        f"{'_' * 85}"
    )

    return preview


########################################################################################################################
################################################## FRONTEND ############################################################
########################################################################################################################

st.set_page_config(
   page_title="🧬 ORF FINDER",
   page_icon="🌐")

Title = st.container()
Intro = st.container()
Input_box = st.form("input")
DNA_prev_box = st.container()
ORF_prev_box = st.container()

st.markdown("<style>.main {background-color: #FFFFFF;color:black;} </style>", unsafe_allow_html=True)

# TITLE
with Title:
    Title.markdown(
        "<h1 style='font-family : century gothic;text-align: center; color: Black;'><u>DNA OPEN READING FRAME FINDER<u></h1>",
        unsafe_allow_html=True)

# INPUT BOX
with Input_box:
    Input_box.markdown(
        "<h4 style='font-family : century gothic;text-align: center; color: Black;'><u>ENTER YOUR SEQUENCE or ACCESSION ID<u></h4>",
        unsafe_allow_html=True)
    input = Input_box.text_area(
        label="Enter NCBI Accession ID of your DNA sequence or Paste the FASTA sequence in the box",
        value="BC005255.1")
    input_type = Input_box.radio(label="| Select Your Input Type |", options=["ACCESSION ID", "FASTA"])
    strand_type = Input_box.radio(label="| Select Strand Type |", options=["Coding strand", "Template strand"])
    input_submit = st.form_submit_button("Submit")

    if input_submit:
        # checking input type
        if input_type == "FASTA":
            chk = FASTA_Parser(input.strip())

            # checking validity of FASTA sequence
            if chk == -1:
                Input_box.error("Please, Enter a Valid FASTA Sequence.")
            else:
                dna_seq = chk.upper()
                seq_validity = dna_seq_validity_checker(dna_seq)    # checking validity of DNA sequence
                if seq_validity == -1:
                    chk = -1
                    Input_box.error("Please, Enter a valid DNA sequence.")
        else:
            chk = fetch_sequence_NCBI(input.strip())
            # checking validity of ACCESSION ID
            if chk == -1:
                Input_box.error("Please, Enter a Valid Accession ID")
            else:
                dna_seq = FASTA_Parser(chk).upper()
                seq_validity = dna_seq_validity_checker(dna_seq)    # checking validity of DNA sequence
                if seq_validity == -1:
                    chk = -1
                    Input_box.error("Please, Enter a valid DNA sequence.")

        if chk != -1:
            
            # generating complementary DBA sequence
            cdna_seq = complement(dna_seq)

            DNA_prev_box.markdown(
                "<h4 style='font-family : century gothic;text-align: center; color: Black;'><u> DNA <u></h4>",
                unsafe_allow_html=True)
            
            # previewing DNA strands
            if strand_type == "Coding strand":
                DNA_prev_box.write("<p style='text-align: center; text-decoration : underline'>5'-CODING STRAND-3'</p>",
                                   unsafe_allow_html=True)
                DNA_prev_box.text(preview_seq(dna_seq))
                DNA_prev_box.write("<p style='text-align: center; text-decoration : underline'>3'-TEMPLATE "
                                   "STRAND-5'</p>",
                                   unsafe_allow_html=True)
                DNA_prev_box.text(preview_seq(cdna_seq))

            elif strand_type == "Template strand":
                DNA_prev_box.write("<p style='text-align: center; text-decoration : underline'>5'-CODING "
                                   "STRAND-3'</p>",
                                   unsafe_allow_html=True)
                DNA_prev_box.text(preview_seq(cdna_seq))
                DNA_prev_box.write(
                    "<p style='text-align: center; text-decoration : underline'>3'-TEMPLATE STRAND-5'</p>",
                    unsafe_allow_html=True)
                DNA_prev_box.text(preview_seq(dna_seq))

            ORF_prev_box.markdown(
                "<h4 style='font-family : century gothic;text-align: center; color: Black;'><u>OPEN READING FRAME(S)<u></h4>",
                unsafe_allow_html=True)
            
            # generating and previewing ORF
            if strand_type == "Template strand":

                template_strand_seq = dna_seq
                coding_orf_span = orf_span_finder(template_strand_seq, strand="T")
                coding_orf = orf_seq_generator(dna=template_strand_seq, span_lst=coding_orf_span)

                for i in range(len(coding_orf)):
                    orf_dna_seq = coding_orf[i]

                    orf_span_start = coding_orf_span[i][0] + 1
                    orf_span_stop = coding_orf_span[i][1] + 3

                    ORF_prev_box.text(preview_ORF(orf_no=(i + 1),
                                                  dna_orf=orf_dna_seq,
                                                  start=orf_span_start,
                                                  stop=orf_span_stop,
                                                  strand="T"))

            elif strand_type == "Coding strand":

                coding_strand_seq = dna_seq
                coding_orf_span = orf_span_finder(coding_strand_seq, strand="C")
                coding_orf = orf_seq_generator(dna=coding_strand_seq, span_lst=coding_orf_span)

                for i in range(len(coding_orf)):
                    orf_dna_seq = coding_orf[i]

                    orf_span_start = coding_orf_span[i][0] + 1
                    orf_span_stop = coding_orf_span[i][1] + 3

                    ORF_prev_box.text(preview_ORF(orf_no=(i + 1),
                                                  dna_orf=orf_dna_seq,
                                                  start=orf_span_start,
                                                  stop=orf_span_stop,
                                                  strand="C"))
