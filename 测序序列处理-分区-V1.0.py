import streamlit as st
import pandas as pd
import re
from io import BytesIO
from Bio import SeqIO
from Bio.SeqIO.AbiIO import AbiIterator

# ===================== 基础工具函数 =====================
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def reverse_complement(dna_seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([complement.get(base, base) for base in reversed(dna_seq)])

def translate_dna_to_protein(dna_seq, start_pos=0):
    protein = []
    for i in range(start_pos, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        protein.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(protein)

def parse_ab1_file(file_obj):
    try:
        record = next(AbiIterator(file_obj))
        dna_seq = str(record.seq).upper()
        return dna_seq
    except:
        return None

def parse_seq_file(file_obj):
    try:
        record = SeqIO.read(file_obj, "fasta")
        dna_seq = str(record.seq).upper()
        return dna_seq
    except:
        return None

def parse_file(file):
    name = file.name
    ext = name.lower().split('.')[-1]
    if ext == 'ab1':
        dna = parse_ab1_file(file)
    elif ext == 'seq':
        dna = parse_seq_file(file)
    else:
        return name, None
    if dna:
        dna = re.sub(r'[^ATCG]', '', dna)
    return name, dna

# ===================== 纳米抗体片段提取 =====================
def extract_nb_fragment(prot_clean):
    starts = ['EVQ', 'QVQ']
    ends = ['TVSS']
    candidates = []
    for start_motif in starts:
        start_idx = prot_clean.find(start_motif)
        if start_idx == -1:
            continue
        for end_motif in ends:
            end_idx = prot_clean.find(end_motif, start_idx)
            if end_idx == -1:
                continue
            fragment = prot_clean[start_idx : end_idx + len(end_motif)]
            if 100 <= len(fragment) <= 150:
                candidates.append(fragment)
    # 优先返回最长的有效片段
    return max(candidates, key=len) if candidates else ""

# ===================== 纳米抗体序列分区 =====================
def nanobody_partition(seq_name, seq):
    seq = seq.strip().upper()
    seq = re.sub(r'[^A-Z]', '', seq)
    original_seq = seq
    total_c = original_seq.count('C')
    error_msg = ""
    
    # FR1 (1-26位)
    if len(seq) < 26:
        return {
            "提取的纳米抗体序列": original_seq,
            "FR1 (1-26位)": "",
            "CDR1 (27-38位)": "",
            "FR2 (39-55位)": "",
            "CDR2 (56-65位)": "",
            "FR3 (66-104位)": "",
            "CDR3 (FR3后-FR4前)": "",
            "FR4 (最后11个)": "",
            "C总数": total_c,
            "合并与原始序列比对结果": "不匹配",
            "错误信息": f"序列长度不足26位（当前{len(seq)}位）"
        }
    fr1 = seq[0:26]

    # W锚点（31-50位找第一个W）
    if len(seq) < 50:
        return {
            "提取的纳米抗体序列": original_seq,
            "FR1 (1-26位)": fr1,
            "CDR1 (27-38位)": "",
            "FR2 (39-55位)": "",
            "CDR2 (56-65位)": "",
            "FR3 (66-104位)": "",
            "CDR3 (FR3后-FR4前)": "",
            "FR4 (最后11个)": "",
            "C总数": total_c,
            "合并与原始序列比对结果": "不匹配",
            "错误信息": f"序列长度不足50位，无法定位W锚点（当前{len(seq)}位）"
        }
    search_w = seq[30:50]
    w_idx = search_w.find('W')
    if w_idx == -1:
        return {
            "提取的纳米抗体序列": original_seq,
            "FR1 (1-26位)": fr1,
            "CDR1 (27-38位)": "",
            "FR2 (39-55位)": "",
            "CDR2 (56-65位)": "",
            "FR3 (66-104位)": "",
            "CDR3 (FR3后-FR4前)": "",
            "FR4 (最后11个)": "",
            "C总数": total_c,
            "合并与原始序列比对结果": "不匹配",
            "错误信息": "31-50位未找到W锚点"
        }
    w_actual_idx = 30 + w_idx
    fr2_start = w_actual_idx - 2
    cdr1_with_dash = seq[26:fr2_start]
    fr2_end = w_actual_idx + 15
    fr2 = seq[fr2_start:fr2_end] if len(seq) >= fr2_end else seq[fr2_start:]

    # C锚点（85-110位找第一个C）
    if len(seq) < 110:
        search_c = seq[84:] if len(seq) > 84 else ""
    else:
        search_c = seq[84:110]
    c_idx = search_c.find('C')
    if c_idx == -1:
        return {
            "提取的纳米抗体序列": original_seq,
            "FR1 (1-26位)": fr1,
            "CDR1 (27-38位)": cdr1_with_dash.replace('-', ''),
            "FR2 (39-55位)": fr2,
            "CDR2 (56-65位)": "",
            "FR3 (66-104位)": "",
            "CDR3 (FR3后-FR4前)": "",
            "FR4 (最后11个)": "",
            "C总数": total_c,
            "合并与原始序列比对结果": "不匹配",
            "错误信息": "85-110位未找到C锚点"
        }
    c_actual_idx = 84 + c_idx
    fr3_start = c_actual_idx - 38
    fr3 = seq[fr3_start:c_actual_idx+1] if fr3_start >= 0 else seq[0:c_actual_idx+1]
    cdr2 = seq[fr2_end:fr3_start] if fr2_end < fr3_start else ""
    fr4 = seq[-11:] if len(seq) >= 11 else seq
    cdr3_start = c_actual_idx + 1
    cdr3 = seq[cdr3_start:-11] if len(seq) > cdr3_start + 11 else ""

    # 合并序列 & 比对
    cdr1_final = cdr1_with_dash.replace('-', '')
    merged_parts = [fr1, cdr1_final, fr2, cdr2, fr3, cdr3, fr4]
    merged_seq = "".join(merged_parts)
    
    if len(merged_seq) != len(original_seq):
        review_result = "不匹配（长度不一致）"
    else:
        mismatch = [i+1 for i in range(len(original_seq)) if original_seq[i] != merged_seq[i]]
        review_result = "匹配" if not mismatch else f"不匹配（差异位置：{mismatch[:10]}...）"

    return {
        "提取的纳米抗体序列": original_seq,
        "FR1 (1-26位)": fr1,
        "CDR1 (27-38位)": cdr1_final,
        "FR2 (39-55位)": fr2,
        "CDR2 (56-65位)": cdr2,
        "FR3 (66-104位)": fr3,
        "CDR3 (FR3后-FR4前)": cdr3,
        "FR4 (最后11个)": fr4,
        "C总数": total_c,
        "合并与原始序列比对结果": review_result,
        "错误信息": ""
    }

# ===================== 主流程 & 界面 =====================
st.set_page_config(page_title="纳米抗体核心结果输出", layout="wide")
st.title("🧬 纳米抗体序列 - 核心结果输出")

# 上传文件
uploaded_files = st.file_uploader(
    "上传 .ab1 或 .seq 文件（可多选）",
    type=["ab1", "seq"],
    accept_multiple_files=True
)

if st.button("开始分析") and uploaded_files:
    final_results = []
    for file in uploaded_files:
        # 1. 解析文件 → 提取DNA → 翻译蛋白 → 提取纳米抗体片段
        name, dna = parse_file(file)
        if not dna:
            final_results.append({
                "文件名": name,
                "提取的纳米抗体序列": "无法提取DNA",
                "FR1 (1-26位)": "",
                "CDR1 (27-38位)": "",
                "FR2 (39-55位)": "",
                "CDR2 (56-65位)": "",
                "FR3 (66-104位)": "",
                "CDR3 (FR3后-FR4前)": "",
                "FR4 (最后11个)": "",
                "C总数": 0,
                "合并与原始序列比对结果": "不匹配"
            })
            continue
        
        # 2. 6框翻译 + 提取最优纳米抗体片段
        best_nb_fragment = ""
        strands = [("正链", dna), ("负链", reverse_complement(dna))]
        for strand_name, strand_seq in strands:
            for shift in [0, 1, 2]:
                prot = translate_dna_to_protein(strand_seq, start_pos=shift)
                prot_clean = prot.replace('*', '')
                fragment = extract_nb_fragment(prot_clean)
                if fragment and len(fragment) > len(best_nb_fragment):
                    best_nb_fragment = fragment
        
        if not best_nb_fragment:
            final_results.append({
                "文件名": name,
                "提取的纳米抗体序列": "未找到EVQ/QVQ+TVSS片段",
                "FR1 (1-26位)": "",
                "CDR1 (27-38位)": "",
                "FR2 (39-55位)": "",
                "CDR2 (56-65位)": "",
                "FR3 (66-104位)": "",
                "CDR3 (FR3后-FR4前)": "",
                "FR4 (最后11个)": "",
                "C总数": 0,
                "合并与原始序列比对结果": "不匹配"
            })
            continue
        
        # 3. 对提取的纳米抗体片段做分区 + 统计C数目 + 比对
        partition_res = nanobody_partition(name, best_nb_fragment)
        final_results.append({
            "文件名": name,
            "提取的纳米抗体序列": partition_res["提取的纳米抗体序列"],
            "FR1 (1-26位)": partition_res["FR1 (1-26位)"],
            "CDR1 (27-38位)": partition_res["CDR1 (27-38位)"],
            "FR2 (39-55位)": partition_res["FR2 (39-55位)"],
            "CDR2 (56-65位)": partition_res["CDR2 (56-65位)"],
            "FR3 (66-104位)": partition_res["FR3 (66-104位)"],
            "CDR3 (FR3后-FR4前)": partition_res["CDR3 (FR3后-FR4前)"],
            "FR4 (最后11个)": partition_res["FR4 (最后11个)"],
            "C总数": partition_res["C总数"],
            "合并与原始序列比对结果": partition_res["合并与原始序列比对结果"]
        })
    
    # 4. 输出核心结果
    df = pd.DataFrame(final_results)
    st.dataframe(df, height=600, use_container_width=True)
    
    # 下载结果
    csv = df.to_csv(index=False, encoding='utf-8-sig')
    st.download_button(
        "下载核心结果CSV",
        data=csv,
        file_name="纳米抗体核心结果.csv",
        mime="text/csv"
    )
