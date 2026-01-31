# %% [markdown]
# # 🧬 ViralSeq-QC: Complete Quality Control Pipeline for Viral Genome Sequences
# 
# [![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
# [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
# [![GitHub](https://img.shields.io/badge/GitHub-ViralSeq--QC-181717?logo=github)](https://github.com/Qasim-Hussain-Code/ViralSeq-QC)
# 
# ---
# 
# ## 📋 Overview
# 
# **ViralSeq-QC** is a zero-dependency, production-ready toolkit for quality control of viral consensus sequences. 
# This notebook demonstrates how to use ViralSeq-QC for rigorous pre-processing of viral genomic data, 
# a critical step before downstream analyses like phylogenetics, variant calling, or outbreak investigation.
# 
# ### Why Quality Control Matters in Viral Genomics
# 
# | Problem | Impact | ViralSeq-QC Solution |
# |---------|--------|---------------------|
# | High N-content | False variants, poor alignments | Filters sequences above threshold |
# | Short fragments | Incomplete genome coverage | Enforces minimum length |
# | Low complexity | Alignment artifacts | Detects repetitive regions |
# | Homopolymer runs | Sequencing errors (esp. Nanopore) | Flags suspicious runs |
# | Terminal Ns | Assembly edge artifacts | Quantifies and filters |
# 
# ### What You'll Learn
# 
# 1. **Installation** - Setting up ViralSeq-QC
# 2. **Basic Usage** - CLI and Python API
# 3. **QC Metrics Deep Dive** - Understanding each metric
# 4. **Real-World Analysis** - Processing viral sequences
# 5. **Visualization** - Creating publication-ready QC plots
# 6. **Best Practices** - Recommendations for viral genomics workflows
# 
# ---

# %% [markdown]
# ## 🔧 Part 1: Installation & Setup

# %%
# Install ViralSeq-QC from GitHub
!pip install git+https://github.com/Qasim-Hussain-Code/ViralSeq-QC.git -q

# Verify installation
!viralseq-qc --version

# %%
# Import required libraries
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Import ViralSeq-QC modules
from src.qc import (
    calculate_gc_content,
    calculate_n_content,
    check_length,
    detect_homopolymer_runs,
    calculate_sequence_complexity,
    check_terminal_ns,
    validate_nucleotides,
    get_sequence_metrics,
    is_high_quality
)
from src.input_output import parse_fasta, write_json_report, write_tsv_report

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (12, 6)
plt.rcParams['font.size'] = 11

print("✅ ViralSeq-QC imported successfully!")
print(f"📦 Python version: {sys.version}")

# %% [markdown]
# ## 🧪 Part 2: Creating Sample Viral Sequences
# 
# Let's create a realistic dataset simulating sequences from a viral outbreak investigation.
# These sequences represent different quality scenarios commonly encountered in real-world viral genomics.

# %%
# Create diverse sample sequences representing real-world scenarios
sample_sequences = {
    # High-quality SARS-CoV-2 spike gene fragment
    "SARS-CoV-2_Spike_HQ": (
        "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAATCTTACAACCAGAACTCAATTAC"
        "CCCCTGCATACACTAATTCTTTCACACGTGGTGTTTATTACCCTGACAAAGTTTTCAGATCCTCAGTTTTACA"
        "TTCAACTCAGGACTTGTTCTTACCTTTCTTTTCCAATGTTACTTGGTTCCATGCTATACATGTCTCTGGGACC"
        "AATGGTACTAAGAGGTTTGATAACCCTGTCCTACCATTTAATGATGGTGTTTATTTTGCTTCCACTGAGAAGT"
        "CTAACATAATAAGAGGCTGGATTTTTGGTACTACTTTAGATTCGAAGACCCAGTCCCTACTTATTGTTAATAA"
        "CGCTACTAATGTTGTTATTAAAGTCTGTGAATTTCAATTTTGTAATGATCCATTTTTGGGTGTTTATTACCAC"
        "AAAAACAACAAAAGTTGGATGGAAAGTGAGTTCAGAGTTTATTCTAGTGCGAATAATTGCACTTTTGAATATG"
        "TCTCTCAGCCTTTTCTTATGGACCTTGAAGGAAAACAGGGTAATTTCAAAAATCTTAGGGAATTTGTGTTTAA"
        "GAATATTGATGGTTATTTTAAAATATATTCTAAGCACACGCCTATTAATTTAGTGCGTGATCTCCCTCAGGGT"
    ),
    
    # HIV-1 reverse transcriptase fragment (different virus type)
    "HIV1_RT_Fragment": (
        "CCCATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGC"
        "CATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAA"
        "AATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAA"
        "TTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCG"
        "CAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGA"
    ),
    
    # Influenza A HA segment with moderate N-content (borderline quality)
    "Influenza_HA_ModerateN": (
        "ATGAAGGCAAAACTACTGGTCCTGNNNATTATATGCACGATCTGCACTGGATAGATGCTCTATTGGATCCTCA"
        "AANCAAGGGAATTTCTGCGATCTGAAATCTAAANNNNNNNNNNNNAATGATAAGGAGCTAAATGTCCAAATAT"
        "GCCTAAATANNNNNNNNNTGAGNATCAGCAATGAGGATCTTTTAGACCCATGGATAAACAATGGAGAAGAAGA"
    ),
    
    # Short fragment (likely PCR artifact or assembly failure)
    "Fragment_TooShort": "ATGCGATCGATCGATCGATCGATCGATCGATCGAT",
    
    # High N-content (poor sequencing quality)
    "LowQuality_HighN": (
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        "ATGCGATCGATCNNNNNNNNNNNNNNNNNNGATCGATCGATCGATNNNNNNNNNNNGATCGATCGATCGATCG"
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    ),
    
    # Low complexity (repetitive region - potential assembly artifact)
    "LowComplexity_Repeat": (
        "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
        "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
        "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGA"
    ),
    
    # Long homopolymer run (common Nanopore sequencing artifact)
    "Homopolymer_Artifact": (
        "ATGCGATCGATCGATCGATCGATCGATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGATCGATCGATCGATC"
        "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
    ),
    
    # Terminal Ns (poor read coverage at contig ends)
    "TerminalN_Assembly": (
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA"
        "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT"
        "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    ),
    
    # Dengue virus with valid IUPAC ambiguity codes
    "Dengue_IUPAC": (
        "ATGAATAACCAACGGAAAAAGGRGCTGGTCATGATCTTAGGRGACTTTAAAGATGTCATCTCCTTGGTCAACA"
        "AATTGAAAGGGGAAGCTGTGCTGGAGATGACTGATAAGGCTTTGGATGCTAAATACAATGCGGAGGATAATGG"
        "TTCTATCCTCATCAAATTGGAGATGGGTGAGGCTATGCAGGCTTTGAATGCTCTGGATGGAACTAAAGGCAAC"
    ),
    
    # Zika virus - perfect quality reference
    "Zika_Reference_HQ": (
        "ATGAAAAACCCAAAAAAGAAATCCGGAGGATTCCGGATTGTCAATATGCTAAAACGCGGAGTAGCCCGTGTGA"
        "GCCCCTTTGGGGGCTTGAAGAGGCTGCCAGCCGGACTTCTGCTGGGTCATGGGCCCATCAGGATGGTCTTGGC"
        "GATTCTAGCCTTTTTGAGATTCACGGCAATCAAGCCATCACTGGGTCTCATCAATAGATGGGGTTCAGTGGGG"
        "AAAAAAGAGGCTATGGAAATAATAAAGAAGTTCAAGAAAGATCTGGCTGCCATGCTGAGAATAATCAATGCTA"
        "GGAAGGAGAAGAAGAGACGAGGCGCAGATACTAGTGTCGGAATTGTTGGCCTCCTGCTGACCACAGCTATGGC"
    ),
}

# Save as FASTA file
fasta_path = "/kaggle/working/viral_sequences.fasta"
with open(fasta_path, 'w') as f:
    for name, seq in sample_sequences.items():
        f.write(f">{name}\n{seq}\n")

print(f"✅ Created FASTA file with {len(sample_sequences)} sequences")
print(f"📁 Saved to: {fasta_path}")

# %% [markdown]
# ## 📊 Part 3: Running ViralSeq-QC Analysis
# 
# ### 3.1 Command Line Interface (CLI)
# 
# ViralSeq-QC provides a powerful CLI for batch processing:

# %%
# Run ViralSeq-QC via CLI with detailed metrics
!viralseq-qc \
    --input /kaggle/working/viral_sequences.fasta \
    --output /kaggle/working/qc_report.tsv \
    --detailed \
    --verbose

# %%
# Also generate JSON report for programmatic analysis
!viralseq-qc \
    --input /kaggle/working/viral_sequences.fasta \
    --output /kaggle/working/qc_report.json \
    --json \
    --detailed

# Load and display results
df_results = pd.read_csv("/kaggle/working/qc_report.tsv", sep='\t')
print("\n📋 QC Results Summary:")
print("=" * 60)
display(df_results)

# %% [markdown]
# ### 3.2 Python API Usage
# 
# For more control, use the Python API directly:

# %%
# Analyze each sequence using the Python API
results = []

for name, sequence in sample_sequences.items():
    # Get all metrics in one call
    metrics = get_sequence_metrics(sequence)
    
    # Determine pass/fail status
    passed = is_high_quality(
        sequence,
        min_length=200,
        max_n_content=5.0,
        min_complexity=0.0,
        max_terminal_n=50
    )
    
    # Add metadata
    metrics['sequence_id'] = name
    metrics['status'] = 'PASS' if passed else 'FAIL'
    metrics['virus_type'] = name.split('_')[0]
    
    results.append(metrics)

# Create DataFrame
df_api = pd.DataFrame(results)
df_api = df_api[['sequence_id', 'virus_type', 'status', 'length', 'gc_content', 
                  'n_content', 'complexity', 'homopolymer_count', '5_prime_ns', '3_prime_ns']]

print("\n📊 Python API Analysis Results:")
display(df_api)

# Summary statistics
print("\n📈 Summary Statistics:")
print(f"   Total sequences:  {len(df_api)}")
print(f"   Passed QC:        {(df_api['status'] == 'PASS').sum()} ({(df_api['status'] == 'PASS').mean()*100:.1f}%)")
print(f"   Failed QC:        {(df_api['status'] == 'FAIL').sum()} ({(df_api['status'] == 'FAIL').mean()*100:.1f}%)")

# %% [markdown]
# ## 📈 Part 4: Visualization & Analysis
# 
# ### 4.1 QC Metrics Distribution

# %%
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('ViralSeq-QC Metrics Distribution', fontsize=16, fontweight='bold')

# Color palette for pass/fail
colors = {'PASS': '#2ecc71', 'FAIL': '#e74c3c'}

# 1. Sequence Length Distribution
ax1 = axes[0, 0]
for status in ['PASS', 'FAIL']:
    subset = df_api[df_api['status'] == status]
    ax1.bar(subset['sequence_id'], subset['length'], 
            color=colors[status], label=status, alpha=0.8)
ax1.axhline(y=200, color='red', linestyle='--', label='Min threshold (200 bp)')
ax1.set_ylabel('Length (bp)')
ax1.set_title('Sequence Length')
ax1.legend()
ax1.tick_params(axis='x', rotation=45)

# 2. GC Content Distribution
ax2 = axes[0, 1]
for status in ['PASS', 'FAIL']:
    subset = df_api[df_api['status'] == status]
    ax2.bar(subset['sequence_id'], subset['gc_content'], 
            color=colors[status], alpha=0.8)
ax2.axhline(y=40, color='gray', linestyle=':', alpha=0.5)
ax2.axhline(y=60, color='gray', linestyle=':', alpha=0.5)
ax2.fill_between(ax2.get_xlim(), 40, 60, alpha=0.1, color='blue', label='Typical range')
ax2.set_ylabel('GC Content (%)')
ax2.set_title('GC Content')
ax2.tick_params(axis='x', rotation=45)

# 3. N-Content (Ambiguity)
ax3 = axes[0, 2]
for status in ['PASS', 'FAIL']:
    subset = df_api[df_api['status'] == status]
    ax3.bar(subset['sequence_id'], subset['n_content'], 
            color=colors[status], alpha=0.8)
ax3.axhline(y=5.0, color='red', linestyle='--', label='Max threshold (5%)')
ax3.set_ylabel('N Content (%)')
ax3.set_title('Ambiguity (N-Content)')
ax3.legend()
ax3.tick_params(axis='x', rotation=45)

# 4. Sequence Complexity
ax4 = axes[1, 0]
for status in ['PASS', 'FAIL']:
    subset = df_api[df_api['status'] == status]
    ax4.bar(subset['sequence_id'], subset['complexity'], 
            color=colors[status], alpha=0.8)
ax4.axhline(y=0.5, color='orange', linestyle='--', alpha=0.7, label='Low complexity threshold')
ax4.set_ylabel('Complexity Score')
ax4.set_title('Sequence Complexity (k-mer diversity)')
ax4.legend()
ax4.tick_params(axis='x', rotation=45)

# 5. Homopolymer Count
ax5 = axes[1, 1]
for status in ['PASS', 'FAIL']:
    subset = df_api[df_api['status'] == status]
    ax5.bar(subset['sequence_id'], subset['homopolymer_count'], 
            color=colors[status], alpha=0.8)
ax5.set_ylabel('Homopolymer Runs (≥6 bp)')
ax5.set_title('Homopolymer Detection')
ax5.tick_params(axis='x', rotation=45)

# 6. Terminal Ns
ax6 = axes[1, 2]
width = 0.35
x = np.arange(len(df_api))
bars1 = ax6.bar(x - width/2, df_api['5_prime_ns'], width, label="5' Ns", color='#3498db', alpha=0.8)
bars2 = ax6.bar(x + width/2, df_api['3_prime_ns'], width, label="3' Ns", color='#9b59b6', alpha=0.8)
ax6.axhline(y=50, color='red', linestyle='--', label='Max threshold (50)')
ax6.set_ylabel('Terminal N Count')
ax6.set_title("Terminal N's (Assembly Artifacts)")
ax6.set_xticks(x)
ax6.set_xticklabels(df_api['sequence_id'], rotation=45, ha='right')
ax6.legend()

plt.tight_layout()
plt.savefig('/kaggle/working/qc_metrics_distribution.png', dpi=300, bbox_inches='tight')
plt.show()

print("📊 Figure saved: qc_metrics_distribution.png")

# %% [markdown]
# ### 4.2 Pass/Fail Analysis Heatmap

# %%
# Create a heatmap showing normalized metrics
fig, ax = plt.subplots(figsize=(14, 8))

# Prepare data for heatmap (normalize each metric 0-1)
metrics_cols = ['length', 'gc_content', 'n_content', 'complexity', 'homopolymer_count']
df_heatmap = df_api[metrics_cols].copy()

# Normalize
for col in metrics_cols:
    df_heatmap[col] = (df_heatmap[col] - df_heatmap[col].min()) / (df_heatmap[col].max() - df_heatmap[col].min() + 1e-10)

df_heatmap.index = df_api['sequence_id']

# Create heatmap
sns.heatmap(df_heatmap.T, annot=True, fmt='.2f', cmap='RdYlGn_r', 
            linewidths=0.5, ax=ax, cbar_kws={'label': 'Normalized Score'})

# Add pass/fail indicators
for i, status in enumerate(df_api['status']):
    color = colors[status]
    ax.add_patch(plt.Rectangle((i, -0.15), 1, 0.15, facecolor=color, 
                                 edgecolor='black', clip_on=False))

ax.set_title('Multi-Dimensional QC Metrics Heatmap\n(Green = Good, Red = Problematic)', 
             fontsize=14, fontweight='bold')
ax.set_xlabel('Sequence ID')
ax.set_ylabel('QC Metric')

plt.tight_layout()
plt.savefig('/kaggle/working/qc_heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

# %% [markdown]
# ### 4.3 Failure Reason Analysis

# %%
# Analyze why sequences failed
failure_reasons = []

for _, row in df_api[df_api['status'] == 'FAIL'].iterrows():
    reasons = []
    if row['length'] < 200:
        reasons.append('Short (<200 bp)')
    if row['n_content'] > 5.0:
        reasons.append('High N-content (>5%)')
    if row['complexity'] < 0.1:
        reasons.append('Low complexity')
    if row['5_prime_ns'] > 50 or row['3_prime_ns'] > 50:
        reasons.append('Terminal Ns (>50)')
    
    failure_reasons.append({
        'sequence_id': row['sequence_id'],
        'reasons': ', '.join(reasons) if reasons else 'Multiple factors'
    })

df_failures = pd.DataFrame(failure_reasons)
print("🔍 Failure Analysis:")
display(df_failures)

# Pie chart of failure reasons
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Pass/Fail pie
status_counts = df_api['status'].value_counts()
ax1.pie(status_counts.values, labels=status_counts.index, autopct='%1.1f%%',
        colors=[colors.get(s, 'gray') for s in status_counts.index],
        explode=[0.05 if s == 'FAIL' else 0 for s in status_counts.index],
        shadow=True, startangle=90)
ax1.set_title('QC Pass/Fail Distribution', fontweight='bold')

# Failure reasons (if any failures)
if len(df_failures) > 0:
    all_reasons = []
    for reasons in df_failures['reasons']:
        all_reasons.extend(reasons.split(', '))
    reason_counts = pd.Series(all_reasons).value_counts()
    
    ax2.barh(reason_counts.index, reason_counts.values, color='#e74c3c', alpha=0.8)
    ax2.set_xlabel('Count')
    ax2.set_title('Failure Reasons Breakdown', fontweight='bold')
    for i, v in enumerate(reason_counts.values):
        ax2.text(v + 0.1, i, str(v), va='center', fontweight='bold')

plt.tight_layout()
plt.savefig('/kaggle/working/failure_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

# %% [markdown]
# ## 🔬 Part 5: Advanced Analysis - Homopolymer Deep Dive
# 
# Homopolymer runs are particularly important for Nanopore/PacBio data quality assessment.

# %%
# Detailed homopolymer analysis
print("🔬 Homopolymer Analysis:")
print("=" * 60)

for name, sequence in sample_sequences.items():
    homopolymers = detect_homopolymer_runs(sequence, min_length=6)
    if homopolymers:
        print(f"\n📍 {name}:")
        for hp in homopolymers:
            print(f"   • Base '{hp['base']}' × {hp['length']} at position {hp['start']}-{hp['end']}")

# %% [markdown]
# ## 🛠️ Part 6: Exporting Results for Downstream Analysis

# %%
# Export passing sequences for downstream analysis
passing_sequences = {name: seq for name, seq in sample_sequences.items() 
                     if df_api[df_api['sequence_id'] == name]['status'].values[0] == 'PASS'}

# Save passing sequences
with open('/kaggle/working/passing_sequences.fasta', 'w') as f:
    for name, seq in passing_sequences.items():
        f.write(f">{name}\n{seq}\n")

print(f"✅ Exported {len(passing_sequences)} passing sequences to passing_sequences.fasta")

# Export detailed report as Excel
df_api.to_excel('/kaggle/working/qc_report_detailed.xlsx', index=False)
print("📊 Detailed report saved to qc_report_detailed.xlsx")

# %% [markdown]
# ## 📚 Part 7: Best Practices & Recommendations
# 
# ### Recommended Thresholds by Application
# 
# | Application | Min Length | Max N% | Min Complexity | Notes |
# |-------------|------------|--------|----------------|-------|
# | **Phylogenetics** | 80% genome | 1% | 0.3 | Strict for accurate trees |
# | **Outbreak Investigation** | 50% genome | 5% | 0.1 | Balance coverage vs quality |
# | **Variant Calling** | 200 bp | 2% | 0.2 | Focus on coding regions |
# | **Surveillance** | 1000 bp | 10% | 0.0 | Maximize sequence recovery |
# | **Metagenomics** | 150 bp | 5% | 0.1 | Short reads acceptable |
# 
# ### Virus-Specific Considerations
# 
# - **SARS-CoV-2**: Expect ~30% GC content, watch for spike gene artifacts
# - **HIV**: High mutation rate, more ambiguity tolerable
# - **Influenza**: Segmented genome, check each segment separately
# - **Dengue/Zika**: ~47% GC, low complexity regions common in UTRs

# %% [markdown]
# ## 🎯 Conclusion
# 
# **ViralSeq-QC** provides a comprehensive, zero-dependency solution for viral sequence quality control:
# 
# ### Key Takeaways
# 
# ✅ **Easy Integration** - Works as CLI or Python library  
# ✅ **Comprehensive Metrics** - Length, GC, N-content, complexity, homopolymers, terminal Ns  
# ✅ **Publication-Ready** - Generates detailed reports in TSV/JSON formats  
# ✅ **Customizable** - Adjustable thresholds for different applications  
# ✅ **Visualization-Friendly** - Output works seamlessly with pandas/matplotlib  
# 
# ### Next Steps
# 
# 1. Apply ViralSeq-QC to your own viral sequence datasets
# 2. Integrate into your bioinformatics pipelines (Snakemake, Nextflow)
# 3. Customize thresholds based on your specific research needs
# 4. Contribute to the project on [GitHub](https://github.com/Qasim-Hussain-Code/ViralSeq-QC)
# 
# ---
# 
# ### 📧 Contact & Citation
# 
# **GitHub**: [Qasim-Hussain-Code/ViralSeq-QC](https://github.com/Qasim-Hussain-Code/ViralSeq-QC)  
# **License**: MIT  
# 
# If you use ViralSeq-QC in your research, please cite:
# 
# ```
# Hussain, Q. (2026). ViralSeq-QC: A zero-dependency toolkit for viral genome quality control. 
# GitHub repository: https://github.com/Qasim-Hussain-Code/ViralSeq-QC
# ```

# %%
# Final summary
print("\n" + "="*60)
print("🎉 ANALYSIS COMPLETE!")
print("="*60)
print(f"\n📁 Output Files Generated:")
print(f"   • qc_report.tsv - TSV format report")
print(f"   • qc_report.json - JSON format report")
print(f"   • qc_report_detailed.xlsx - Excel workbook")
print(f"   • passing_sequences.fasta - Filtered high-quality sequences")
print(f"   • qc_metrics_distribution.png - Metrics visualization")
print(f"   • qc_heatmap.png - Multi-dimensional heatmap")
print(f"   • failure_analysis.png - Failure analysis charts")
print("\n🧬 Happy sequencing!")
