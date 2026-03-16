#!/bin/bash
#
# Pipeline Status Checker
# Shows progress of genome annotations and consolidation readiness
#

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║         BV-BRC Pipeline Status                                 ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Check if logged in
echo "🔐 Authentication Status:"
if p3-whoami &>/dev/null; then
    USER=$(p3-whoami 2>/dev/null | grep -o 'user .*' | awk '{print $2}')
    echo "   ✓ Logged in as: $USER"
else
    echo "   ✗ Not logged in. Run: p3-login <username>"
fi
echo ""

# Count input genomes
echo "📥 Input Genomes:"
input_count=$(find input -name "*.fna" -o -name "*.fasta" | wc -l | tr -d ' ')
if [ "$input_count" -gt 0 ]; then
    echo "   Found $input_count FASTA file(s):"
    find input -name "*.fna" -o -name "*.fasta" | sed 's/^/     • /'
else
    echo "   No FASTA files found in input/"
fi
echo ""

# Check annotated genomes
echo "🧬 Annotated Genomes:"
if [ -d "output/rasttk" ]; then
    annotated=0
    total_features=0
    
    for genome_dir in output/rasttk/*/; do
        if [ -d "$genome_dir" ]; then
            genome=$(basename "$genome_dir")
            rast_file="${genome_dir}rast.tsv"
            
            if [ -f "$rast_file" ]; then
                feature_count=$(tail -n +2 "$rast_file" 2>/dev/null | wc -l | tr -d ' ')
                total_features=$((total_features + feature_count))
                echo "   ✓ $genome ($feature_count features)"
                annotated=$((annotated + 1))
            else
                echo "   ⏳ $genome (in progress or failed)"
            fi
        fi
    done
    
    echo ""
    echo "   Summary: $annotated genome(s) completed, $total_features total features"
else
    echo "   No output directory found yet"
    annotated=0
fi
echo ""

# Check consolidated output
echo "📊 Consolidated Output:"
if [ -f "output/consolidated_all_genomes.tsv" ]; then
    consolidated_lines=$(wc -l < output/consolidated_all_genomes.tsv | tr -d ' ')
    consolidated_features=$((consolidated_lines - 1))
    file_size=$(du -h output/consolidated_all_genomes.tsv | cut -f1)
    echo "   ✓ output/consolidated_all_genomes.tsv exists"
    echo "     • Size: $file_size"
    echo "     • Features: $consolidated_features"
    
    # Check if it's up to date
    newest_rast=$(find output/rasttk -name "rast.tsv" -type f -print0 2>/dev/null | xargs -0 ls -t 2>/dev/null | head -1)
    if [ -n "$newest_rast" ] && [ "$newest_rast" -nt "output/consolidated_all_genomes.tsv" ]; then
        echo "     ⚠️  WARNING: Some rast.tsv files are newer. Consider re-running consolidate_all"
    fi
else
    echo "   ✗ Not created yet"
    if [ "$annotated" -gt 0 ]; then
        echo "     → Ready to run: consolidate_all"
    fi
fi
echo ""

# Database status
echo "💾 SEED Variant Database:"
if [ -d "/database/variant_definitions" ]; then
    json_count=$(find /database/variant_definitions/structured_database/json_files -name "*.json" 2>/dev/null | wc -l | tr -d ' ')
    metadata_file="/database/variant_definitions/local_references/subsystem_metadata.tsv"
    
    echo "   ✓ Variant definitions loaded"
    echo "     • JSON files: $json_count subsystems"
    
    if [ -f "$metadata_file" ]; then
        metadata_count=$(tail -n +2 "$metadata_file" 2>/dev/null | wc -l | tr -d ' ')
        echo "     • Metadata entries: $metadata_count"
    fi
else
    echo "   ✗ Database not found at /database/variant_definitions"
fi
echo ""

# Next steps
echo "📋 Next Steps:"
if ! p3-whoami &>/dev/null; then
    echo "   1. Login to BV-BRC: p3-login <username>"
elif [ "$annotated" -eq 0 ]; then
    echo "   1. Annotate genomes: run_genome <name> <fasta_file> [scientific_name]"
    echo "      Example: run_genome ref_colik12 input/ref_colik12/ref_colik12.fna \"E. coli\""
elif [ ! -f "output/consolidated_all_genomes.tsv" ]; then
    echo "   1. Consolidate results: consolidate_all"
else
    echo "   ✓ Pipeline complete! Results in: output/consolidated_all_genomes.tsv"
    echo "   • View results: head output/consolidated_all_genomes.tsv"
    echo "   • Copy to host: Results are already on your mounted volume"
fi
echo ""

echo "════════════════════════════════════════════════════════════════"
