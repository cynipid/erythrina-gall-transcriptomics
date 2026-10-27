#!/usr/bin/env python3
"""Audit retained featureCounts annotation/count records without recounting BAMs.

Python standard library only. The records directory must contain the seven files
in EGW_featureCounts_records_comment8.tar.gz, preserving their relative paths.
"""
import argparse
import csv
import hashlib
import json
import re
import shlex
from collections import Counter, defaultdict
from pathlib import Path


def union_length(features):
    groups = defaultdict(list)
    for chrom, start, end, strand in features:
        groups[(chrom, strand)].append((start, end))
    total = 0
    for spans in groups.values():
        merged = []
        for start, end in sorted(spans):
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        total += sum(end - start + 1 for start, end in merged)
    return total


def run(records, output):
    csv.field_size_limit(10_000_000)
    names = {
        'annotation': '01_reference/complete.genomic.gff',
        'stdout': '03_truba_scripts_logs/logs/featurecounts_strand_test_5981196.out',
        'stderr': '03_truba_scripts_logs/logs/featurecounts_strand_test_5981196.err',
        'raw_counts': '06_counts/gene_counts_s2.txt',
        'count_summary': '06_counts/gene_counts_s2.txt.summary',
        'final_matrix': '06_counts/gene_counts_matrix.tsv',
        'sample_metadata': '06_counts/sample_metadata.tsv',
    }
    for rel in names.values():
        if not (records / rel).is_file():
            raise FileNotFoundError(records / rel)
    groups = defaultdict(list)
    exon_rows = 0
    missing_or_ambiguous_tags = []
    with (records / names['annotation']).open() as fh:
        for line_no, line in enumerate(fh, 1):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'GFF line {line_no}: not nine fields')
            if fields[2] != 'exon':
                continue
            exon_rows += 1
            tags = [v for term in fields[8].split(';')
                    for k, sep, v in [term.partition('=')]
                    if sep and k == 'locus_tag' and v]
            if len(tags) != 1:
                missing_or_ambiguous_tags.append(line_no)
                continue
            groups[tags[0]].append((fields[0], int(fields[3]), int(fields[4]), fields[6]))

    raw = {}
    duplicates = []
    with (records / names['raw_counts']).open() as fh:
        command_header = fh.readline().rstrip('\n')
        reader = csv.DictReader(fh, delimiter='\t')
        raw_columns = reader.fieldnames[6:]
        original_ids = [Path(col).parent.name for col in raw_columns]
        for row in reader:
            gene = row['Geneid']
            if gene in raw:
                duplicates.append(gene)
            coordinates = list(zip(row['Chr'].split(';'),
                                   map(int, row['Start'].split(';')),
                                   map(int, row['End'].split(';')),
                                   row['Strand'].split(';'), strict=True))
            raw[gene] = {'features': coordinates, 'length': int(row['Length']),
                         'counts': tuple(int(row[col]) for col in raw_columns)}
    tokens = shlex.split(command_header.split('Command:', 1)[1])
    requested_options = {'-F': 'GFF', '-t': 'exon', '-g': 'locus_tag', '-s': '2'}
    option_checks = {key: tokens[tokens.index(key) + 1] == value
                     for key, value in requested_options.items()}

    matrix = {}
    matrix_duplicates = []
    with (records / names['final_matrix']).open() as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        matrix_columns = reader.fieldnames[1:]
        for row in reader:
            gene = row['Geneid']
            if gene in matrix:
                matrix_duplicates.append(gene)
            matrix[gene] = tuple(int(row[s]) for s in matrix_columns)
    with (records / names['sample_metadata']).open() as fh:
        metadata = list(csv.DictReader(fh, delimiter='\t'))
    with (records / names['count_summary']).open() as fh:
        reader = csv.reader(fh, delimiter='\t')
        summary_header = next(reader)
        summary = {row[0]: tuple(map(int, row[1:])) for row in reader}

    errors = defaultdict(list)
    gene_records = []
    raw_sums = [0] * len(original_ids)
    mismatch_cells = 0
    for gene, data in raw.items():
        if gene not in groups:
            errors['raw_ids_without_exons'].append(gene)
            continue
        coordinate_match = Counter(groups[gene]) == Counter(data['features'])
        length = union_length(groups[gene])
        length_match = length == data['length']
        count_match = gene in matrix and data['counts'] == matrix[gene]
        if not coordinate_match:
            errors['coordinate_or_strand_mismatch'].append(gene)
        if not length_match:
            errors['union_length_mismatch'].append(gene)
        if not count_match:
            errors['matrix_mismatch_genes'].append(gene)
            mismatch_cells += sum(a != b for a, b in zip(data['counts'], matrix.get(gene, (None,) * len(original_ids))))
        raw_sums = [a + b for a, b in zip(raw_sums, data['counts'])]
        gene_records.append({'locus_tag': gene, 'exon_records': len(groups[gene]),
                             'coordinate_and_strand_multiset_matches': coordinate_match,
                             'union_exon_length': length, 'output_length': data['length'],
                             'union_exon_length_matches': length_match,
                             'all_sample_counts_match_final_matrix': count_match})
    errors['annotation_ids_without_raw_counts'] = sorted(set(groups) - set(raw))
    errors['matrix_ids_without_raw_counts'] = sorted(set(matrix) - set(raw))
    errors['raw_ids_without_matrix'] = sorted(set(raw) - set(matrix))
    stderr = (records / names['stderr']).read_text()
    warning = 'Warning: Unknown annotation format: GFF. GTF format is used.'
    reverse = re.split(r'Output file\s*:\s*gene_counts_s2\.txt', stderr, maxsplit=1)[1]
    log_features = int(re.search(r'Features\s*:\s*(\d+)', reverse).group(1))
    log_groups = int(re.search(r'Meta-features\s*:\s*(\d+)', reverse).group(1))
    log_totals = tuple(map(int, re.findall(r'Total alignments\s*:\s*(\d+)', reverse)))
    log_assigned = tuple(map(int, re.findall(r'Successfully assigned alignments\s*:\s*(\d+)', reverse)))
    summary_totals = tuple(sum(v[i] for v in summary.values()) for i in range(len(original_ids)))
    feature_checks = {
        'every_exon_has_one_nonempty_locus_tag': not missing_or_ambiguous_tags,
        'no_duplicate_raw_gene_ids': not duplicates,
        'no_duplicate_matrix_gene_ids': not matrix_duplicates,
        'all_annotation_ids_and_raw_count_ids_agree': set(groups) == set(raw),
        'all_matrix_ids_and_raw_count_ids_agree': set(matrix) == set(raw),
        'all_exon_coordinate_and_strand_multisets_agree': not errors['coordinate_or_strand_mismatch'],
        'all_union_exon_lengths_agree': not errors['union_length_mismatch'],
        'all_raw_and_final_count_cells_agree': mismatch_cells == 0 and not errors['matrix_mismatch_genes'],
        'raw_and_matrix_sample_order_agree': original_ids == matrix_columns,
        'metadata_sample_order_agrees': [r['sample'] for r in metadata] == original_ids,
        'log_feature_total_agrees': log_features == exon_rows,
        'log_group_total_agrees': log_groups == len(groups),
        'column_sums_match_summary_assigned': tuple(raw_sums) == summary['Assigned'],
        'column_sums_match_log_assigned': tuple(raw_sums) == log_assigned,
        'summary_totals_match_log_total_alignments': summary_totals == log_totals,
        'warning_retained_for_three_strand_runs': stderr.count(warning) == 3,
        'software_version_recorded': 'v2.1.1' in stderr and 'featureCounts v2.1.1' in command_header,
        'executed_options_match_reported_options': all(option_checks.values()),
    }
    display = dict(zip(['L1', 'L3', 'L4', 'P3', 'P4', 'P5'], ['L1', 'L2', 'L3', 'P1', 'P2', 'P3']))
    samples = [{'original_id': sample, 'manuscript_label': display[sample],
                'assigned_count': raw_sums[i], 'summary_total_alignments': summary_totals[i],
                'assignment_rate': raw_sums[i] / summary_totals[i]}
               for i, sample in enumerate(original_ids)]
    hashes = {rel: hashlib.sha256((records / rel).read_bytes()).hexdigest() for rel in names.values()}
    result = {'software_version': 'featureCounts 2.1.1', 'exon_records': exon_rows,
              'locus_tag_groups': len(groups), 'raw_count_rows': len(raw), 'final_matrix_rows': len(matrix),
              'sample_columns': len(original_ids), 'count_cells_compared': len(raw) * len(original_ids),
              'mismatched_count_cells': mismatch_cells, 'checks': feature_checks,
              'all_checks_passed': all(feature_checks.values()), 'samples': samples,
              'warning': warning, 'command_header': command_header, 'source_sha256': hashes,
              'missing_or_ambiguous_locus_tag_lines': missing_or_ambiguous_tags,
              'duplicate_raw_ids': duplicates, 'duplicate_matrix_ids': matrix_duplicates,
              'gene_mismatches': dict(errors),
              'scope': 'Existing annotation/count records only. No BAM recount or downstream analysis rerun.'}
    output.mkdir(parents=True, exist_ok=True)
    (output / 'featureCounts_audit_summary.json').write_text(json.dumps(result, indent=2) + '\n')
    with (output / 'gene_annotation_audit.jsonl').open('w') as fh:
        for row in gene_records:
            fh.write(json.dumps(row) + '\n')
    (output / 'featureCounts_exact_command.txt').write_text(command_header + '\n\n' + shlex.join(tokens) + '\n')
    (output / 'SOURCE_SHA256SUMS.txt').write_text('\n'.join(f'{h}  {rel}' for rel, h in hashes.items()) + '\n')
    report = [
        'featureCounts annotation and matrix audit for review comment 8',
        f'Exon records: {exon_rows:,}', f'Unique locus_tag groups: {len(groups):,}',
        f'Count cells compared: {len(raw) * len(original_ids):,}',
        f'Count-cell mismatches: {mismatch_cells}',
        'All recorded checks passed: ' + str(result['all_checks_passed']),
        '', *[f'{k}: {v}' for k, v in feature_checks.items()], '',
        'The historical warning and exact executed command are retained.',
        'Annotation comparison includes repeated exon records, chromosomes, coordinates, strands and union lengths.',
        'Assignment-rate denominators are the featureCounts-reported Total alignments.',
        'No BAM recount, DESeq2 fit, enrichment analysis or candidate selection was rerun.',
    ]
    (output / 'featureCounts_audit_report.txt').write_text('\n'.join(report) + '\n')
    print(json.dumps({k: result[k] for k in ['exon_records', 'locus_tag_groups', 'count_cells_compared', 'mismatched_count_cells', 'all_checks_passed']}, indent=2))
    if not result['all_checks_passed']:
        raise SystemExit('Audit found a discrepancy; inspect the JSON summary before using the results.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.records, args.output)
