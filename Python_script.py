import pandas as pd
import re
import os
from flask import Flask, request, jsonify
from flask_cors import CORS
from tempfile import NamedTemporaryFile
from Bio import SeqIO


# Initialize Flask App
app = Flask(__name__)
CORS(app)


REGULON_STATIC_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'static_dataset_regulonDB.csv')
AMR_GENELIST_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'AMR_genelist_ncbi.txt')



# ID NORMALIZATION

def normalize_id(raw_id):
    # accn|CP183138 → CP183138
    # accn_CP183138 → CP183138
    return re.split(r'[|_]', raw_id)[-1]


# Load AMR Genes

def load_amr_genes(filepath):
    try:
        with open(filepath, 'r') as f:
            return {line.strip().lower() for line in f if line.strip()}
    except FileNotFoundError:
        return set()


AMR_GENE_SET = load_amr_genes(AMR_GENELIST_PATH)


# FASTA Parsing
def parse_genome_fasta(genome_path):
    chromosome_ids, plasmid_ids = [], []

    for record in SeqIO.parse(genome_path, "fasta"):
        header = record.description.lower()
        norm_id = normalize_id(record.id)

        if 'plasmid' in header:
            plasmid_ids.append(norm_id)
        else:
            chromosome_ids.append(norm_id)

    return chromosome_ids, plasmid_ids


# GFF UNIQUE GENE EXTRACTION
def extract_unique_genes_from_gff(gff_path):
    gene_regex = re.compile(r'gene=([^;]+)')
    genes = set()

    with open(gff_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            m = gene_regex.search(parts[8])
            if not m:
                continue

            g = re.match(r'^[A-Za-z]+', m.group(1))
            if g:
                genes.add(g.group(0).lower())

    return sorted(genes)


# GFF gene extraction by replicon
def extract_genes_by_replicon(gff_path, chromosome_ids, plasmid_ids):
    gene_regex = re.compile(r'gene=([^;]+)')
    chromosome_genes, plasmid_genes = set(), set()

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            seqid = normalize_id(parts[0])
            m = gene_regex.search(parts[8])
            if not m:
                continue

            g = re.match(r'^[A-Za-z]+', m.group(1))
            if not g:
                continue

            gene = g.group(0).lower()

            if seqid in chromosome_ids:
                chromosome_genes.add(gene)
            elif seqid in plasmid_ids:
                plasmid_genes.add(gene)

    return sorted(chromosome_genes), sorted(plasmid_genes)


# NETWORK (chromosome-only construction)
def compare_and_extract_network(regulon_path, gff_path, allowed_genes=None):
    try:
        gff_gene_list = extract_unique_genes_from_gff(gff_path)
        gff_genes = set(gff_gene_list)

        if allowed_genes is not None:
            gff_genes = set(allowed_genes)
            gff_gene_list = sorted(gff_genes)

        if not gff_genes:
            return [], [], gff_gene_list, [], 0, 0, len(gff_gene_list), 0, 0, 0

        df = pd.read_csv(regulon_path)[['regulatorName', 'riFunction', 'tuGenes']].dropna().drop_duplicates()

        def clean(x):
            m = re.match(r'^[A-Za-z]+', str(x))
            return m.group(0).lower() if m else None

        df['reg'] = df['regulatorName'].apply(clean)
        df = df[df['reg'].isin(gff_genes)]

        edges, nodes, tugenes, tu_list = [], set(), set(), []

        for _, r in df.iterrows():
            genes = [clean(g) for g in r['tuGenes'].split(';') if clean(g) in gff_genes]
            if not genes:
                continue

            nodes.add(r['reg'])
            edges.append({'source': r['reg'], 'target': genes[0], 'type': r['riFunction']})
            nodes.add(genes[0])
            tugenes.add(genes[0])

            for i in range(len(genes) - 1):
                edges.append({'source': genes[i], 'target': genes[i + 1], 'type': 'tu_link'})
                nodes.add(genes[i + 1])
                tugenes.add(genes[i + 1])

            tu_list.append({
                'Regulator': r['reg'],
                'Regulation_Type': r['riFunction'],
                'Transcription_Unit_Genes': ";".join(genes)
            })

        node_objs = [{'id': n, 'type': 'regulator' if n in df['reg'].values else 'tugene'} for n in nodes]

        return (
            node_objs, edges, gff_gene_list, tu_list,
            df['reg'].nunique(), df['tuGenes'].nunique(),
            len(gff_gene_list), len(tugenes),
            len(edges), len(node_objs)
        )

    except Exception as e:
        print(e)
        return [], [], [], [], 0, 0, 0, 0, 0, 0


# Flask Route
@app.route('/process_files', methods=['POST'])
def process_files():

    gff_tmp = NamedTemporaryFile(delete=False)
    genome_tmp = NamedTemporaryFile(delete=False)

    try:
        request.files['file2'].save(gff_tmp.name)
        request.files['genome'].save(genome_tmp.name)

        gff_tmp.close()
        genome_tmp.close()

        chromosome_ids, plasmid_ids = parse_genome_fasta(genome_tmp.name)

        chromosome_genes, plasmid_genes = extract_genes_by_replicon(
            gff_tmp.name, chromosome_ids, plasmid_ids
        )
        total_gff_gene_count = len(chromosome_genes) + len(plasmid_genes)


        nodes, edges, gene_list, regulator_tu_list, regulator_count, tu_count, \
        gff_gene_count, unique_tugene_nodes_count, total_edges, total_nodes = \
            compare_and_extract_network(
                REGULON_STATIC_PATH,
                gff_tmp.name,
                allowed_genes=chromosome_genes
            )

        plasmid_overlap_genes = list(
            set(plasmid_genes).intersection(n['id'] for n in nodes)
        )

        return jsonify({
            'success': True,
            'nodes': nodes,
            'edges': edges,

            # REQUIRED BY HTML
            'chromosome_gene_list': chromosome_genes,
            'plasmid_gene_list': plasmid_genes,

            'regulator_tu_list': regulator_tu_list,
            'amr_gene_list': [],

            'chromosome_genes': chromosome_genes,
            'plasmid_genes': plasmid_genes,
            'plasmid_overlap_genes': plasmid_overlap_genes,

            'stats': {
                'gff_genes': total_gff_gene_count,
                'regulators': regulator_count,
                'tus': tu_count,
                'tugene_targets': unique_tugene_nodes_count,
                'edges': total_edges,
                'nodes': total_nodes,

        
                'plasmid_count': len(plasmid_ids),                 # No. of plasmids
                'chromosome_gene_count': len(chromosome_genes),   # Genes in chromosome
                'plasmid_gene_count': len(plasmid_genes),         # Genes in plasmid
                'overlap_gene_count': len(plasmid_overlap_genes)  # Overlapping genes
            }
        }), 200

    except Exception as e:
        print("SERVER ERROR:", e)
        return jsonify({'success': False, 'error': str(e)}), 500

    finally:
        try:
            if os.path.exists(gff_tmp.name):
                os.remove(gff_tmp.name)
            if os.path.exists(genome_tmp.name):
                os.remove(genome_tmp.name)
        except PermissionError:
            pass


if __name__ == '__main__':
    app.run(debug=True, port=5000)
