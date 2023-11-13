import sys, os, gtfparse, optparse, warnings, pandas

warnings.simplefilter(action='ignore', category=FutureWarning)

optParser = optparse.OptionParser( 
   
   usage = "python %prog <in.gtf> <out.gff>",
   
   description=
      "Script to prepare annotation for DEXSeq." +
      "This script takes the splice_lib_events.gtf annotation file" +
      "and outputs a 'flattened' annotation file suitable for" +
      "downstream exonic part counting followed by DEXSeq.",
      
   epilog = 
      "Written by Alexander J Ritter (achoi18@ucsc.edu), Sanford Lab at UCSC.")

(opts, args) = optParser.parse_args()

if len( args ) != 2:
   sys.stderr.write( "Script to prepare annotation for DEXSeq.\n\n" )
   sys.stderr.write( "Usage: python %s <in.gtf> <out.gff>\n\n" % os.path.basename(sys.argv[0]) )
   sys.stderr.write( "This script takes the splice_lib_events.gtf annotation file\n" )
   sys.stderr.write( "and outputs a 'flattened' annotation file suitable for\n" )
   sys.stderr.write( "downstream exonic part counting followed by DEXSeq.\n" )
   sys.exit(1)

# Prepare the annotation file.
gtf_file = args[0]
out_file = args[1]

df = gtfparse.read_gtf(gtf_file)
dfDict = df.groupby('gene_id')[['transcript_id', 'seqname', 'source',
                                'feature', 'start', 'end', 'score',
                                'strand', 'frame']].apply(lambda g: g.values.tolist()).to_dict()
exonicParts = []
for gene, f in dfDict.items():
    seqname, source, score, strand, frame = f[0][1], f[0][2], '.', f[0][7], '.'
    geneID = gene
    exons = sorted([(exon[0], exon[4], exon[5]) for exon in f], key=lambda x: (x[1], x[2]))
    gene_start, gene_end = exons[0][1], exons[-1][2]
    exonicParts.append([seqname, source, 'aggregate_gene', gene_start, gene_end, '.',
                        strand, '.', 'gene_id ' + '"{}"'.format(geneID)])
    included, excluded = geneID + '_included', geneID + '_excluded'
    counter = 1
    for transcript in [included, excluded]:
        info = ['transcripts ', '"{}"'.format(transcript), '; exonic_part_number ',
                    '"{:03d}"'.format(counter), '; gene_id ', '"{}"'.format(geneID)]
        exonicParts.append([seqname, source, 'exonic_part', gene_start, gene_end, '.',
                            strand, '.', ''.join(info)])
        counter += 1

with open(out_file, 'w') as file:
   for x in exonicParts:
      line = [str(el) for el in x]
      file.write('\t'.join(line) + '\n')

file.close()      
