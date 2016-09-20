import os
import subprocess
import datetime
import re
import uuid
import json

from chunkypipes.components import Software, Parameter, Redirect, BasePipeline

FIRST_READS_PAIR = 0
FIRST_CHAR = 0
PERCENT_DUPLICATION = 7
MAPPED_READS_COUNT = 0

JAVA_DEFAULT_HEAP_SIZE = '6'


class Pipeline(BasePipeline):
    def description(self):
        return """RNAseq pipeline used at the University of Chicago."""

    def configure(self):
        return {
            'cutadapt': {
                'path': 'Full path to cutadapt',
                'quality-base': 'Phred quality scale of raw reads [33|64]'
            },
            'fastqc': {
                'path': 'Full path to fastqc'
            },
            'STAR': {
                'path': 'Full path to STAR',
                'threads': 'Number of threads to run STAR',
                'genome-dir': 'Directory containing a STAR genome index'
            },
            'novosort': {
                'path': 'Full path to novosort'
            },
            'samtools': {
                'path': 'Full path to samtools (must be >= version 1.2)'
            },
            'picard': {
                'path': 'Full path to the picard jar [Ex. /path/to/picard.jar]',
                'heap_size': 'Java heap size in GB [Ex. 6]'
            },
            'cufflinks': {
                'path': 'Full path to cufflinks',
                'transcriptome-gtf': 'Annotation GTF file for cufflinks',
                'threads': 'Number of threads to run cufflinks'
            },
            'htseq': {
                'path': 'Full path to htseq-count',
                'transcriptome-gtf': 'Annotation GTF file for HTSeq'
            },
            'bedtools': {
                'path': 'Full path to bedtools'
            },
            'qc': {
                'rRNA-bed': 'Full path to BED6 file containing rRNA transcript regions',
                'genome-fa': 'Full path to genome fasta'
            },
            'RNAseQC': {
                'path': 'Full path to RNAseQC jar [Ex. /path/to/RNAseQC.jar]'
            }
        }

    def add_pipeline_args(self, parser):
        parser.add_argument('--reads', required=True, action='append',
                            help=('Reads to process with this pipeline. Denote paired-end reads with ' +
                                  'a colon (Ex. read1.fastq:read2.fastq). Specify multiple times to ' +
                                  'align multiple libraries (or pairs).'))
        parser.add_argument('--output', required=True,
                            help='Full path to output directory.')
        parser.add_argument('--lib', default=datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S'),
                            help=('Name of the library, prepended to output file names. Defaults to ' +
                                  'a date string (YYYY-MM-DD-hh-mm-ss).'))
        parser.add_argument('--step', default=0,
                            help='Pipeline step to start at.')
        parser.add_argument('--forward-adapter', default='ZZZ',
                            help='Adapter sequence for the forward strand.')
        parser.add_argument('--reverse-adapter', default='ZZZ',
                            help='Adapter sequnce for the reverse strand.')
        parser.add_argument('--is-stranded', action='store_true',
                            help='Provide this argument if library is stranded.')
        parser.add_argument('--cufflinks-lib-type', default='fr-firststrand',
                            choices=['ff-firststrand',
                                     'ff-secondstrand',
                                     'ff-unstranded',
                                     'fr-firststrand',
                                     'fr-secondstrand',
                                     'fr-unstranded',
                                     'transfrags'],
                            help='Library type for cufflinks. Defaults to fr-firststrand.')
        parser.add_argument('--htseq-stranded', default='yes',
                            choices=['yes', 'no', 'reverse'],
                            help='Strandedness for HTSeq. Defaults to yes.')
        return parser

    def count_gzipped_lines(self, filepath):
        zcat = subprocess.Popen(['zcat', filepath], stdout=subprocess.PIPE)
        num_lines = subprocess.check_output(['wc', '-l'], stdin=zcat.stdout)
        return num_lines.strip()

    def run_pipeline(self, pipeline_args, pipeline_config):
        # Instantiate options
        reads = pipeline_args['reads']
        output_dir = pipeline_args['output']
        logs_dir = os.path.join(output_dir, 'logs')
        lib_prefix = pipeline_args['lib']
        step = pipeline_args['step']
        forward_adapter = pipeline_args['forward_adapter']
        reverse_adapter = pipeline_args['reverse_adapter']
        run_is_stranded = pipeline_args['is_stranded']
        cufflinks_lib_type = pipeline_args['cufflinks_lib_type']
        htseq_stranded = pipeline_args['htseq_stranded']

        # Determine if run is paired-end from input
        run_is_paired_end = len(reads[FIRST_READS_PAIR].split(':')) > 1

        # Create output, tmp, and logs directories
        tmp_dir = os.path.join(output_dir, 'tmp')
        subprocess.call(['mkdir', '-p', output_dir, logs_dir, tmp_dir])

        # Keep list of items to delete
        staging_delete = [os.path.join(output_dir, 'tmp')]

        qc_metrics = {
            'total_raw_reads': [],
            'total_trimmed_reads': [],
            'percent_num_reads_mapped_genome': [],
            'percent_num_reads_mapped_transcriptome': [],
            'percent_duplicate_reads': '0',
            'num_reads_multimapped': [],
            'percent_num_reads_rrna': '',
            'viral_rna': []

        }

        synapse_metadata = {
            'Assay': 'RNAseq',
            'Individual_ID': '',
            'Sample_ID': '',
            'File_Name': '',
            'BrodmannArea': '',
            'BrainRegion': '',
            'Hemisphere': '',
            'CellType': 'NA',
            'TissueState': '',
            'RNAIsolationBatch': '',
            'RIN': '',
            'LibraryBatch': '',
            'LibraryPrep': 'stranded, rRNA depletion',
            'LibraryKit': 'Illumina RS-122-2301',
            'ERCC_Added': '',
            'RunType': 'paired-end',
            'ReadLength': '100bp',
            'FlowcellBatch': '',
            'SequencingPlatform': '',
            'TotalReads': '',
            'MappedReads_Primary': '0',
            'MappedReads_Multimapped': '0',
            'rRNARate': '0',
            'Notes': ''
        }

        # Establish Software instances
        cutadapt = Software('Cutadapt', pipeline_config['cutadapt']['path'])
        fastqc = Software('FastQC', pipeline_config['fastqc']['path'])
        star = Software('STAR Two-Pass', pipeline_config['STAR']['path'])
        novosort = Software('Novosort', pipeline_config['novosort']['path'])
        samtools_flagstat = Software('Samtools Flagstat', pipeline_config['samtools']['path'] + ' flagstat')
        samtools_index = Software('Samtools Index', pipeline_config['samtools']['path'] + ' index')
        samtools_faidx = Software('Samtools Faidx', pipeline_config['samtools']['path'] + ' faidx')
        picard_markduplicates = Software('Picard MarkDuplicates',
                                         'java -Xmx{heap_size}g -jar {path} MarkDuplicates'.format(
                                             heap_size=pipeline_config['picard'].get('heap_size',
                                                                                     JAVA_DEFAULT_HEAP_SIZE),
                                             path=pipeline_config['picard']['path']
                                         ))
        picard_create_seq_dict = Software('Picard CreateSequenceDictionary',
                                          'java -Xmx{heap_size}g -jar {path} CreateSequenceDictionary'.format(
                                             heap_size=pipeline_config['picard'].get('heap_size',
                                                                                     JAVA_DEFAULT_HEAP_SIZE),
                                             path=pipeline_config['picard']['path']
                                          ))
        rnaseqc = Software('RNAseQC', 'java -Xmx{heap_size}g -jar {path}'.format(
                                             heap_size=pipeline_config['picard'].get('heap_size',
                                                                                     JAVA_DEFAULT_HEAP_SIZE),
                                             path=pipeline_config['RNAseQC']['path']
                                      ))
        picard_add_read_groups = Software('Picard AddOrReplaceReadGroups',
                                          'java -Xmx{heap_size}g -jar {path} AddOrReplaceReadGroups'.format(
                                             heap_size=pipeline_config['picard'].get('heap_size',
                                                                                     JAVA_DEFAULT_HEAP_SIZE),
                                             path=pipeline_config['picard']['path']
                                          ))
        cufflinks = Software('Cufflinks', pipeline_config['cufflinks']['path'])
        htseq = Software('HTSeq', pipeline_config['htseq']['path'])
        bedtools_coverage = Software('Bedtools Coverage', pipeline_config['bedtools']['path'] + ' coverage')
        bedtools_bamtobed = Software('Bedtools Bamtobed', pipeline_config['bedtools']['path'] + ' bamtobed')

        # Housekeeping
        star_output = []
        novosort_outfile = ''

        # Step 1: Trimming | Cutadapt
        if step <= 1:
            for i, read in enumerate(reads):
                if run_is_paired_end:
                    # Get paired-end reads, construct new filenames
                    read1, read2 = read.split(':')

                    # QC: Get raw fastq read counts
                    qc_metrics['total_raw_reads'].append([
                        str(int(self.count_gzipped_lines(read1))/4),
                        str(int(self.count_gzipped_lines(read2))/4)
                    ])

                    trimmed_read1_filename = os.path.join(output_dir,
                                                          lib_prefix + '_{}_read1.trimmed.fastq.gz'.format(i))
                    trimmed_read2_filename = os.path.join(output_dir,
                                                          lib_prefix + '_{}_read2.trimmed.fastq.gz'.format(i))
                    staging_delete.extend([
                        trimmed_read1_filename,
                        trimmed_read2_filename
                    ])

                    # Run cutadapt
                    cutadapt.run(
                        Parameter('--quality-base={}'.format(pipeline_config['cutadapt']['quality-base'])),
                        Parameter('--minimum-length=5'),
                        Parameter('--output={}'.format(trimmed_read1_filename)),
                        Parameter('--paired-output={}'.format(trimmed_read2_filename)),
                        Parameter('-a', forward_adapter),
                        Parameter('-A', reverse_adapter),
                        Parameter('-q', '30'),
                        Parameter(read1),
                        Parameter(read2),
                        Redirect(stream=Redirect.STDOUT, dest=os.path.join(logs_dir, 'cutadapt.summary'))
                    )

                    # QC: Get trimmed fastq read counts
                    qc_metrics['total_trimmed_reads'].append([
                        str(int(self.count_gzipped_lines(trimmed_read1_filename))/4),
                        str(int(self.count_gzipped_lines(trimmed_read2_filename))/4)
                    ])

                    # Update reads list
                    reads[i] = ':'.join([trimmed_read1_filename, trimmed_read2_filename])
                else:
                    # QC: Get raw fastq read counts
                    qc_metrics['total_raw_reads'].append([
                        str(int(self.count_gzipped_lines(read))/4)
                    ])

                    # Construct new filename
                    trimmed_read_filename = os.path.join(output_dir,
                                                         lib_prefix + '_{}.trimmed.fastq.gz'.format(i))
                    staging_delete.append(trimmed_read_filename)

                    # Run cutadapt
                    cutadapt.run(
                        Parameter('--quality-base={}'.format(pipeline_config['cutadapt']['quality-base'])),
                        Parameter('--minimum-length=5'),
                        Parameter('--output={}'.format(trimmed_read_filename)),
                        Parameter('-a', forward_adapter),
                        Parameter('-q', '30'),
                        Parameter(read),
                        Redirect(stream=Redirect.STDOUT, dest=os.path.join(logs_dir, 'cutadapt.chicago.summary'))
                    )

                    # QC: Get trimmed fastq read counts
                    qc_metrics['total_trimmed_reads'].append([
                        str(int(self.count_gzipped_lines(trimmed_read_filename))/4)
                    ])

                    # Update reads list
                    reads[i] = trimmed_read_filename

        # Step 2: FastQC
        if step <= 2:
            # Make FastQC directory
            fastqc_output_dir = os.path.join(output_dir, 'fastqc')
            subprocess.call(['mkdir', '-p', fastqc_output_dir])

            all_fastqs = []

            if run_is_paired_end:
                for read in reads:
                    all_fastqs.extend(read.split(':'))
            else:
                all_fastqs.extend(reads)

            for fastq in all_fastqs:
                fastqc.run(
                    Parameter('--outdir={}'.format(fastqc_output_dir)),
                    Parameter(fastq)
                )

        # Step 3: Alignment | STAR 2-pass, Alignment Stats | samtools flagstat
        if step <= 3:
            # Set up common STAR parameters
            star_common = [
                Parameter('--runMode', 'alignReads'),
                Parameter('--twopassMode', 'Basic'),
                Parameter('--runThreadN', pipeline_config['STAR']['threads']),
                Parameter('--genomeDir', pipeline_config['STAR']['genome-dir']),
                Parameter('--readFilesCommand', 'zcat'),
                Parameter('--quantMode', 'TranscriptomeSAM', 'GeneCounts'),
                Parameter('--outSAMtype', 'BAM', 'Unsorted'),
                Parameter('--outFilterType', 'BySJout'),
                Parameter('--outFilterMultimapNmax', '20'),
                Parameter('--alignSJoverhangMin', '8'),
                Parameter('--alignSJDBoverhangMin', '1'),
                Parameter('--outFilterMismatchNmax', '2'),
                Parameter('--alignIntronMin', '20'),
                Parameter('--alignIntronMax', '1000000'),
                Parameter('--alignMatesGapMax', '1000000'),
                (
                    Parameter('--outFilterIntronMotifs', 'RemoveNoncanonical') if run_is_stranded
                    else Parameter('--outSAMstrandField', 'intronMotif')
                )
            ]

            # Get STAR output file prefix
            star_outfile_prefix = os.path.join(output_dir,
                                               lib_prefix + ('_' if lib_prefix[-1] != '.' else '') + '{}.')

            # Align each read or read pair
            for i, read in enumerate(reads):
                star_output_bam = star_outfile_prefix.format(i) + 'Aligned.out.bam'
                star_output_transcriptome_bam = star_outfile_prefix.format(i) + 'Aligned.toTranscriptome.out.bam'
                star_output.append(star_output_bam)

                if run_is_paired_end:
                    read1, read2 = read.split(':')

                    star_paired_end = [
                        Parameter('--readFilesIn', read1, read2),
                        Parameter('--outFileNamePrefix', star_outfile_prefix.format(i))
                    ]

                    star.run(*(star_common + star_paired_end))
                else:
                    star_single_end = [
                        Parameter('--readFilesIn', read),
                        Parameter('--outFileNamePrefix', star_outfile_prefix.format(i))
                    ]

                    star.run(*(star_common + star_single_end))

                # Get flagstats for both alignments
                samtools_flagstat.run(
                    Parameter(star_output_bam),
                    Redirect(stream=Redirect.STDOUT, dest=star_output_bam + '.flagstat')
                )
                samtools_flagstat.run(
                    Parameter(star_output_transcriptome_bam),
                    Redirect(stream=Redirect.STDOUT, dest=star_output_transcriptome_bam + '.flagstat')
                )

                # QC: Get number of mapped reads to the genome from this BAM
                try:
                    with open(star_output_bam + '.flagstat') as flagstats:
                        flagstats_contents = flagstats.read()

                        # Pull out mapped reads
                        target_line = re.search(r'(\d+) \+ \d+ mapped \(([0-9\.]+)%', flagstats_contents)
                        if target_line is not None:
                            num_mapped = int(target_line.group(1))
                            qc_metrics['percent_num_reads_mapped_genome'].append(
                                [str(num_mapped/2), '{}%'.format(target_line.group(2))]
                            )

                            num_secondary = int(re.search(r'(\d+) \+ \d+ secondary', flagstats_contents)
                                                .group(1)
                                                )
                            num_supplementary = int(re.search(r'(\d+) \+ \d+ supplementary', flagstats_contents)
                                                    .group(1)
                                                    )

                            synapse_metadata['MappedReads_Primary'] = str(
                                int(synapse_metadata['MappedReads_Primary']) +
                                num_mapped - num_secondary - num_supplementary
                            )
                            synapse_metadata['MappedReads_Multimapped'] = str(
                                int(synapse_metadata['MappedReads_Multimapped']) + num_secondary
                            )
                        else:
                            qc_metrics['percent_num_reads_mapped_genome'].append('0')

                        # Pull out multimapped reads
                        target_line = re.search(r'(\d+) \+ \d+ secondary', flagstats_contents)
                        if target_line is not None:
                            qc_metrics['num_reads_multimapped'].append(
                                str(int(target_line.group(1))/2)
                            )
                        else:
                            qc_metrics['num_reads_multimapped'].append('0')
                except:
                    qc_metrics['percent_num_reads_mapped_genome'].append(
                        'Could not open flagstats for {}'.format(star_output_bam)
                    )
                    qc_metrics['num_reads_multimapped'].append(
                        'Could not open flagstats for {}'.format(star_output_bam)
                    )

                # QC: Get number of mapped reads to the transcriptome from this BAM
                try:
                    with open(star_output_transcriptome_bam + '.flagstat') as flagstats:
                        flagstats_contents = flagstats.read()
                        target_line = re.search(r'(\d+) \+ \d+ mapped \(([0-9\.]+)%', flagstats_contents)
                        if target_line is not None:
                            qc_metrics['percent_num_reads_mapped_transcriptome'].append(
                                [str(int(target_line.group(1))/2), '{}%'.format(target_line.group(2))]
                            )
                        else:
                            qc_metrics['percent_num_reads_mapped_transcriptome'].append('0')
                except:
                    qc_metrics['percent_num_reads_mapped_transcriptome'].append(
                        'Could not open flagstats for {}'.format(star_output_bam)
                    )

        # Step 4: BAM Merge | Novosort
        if step <= 4:
            # Novosort to sort and merge BAM files
            novosort_outfile = os.path.join(output_dir,
                                            lib_prefix + ('.' if lib_prefix[-1] != '.' else '') +
                                            'merged.Aligned.out.bam')
            novosort.run(
                Parameter('--tmpdir', os.path.join(output_dir, 'tmp')),
                Parameter(*[bam for bam in star_output]),
                Redirect(stream=Redirect.STDOUT, dest=novosort_outfile)
            )

            # QC: Get number of reads mapped to rRNA regions
            aligned_bed_file = os.path.join(output_dir, str(uuid.uuid4()) + '.bed')
            coverage_file = os.path.join(output_dir, str(uuid.uuid4()) + '.coverage.bed')
            staging_delete.extend([aligned_bed_file, coverage_file])

            bedtools_bamtobed.run(
                Parameter('-i', novosort_outfile),
                Redirect(stream=Redirect.STDOUT, dest=aligned_bed_file)
            )
            bedtools_coverage.run(
                Parameter('-s'),
                Parameter('-counts'),
                Parameter('-a', pipeline_config['qc']['rRNA-bed']),
                Parameter('-b', aligned_bed_file),
                Redirect(stream=Redirect.STDOUT, dest=coverage_file)
            )
            try:
                rRNA_count = 0
                with open(coverage_file) as coverage:
                    for line in coverage:
                        rRNA_count += int(line.strip().split('\t')[6])
                percent_rRNA = (rRNA_count /
                                float(sum([int(aln[MAPPED_READS_COUNT])
                                           for aln
                                           in qc_metrics['percent_num_reads_mapped_transcriptome']]))
                                )
                qc_metrics['percent_num_reads_rrna'] = [str(rRNA_count), str(percent_rRNA)]
                synapse_metadata['rRNARate'] = str(percent_rRNA)
            except Exception as e:
                qc_metrics['percent_num_reads_rrna'] = ['error', 'error', e.message]

            # Prepare genome fasta for RNAseQC
            genome_fa = pipeline_config['qc']['genome-fa']
            genome_fai = genome_fa + '.fai'
            genome_dict = os.path.splitext(genome_fa)[0] + '.dict'

            if not os.path.isfile(genome_fai):
                samtools_faidx.run(
                    Parameter(genome_fa)
                )
            if not os.path.isfile(genome_dict):
                picard_create_seq_dict.run(
                    Parameter('REFERENCE={}'.format(genome_fa)),
                    Parameter('OUTPUT={}'.format(genome_dict))
                )

            # Add read group to alignment file
            read_group_bam = os.path.join(output_dir, 'readgroup.bam')
            staging_delete.append(read_group_bam)
            picard_add_read_groups.run(
                Parameter('INPUT={}'.format(novosort_outfile)),
                Parameter('OUTPUT={}'.format(read_group_bam)),
                Parameter('RGLB={}'.format(lib_prefix)),
                Parameter('RGPL=Illumina'),
                Parameter('RGPU=1'),
                Parameter('RGSM=Sample')
            )

            # Generate BAM index for RNAseQC
            samtools_index.run(
                Parameter(read_group_bam)
            )
            staging_delete.append(read_group_bam + '.bai')

            # QC: Get RNAseQC output
            rnaseqc_output_dir = os.path.join(output_dir, 'RNAseQC')
            subprocess.call(['mkdir', '-p', rnaseqc_output_dir])
            rnaseqc.run(
                Parameter('-o', rnaseqc_output_dir),
                Parameter('-r', genome_fa),
                Parameter('-t', pipeline_config['cufflinks']['transcriptome-gtf']),
                Parameter('-s', '"{sample_id}|{bam_file}|{notes}"'.format(
                    sample_id=lib_prefix,
                    bam_file=read_group_bam,
                    notes='None'
                )),
                Parameter('-singleEnd') if not run_is_paired_end else Parameter()
            )

            # Picard MarkDuplicates to get duplicates metrics
            markduplicates_metrics_filepath = os.path.join(logs_dir, 'mark_dup.metrics')
            picard_markduplicates.run(
                Parameter('INPUT={}'.format(novosort_outfile)),
                Parameter('OUTPUT={}.duprm.bam'.format(lib_prefix)),
                Parameter('TMP_DIR={}'.format(tmp_dir)),
                Parameter('METRICS_FILE={}'.format(markduplicates_metrics_filepath)),
                Redirect(stream=Redirect.BOTH, dest=os.path.join(logs_dir, 'mark_dup.log'))
            )

            # QC: Get percent duplicates
            try:
                with open(markduplicates_metrics_filepath) as markdup_metrics:
                    for line in markdup_metrics:
                        if line[FIRST_CHAR] == '#':
                            continue
                        record = line.strip().split('\t')
                        if len(record) == 9:
                            if re.match(r'\d\.\d+', record[7]) is not None:
                                qc_metrics['percent_duplicate_reads'] = record[7]
            except Exception as e:
                qc_metrics['percent_duplicate_reads'] = ['Could not open MarkDuplicates metrics', e.message]

        # Step 5a: Quantification | Cufflinks
        if step <= 5:
            cufflinks_output_dir = os.path.join(output_dir, 'cufflinks')
            subprocess.call(['mkdir', '-p', cufflinks_output_dir])
            cufflinks.run(
                Parameter('--GTF', pipeline_config['cufflinks']['transcriptome-gtf']),
                Parameter('-p', pipeline_config['cufflinks']['threads']),
                Parameter('--library-type', cufflinks_lib_type),
                Parameter('--upper-quartile-norm'),
                Parameter('-o', cufflinks_output_dir),
                Parameter('--max-bundle-frags', '1000000000'),
                Parameter(novosort_outfile)
            )

        # Step 5b: Quantification | HTSeq
        if step <= 6:
            htseq_output_dir = os.path.join(output_dir, 'htseq')
            subprocess.call(['mkdir', '-p', htseq_output_dir])
            for id_attr in ['gene_id', 'gene_name']:
                for feature_type in ['gene', 'transcript', 'exon']:
                    htseq.run(
                        Parameter('-f', 'bam'),
                        Parameter('-r', 'name'),
                        Parameter('-s', htseq_stranded),
                        Parameter('-t', feature_type),
                        Parameter('-i', id_attr),
                        Parameter(novosort_outfile),
                        Parameter(pipeline_config['htseq']['transcriptome-gtf']),
                        Redirect(stream=Redirect.STDOUT, dest=os.path.join(htseq_output_dir,
                                                                           '{}.{}.counts'.format(feature_type,
                                                                                                 id_attr)))
                    )

        with open(os.path.join(logs_dir, 'qc_metrics.txt'), 'w') as qc_data_file:
            qc_data_file.write(json.dumps(qc_metrics, indent=4) + '\n')

        # Populate Synapse QC matrix
        if re.match(r'\d{4}-\d{4}', lib_prefix.strip()) is not None:
            synapse_metadata['Individual_ID'] = lib_prefix
            synapse_metadata['File_Name'] = 'PEC_BrainGVEX_UIC-UChicago_FC_mRNA_HiSeq2000_{}'.format(lib_prefix)

        re_raw_filename = re.match(r'\d{4}-\d{4}_.+_(.+)_.+_(.+_\d)_\d_sequence\.txt\.gz',
                                   os.path.basename(pipeline_args['reads'][0].split(':')[0]))
        if re_raw_filename is not None:
            sequencing_inst_name = re_raw_filename.group(1)
            if '673' in sequencing_inst_name or '484' in sequencing_inst_name:
                synapse_metadata['SequencingPlatform'] = 'HiSeq2000'
            elif '1070' in sequencing_inst_name:
                synapse_metadata['SequencingPlatform'] = 'HiSeq2500'
            flowcell_batch = re_raw_filename.group(2)
            synapse_metadata['FlowcellBatch'] = flowcell_batch

        total_raw_reads_end1 = sum([int(count[0]) for count in qc_metrics['total_raw_reads']])/4
        synapse_metadata['TotalReads'] = str(total_raw_reads_end1)

        with open(os.path.join(logs_dir, 'synapse_metadata.txt'), 'w') as synapse_metadata_file:
            synapse_metadata_file.write(json.dumps(synapse_metadata, indent=4) + '\n')

        # Delete temporary files
        for delete_file in staging_delete:
            subprocess.call(['rm', '-rf', delete_file])
