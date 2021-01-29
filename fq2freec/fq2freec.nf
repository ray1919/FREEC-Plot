EXAC_FILE = file(params.small_exac)
REF_FILE = ["hg38":file(params.hg38_ref), "hg19":file(params.liftover_ref)]
PON_FILE = file(params.pon_file)
AF_FILE = file(params.af_file)
EXOME_REGIONS = file(params.exome_bed)
LIFTOVER_FILE = file(params.liftover_chain)
FUNCO_DATA = file(params.funco_data_source)

inputSample = extractFastq(file(params.infile))

(genderMap, statusMap, inputPairReads) = extractInfos(inputSample)

process saveInputFile {
  publishDir "${params.outdir}/data"

  input:
  file 'samples.tsv' from Channel.fromPath(params.infile)

  output:
  file 'samples.tsv'

    """
    """
}

// preprocess fastq
// align to hg38
process runFastpAndBwa {
  tag "${idSample}"
  publishDir "${params.outdir}/data/fastp", pattern: "*.{html,json}"
  publishDir "${params.outdir}/temp/rawsam", pattern: "*.sam"
  executor 'lsf'
  cpus 28
  clusterOptions "-q ${params.lsf_queue} -P twes -J  ${idSample}_FastpAndBwa -oo  ${idSample}_FastpAndBwa.log"
  
  input:
  set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from inputPairReads
  
  output:
    tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_aligned.sam") into samMapped
  // tuple val(sampleID), file(samfile) into runFastpAndBwaOutput
  file "${idSample}.html"
  file "${idSample}.json"

  script:
  // samfile = idSample + "_aligned.sam"
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    // readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    readGroup = "@RG\\tID:${idSample}\\t${CN}PU:${idSample}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    // adjust mismatch penalty for tumor samples
    status = statusMap[idPatient, idSample]
    extra = status == 1 ? "-B 3" : ""
    test_args = params.test_mode == true ? "--reads_to_process 2000000" : ""

    """
dt=`date +%Y-%m-%d`
fastp --in1 ${inputFile1} --in2 ${inputFile2} \
  --stdout \
  ${test_args} \
  --html ${idSample}.html \
  --json ${idSample}.json \
  --thread 16 \
  --length_required 30 \
  --n_base_limit 5 \
  --cut_front --cut_front_window_size 3 \
  --cut_tail --cut_tail_window_size 3 \
  --cut_right --cut_right_window_size 10 \
  --cut_mean_quality 25 | \
  bwa mem -t 28 \
    -v 1 -Y \
    ${extra} \
    -o ${idSample}_${idRun}_aligned.sam \
    -R "${readGroup}" \
    -p ${REF_FILE["hg38"]} -
    """
}

// clean sam tags and map qualities
// This filter fixes alignments in two ways:
// - it soft clips alignments that hang off the end of its reference sequence
// - it sets the mapping quality to 0 if the alignment is unmapped
process runCleanSam {
  tag "${idSample}"
  publishDir "${params.outdir}/temp/cleansam"
  executor 'lsf'
  cpus 2
  clusterOptions "-q ${params.lsf_queue} -P twes -J  ${idSample}_CleanSam -oo  ${idSample}_CleanSam.log"
  
  input:
  tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_aligned.sam") from samMapped
  
  output:
  tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_cleaned.sam") into CleanSam

  """
gatk CleanSam \
  --VALIDATION_STRINGENCY LENIENT \
  --INPUT "${idSample}_${idRun}_aligned.sam" \
  --OUTPUT "${idSample}_${idRun}_cleaned.sam"
  """

}

// sort sam file and convert to bam
process runSortSam {
  tag "${idSample}"
  publishDir "${params.outdir}/temp/sortsam"
  executor 'lsf'
  cpus 5
  clusterOptions "-q ${params.lsf_queue} -P twes -J  ${idSample}_SortSam -oo  ${idSample}_SortSam.log"
  
  input:
  tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_cleaned.sam") from CleanSam
  
  output:
  tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into BamMapped

  """
samtools sort --threads 5 -m 2G \
  -T ${idSample} \
  "${idSample}_${idRun}_cleaned.sam" \
  -o "${idSample}_${idRun}.bam"
  """
}

// Mark alignment duplicates
process runMarkDuplicates {
  tag "${idSample}"
  publishDir "${params.outdir}/data/MarkDuplicates", pattern: "${idSample}.bam.metrics"
  publishDir "${params.outdir}/temp/markdup", pattern: "${idSample}.md.ba*"
  executor 'lsf'
  cpus 26
  clusterOptions "-q ${params.lsf_queue} -P twes -J  ${idSample}_MarkDuplicates -oo  ${idSample}_MarkDuplicates.log"
  
  input:
  tuple idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") from BamMapped

  output:
    tuple idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bam.bai") into bam_duplicates_marked1, bam_duplicates_marked2
    tuple idPatient, idSample into tsv_bam_duplicates_marked
    file ("${idSample}.bam.metrics") optional true into duplicates_marked_report

  script:
  // --remove-sequencing-duplicates \

  """
gatk MarkDuplicatesSpark \
  --create-output-bam-index \
  --input ${idSample}_${idRun}.bam \
  --output ${idSample}.md.bam \
  --metrics-file ${idSample}.bam.metrics \
  --spark-master local[*]
  """
}

// do alignment statistics
process runBamdst {
  tag "${idSample}"
  publishDir "${params.outdir}/data/bamdst/${idSample}"
  executor 'lsf'
  cpus 1
  clusterOptions "-q ${params.lsf_queue} -P twes -J  ${idSample}_Bamdst -oo  ${idSample}_Bamdst.log"
  
  input:
  tuple idPatient, idSample, file(bamfile), file(baifile) from bam_duplicates_marked1

  output:
  file "coverage.report"
  file 'insertsize.plot'

  """
  bamdst -p ${EXOME_REGIONS} \
      --flank 200 --cutoffdepth 20 \
      --isize 1000 \
      -o ./ ${bamfile}
  """
}

// do bam pileup
process runPileup {
  tag "${idSample}"
  publishDir "${params.outdir}/temp/pileup"
  executor 'lsf'
  cpus 1
  clusterOptions "-q ${params.lsf_queue} -P twes -J  ${idSample}_pileup -oo  ${idSample}_pileup.log"
  
  input:
  tuple idPatient, idSample, file(bamfile), file(baifile) from bam_duplicates_marked2

  output:
  tuple idPatient, idSample, file("${idSample}.pileup") into mpileupOut

  """
  # Control-FREEC reads uncompresses the zipped file TWICE in single-threaded mode.
  # we are therefore not using compressed pileups here
  samtools mpileup \
      -f ${REF_FILE["hg38"]} \
      --positions ${EXOME_REGIONS} \
      ${bamfile} > ${idSample}.pileup
  """
}

// STEP CONTROLFREEC.1 - CONTROLFREEC

mpileupOutNormal = Channel.create()
mpileupOutTumor = Channel.create()

mpileupOut
    .choice(mpileupOutTumor, mpileupOutNormal) {statusMap[it[0], it[1]] == 0 ? 1 : 0}

mpileupOut = mpileupOutNormal.combine(mpileupOutTumor, by:0)

mpileupOut = mpileupOut.map {
    idPatientNormal, idSampleNormal, mpileupOutNormal,
    idSampleTumor, mpileupOutTumor ->
    [idPatientNormal, idSampleNormal, idSampleTumor, mpileupOutNormal, mpileupOutTumor]
}


process runControlFREEC {
  tag "${idPatient}"
  publishDir "${params.outdir}/temp/freec/${idPatient}"
  executor 'lsf'
  cpus 4
  clusterOptions "-q ${params.lsf_queue} -R 'rusage[mem=45G]' -P twes -J  ${idPatient}_freec -oo  ${idPatient}_freec.log"
  conda '/Business/gene_company/zhaorui/miniconda3/envs/freec'
  // use dbSNP154 will require about 42G memory per run

  tag "${idSampleTumor}_vs_${idSampleNormal}"

  input:
  tuple idPatient, idSampleNormal, idSampleTumor, file(mpileupNormal), "${idPatient}" from mpileupOut

  output:
  tuple idPatient, idSampleNormal, idSampleTumor, file("${idPatient}_CNVs"), file("${idPatient}_ratio.txt"), file("${idPatient}_BAF.txt") into controlFreecViz
  tuple file("${idPatient}_ratio.BedGraph"), file("${idSampleTumor}_vs_${idSampleNormal}.config.txt") into controlFreecOut
  tuple idPatient, idSampleNormal, idSampleTumor, file("${idPatient}_sample.cpn"), file("${idSampleNormal}.pileup_control.cpn") into controlFreecCPN
  tuple idPatient, idSampleNormal, idSampleTumor, file("${idPatient}_subclones.txt") into controlFreecSubclones
  file ("GC_profile.targetedRegions.cnp")

  // when: 'controlfreec' in tools

  script:
  config = "${idSampleTumor}_vs_${idSampleNormal}.config.txt"
  gender = genderMap[idPatient]

  """
  touch ${config}
  echo "[general]" >> ${config}
  echo "BedGraphOutput = TRUE" >> ${config}
  echo "chrFiles = ${params.freec_chrFiles}" >> ${config}
  echo "chrLenFile = ${params.freec_chrLenFile}" >> ${config}
  echo "window = 0" >> ${config}
  echo "contaminationAdjustment = TRUE" >> ${config}
  echo "degree = 1" >> ${config}
  echo "forceGCcontentNormalization = 1" >> ${config}
  echo "maxThreads = 7" >> ${config}
  echo "minimalSubclonePresence = 30" >> ${config}
  echo "ploidy = 2" >> ${config}
  echo "intercept = 1" >> ${config}
  echo "minCNAlength = 3" >> ${config}
  echo "noisyData = TRUE" >> ${config}
  echo "printNA = FALSE" >> ${config}
  echo "readCountThreshold = 20" >> ${config}
  echo "sex = ${gender}" >> ${config}
  echo "samtools = /Business/gene_company/zhaorui/bin/samtools" >> ${config}
  echo "sambamba = /Business/gene_company/zhaorui/miniconda3/envs/freec/bin/sambamba " >> ${config}
  echo "SambambaThreads = 2" >> ${config}
  echo "" >> ${config}
  echo "[control]" >> ${config}
  echo "inputFormat = pileup" >> ${config}
  echo "mateFile = ${mpileupNormal}" >> ${config}
  echo "mateOrientation = FR" >> ${config}
  echo "" >> ${config}
  echo "[sample]" >> ${config}
  echo "inputFormat = pileup" >> ${config}
  echo "mateFile = ${idPatient}" >> ${config}
  echo "mateOrientation = FR" >> ${config}
  echo "" >> ${config}
  echo "[BAF]" >> ${config}
  echo "shiftInQuality = 33" >> ${config}
  echo "SNPfile = ${params.dbsnp_vcf}" >> ${config}
  echo "" >> ${config}
  echo "[target]" >> ${config}
  echo "captureRegions = /Business/gene_company/zhaorui/data/bed/SureSelect_Human_All_Exon_V6_r2_hs_hg38/S07604514_Covered.bed" >> ${config}
  freec -conf ${config}
  """
}



// Extract gender and status from Channel
def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [genderMap, statusMap, channel]
}

// Channeling the TSV file containing FASTQ
// Format is: "subject gender status sample lane fastq1 fastq2"
def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def idRun      = row[4]
            def file1      = returnFile(row[5])
            def file2      = "null"
            if (hasExtension(file1, "fastq.gz") || hasExtension(file1, "fq.gz") || hasExtension(file1, "fastq") || hasExtension(file1, "fq")) {
                checkNumberOfItem(row, 7)
                file2 = returnFile(row[6])
            if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")  && !hasExtension(file2, "fastq") && !hasExtension(file2, "fq")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
            if (hasExtension(file1, "fastq") || hasExtension(file1, "fq") || hasExtension(file2, "fastq") || hasExtension(file2, "fq")) {
                exit 1, "We do recommend to use gziped fastq file to help you reduce your data footprint."
            }
        }
        else if (hasExtension(file1, "bam")) checkNumberOfItem(row, 6)
        else "No recognisable extention for input file: ${file1}"

        [idPatient, gender, status, idSample, idRun, file1, file2]
    }
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}"
    return file(it)
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}"
    return it
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}"
    return true
}
