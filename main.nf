$HOSTNAME = ""
params.outdir = 'results'  

//* params.species =  "No selection"  //* @dropdown @description:"If selected, will use IgBLAST ndm and aux files for the species and chain unless overridden by other choices" @options:"No selection",human","mouse","rabbit","rat","rhesus_monkey"
//* params.chain =  ""  //* @dropdown @description:"Chain" @options:"IGH","IGK", "IGL", "TRA", "TRB", "TRD", "TRG"
//* params.igblast_ndm =  ""  //* @input @optional @file @description:"IgBLAST NDM file"
//* params.igblast_aux =  ""  //* @input @optional @file @description:"IgBLAST AUX file"
//* params.cpu =  1  //* @input @description"Number of CPUs to use in this workflow"
detailed_log = params.pipeline.detailed_log
infer_novel_v_alleles = params.pipeline.infer_novel_v_alleles


if (!params.input_reads){params.input_reads = ""} 
if (!params.v_germline_set){params.v_germline_set = ""} 
if (!params.v_germline_set_gapped){params.v_germline_set_gapped = ""} 
if (!params.d_germline_set){params.d_germline_set = ""} 
if (!params.j_germline_set){params.j_germline_set = ""} 
if (!params.infer_novels){params.infer_novels = ""} 

Channel.fromPath(params.input_reads, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_1_reads_g_21}
g_22_germline_set_g_21 = file(params.v_germline_set, type: 'any')
g_23_germline_set_g_21 = file(params.v_germline_set_gapped, type: 'any')
g_24_germline_set_g_21 = file(params.d_germline_set, type: 'any')
g_25_germline_set_g_21 = file(params.j_germline_set, type: 'any')


process init_igblast {

input:
 set val(name),file(reads) from g_1_reads_g_21
 file g_v from g_22_germline_set_g_21
 file g_v_gapped from g_23_germline_set_g_21
 file g_d from g_24_germline_set_g_21
 file g_j from g_25_germline_set_g_21

output:
 set val(name),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j)  into g_21_igblast_input0_g_9

script:
println "basedir: ${baseDir}"
println "init_igblast g_v: ${g_v}"
println "init_igblast g_v realpath: ${g_v.toRealPath()}"
println "init_igblast g_d: ${g_d}"
println "init_igblast g_d realpath: ${g_d.toRealPath()}"
println "init_igblast g_j: ${g_j}"
println "init_igblast g_j realpath: ${g_j.toRealPath()}"

cmds = []

if (g_v.toRealPath().startsWith(baseDir)) {
        fp = baseDir.relativize(g_v.toRealPath())

        println "g_v fp -> ${fp}"

        if (fp.startsWith("database/")) {
                cmds.add("blastdbcmd -entry all -db \$IGDATA/${fp} -out ${g_v}.fasta")
                g_v = g_v.fileName + ".fasta"
        }
}

if (g_d.toRealPath().startsWith(baseDir)) {
        fp = baseDir.relativize(g_d.toRealPath())

        println "g_d fp -> ${fp}"
        if (fp.startsWith("database/")) {
                cmds.add("blastdbcmd -entry all -db \$IGDATA/${fp} -out ${g_d}.fasta")
                g_d = g_d.fileName + ".fasta"
        }
}

if (g_j.toRealPath().startsWith(baseDir)) {
        fp = baseDir.relativize(g_j.toRealPath())

        if (fp.startsWith("database/")) {
                cmds.add("blastdbcmd -entry all -db \$IGDATA/${fp} -out ${g_j}.fasta")
                g_j = g_j.fileName + ".fasta"
        }
}

cmds = cmds.join("\n")
println "${cmds}"

"""
export IGDATA=/usr/local/share/igblast
${cmds}
"""

}


process igblast {

input:
 set val(name),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j) from g_21_igblast_input0_g_9

output:
 set val(name),file("${name}.txt"),file("${name}.fasta"),file(g_v),file(g_v_gapped),file(g_d),file(g_j)  into g_9_igblast_output0_g_8

script:
println "igblast g_v: ${g_v}"
println "igblast g_d: ${g_d}"
println "igblast g_j: ${g_j}"

seq_tpe = 'Tr'

if (params.chain.contains('IG')) {
        seq_type = 'Ig'
}

optional_args = ""

if (params.species != 'No selection') {
        optional_args += " -organism ${params.species}"
}

if (params.igblast_ndm) {
        optional_args += " -custom_internal_data ${params.igblast_ndm}"
}

if (params.igblast_aux) {
        optional_args += " -auxiliary_data ${params.igblast_aux}"
} else {
        optional_args += " -auxiliary_data optional_file/${params.species}_gl.aux"
}

cmds = []

cmds.add("makeblastdb -parse_seqids -dbtype nucl -in ${g_v}")
cmds.add("makeblastdb -parse_seqids -dbtype nucl -in ${g_d}")
cmds.add("makeblastdb -parse_seqids -dbtype nucl -in ${g_j}")

if(reads.toString().contains('.fastq')) {
        cmds.add("cp -n ${reads} ${name}.fasta")
}


cmds = cmds.join("\n")

println "${cmds}"

"""
export IGDATA=/usr/local/share/igblast
${cmds}
igblastn -germline_db_V ${g_v} -germline_db_D ${g_d} -germline_db_J ${g_j} -num_threads ${params.cpu} ${optional_args} -show_translation -ig_seqtype ${seq_type} -outfmt "7 std qseq sseq btop" -query ${name}.fasta >${name}.txt
"""

}


process makedb_igblast {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /makedb_${name}_log.txt$/) "ungenotyped_igblast_log_txt/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /makedb_${name}_detailed_log.txt$/) "detailed_log/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}.tsv$/) "ungenotyped_annotations/$filename"}
input:
 set val(name),file(recs),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j) from g_9_igblast_output0_g_8

output:
 set val(name),file("${name}.tsv"),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j)  into g_8_changeo_plus_reads0_g_11
 file "makedb_${name}_log.txt"  into g_8_outputFileTxt11
 file "makedb_${name}_detailed_log.txt" optional true  into g_8_outputFileTxt22
 file "${name}.tsv"  into g_8_outputFileTSV33
 file "${g_v}"  into g_8_germline_set44
 file "${g_v_gapped}"  into g_8_germline_set55
 file "${g_d}"  into g_8_germline_set66
 file "${g_j}"  into g_8_germline_set77

script:

if (!detailed_log)
{
"""
MakeDb.py igblast -i ${recs} -o ${name}.tsv -s ${reads} -r ${g_v_gapped} ${g_d} ${g_j} --extended >makedb_${name}_log.txt
"""
} else {
"""
MakeDb.py igblast -i ${recs} -o ${name}.tsv -s ${reads} -r ${g_v_gapped} ${g_d} ${g_j} --extended --log makedb_${name}_detailed_log.txt >makedb_${name}_log.txt
"""
}
}

sample_limit = params.downsample.sample_limit

process downsample {

input:
 set val(name),file(recs),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j) from g_8_changeo_plus_reads0_g_11

output:
 set val(name),file("${name}_downsampled.tsv"),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j)  into g_11_changeo_plus_reads0_g_13

"""
Rscript /software/downsample.R ${recs} ${name}_downsampled.tsv ${g_v} ${params.chain} ${sample_limit}
"""

}


process define_clones {

input:
 set val(name), file(recs),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j) from g_11_changeo_plus_reads0_g_13

output:
 set val(name),file("${name}_clones.tsv"),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j)  into g_13_changeo_plus_reads0_g_14

"""
DefineClones.py -d ${recs} -o ${name}_clones.tsv --model hh_s5f --dist 0.1 --mode gene --norm len --act set --sym min --nproc 2
"""

}


process tigger_genotype {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_genotype.tsv$/) "genotype/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_genotype.pdf$/) "genotype_plot/$filename"}
input:
 set val(name), file(recs),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j) from g_13_changeo_plus_reads0_g_14

output:
 set val(name),file(reads),file("${name}_v_germline.fasta"),file("${name}_v_germline_gapped.fasta"),file("${name}_d_germline.fasta"),file("${name}_j_germline.fasta")  into g_14_igblast_input0_g_17
 file "${name}_genotype.tsv"  into g_14_outputFileTSV11
 file "${name}_genotype.pdf"  into g_14_outputFilePdf22

"""
echo ${infer_novel_v_alleles}
Rscript /software/to_tigger.R ${recs} ${infer_novel_v_alleles} ${g_v_gapped} ${g_d} ${g_j} ${name}_
"""

}


process igblast_genotyped {

input:
 set val(name),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j) from g_14_igblast_input0_g_17

output:
 set val(name),file("${name}.txt"),file("${name}.fasta"),file(g_v),file(g_v_gapped),file(g_d),file(g_j)  into g_17_igblast_output0_g_18

script:
println "igblast g_v: ${g_v}"
println "igblast g_d: ${g_d}"
println "igblast g_j: ${g_j}"

seq_tpe = 'Tr'

if (params.chain.contains('IG')) {
        seq_type = 'Ig'
}

optional_args = ""

if (params.species != 'No selection') {
        optional_args += " -organism ${params.species}"
}

if (params.igblast_ndm) {
        optional_args += " -custom_internal_data ${params.igblast_ndm}"
}

if (params.igblast_aux) {
        optional_args += " -auxiliary_data ${params.igblast_aux}"
} else {
        optional_args += " -auxiliary_data optional_file/${params.species}_gl.aux"
}

cmds = []

cmds.add("makeblastdb -parse_seqids -dbtype nucl -in ${g_v}")
cmds.add("makeblastdb -parse_seqids -dbtype nucl -in ${g_d}")
cmds.add("makeblastdb -parse_seqids -dbtype nucl -in ${g_j}")

if(reads.toString().contains('.fastq')) {
        cmds.add("cp -n ${reads} ${name}.fasta")
}


cmds = cmds.join("\n")

println "${cmds}"

"""
export IGDATA=/usr/local/share/igblast
${cmds}
igblastn -germline_db_V ${g_v} -germline_db_D ${g_d} -germline_db_J ${g_j} -num_threads ${params.cpu} ${optional_args} -show_translation -ig_seqtype ${seq_type} -outfmt "7 std qseq sseq btop" -query ${name}.fasta >${name}.txt
"""

}


process makedb_igblast_genotyped {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /makedb_${name}_log.txt$/) "genotyped_igblast_log_txt/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /makedb_${name}_detailed_log.txt$/) "detailed_log/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}.tsv$/) "genotyped_annotations/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${g_v}$/) "genotyped_germline__sets/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${g_v_gapped}$/) "genotyped_germline__sets/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${g_d}$/) "genotyped_germline__sets/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${g_j}$/) "genotyped_germline__sets/$filename"}
input:
 set val(name),file(recs),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j) from g_17_igblast_output0_g_18

output:
 set val(name),file("${name}.tsv"),file(reads),file(g_v),file(g_v_gapped),file(g_d),file(g_j)  into g_18_changeo_plus_reads0_g_19
 file "makedb_${name}_log.txt"  into g_18_outputFileTxt11
 file "makedb_${name}_detailed_log.txt" optional true  into g_18_outputFileTxt22
 file "${name}.tsv"  into g_18_outputFileTSV33
 file "${g_v}"  into g_18_germline_set44
 file "${g_v_gapped}"  into g_18_germline_set55
 file "${g_d}"  into g_18_germline_set66
 file "${g_j}"  into g_18_germline_set77

script:

if (!detailed_log)
{
"""
MakeDb.py igblast -i ${recs} -o ${name}.tsv -s ${reads} -r ${g_v_gapped} ${g_d} ${g_j} --extended >makedb_${name}_log.txt
"""
} else {
"""
MakeDb.py igblast -i ${recs} -o ${name}.tsv -s ${reads} -r ${g_v_gapped} ${g_d} ${g_j} --extended --log makedb_${name}_detailed_log.txt >makedb_${name}_log.txt
"""
}
}


process rabhit_j_haplotype {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.tab$/) "haplotypes/$filename"}
input:
 set val(name),file(recs),file(reads),val(g_v),val(g_v_gapped),val(g_d),val(g_j) from g_18_changeo_plus_reads0_g_19

output:
 file "*.tab" optional true  into g_19_outputFileTSV00


"""
ls /usr/bin
/usr/bin/python3 /software/process_haplotypes.py ${recs} ${name} ${params.chain} -v ${g_v_gapped} -d ${g_d}
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
