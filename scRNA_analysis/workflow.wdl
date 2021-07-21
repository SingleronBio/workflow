version 1.0

struct RuntimeAttr {
  Int? cpu
  Int? memory_gb
  String? docker
  String? queue
}

workflow scRNA_analysis {
  input {
    String sample_name
    String datatype
    String species
    File exp_matrix
    File bam
    String KEGG
    String TF
    String GO
    String ppi
    String geneinf
    String keggdir
    String keggpathway
    String stringdb
    String abbr
    String ppinum
    String Veloref1
    String Veloref2
    String projectID
    String language
    String company
    File mtx_report

    RuntimeAttr runtime_attr_seurat
    RuntimeAttr runtime_attr_generate_anno
    RuntimeAttr runtime_attr_monocle
    RuntimeAttr runtime_attr_velocity
    RuntimeAttr runtime_attr_enrichment
    RuntimeAttr runtime_attr_report

  }
  call scrna_seurat {
    input:
      sample_name = sample_name,
      exp_matrix = exp_matrix,
      datatype = datatype,

      runtime_attr_override = runtime_attr_seurat
  }
  call generate_anno {
    input:
      exp_rds = scrna_seurat.RDSout,
      sample_name = sample_name,

      runtime_attr_override = runtime_attr_generate_anno
  }
  call scrna_monocle {
    input:
      exp_rds = scrna_seurat.RDSout,
      sample_name = sample_name,

      runtime_attr_override = runtime_attr_monocle
  }
  call scrna_velocity {
    input:
      exp_rds = scrna_seurat.RDSout,
      bam = bam,
      species = species,
      sample_name = sample_name,
      ref_anno = generate_anno.rds_anno,
      ref_umap = scrna_seurat.refUMAP,
      exp_matrix = exp_matrix,
      Veloref1 = Veloref1,
      Veloref2 = Veloref2,

      runtime_attr_override = runtime_attr_velocity
  }
  call scrna_enrichment {
    input:
      sample_name = sample_name,
      species = species,
      diffgenes = scrna_seurat.diffgenes,
      TF = TF,
      KEGG = KEGG,
      GO = GO,
      geneinf = geneinf,
      abbr = abbr,
      ppi = ppi,
      stringdb = stringdb,
      keggdir = keggdir,
      ppinum = ppinum,
      keggpathway = keggpathway,

      runtime_attr_override = runtime_attr_enrichment
  }
  call scrna_report {
   input:
     sample_name = sample_name,
     company = company,
     language = language,
     species = species,
     projectID = projectID,
     seurat_result = scrna_seurat.resultFolder,
     monocle_result = scrna_monocle.resultFolder,
     velocity_result = scrna_velocity.resultFolder,
     enrichment_result = scrna_enrichment.resultFolder,
     report_temp = mtx_report,

     runtime_attr_override = runtime_attr_report
  }
  output {
    File resultFolder1 = scrna_report.resultFolder1
    File resultFolder2 = scrna_report.resultFolder2
  }
}

task scrna_seurat {
  input {
    String datatype
    String sample_name
    File exp_matrix

    RuntimeAttr runtime_attr_override
    String targetFolder = "Seurat"
  }
  command {
    set -euo pipefail
    mkdir ~{targetFolder}
    Rscript /opt/scRNA_seq/Seurat/Seruat_analysis.R \
    --datatype  ~{datatype} \
    --expM ~{exp_matrix} \
    --spname ~{sample_name} \
    --prefix ~{sample_name} \
    --outdir ~{targetFolder}
    cp  /opt/scRNA_seq/Seurat/1.QC_mito.readme.txt   "~{targetFolder}/1.QC_mito"
    cp  /opt/scRNA_seq/Seurat/2.Cell_Cluster.readme.txt    "~{targetFolder}/2.Cell_Cluster"
    cp  /opt/scRNA_seq/Seurat/3.Cell_Dimensional_Reduction.readme.txt   "~{targetFolder}/3.Cell_Dimensional_Reduction"
    cp  /opt/scRNA_seq/Seurat/4.Differential_Gene.readme.txt      "~{targetFolder}/4.Differential_Gene"
    tar -zcf  "Differential_Gene.tar.gz"     "~{targetFolder}/4.Differential_Gene"
    tar -zcf ~{targetFolder}.tar.gz ~{targetFolder}
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 8
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  output {
    File RDSout = "~{targetFolder}/Rds_data/~{sample_name}.PRO.rds"
    File refUMAP = "~{targetFolder}/3.Cell_Dimensional_Reduction/~{sample_name}.umap_file.xls"
    File resultFolder = "~{targetFolder}.tar.gz"
    File diffgenes = "Differential_Gene.tar.gz"
  }
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
}

task generate_anno {
  input {
    File exp_rds
    String sample_name

    RuntimeAttr runtime_attr_override
    String outdir = "result"
  }
  command {
    set -euo pipefail
    mkdir ~{outdir}
    Rscript /opt/scRNA_seq/Seurat/Pre_Anno.R \
    --rds ~{exp_rds} \
    --prefix ~{sample_name} \
    --outdir ~{outdir}
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  output {
    File rds_anno = "~{outdir}/~{sample_name}_anno.txt"
  }
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
}

task scrna_monocle {
  input {
    String sample_name
    File exp_rds

    RuntimeAttr runtime_attr_override
    String targetFolder = "monocle"
  }
  command {
    set -euo pipefail
    mkdir ~{targetFolder}
    Rscript  /opt/scRNA_seq/Monocle/Monocle_analysis.R \
    --loadrds ~{exp_rds} \
    --prefix ~{sample_name} \
    --outdir ~{targetFolder}
    cp /opt/scRNA_seq/Monocle/monocle_readme.txt   ~{targetFolder}
    tar -zcf ~{targetFolder}.tar.gz ~{targetFolder}
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 8
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  output {
    File resultFolder = "~{targetFolder}.tar.gz"
  }
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
}

task scrna_velocity {
  input {
    File exp_rds
    File exp_matrix
    File ref_anno
    File ref_umap
    File bam
    String Veloref1
    String Veloref2
    String sample_name
    String species

    RuntimeAttr runtime_attr_override
    String targetFolder = "RNAvelocity"
  }
  command {
    set -euo pipefail
    mkdir -p "~{targetFolder}/~{sample_name}/outs/filtered_gene_bc_matrices"
    Rscript /opt/scRNA_seq/RNAvelocity/MatrixTo10X.R \
    --inputfile ~{exp_matrix} \
    --outdir "~{targetFolder}/~{sample_name}/outs/filtered_gene_bc_matrices/filtered_feature_bc_matrix/"
    python /opt/scRNA_seq/RNAvelocity/AddTagtoBAM.py \
    -l  ~{bam} \
    -o "~{targetFolder}/~{sample_name}/outs" \
    -v    scopeV2
    if [ "~{species}" = "human" ]; then
    velocyto run10x -m  ~{Veloref1}  \
    "~{targetFolder}/~{sample_name}"   ~{Veloref2}
    elif [ "~{species}" = "mouse" ]; then
    velocyto run10x -m  ~{Veloref1} \
    "~{targetFolder}/~{sample_name}"    ~{Veloref2}
    fi
    python /opt/scRNA_seq/RNAvelocity/RNAvelocity_analysis.py \
    -l  "~{targetFolder}/~{sample_name}/velocyto/~{sample_name}.loom" \
    -r   ~{ref_anno} \
    -u   ~{ref_umap} \
    -p   ~{sample_name} \
    -o   ~{targetFolder}
    mv "~{targetFolder}/~{sample_name}/velocyto/~{sample_name}.loom"     ./
    rm -rf ~{targetFolder}/~{sample_name}
    tar -zcf ~{targetFolder}.tar.gz ~{targetFolder}
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 2,
    memory_gb: 15
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  output {
    File resultFolder = "~{targetFolder}.tar.gz"
    File resultloom = "~{sample_name}.loom"
  }
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
}

task scrna_enrichment {
  input {
    File diffgenes
    String KEGG
    String TF
    String GO
    String ppi
    String geneinf
    String stringdb
    String abbr
    String sample_name
    String species
    String keggdir
    String keggpathway
    String ppinum

    RuntimeAttr runtime_attr_override
    String targetFolder = "enrichment"
  }
  command {
    set -euo pipefail
    mkdir ~{targetFolder}
    cp ~{diffgenes}    ~{targetFolder}
    tar -xzf  "~{targetFolder}/Differential_Gene.tar.gz"
    cp -r "Seurat/4.Differential_Gene"     "~{targetFolder}/Differential_Gene"
    cp -r "~{targetFolder}/Differential_Gene"     "~{targetFolder}/Differential_Gene.temp"
    perl  /opt/scRNA_seq/Enrich_anno/enrich_preparation.pl \
    -indir  "~{targetFolder}/Differential_Gene.temp" \
    -outdir   ~{targetFolder} \
    -geneinf   ~{geneinf} \
    -samplename  ~{sample_name}
    python /opt/scRNA_seq/Enrich_anno/enrich_analysis.py \
    --indir "~{targetFolder}" \
    --go ~{GO} \
    --kegg ~{KEGG} \
    --ppi  ~{ppi} \
    --tf ~{TF} \
    --abbr ~{abbr} \
    --samplename   ~{sample_name} \
    --keggid  "~{keggdir}/keggid.txt"  \
    --keggdir  ~{keggdir}  \
    --stringdb ~{stringdb} \
    --ppinum   ~{ppinum} \
    --geneinf  ~{geneinf} \
    --keggpathway ~{keggpathway}
    rm -rf  "~{targetFolder}/Differential_Gene.temp"  "~{targetFolder}/Differential_Gene"  "~{targetFolder}/Differential_Gene.tar.gz"
    tar -zcf ~{targetFolder}.tar.gz ~{targetFolder}
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 10
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  output {
    File resultFolder = "~{targetFolder}.tar.gz"
  }
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
}

task scrna_report {
  input {
    String company
    String species
    String language
    String sample_name
    File seurat_result
    File monocle_result
    File velocity_result
    File enrichment_result
    File report_temp
    String projectID
    String targetFolder = "Report"
    RuntimeAttr runtime_attr_override
  }
  command {
    set -euo pipefail
    mkdir -p "~{targetFolder}/~{sample_name}_RNA_result"
    mkdir -p "~{targetFolder}/~{sample_name}_RNA_report"
    cp ~{seurat_result}    "~{targetFolder}/~{sample_name}_RNA_result"
    cp ~{monocle_result}  "~{targetFolder}/~{sample_name}_RNA_result"
    cp ~{velocity_result}   "~{targetFolder}/~{sample_name}_RNA_result"
    cp ~{enrichment_result}  "~{targetFolder}/~{sample_name}_RNA_result"
    cd  "~{targetFolder}/~{sample_name}_RNA_result"
    tar -xzf  Seurat.tar.gz
    tar -xzf  monocle.tar.gz  
    tar -xzf  RNAvelocity.tar.gz
    tar -xzf  enrichment.tar.gz
    rm -rf Seurat.tar.gz monocle.tar.gz RNAvelocity.tar.gz enrichment.tar.gz
    cd  ../..
    mkdir -p "~{targetFolder}/~{sample_name}_RNA_report/src"
    if [ "~{language}" = "Chinese" ]; then
       if [ "~{company}" = "Singleron" ]; then   
           cp ~{report_temp}  "~{targetFolder}/~{sample_name}_RNA_report/src/web_summary.html"
           perl /opt/scRNA_Report/Report/Singleron/report.V2.0.1.pl \
            -indir  ~{targetFolder} \
            -projectID ~{projectID} \
            -type SCOPEv1 \
            -prefix "~{sample_name}_RNA" \
            -sample ~{sample_name} \
            -rmd  /opt/scRNA_Report/Report/Singleron/report_summary.rmd \
            -pdf /opt/conda/bin/wkhtmltopdf \
            -species ~{species} \
            -rmarkdown  /opt/scRNA_Report/Report/rmarkdown
       elif [ "~{company}" = "Tongyuan" ]; then
           python /opt/scRNA_Report/Report/Tongyuan/logo.py \
            -r ~{report_temp} \
            -n "~{targetFolder}/~{sample_name}_RNA_report/src/web_summary.html"
           perl /opt/scRNA_Report/Report/Tongyuan/tongyuan.report.V2.0.1.pl \
            -indir  ~{targetFolder} \
            -projectID ~{projectID} \
            -type SCOPEv1 \
            -prefix "~{sample_name}_RNA" \
            -sample ~{sample_name} \
            -rmd  /opt/scRNA_Report/Report/Tongyuan/report_summary.rmd \
            -pdf /opt/conda/bin/wkhtmltopdf \
            -species ~{species} \
            -rmarkdown  /opt/scRNA_Report/Report/rmarkdown
       fi
    elif [ "~{language}" = "English" ]; then
       cp ~{report_temp}  "~{targetFolder}/~{sample_name}_RNA_report/src/web_summary.html"
       perl /opt/scRNA_Report/Report/Eng.Singleron/report.V2.0.1.pl \
            -indir  ~{targetFolder} \
            -projectID ~{projectID} \
            -type SCOPEv1 \
            -prefix "~{sample_name}_RNA" \
            -sample ~{sample_name} \
            -rmd  /opt/scRNA_Report/Report/Eng.Singleron/report_summary.rmd \
            -pdf /opt/conda/bin/wkhtmltopdf \
            -species ~{species} \
            -rmarkdown  /opt/scRNA_Report/Report/rmarkdown
    fi
    cd  ~{targetFolder}
    tar -zcf "~{sample_name}_RNA_report.tar.gz" "~{sample_name}_RNA_report"
    tar -zcf "~{sample_name}_RNA_result.tar.gz" "~{sample_name}_RNA_result"
    }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1
  }
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
  output{
    File resultFolder1 = "~{targetFolder}/~{sample_name}_RNA_report.tar.gz"
    File resultFolder2 = "~{targetFolder}/~{sample_name}_RNA_result.tar.gz"
  }
}
