version 1.0

struct RuntimeAttr {
  Int? cpu
  Int? memory_gb
  String? docker
  String? queue
}

workflow rawdata {
  input {
    String project_id
    String sample_id
    String library_id
    String batch_id
    String type
    String bucket
    String endpoint
    String aki
    String aks
    Array[Map[String, String]] files

    RuntimeAttr runtime_attr_download
  }
  scatter(file in files) {
    if (type == "obs") {
      call obs {
        input:
          type = type,
          bucket = bucket,
          endpoint = endpoint,
          aki = aki,
          aks = aks,
          file = file,

          runtime_attr_override = runtime_attr_download
      }
    }
    if (type == "oss") {
      call oss {
        input:
          type = type,
          bucket = bucket,
          endpoint = endpoint,
          aki = aki,
          aks = aks,
          file = file,

          runtime_attr_override = runtime_attr_download
      }
    }
    if (type == "s3") {
      call s3 {
        input:
          type = type,
          bucket = bucket,
          endpoint = endpoint,
          aki = aki,
          aks = aks,
          file = file,

          runtime_attr_override = runtime_attr_download
      }
    }
  }
  output {
    Array[File?]? fastqs = select_all(flatten([oss.fastq, obs.fastq, s3.fastq]))
  }
}

task obs {
  input {
    String type
    String bucket
    String endpoint
    String aki
    String aks
    Map[String, String] file

    RuntimeAttr runtime_attr_override
  }
  command {
    set -euo pipefail
    obsutil cp -p 100 -e "~{endpoint}" -i "~{aki}" -k "~{aks}" "~{type}://~{bucket}/~{file.path}" . 1>/dev/null
    echo "~{file.md5}" "$(basename "~{file.path}")" | md5sum --check
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
  output {
    File fastq = basename(file.path)
  }
}

task oss {
  input {
    String type
    String bucket
    String endpoint
    String aki
    String aks
    Map[String, String] file

    RuntimeAttr runtime_attr_override
  }
  command {
    set -euo pipefail
    ossutil64 cp -e "~{endpoint}" -i "~{aki}" -k "~{aks}" "~{type}://~{bucket}/~{file.path}" . 1>/dev/null
    echo "~{file.md5}" "$(basename "~{file.path}")" | md5sum --check
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
  output {
    File fastq = basename(file.path)
  }
}

task s3 {
  input {
    String type
    String bucket
    String endpoint
    String aki
    String aks
    Map[String, String] file

    RuntimeAttr runtime_attr_override
  }
  command {
    set -euo pipefail
    export AWS_ACCESS_KEY_ID="~{aki}"
    export AWS_SECRET_ACCESS_KEY="~{aks}"
    export AWS_DEFAULT_REGION="~{endpoint}"
    /usr/local/aws-cli/v2/current/bin/aws s3 cp --region "~{endpoint}" "~{type}://~{bucket}/~{file.path}" . 1>/dev/null
    echo "~{file.md5}" "$(basename "~{file.path}")" | md5sum --check
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1
  }
  Int cpu = select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
  Int memory_gb = select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])
  runtime {
    cpu: select_first([runtime_attr_override.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr_override.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr_override.docker, runtime_attr_default.docker])
    queueArn: select_first([runtime_attr_override.queue, runtime_attr_default.queue])
  }
  output {
    File fastq = basename(file.path)
  }
}
