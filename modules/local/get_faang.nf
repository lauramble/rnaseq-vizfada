process GET_FAANG {
  publishDir "${params.outdir}/metadata",
      mode: params.publish_dir_mode,
      pattern: "*.tsv"

  conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.8.3"
  } else {
      container "quay.io/biocontainers/python:3.8.3"
  }
  
  output:
  path "input_*.txt", emit: ids
  path "metadata.tsv", emit: metadata
  
  script:
  def species = "${params.species}".capitalize().replace("_", " ")
  """
  get_faang_data.py "$species" ${params.n_exp}
  """
}
