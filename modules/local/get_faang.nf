process GET_FAANG {
  publishDir "${params.outdir}/metadata",
      mode: params.publish_dir_mode,
      pattern: "*.tsv"

  container "lauramble/python-vizfada"
  
  output:
  path "input_*.txt", emit: ids
  path "metadata.tsv", emit: metadata
  
  script:
  def species = "${params.species}".capitalize().replace("_", " ")
  if (params.ids) {
    def idsFile = file(params.ids)
    """
    get_faang_data.py "$species" ${params.n_exp} $idsFile
    """
  } else {
    """
    get_faang_data.py "$species" ${params.n_exp}
    """
  }
}
