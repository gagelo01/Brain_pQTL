run_all_pqtl_analyses_rosmap_banner <- function(uniprot, 
                                                out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf/",
                                                vcffile_out,
                                                dt_gene_region,
                                                study_toinclude = unique(dt_gene_region$study),
                                                method_list = list("get_uni_cis", "get_coloc"),
                                                parameters = default_param()
) {
  
  stopifnot(all(c("trait", "id",  "gene_region", "hgnc", "study", "UniProt") %in% colnames(dt_gene_region)))
  stopifnot(!(dt_gene_region[, sapply(.SD, function(x) any(is.na(x)))] %>% any))
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  
  gene_region <- dt_gene_region[study %in% study_toinclude & UniProt == uniprot, ]$gene_region %>% unique
  ID <- dt_gene_region[study %in% study_toinclude & UniProt == uniprot, ]$id %>% unique
  vcffile_exp <- paste0(out_wd, ID, "/", ID, ".vcf.gz")
  

  res_all<- GagnonMR::run_all_pqtl_analyses(vcffile_exp = vcffile_exp, 
                                            vcffile_out = vcffile_out, 
                                            chrompos = gene_region, 
                                            method_list = method_list, 
                                            parameters = parameters)
  setDT(res_all)
  res_all <- merge(res_all, dt_gene_region[,.(id,study)], by.x = "id.exposure", by.y = "id")
  return(res_all)
}
