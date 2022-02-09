ryead_all_sheets = function(xlsxFile) {
  sheet_names = openxlsx::getSheetNames(xlsxFile)
  sheet_list = as.list(rep(NA, length(sheet_names)))
  names(sheet_list) = sheet_names
  for (sn in sheet_names) {
    sheet_list[[sn]] = openxlsx::read.xlsx(xlsxFile, sheet=sn)
  }
  return(sheet_list)
}

write_to_xlsx_multitab <- function(data, name, file){
  wb <- createWorkbook()
  
  for( i in seq_along(name)){
    addWorksheet(wb, name[i])
    writeData(wb, name[i], data[[i]])
  }
  saveWorkbook(wb, file = file, overwrite = TRUE)
}
