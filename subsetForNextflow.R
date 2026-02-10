transcript_matrix <- t_matrix_batch_effect_correction_R30@transcriptogramS2

col_names <- colnames(transcript_matrix)

selected_cols <- c(
  1, # Mantém a coluna Position
  grep("notreated-batch1_", col_names)[1], # 1ª célula dia 0
  grep("notreated-batch2_", col_names)[1], # 1ª célula dia 0b
  grep("TGFbeta1-1day-batch2_", col_names)[1], # 1ª célula dia 1
  grep("TGFbeta1-2day-batch2_", col_names)[1], # 1ª célula dia 2
  grep("TGFbeta1-3day-batch2_", col_names)[1], # 1ª célula dia 3
  grep("TGFbeta1-4day-batch1_", col_names)[1], # 1ª célula dia 4
  grep("TGFbeta1-8day-batch1_", col_names)[1]  # 1ª célula dia 8
) %>% na.omit() # Remove dias faltantes


single_cell_transcript <- transcript_matrix[, selected_cols]


saveRDS(single_cell_transcript, "single_cell_transcriptogram_R30.rds")




# Conferir as células selecionadas
colnames(single_cell_transcript)

