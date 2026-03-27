# ==============================================================================
# Utiliza o modelo Hurdle para analisar PctWater (Logit + Gamma GAM).
# ==============================================================================


# ==============================================================================
# PARTE 1 - CONFIGURAÇÃO DO SISTEMA E BIBLIOTECAS
# ==============================================================================
# FINALIDADE: Preparar o ambiente de trabalho do R, garantindo que todos os 
# pacotes necessários para manipulação, modelagem e visualização gráfica estejam disponíveis.
#
# ETAPAS DE PROCESSAMENTO:
# 1. Definição do diretório de trabalho principal.
# 2. Listagem dos pacotes (bibliotecas) exigidos pelo script.
# 3. Verificação e instalação automática de pacotes ausentes.
# 4. Carregamento de todas as bibliotecas no ambiente.
# ==============================================================================
# Definir diretório de trabalho onde os dados estão armazenados e onde os resultados serão salvos
setwd("Caminho/para/seu/diretorio_de_trabalho") # Exemplo: "C:/Projetos/Meu_Projeto"

# --- Bibliotecas ---
pkgs <- c(
  "tidyverse","mgcv","car","corrplot","performance","patchwork","ggpubr","caret",
  "sf","spatialreg","spdep","sjPlot","ggeffects","lmtest","writexl","tidyr"
)

# Instalar pacotes ausentes
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# Carregar bibliotecas
library(tidyverse)
library(mgcv)
library(car)
library(corrplot)
library(performance)
library(patchwork)
library(ggpubr)
library(caret)
library(sf)
library(spatialreg)
library(spdep)
library(sjPlot)
library(ggeffects)
library(lmtest)
library(writexl)
library(tidyr)

message("Bibliotecas carregadas com sucesso.")


# ==============================================================================
# PARTE 2 - CARREGAMENTO E PREPARAÇÃO DOS DADOS
# ==============================================================================
# FINALIDADE: Importar os dados brutos e realizar a limpeza e formatação inicial 
# necessária para as análises estatísticas subsequentes.
#
# ETAPAS DE PROCESSAMENTO:
# 1. Leitura do arquivo de dados tabulares (CSV).
# 2. Conversão de colunas chave para formatos numéricos adequados.
# 3. Tratamento de valores ausentes (NAs) em variáveis específicas.
# 4. Remoção de outliers conhecidos e filtragem de dados inválidos.
# 5. Criação de novas variáveis (ex: classes de área) e ajuste de unidades.
# ==============================================================================
# Definir os nomes dos arquivos de entrada
file_path_csv <- "nome_do_arquivo_de_dados.csv" # Arquivo em formato CSV com todas as variáveis de cada unidade amostral (No presente projeto foram consideradas bacias hidrográficas de nível 9 da base HydroSHEDS)
file_path_shp <- "nome_do_arquivo_espacial.shp" # Arquivo shapefile correspondente às unidades amostrais para análises espaciais (se aplicável)

if (!file.exists(file_path_csv)) stop("Arquivo CSV não encontrado: ", file_path_csv)

bhs_metrics <- read.csv(file_path_csv, sep = ";", header = TRUE, dec = ",")

# Converter colunas específicas para formato numérico (Ajuste os nomes das colunas conforme sua base de dados)
cols_to_convert <- c(
  "DCI_index",           # Índice de Conectividade Dendrítica (DCI)
  "SUB_AREA",            # Área total da bacia hidrográfica
  "ratio_DEFOREST",      # Proporção de cobertura NÃO floresta em cada bacia hidrográfica
  "RdLength",            # Comprimento total de estradas na bacia (km)
  "RdDensity",           # Densidade de estradas na bacia 
  "HyLength",            # Comprimento total da rede hidrográfica na bacia (km)
  "HyDensity",           # Densidade da rede de drenagem na bacia
  "COUNT_area",          # Contagem de áreas de superfície de água exposta (SWA) 
  "MEAN_area",           # Área média das superfície de água exposta (km²)
  "SUM_area_swa",        # Soma total da área de superfície de água exposta (km²)
  "PctWater_swa",        # Proporção de superfície de água exposta na bacia (%)
  "Centroid_X",          # Coordenada X do centroide da bacia
  "Centroid_Y",          # Coordenada Y do centroide da bacia
  "highway_DIST_meters"  # Distância da bacia hidrográfica até rodovias (km)
)
bhs_metrics <- bhs_metrics %>% mutate(across(all_of(cols_to_convert), as.numeric))

# Substituir NAs em SUM_area_swa por 0
bhs_metrics$SUM_area_swa[is.na(bhs_metrics$SUM_area_swa)] <- 0

# Remoção opcional de outliers (Ajuste os IDs conforme a necessidade do seu projeto)
remove_ids <- c("PGM_02", "PGM_07", "PGM_32", "STM_10") # Exemplo: PGM_32 é outlier de DCI vs Desmatamento

bhs_metrics <- bhs_metrics %>%
  filter(!(ID_bh_regi %in% remove_ids),
         CLASS_analyses != 0,
         !is.na(DCI_index),
         DCI_index != 999,
         DCI_index > 0) %>%
  mutate(
    REGION = as.factor(REGION),
    ratio_DEFOREST_adjusted = ratio_DEFOREST + 0.0001,
    RdDensity_div_ratio_DEFOREST = RdDensity / ratio_DEFOREST_adjusted
  )

message("Base de dados carregada, outliers removidos e pré-processamento concluído.")


# ==============================================================================
# PARTE 3 - TRANSFORMAÇÕES DE VARIÁVEIS (YEO-JOHNSON) E ANÁLISES EXPLORATÓRIAS
# ==============================================================================
# FINALIDADE: Padronizar e normalizar as variáveis preditoras e resposta para 
# melhorar a estabilidade e a acurácia dos modelos estatísticos.
#
# ETAPAS DE PROCESSAMENTO:
# 1. Identificação de colunas constantes ou contendo apenas valores ausentes (NAs).
# 2. Aplicação da transformação de Yeo-Johnson nAs variáveis preditoras e resposta.
# 3. Centralização e escalonamento (Z-score) das variáveis transformadas.
# 4. Criação de variáveis logarítmicas e binárias (0/1) auxiliares.
# 5. (Opcional) Geração de matriz de correlação para exploração visual.
# ==============================================================================
# Função auxiliar: detectar colunas constantes ou contendo apenas NA
detect_constant_or_all_na <- function(df) {
  sapply(df, function(x) {
    all(is.na(x)) || (length(na.omit(x)) > 0 && var(x, na.rm = TRUE) == 0)
  })
}

# Preditores para transformar (Yeo-Johnson + center + scale)
predictors_raw <- c("HyDensity","ratio_DEFOREST","highway_DIST_meters", "RSC_density")
bhs_metrics <- bhs_metrics %>% mutate(across(all_of(predictors_raw), as.numeric))

pred_df <- bhs_metrics %>% select(all_of(predictors_raw))
pred_flags <- detect_constant_or_all_na(pred_df)
pred_const_cols <- names(pred_flags)[pred_flags]
preds_to_do <- setdiff(predictors_raw, pred_const_cols)

if (length(preds_to_do) > 0) {
  pred_input <- pred_df %>% select(all_of(preds_to_do))
  pred_medians <- sapply(pred_input, function(x) median(x, na.rm = TRUE))
  for (col in names(pred_input)) pred_input[[col]][is.na(pred_input[[col]])] <- pred_medians[[col]]
  pred_preproc <- caret::preProcess(pred_input, method = c("YeoJohnson","center","scale"))
  pred_trans <- predict(pred_preproc, newdata = pred_input)
  names(pred_trans) <- paste0(names(pred_trans), "_yj")
  bhs_metrics <- bind_cols(bhs_metrics, as_tibble(pred_trans))
}

# Copiar constantes como estão (sufixo _yj)
if (length(pred_const_cols) > 0) {
  for (col in pred_const_cols) bhs_metrics[[paste0(col, "_yj")]] <- bhs_metrics[[col]]
  warning("Preditor(es) constante(s) ou totalmente NA (mantidos sem alteração): ", paste(pred_const_cols, collapse = ", "))
}

message("Transformação de preditores concluída: ", paste(predictors_raw, collapse = ", "))

# Variáveis de resposta: manter originais + log1p + YJ
responses_raw <- c("DCI_index","PctWater_swa")
resp_df <- bhs_metrics %>% select(all_of(responses_raw))
resp_flags <- detect_constant_or_all_na(resp_df)
resp_const_cols <- names(resp_flags)[resp_flags]
resps_to_do <- setdiff(responses_raw, resp_const_cols)

if (length(resps_to_do) > 0) {
  resp_input <- resp_df %>% select(all_of(resps_to_do))
  resp_medians <- sapply(resp_input, function(x) median(x, na.rm = TRUE))
  for (col in names(resp_input)) resp_input[[col]][is.na(resp_input[[col]])] <- resp_medians[[col]]
  resp_preproc <- caret::preProcess(resp_input, method = c("YeoJohnson","center","scale"))
  resp_trans <- predict(resp_preproc, newdata = resp_input)
  names(resp_trans) <- paste0(names(resp_trans), "_yj")
  bhs_metrics <- bind_cols(bhs_metrics, as_tibble(resp_trans))
}

bhs_metrics <- bhs_metrics %>%
  mutate(
    DCI_index_log = log1p(DCI_index),
    PctWater_swa_log = log1p(PctWater_swa),
    Pct_pos = as.integer(PctWater_swa > 0)  # Utilizado para o modelo hurdle
  )

# EXPLORAÇÃO OPCIONAL 
cor_matrix <- cor(
  bhs_metrics %>% select(ends_with("_yj")) %>% select(matches("HyDensity|RdDensity|ratio_DEFOREST|highway_DIST_meters")),
  use = "pairwise.complete.obs"
)
corrplot::corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black")


# ==============================================================================
# PARTE 4 - MODELAGEM E SELEÇÃO DOS MODELOS DE MELHOR PERFORMANCE (GAM E HURDLE)
# ==============================================================================
# FINALIDADE: Construir os modelos com diferentes combinações das variáveis
# preditoras e resposta, a fim de identificar o modelo com melhor perfomance 
# e, assim, possibilitar a análise dos preditores mais significativos.
#
# ETAPAS DE PROCESSAMENTO:
# 1. Criação de modelos com todas as combinações possíveis das variáveis.
# 2. Definição de funções para extração automatizada de métricas (AIC, BIC, R2).
# 3. Ajuste de modelos GAM para a variável DCI (Índice de Conectividade).
# 4. Ajuste de modelos Hurdle (Logit para presença + Gamma para proporção) para PctWater.
# 5. Exportar os resultados brutos e seleção do melhor modelo de cada tipo via AIC.
# ==============================================================================
build_predictor_sets <- function(base_always = "ratio_DEFOREST_yj",
                                 others = c("HyDensity_yj","highway_DIST_meters_yj")) {
  combos <- list()
  combos[[length(combos) + 1]] <- c(base_always)
  for (o in others) combos[[length(combos) + 1]] <- c(base_always, o)
  if (length(others) >= 2) {
    for (pair in combn(others, 2, simplify = FALSE)) combos[[length(combos) + 1]] <- c(base_always, pair)
  }
  combos[[length(combos) + 1]] <- c(base_always, others)
  unique(lapply(combos, function(x) unique(x)))
}

predictor_combos <- build_predictor_sets()

# Função de sumário de métricas
extract_model_info <- function(model, model_name) {
  if (is.null(model)) return(data.frame(Model = model_name, Class = "Failed", Formula = NA, AIC = NA, BIC = NA, R2 = NA, stringsAsFactors = FALSE))
  form_str <- tryCatch({ f <- formula(model); paste(deparse(f), collapse = " ") }, error = function(e) NA)
  aic_val <- tryCatch(AIC(model), error = function(e) NA)
  bic_val <- tryCatch(BIC(model), error = function(e) NA)
  r2_val <- tryCatch({
    if (inherits(model, "lm")) summary(model)$adj.r.squared 
    else if (inherits(model, "betareg")) summary(model)$pseudo.r.squared 
    else performance::r2_nagelkerke(model)
  }, error = function(e) NA)
  
  if (length(aic_val) != 1) aic_val <- NA
  if (length(bic_val) != 1) bic_val <- NA
  if (length(r2_val) != 1) r2_val <- NA
  
  data.frame(Model = model_name, Class = class(model)[1], Formula = form_str, AIC = aic_val, BIC = bic_val, R2 = r2_val, stringsAsFactors = FALSE)
}

# --- DCI: Modelos GAM ---
run_dci_gams <- function(data, predictor_combos, response = "DCI_index_log") {
  results <- list()
  for (preds in predictor_combos) {
    other_preds <- setdiff(preds, "ratio_DEFOREST_yj")
    rhs_terms <- c("s(ratio_DEFOREST_yj)", "REGION")
    if (length(other_preds) > 0) rhs_terms <- c(rhs_terms, paste0("s(", other_preds, ")"))
    
    rhs <- paste(rhs_terms, collapse = " + ")
    form <- as.formula(paste(response, "~", rhs))
    name <- paste0("GAM_DCI | ratio_DEFOREST_yj + REGION | ", paste(preds, collapse = "_"))
    fit <- tryCatch(gam(form, data = data, method = "REML"), error = function(e) NULL)
    results[[name]] <- extract_model_info(fit, name)
  }
  bind_rows(results)
}

dci_summary <- run_dci_gams(bhs_metrics, predictor_combos, response = "DCI_index_log")
# Salva tabela comparativa com todas as combinações do modelo GAM para DCI
write.csv(dci_summary, "resultados_modelos_gam_dci.csv", row.names = FALSE)

# --- PctWater: Modelos Hurdle (Logit + Gamma GAM) ---
run_pctwater_hurdle <- function(data, predictor_combos, response_pos = "Pct_pos", response_gamma = "PctWater_swa") {
  summaries <- list()
  data_pos <- data %>% filter(.data[[response_gamma]] > 0)
  
  for (preds in predictor_combos) {
    other_preds <- setdiff(preds, "ratio_DEFOREST_yj")
    
    # Parte A: Logit
    rhs_logit <- c("ratio_DEFOREST_yj", "REGION")
    if (length(other_preds) > 0) rhs_logit <- c(rhs_logit, other_preds)
    form_logit <- as.formula(paste(response_pos, "~", paste(rhs_logit, collapse = " + ")))
    name_logit <- paste0("Hurdle_Logit | ratio_DEFOREST_yj + REGION | ", paste(preds, collapse = "_"))
    fit_logit <- tryCatch(glm(form_logit, data = data, family = binomial()), error = function(e) NULL)
    summaries[[name_logit]] <- extract_model_info(fit_logit, name_logit)
    
    # Parte B: Gamma
    rhs_gamma <- c("s(ratio_DEFOREST_yj)", "REGION")
    if (length(other_preds) > 0) rhs_gamma <- c(rhs_gamma, paste0("s(", other_preds, ")"))
    form_gamma <- as.formula(paste(response_gamma, "~", paste(rhs_gamma, collapse = " + ")))
    name_gamma <- paste0("Hurdle_PosGamma | ratio_DEFOREST_yj + REGION | ", paste(preds, collapse = "_"))
    fit_gamma <- tryCatch(gam(form_gamma, data = data_pos, family = Gamma(link = "log"), method = "REML"), error = function(e) NULL)
    summaries[[name_gamma]] <- extract_model_info(fit_gamma, name_gamma)
  }
  bind_rows(summaries)
}

pct_summary <- run_pctwater_hurdle(bhs_metrics, predictor_combos)
# Salva tabela comparativa com todas as combinações do modelo Hurdle (Logit + Gamma)
write.csv(pct_summary, "resultados_modelos_hurdle_pctwater.csv", row.names = FALSE)

# --- Seleção dos Melhores Modelos por AIC ---
all_results <- bind_rows(dci_summary, pct_summary)
# Salva um compilado com todos os modelos executados
write.csv(all_results, "comparacao_todos_modelos_aic.csv", row.names = FALSE)

safe_as_num <- function(x) { x2 <- suppressWarnings(as.numeric(as.character(x))); ifelse(is.na(x2), NA_real_, x2) }
extract_response_from_formula <- function(formula_str) { parts <- unlist(strsplit(formula_str, "~", fixed = TRUE)); trimws(parts[1]) }

select_best_by_aic <- function(summary_df, top_n = 1) {
  df <- summary_df
  if (!"Response" %in% names(df)) df$Response <- sapply(df$Formula, extract_response_from_formula)
  if (!"AIC_num" %in% names(df))  df$AIC_num  <- safe_as_num(df$AIC)
  df <- df[!is.na(df$AIC_num), , drop = FALSE]
  df <- df[order(df$Response, df$AIC_num), , drop = FALSE]
  do.call(rbind, lapply(split(df, df$Response), function(sub) head(sub, n = top_n)))
}

best_dci <- select_best_by_aic(dci_summary, top_n = 1)
best_pct <- select_best_by_aic(pct_summary, top_n = 1)

cat("\n--- Top 5 Modelos DCI (por AIC) ---\n")
print(head(dci_summary[order(dci_summary$AIC), ], 5))


# ==============================================================================
# PARTE 5 - TABELA FINAL COM RESULTADOS DOS MODELOS (PADRÃO ABNT/APA)
# ==============================================================================
# FINALIDADE: Sintetizar os resultados dos melhores modelos selecionados em uma 
# tabela formatada e pronta para inclusão em relatórios ou artigos científicos.
#
# ETAPAS DE PROCESSAMENTO:
# 1. Re-ajuste dos melhores modelos (DCI, Logit e Gamma) com toda a base de dados.
# 2. Extração dos coeficientes paramétricos (estimativas/erros) e de suavização.
# 3. Substituição dos nomes técnicos por termos claros e padronizados.
# 4. Formatação de valores de p e inclusão de asteriscos para identificar o grau de significância (*, **, ***).
# 5. Construção e exportação de uma tabela para documento MS Word (.docx).
# ==============================================================================
if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
library(broom); library(flextable); library(officer)

# Refitar os melhores modelos
best_dci_model <- gam(as.formula(best_dci$Formula[1]), data = bhs_metrics, method = "REML")
best_logit_model <- glm(as.formula((best_pct %>% filter(grepl("Logit", Model)))$Formula[1]), data = bhs_metrics, family = binomial())
data_pos <- bhs_metrics %>% filter(PctWater_swa > 0)
best_gamma_model <- gam(as.formula((best_pct %>% filter(grepl("PosGamma", Model)))$Formula[1]), data = data_pos, family = Gamma(link="log"), method = "REML")

# Extração de resultados
extract_model_results <- function(model, model_label) {
  if (is.null(model)) return(NULL)
  s <- summary(model)
  
  if (inherits(model, "gam")) {
    df_param <- data.frame(Termo = rownames(s$p.table), Estimate = s$p.table[, "Estimate"], StdError = s$p.table[, "Std. Error"], Statistic = s$p.table[, "t value"], P_Value = s$p.table[, "Pr(>|t|)"], Type = "Parametric")
  } else {
    df_param <- data.frame(Termo = rownames(coef(s)), Estimate = coef(s)[, "Estimate"], StdError = coef(s)[, "Std. Error"], Statistic = coef(s)[, "z value"], P_Value = coef(s)[, "Pr(>|z|)"], Type = "Parametric")
  }
  
  df_smooth <- NULL
  if (inherits(model, "gam") && !is.null(s$s.table) && nrow(s$s.table) > 0) {
    df_smooth <- data.frame(Termo = rownames(s$s.table), Estimate = s$s.table[, "edf"], StdError = NA, Statistic = s$s.table[, "F"], P_Value = s$s.table[, "p-value"], Type = "Smooth (edf/F)")
  }
  bind_rows(df_param, df_smooth) %>% mutate(Model = model_label)
}

df_all <- bind_rows(
  extract_model_results(best_dci_model, "DCI (Conectividade)"),
  extract_model_results(best_logit_model, "PctWater (Presença/Logit)"),
  extract_model_results(best_gamma_model, "PctWater (Extensão/Gamma)")
)

clean_names <- c(
  "(Intercept)" = "Intercepto", "ratio_DEFOREST_yj" = "Desmatamento (YJ)", "s(ratio_DEFOREST_yj)" = "Desmatamento (Smooth)",
  "REGIONSTM" = "Região [STM]", "HyDensity_yj" = "Densidade Hidrográfica (YJ)", "s(HyDensity_yj)" = "Densidade Hidrográfica (Smooth)",
  "highway_DIST_meters_yj" = "Dist. Rodovias (YJ)", "s(highway_DIST_meters_yj)" = "Dist. Rodovias (Smooth)"
)

table_wide <- df_all %>%
  mutate(
    Termo_Clean = ifelse(Termo %in% names(clean_names), clean_names[Termo], Termo),
    P_Value_Fmt = ifelse(P_Value < 0.001, "< 0.001", sprintf("%.3f", P_Value)),
    Signif = case_when(P_Value < 0.001 ~ "***", P_Value < 0.01 ~ "**", P_Value < 0.05 ~ "*", P_Value < 0.1 ~ ".", TRUE ~ ""),
    Estimate_Fmt = ifelse(Type == "Parametric", sprintf("%.2f (%.2f)", Estimate, StdError), sprintf("%.2f (F=%.2f)", Estimate, Statistic)),
    Value_Cell = paste0(Estimate_Fmt, " ", Signif)
  ) %>%
  select(Termo_Clean, Model, Value_Cell) %>%
  pivot_wider(names_from = Model, values_from = Value_Cell, values_fill = "-")

# Exportar tabela formatada para Word
ft <- flextable(table_wide) %>%
  set_caption(caption = "Tabela 1. Resultados dos modelos aditivos generalizados (GAM) e modelo Hurdle para conectividade (DCI) e disponibilidade hídrica.") %>%
  theme_vanilla() %>% autofit() %>%
  add_footer_lines("Nota: Termos paramétricos = Estimativa (Erro Padrão). Termos de suavização = Graus de Liberdade Efetivos (Estatística F). *** p<0.001, ** p<0.01, * p<0.05.") %>%
  color(part = "footer", color = "#666666") %>% fontsize(part = "footer", size = 9) %>% bold(part = "header")

# Salva a tabela final sumarizada formatada para MS Word
save_as_docx(ft, path = "Tabela_Final_Resultados_Modelos.docx")


# ==============================================================================
# PARTE 6 - DIAGNÓSTICOS DOS MODELOS (DHARMa, Concurvity, ROC)
# ==============================================================================
# FINALIDADE: Avaliar a validade e a qualidade do ajuste dos modelos selecionados, 
# verificando o atendimento às premissas estatísticas e sua capacidade preditiva.
#
# ETAPAS DE PROCESSAMENTO:
# 1. Plotagem de resíduos simulados (DHARMa) para avaliar dispersão e distribuição.
# 2. Cálculo de 'concurvity' para verificar redundância/colinearidade não linear nos GAMs.
# 3. Construção e plotagem da Curva ROC para validar a acurácia do modelo Logit (Presença/Ausência).
# ==============================================================================
if (!requireNamespace("DHARMa", quietly = TRUE)) install.packages("DHARMa")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(DHARMa); library(pROC)

# Diagnósticos DHARMa para resíduos
sim_dci <- simulateResiduals(fittedModel = best_dci_model, n = 200)
plot(sim_dci)
sim_logit <- simulateResiduals(fittedModel = best_logit_model, n = 200)
plot(sim_logit)
sim_gamma <- simulateResiduals(fittedModel = best_gamma_model, n = 200)
plot(sim_gamma)

# Concurvity para GAMs
print(mgcv::concurvity(best_dci_model, full = TRUE))
print(mgcv::concurvity(best_gamma_model, full = TRUE))

# Curva ROC para PctWater (Logit)
probs <- predict(best_logit_model, type = "response")
actual <- bhs_metrics$Pct_pos[complete.cases(bhs_metrics[, all.vars(formula(best_logit_model))])]
roc_obj <- pROC::roc(actual, probs)
plot.roc(roc_obj, main = "ROC - Logit (Presença de Água)")


# ==============================================================================
# PARTE 7 - GRÁFICOS DCI + Ponto de máxima curvatura
# ==============================================================================
# FINALIDADE: Criar gráficos de alta qualidade para o modelo de conectividade (DCI), 
# com identificação e destaque analítico de limiares críticos ecológicos.
#
# ETAPAS DE PROCESSAMENTO:
# 1. Cálculo e identificação do ponto de máxima curvatura (elbow point) na relação DCI vs Desmatamento.
# 2. Extração do R² ajustado do modelo principal para inclusão no gráfico.
# 3. Construção do gráfico de dispersão com curva de tendência e destaque do ponto de máxima curvatura.
# 4. Gráfico finalizado salvo em alta resolução (formato PNG).
# ==============================================================================
library(gratia)

# --- CÁLCULO DO PONTO DE MÁXIMA CURVATURA (DCI vs DEFOREST) ---
get_elbow_point <- function(data, x_col, y_col) {
  df_sub <- na.omit(data.frame(x = data[[x_col]], y = data[[y_col]]))
  model <- lm(y ~ log(x + 1), data = df_sub)
  b1 <- coef(model)[2]
  x_seq <- seq(min(df_sub$x), max(df_sub$x), length.out = 5000)
  
  # Fórmula da Curvatura Numérica
  y_prime <- b1 / (x_seq + 1)
  y_dbl_prime <- -b1 / (x_seq + 1)^2
  kappa <- abs(y_dbl_prime) / (1 + y_prime^2)^(1.5)
  
  max_k_idx <- which.max(kappa)
  return(data.frame(x = x_seq[max_k_idx], y = predict(model, newdata = data.frame(x = x_seq[max_k_idx]))))
}

elbow_pt_dci <- get_elbow_point(bhs_metrics, "ratio_DEFOREST", "DCI_index")

# R2 Ajustado do GAM
gam_fit <- gam(DCI_index_log ~ s(ratio_DEFOREST_yj), data = bhs_metrics, method = "REML")
r_label_dci <- paste0("bold(italic(R)[adj]^2 == '", round(summary(gam_fit)$r.sq, 3), "')")

# Plot DCI Final com Ponto de Inflexão
p_dci_final <- ggplot(bhs_metrics, aes(x = ratio_DEFOREST, y = DCI_index)) +
  geom_point(aes(color = REGION), alpha = 0.8, size = 3) +
  scale_color_manual(values = c("STM" = "#FF0F80", "PGM" = "#FFB400"), name = "Região:") +
  geom_smooth(method = "lm", formula = y ~ log(x + 1), se = TRUE, color = "black", linewidth = 0.8) +
  geom_point(data = elbow_pt_dci, aes(x = x, y = y, fill = "Ponto crítico"), color = "black", size = 5, shape = 21, stroke = 1.5) +
  scale_fill_manual(name = NULL, values = c("Ponto crítico" = "red")) +
  geom_segment(data = elbow_pt_dci, aes(x = x, xend = x, y = 0, yend = y), color = "red", linetype = "dashed", linewidth = 0.6) +
  geom_segment(data = elbow_pt_dci, aes(x = 0, xend = x, y = y, yend = y), color = "red", linetype = "dashed", linewidth = 0.6) +
  ggrepel::geom_label_repel(data = elbow_pt_dci, aes(x = x, y = y, label = paste0("X = ", round(x, 1), "%\nY = ", round(y, 1))), box.padding = 1.5, point.padding = 0.5, nudge_x = 15, nudge_y = 10, color = "red", size = 4.5, fontface = "bold") +
  annotate("text", x = max(bhs_metrics$ratio_DEFOREST), y = max(bhs_metrics$DCI_index), label = r_label_dci, parse = TRUE, hjust = 1, vjust = 1, size = 5) +
  scale_x_continuous(breaks = c(0, 25, 50, 75)) + scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  labs(x = "Desmatamento (%)", y = "Índice de Conectividade Dendrítica (DCI)") +
  theme_classic()

# Salva o gráfico do modelo DCI apontando o ponto crítico (inflexão)
ggsave("grafico_dci_ponto_critico.png", plot = p_dci_final, width = 28, height = 18, units = "cm", dpi = 300)


# ==============================================================================
# PARTE 8 - GRÁFICOS SWA (PctWater LOGIT e GAMMA)
# ==============================================================================
# FINALIDADE: Produzir gráficos detalhados para o modelo Hurdle, demonstrando as 
# tendências na presença e na quantidade de superfície de água exposta.
#
# ETAPAS DE PROCESSAMENTO:
# 1. Cálculo do R² (Nagelkerke) e do ponto de inflexão (probabilidade de 50%) para o modelo Logit.
# 2. Gráfico da curva logística com destaque do ponto de inflexão.
# 3. Cálculo do pseudo-R² e do ponto de máxima curvatura (elbow point) para o modelo Gamma.
# 4. Construção do gráfico de dispersão para o modelo GAMMA com curva de tendência e destaque do ponto de máxima curvatura.
# 5. Gráfico finalizado salvo em alta resolução (formato PNG).
                        
# ==============================================================================
# --- PLOT PctWATER LOGIT ---
r2_logit <- as.numeric(performance::r2_nagelkerke(best_logit_model))
r_label_logit <- paste0("R² (Nagelkerke) = ", round(r2_logit, 3))

# Encontrar o X onde a probabilidade atinge 50%
logit_inflection <- data.frame(x = as.numeric(-coef(best_logit_model)[1] / coef(best_logit_model)[2]), y = 0.5)

p_logit_final <- ggplot(bhs_metrics, aes(x = ratio_DEFOREST, y = Pct_pos)) +
  geom_jitter(aes(color = REGION), height = 0.05, width = 0, alpha = 0.8, size = 3) +
  scale_color_manual(values = c("STM" = "#FF0F80", "PGM" = "#FFB400"), name = "Região:") +
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, color = "black", linewidth = 0.8) +
  geom_point(data = logit_inflection, aes(x = x, y = y, fill = "Ponto de Inflexão"), color = "black", size = 5, shape = 21, stroke = 1.5) +
  scale_fill_manual(name = NULL, values = c("Ponto de Inflexão" = "red")) +
  geom_segment(data = logit_inflection, aes(x = x, xend = x, y = -Inf, yend = y), color = "red", linetype = "dashed", linewidth = 0.6) +
  geom_segment(data = logit_inflection, aes(x = -Inf, xend = x, y = y, yend = y), color = "red", linetype = "dashed", linewidth = 0.6) +
  ggrepel::geom_label_repel(data = logit_inflection, aes(x = x, y = y, label = paste0("X = ", round(x, 1), "%\nY = 50%")), box.padding = 1.5, nudge_x = 15, nudge_y = -0.15, color = "red", fontface = "bold") +
  annotate("text", x = max(bhs_metrics$ratio_DEFOREST), y = 0.7, label = r_label_logit, hjust = 1, vjust = 1, size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(0, 25, 50, 75)) + scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.5, 1)) +
  labs(x = "Desmatamento (%)", y = "Presença de superfície de água (0/1)") +
  theme_classic()

# Salva o gráfico do modelo Hurdle Logit indicando onde a chance de presença de água é de 50%
ggsave("grafico_pctwater_logit_inflexao.png", plot = p_logit_final, width = 28, height = 18, units = "cm", dpi = 300)

# --- PLOT PctWATER GAMMA ---
r2_gamma <- 1 - (best_gamma_model$deviance / best_gamma_model$null.deviance)
r_label_gamma <- paste0("R² (pseudo) = ", round(r2_gamma, 2))

get_gamma_elbow <- function(model, data, x_col) {
  x_seq <- seq(min(data[[x_col]]), max(data[[x_col]]), length.out = 5000)
  preds <- predict(model, newdata = data.frame(setNames(list(x_seq), x_col)), type = "response")
  b1 <- coef(model)[2]
  y_prime <- b1 * preds
  y_dbl_prime <- (b1^2) * preds
  kappa <- abs(y_dbl_prime) / (1 + y_prime^2)^(1.5)
  return(data.frame(x = x_seq[which.max(kappa)], y = preds[which.max(kappa)]))
}

elbow_gamma <- get_gamma_elbow(glm(PctWater_swa ~ ratio_DEFOREST, data=data_pos, family=Gamma(link="log")), data_pos, "ratio_DEFOREST")

p_pct_gamma <- ggplot(data_pos, aes(x = ratio_DEFOREST, y = PctWater_swa)) +
  geom_point(aes(color = REGION), alpha = 0.8, size = 3) +
  scale_color_manual(values = c("STM" = "#FF0F80", "PGM" = "#FFB400"), name = "Região:") +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")), se = TRUE, color = "black", linewidth = 0.8) +
  geom_point(data = elbow_gamma, aes(x = x, y = y, fill = "Ponto crítico"), color = "black", size = 5, shape = 21, stroke = 1.5) +
  scale_fill_manual(name = NULL, values = c("Ponto crítico" = "red")) +
  geom_segment(data = elbow_gamma, aes(x = x, xend = x, y = 0, yend = y), color = "red", linetype = "dashed", linewidth = 0.6) +
  geom_segment(data = elbow_gamma, aes(x = 0, xend = x, y = y, yend = y), color = "red", linetype = "dashed", linewidth = 0.6) +
  ggrepel::geom_label_repel(data = elbow_gamma, aes(x = x, y = y, label = paste0("X = ", round(x, 1), "%\nY = ", round(y, 1), "%")), box.padding = 1.5, nudge_x = 15, nudge_y = 0.5, color = "red", fontface = "bold") +
  annotate("text", x = 65, y = 1.2, label = r_label_gamma, hjust = 1, vjust = 1, size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(0, 25, 50, 75)) +
  labs(x = "Desmatamento (%)", y = "Proporção de superfície de água (%)") +
  theme_classic()

# Salva o gráfico do modelo Hurdle Gamma mostrando o ponto crítico de proporção de água
ggsave("grafico_pctwater_gamma_ponto_critico.png", plot = p_pct_gamma, width = 28, height = 18, units = "cm", dpi = 300)
