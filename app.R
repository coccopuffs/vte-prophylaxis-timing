# app.R: Shiny app for patient-specific VTE partial-dependence (PD) plot
# (with peak CI reporting)

library(shiny)
library(randomForestSRC)
library(ggplot2)
library(rsconnect)
library(pins)

op <- options(timeout = 2400)                   # 10 min
rsf_wt <- readRDS(url(
  "https://vte-prophylaxis-timing.s3.us-east-2.amazonaws.com/rsf_wt.RDS",
  open = "rb"                                  # binary read
))
options(op)                    # reset to original timeout

ui <- fluidPage(
  titlePanel("Patient-specific VTE PD Plot"),
  sidebarLayout(
    sidebarPanel(
      numericInput("gcs", "TBI Highest Total GCS:", 15, min = 3, max = 15),
      numericInput("iss", "Injury Severity Score (ISS):", 10, min = 0),
      numericInput("age", "Age (years):", 50, min = 0),
      checkboxInput("preAnticoag", "Pre-hospital Anticoagulant Therapy", FALSE),
      checkboxInput("cva", "History of CVA", FALSE),
      checkboxInput("mi", "History of MI", FALSE),
      checkboxInput("cirr", "History of Cirrhosis", FALSE),
      checkboxInput("bleed", "Bleeding Disorder", FALSE),
      checkboxInput("htn", "Hypertension", FALSE),
      selectInput("vteType", "VTE Prophylaxis Type:", c("LMWH", "Unfractionated Heparin")),
      checkboxInput("icu", "ICU Stay", TRUE),
      selectInput("sex", "Sex:", c("Male" = 0, "Female" = 1)),
      actionButton("plotPD", "Generate PD Plot")
    ),
    mainPanel(
      plotOutput("pdPlot"),
      verbatimTextOutput("peakInfo")
    )
  )
)

server <- function(input, output) {
  base_covariates <- reactive({
    list(
      TBIHIGHESTTOTALGCS            = input$gcs,
      ISS                           = input$iss,
      AGEyears                      = input$age,
      PreHospital_Anticoagulant_Therapy = as.numeric(input$preAnticoag),
      History_of_CVA                = as.numeric(input$cva),
      History_of_MI                 = as.numeric(input$mi),
      History_of_Cirrhosis          = as.numeric(input$cirr),
      Bleeding_Disorder             = as.numeric(input$bleed),
      Hypertension                  = as.numeric(input$htn),
      VTEPROPHYLAXISTYPE            = ifelse(input$vteType == "Unfractionated Heparin", 1, 0),
      ICU_Stay                      = as.numeric(input$icu),
      SEX                           = as.numeric(input$sex)
    )
  })
  
  pd_data <- eventReactive(input$plotPD, {
    vte_seq <- seq(0, 168, length.out = 300)
    bcov <- base_covariates()
    newdata <- do.call(rbind, lapply(vte_seq, function(x) as.data.frame(c(bcov, VTEPROPHYLAXISHRS = x))))
    # Ensemble predictions
    pred <- predict(rsf_wt, newdata = newdata)
    surv <- pred$survival[, ncol(pred$survival)]
    best_idx <- which.max(surv)
    list(seq = vte_seq, surv = surv, best_idx = best_idx, newdata = newdata)
  })
  
  output$pdPlot <- renderPlot({
    res <- pd_data()
    df <- data.frame(VTEHR = res$seq, Surv = res$surv)
    best_hr <- df$VTEHR[res$best_idx]
    best_surv <- df$Surv[res$best_idx]
    
    ggplot(df, aes(x = VTEHR, y = Surv)) +
      geom_line(linetype = "dashed") +
      geom_vline(xintercept = best_hr, linetype = "dotted") +
      geom_point(aes(x = best_hr, y = best_surv), size = 3) +
      labs(
        title = paste0("Patient-specific PD: Survival to ",
                       round(max(rsf_wt$time.interest)), " hrs"),
        x = "Time to VTE prophylaxis (hrs)",
        y = "Predicted survival probability"
      ) +
      theme_minimal()
  })
  
  output$peakInfo <- renderPrint({
    res <- pd_data()
    try({
      allpred <- predict(rsf_wt, newdata = res$newdata, predict.all = TRUE)
      surv_all <- allpred$survival
      if (is.list(surv_all)) surv_arr <- simplify2array(surv_all) else surv_arr <- surv_all
      tp <- if (length(dim(surv_arr)) == 3) dim(surv_arr)[2] else 1
      mat <- if (length(dim(surv_arr)) == 3) surv_arr[, tp, ] else matrix(surv_arr, ncol = 1)
      vals <- mat[res$best_idx, ]
      ci <- quantile(vals, c(0.025, 0.975), na.rm = TRUE)
      cat(sprintf(
        "Peak at %.1f hrs: Survival = %.3f (95%% CI: %.3f - %.3f)",
        res$seq[res$best_idx], res$surv[res$best_idx], ci[1], ci[2]
      ))
    }, silent = TRUE)
  })
}

shinyApp(ui, server)

rsconnect::writeManifest(".", force = TRUE)
