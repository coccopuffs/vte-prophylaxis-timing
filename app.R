# app.R: Shiny app for patient-specific VTE partial-dependence (PD) plot
# (with peak CI reporting)

# ── libraries ──────────────────────────────────────────────────────────────
library(shiny)
library(future); plan(multisession, workers = 1)   # one background worker
library(promises)
library(randomForestSRC)
library(ggplot2)

MODEL_URL <- "https://vte-prophylaxis-timing.s3.us-east-2.amazonaws.com/rsf_wt.RDS"

# ── UI ─────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Patient-specific VTE PD Plot"),
  sidebarLayout(
    sidebarPanel(
      numericInput("gcs",  "TBI Highest Total GCS:",        15, min = 3,  max = 15),
      numericInput("iss",  "Injury Severity Score (ISS):",  10, min = 0),
      numericInput("age",  "Age (years):",                  50, min = 0),
      checkboxInput("preAnticoag", "Pre-hospital Anticoagulant Therapy", FALSE),
      checkboxInput("cva",  "History of CVA",        FALSE),
      checkboxInput("mi",   "History of MI",         FALSE),
      checkboxInput("cirr", "History of Cirrhosis",  FALSE),
      checkboxInput("bleed","Bleeding Disorder",     FALSE),
      checkboxInput("htn",  "Hypertension",          FALSE),
      selectInput("vteType", "VTE Prophylaxis Type:",
                  c("LMWH", "Unfractionated Heparin")),
      checkboxInput("icu", "ICU Stay", TRUE),
      selectInput("sex", "Sex:", c("Male" = 0, "Female" = 1)),
      actionButton("plotPD", "Generate PD Plot", class = "btn-primary"),
      hr(),
      textOutput("status")
    ),
    mainPanel(
      plotOutput("pdPlot"),
      verbatimTextOutput("peakInfo")
    )
  )
)

# ── server ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  ## --- 1.  kick off background download -----------------------------------
  output$status <- renderText("Downloading RSF model …")
  enable_btn          <- reactiveVal(FALSE)    # flag when model is ready
  rsf_model           <- reactiveVal(NULL)     # will hold the fitted object
  
  fut <- future({
    options(timeout = 1200)                    # 20-min window just in case
    readRDS(url(MODEL_URL, "rb"))
  })
  
  then(
    fut,
    onFulfilled = function(mod) {
      rsf_model(mod)
      enable_btn(TRUE)
      output$status <- renderText("Model ready ✔")
    },
    onRejected  = function(err) {
      msg <- paste("❌ Model failed to download:", conditionMessage(err))
      output$status <- renderText(msg)
    }
  )
  
  ## --- 2.  observe to enable / disable button -----------------------------
  observe({
    shinyjs::toggleState("plotPD", condition = enable_btn())
  })
  
  ## --- 3.  base covariate reactive ---------------------------------------
  base_covariates <- reactive({
    list(
      TBIHIGHESTTOTALGCS = input$gcs,
      ISS                = input$iss,
      AGEyears           = input$age,
      PreHospital_Anticoagulant_Therapy = as.numeric(input$preAnticoag),
      History_of_CVA     = as.numeric(input$cva),
      History_of_MI      = as.numeric(input$mi),
      History_of_Cirrhosis = as.numeric(input$cirr),
      Bleeding_Disorder  = as.numeric(input$bleed),
      Hypertension       = as.numeric(input$htn),
      VTEPROPHYLAXISTYPE = ifelse(input$vteType == "Unfractionated Heparin", 1, 0),
      ICU_Stay           = as.numeric(input$icu),
      SEX                = as.numeric(input$sex)
    )
  })
  
  ## --- 4.  partial-dependence data ---------------------------------------
  pd_data <- eventReactive(input$plotPD, {
    req(enable_btn())                        # ensure model is loaded
    mod  <- rsf_model()
    
    vte_seq <- seq(0, 168, length.out = 300)
    bcov    <- base_covariates()
    
    newdata <- do.call(
      rbind,
      lapply(vte_seq, \(x) as.data.frame(c(bcov, VTEPROPHYLAXISHRS = x)))
    )
    
    pred <- predict(mod, newdata = newdata)
    surv <- pred$survival[, ncol(pred$survival)]
    best_idx <- which.max(surv)
    
    list(seq = vte_seq, surv = surv, best_idx = best_idx, newdata = newdata, mod = mod)
  })
  
  ## --- 5.  plot -----------------------------------------------------------
  output$pdPlot <- renderPlot({
    res <- pd_data()
    df  <- data.frame(VTEHR = res$seq, Surv = res$surv)
    best_hr   <- df$VTEHR[res$best_idx]
    best_surv <- df$Surv[res$best_idx]
    
    ggplot(df, aes(VTEHR, Surv)) +
      geom_line(linetype = "dashed") +
      geom_vline(xintercept = best_hr, linetype = "dotted") +
      geom_point(aes(best_hr, best_surv), size = 3) +
      labs(
        title = sprintf(
          "Patient-specific PD: Survival to %.0f hrs",
          max(res$mod$time.interest)
        ),
        x = "Time to first VTE prophylaxis (hrs)",
        y = "Predicted survival probability"
      ) +
      theme_minimal()
  })
  
  ## --- 6.  95 % interval at peak -----------------------------------------
  output$peakInfo <- renderPrint({
    res <- pd_data()
    allpred <- predict(res$mod, newdata = res$newdata, predict.all = TRUE)
    surv_all <- if (is.list(allpred$survival)) simplify2array(allpred$survival)
    else allpred$survival
    tp_idx <- if (length(dim(surv_all)) == 3) dim(surv_all)[2] else 1
    mat    <- if (length(dim(surv_all)) == 3) surv_all[, tp_idx, ] else
      matrix(surv_all, ncol = 1)
    vals <- mat[res$best_idx, ]
    ci   <- quantile(vals, c(.025, .975), na.rm = TRUE)
    cat(sprintf(
      "Peak at %.1f hrs: Survival = %.3f  (95%% CI: %.3f – %.3f)",
      res$seq[res$best_idx], res$surv[res$best_idx], ci[1], ci[2]
    ))
  })
}

# ── launch app ─────────────────────────────────────────────────────────────
shinyApp(ui, server)

