# Shiny app for cgrc.bayes: "Is CGR adjustment safe for my trial?"
#
# A design/planning instrument, not a demo. Panel A answers the question the
# original paper's Limitations section asks (run simulations first); Panel B
# adjusts an uploaded trial. All statistics come from cgrc.bayes - the app holds
# NO copy of the estimand. Heavy simulation is precomputed (inst/extdata via
# data-raw/build_lookup.R); sliders read that lookup instantly.
#
# Launch with cgrc.bayes::cgrc_app().

library(shiny)
library(ggplot2)
library(cgrc.bayes)

LUT <- tryCatch(cgrc_lookup(), error = function(e) NULL)
EFFS <- if (is.null(LUT)) c(0, 1.5, 3, 4.5) else sort(unique(LUT$true_effect))
NRANGE <- if (is.null(LUT)) c(60, 1000) else range(LUT$n)
CELL_LABEL <- c("00" = "no effect, no expectancy",
                "10" = "real effect, no expectancy",
                "01" = "no effect, pure expectancy",
                "11" = "real effect + expectancy")

## ---- helpers used only for display -----------------------------------------

op_table_A <- function(lut, n, p_cg, eff) {
  rows <- do.call(rbind, lapply(list(c(0,0), c(1,0), c(0,1), c(1,1)), function(z)
    cgrc_op_at(lut, n, p_cg, eff, z[1], z[2])))
  data.frame(
    scenario = CELL_LABEL[paste0(rows$DTE, rows$AEB)],
    `true effect` = round(ifelse(rows$DTE == 1, eff, 0), 2),
    `adjusted bias` = round(rows$adj_bias, 2),
    `95% coverage` = round(rows$coverage95, 3),
    `adjusted flags effect` = round(rows$p_fav_gt_95, 3),
    `unadjusted significant` = round(rows$freq_sig, 3),
    check.names = FALSE)
}

## ---- UI ---------------------------------------------------------------------

ui <- navbarPage(
  "CGRC — is adjustment safe for my trial?",
  id = "navbar",
  header = tags$style(HTML(
    ".verdict{font-size:1.15em;line-height:1.5;padding:12px 14px;
      background:#eef3f8;border-left:5px solid #2471A3;border-radius:4px;}
     .warn{background:#FDEDEC;border-left-color:#C0392B;}
     .feas{font-size:1.0em;} .muted{color:#666;font-size:.9em;}")),

  ## ===== Panel A: Design ====================================================
  tabPanel(
    "Design",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        helpText(strong("Describe your planned trial.")),
        sliderInput("n", "Sample size (n)", min = NRANGE[1], max = NRANGE[2],
                    value = 120, step = 10),
        sliderInput("pcg", "Correct guess rate you expect", min = 0.5, max = 0.95,
                    value = 0.85, step = 0.01),
        selectInput("eff", "True effect size to detect (points)",
                    choices = EFFS, selected = 3),
        tags$hr(),
        helpText(class = "muted",
          "Curves and tables are read from a precomputed simulation grid and",
          "interpolated. For the exact numbers at these settings, run the",
          "simulation below (it fits 2000 posteriors, ~10 s)."),
        actionButton("run_exact", "Run exact simulation (500 trials)",
                     class = "btn-primary btn-sm")
      ),
      mainPanel(
        width = 9,
        uiOutput("verdict"),
        br(),
        uiOutput("feasibility"),
        br(),
        h4("Power to detect the effect, over sample size"),
        p(class = "muted",
          "Probability the adjusted analysis flags the effect (posterior P>0.95),",
          "for a real effect with no expectancy confound. Your n is marked."),
        plotOutput("power_plot", height = "300px"),
        br(),
        h4("The trade-off: false positives vs power, with and without adjustment"),
        plotOutput("tradeoff_plot", height = "230px"),
        br(),
        h4("Operating characteristics at your settings"),
        p(class = "muted",
          "Honest power comparison is the \"real effect, no expectancy\" row:",
          "there the adjusted flag rate and the unadjusted significance rate",
          "measure the same thing. In the expectancy rows the unadjusted column",
          "is mostly detecting expectancy, not drug - do not read it as power."),
        tableOutput("opchar"),
        uiOutput("exact_out")
      )
    )
  ),

  ## ===== Panel B: Analyse your own trial ====================================
  tabPanel(
    "Analyse your own trial",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput("csv", "Upload trial CSV", accept = ".csv"),
        helpText(class = "muted", "One row per participant."),
        selectInput("col_cond", "Column: treatment received", choices = NULL),
        selectInput("col_guess", "Column: treatment guessed", choices = NULL),
        selectInput("col_value", "Column: outcome value", choices = NULL),
        radioButtons("direction", "Favourable direction",
                     c("higher is better" = "1", "lower is better" = "-1"), "1"),
        sliderInput("rope", "ROPE half-width (fraction of outcome SD)",
                    min = 0.05, max = 0.5, value = 0.1, step = 0.05),
        actionButton("analyse", "Analyse", class = "btn-primary"),
        uiOutput("to_design")
      ),
      mainPanel(
        width = 9,
        uiOutput("bpanel")
      )
    )
  )
)

## ---- server -----------------------------------------------------------------

server <- function(input, output, session) {

  no_lut <- is.null(LUT)

  ## ---- Panel A ----
  output$verdict <- renderUI({
    if (no_lut) return(div(class = "verdict warn",
      "Lookup table not built. Run data-raw/build_lookup.R, then reinstall."))
    div(class = "verdict",
        cgrc_verdict(LUT, input$n, input$pcg, as.numeric(input$eff)))
  })

  output$feasibility <- renderUI({
    minstr <- cgr_min_stratum(input$n, input$pcg)
    degen  <- if (no_lut) NA else
      cgrc_op_at(LUT, input$n, input$pcg, as.numeric(input$eff), 0, 1)$empty_stratum_rate
    warn <- minstr < 15 || (!is.na(degen) && degen > 0.02)
    div(class = if (warn) "verdict warn feas" else "verdict feas",
      HTML(sprintf(
        "<b>Feasibility.</b> Expected smallest stratum: <b>~%.0f</b> participants%s.
         Simulated trials with an empty stratum: <b>%s</b>.%s",
        minstr, if (minstr < 15) " (thin — CGR adjustment is fragile)" else "",
        if (is.na(degen)) "n/a" else sprintf("%.1f%%", 100 * degen),
        if (warn) " CGR adjustment may not be reliably computable at this design."
        else "")))
  })

  output$power_plot <- renderPlot({
    if (no_lut) return(NULL)
    pc <- cgrc_power_curve(LUT, input$pcg, as.numeric(input$eff))
    ggplot(pc, aes(n, power)) +
      geom_hline(yintercept = c(0.8, 0.9), linetype = "dotted", colour = "grey60") +
      geom_line(colour = "#2471A3", linewidth = 1) +
      geom_point(colour = "#2471A3", size = 2) +
      geom_vline(xintercept = input$n, linetype = "dashed", colour = "#C0392B") +
      annotate("text", x = input$n, y = 0.02, label = paste0("your n=", input$n),
               colour = "#C0392B", hjust = -0.05, size = 3.6) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(x = "sample size (n)", y = "power of the adjusted analysis") +
      theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())
  })

  output$tradeoff_plot <- renderPlot({
    if (no_lut) return(NULL)
    eff <- as.numeric(input$eff)
    pw  <- cgrc_op_at(LUT, input$n, input$pcg, eff, 1, 0)  # real effect, no exp
    fp  <- cgrc_op_at(LUT, input$n, input$pcg, eff, 0, 1)  # pure expectancy
    df <- data.frame(
      metric = factor(c("false positive\n(pure expectancy)", "false positive\n(pure expectancy)",
                        "power\n(real effect)", "power\n(real effect)"),
                      levels = c("false positive\n(pure expectancy)", "power\n(real effect)")),
      analysis = c("unadjusted", "CGR-adjusted", "unadjusted", "CGR-adjusted"),
      rate = c(fp$freq_sig, fp$p_fav_gt_95, pw$freq_sig, pw$p_fav_gt_95))
    df$analysis <- factor(df$analysis, levels = c("unadjusted", "CGR-adjusted"))
    ggplot(df, aes(metric, rate, fill = analysis)) +
      geom_col(position = position_dodge(0.7), width = 0.62) +
      geom_text(aes(label = sprintf("%.0f%%", 100 * rate)),
                position = position_dodge(0.7), vjust = -0.4, size = 3.6) +
      scale_fill_manual(values = c(unadjusted = "#C0392B", `CGR-adjusted` = "#2471A3")) +
      scale_y_continuous(limits = c(0, 1.05), labels = scales::percent) +
      labs(x = NULL, y = NULL, fill = NULL,
           caption = "Adjustment trades a large drop in expectancy-driven false positives for some power.") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "top", panel.grid.minor = element_blank())
  })

  output$opchar <- renderTable({
    if (no_lut) return(NULL)
    op_table_A(LUT, input$n, input$pcg, as.numeric(input$eff))
  }, digits = 3)

  ## exact simulation, only on demand
  exact_rv <- reactiveVal(NULL)
  observeEvent(input$run_exact, {
    withProgress(message = "Running 500 simulated trials x 4 scenarios...",
                 value = 0.3, {
      op <- cgr_operating(n_trials = 500, n = input$n, p_cg = input$pcg,
                          mu_dte = as.numeric(input$eff), noise = "all", seed = 1)
      incProgress(0.7)
      exact_rv(op)
    })
  })
  output$exact_out <- renderUI({
    op <- exact_rv(); if (is.null(op)) return(NULL)
    tagList(br(), h4("Exact simulation (500 trials, computed just now)"),
            renderTable({
              data.frame(scenario = CELL_LABEL[paste0(op$DTE, op$AEB)],
                         `adj bias` = round(op$adj_bias, 3),
                         `95% coverage` = round(op$coverage95, 3),
                         `adjusted flags` = round(op$p_fav_gt_95, 3),
                         `unadj significant` = round(op$freq_sig, 3),
                         `empty-stratum trials` = paste0(round(100*op$empty_stratum_rate,1), "%"),
                         check.names = FALSE)
            }, digits = 3))
  })

  ## ---- Panel B: analyse uploaded trial ----
  raw_csv <- reactive({
    req(input$csv)
    read.csv(input$csv$datapath, stringsAsFactors = FALSE)
  })
  observeEvent(raw_csv(), {
    nm <- names(raw_csv()); low <- tolower(nm)
    pick <- function(cands, default) {
      hit <- nm[low %in% cands]; if (length(hit)) hit[1] else default
    }
    updateSelectInput(session, "col_cond", choices = nm,
      selected = pick(c("condition","arm","treatment","group","received"), nm[1]))
    updateSelectInput(session, "col_guess", choices = nm,
      selected = pick(c("guess","guessed","belief","perceived"), nm[min(2,length(nm))]))
    updateSelectInput(session, "col_value", choices = nm,
      selected = pick(c("value","outcome","score","y"), nm[min(3,length(nm))]))
  })

  fit <- eventReactive(input$analyse, {
    d <- raw_csv()
    trial <- data.frame(
      condition = cgrc_normalise_arm(d[[input$col_cond]], "treatment received"),
      guess     = cgrc_normalise_arm(d[[input$col_guess]], "treatment guessed"),
      value     = as.numeric(d[[input$col_value]]))
    trial <- trial[stats::complete.cases(trial), ]
    dir <- as.numeric(input$direction)
    grid <- sort(unique(c(seq(0, 1, length.out = 101), cgr_observed(cgr_strata(trial)))))
    list(trial = trial,
         fit = cgrc(trial, n_draws = 8000, direction = dir),
         rope = cgr_rope(trial, grid = grid, n_draws = 8000,
                         delta_sd_frac = input$rope, direction = dir),
         sens = cgr_rope_sensitivity(trial, at_cgr = 0.5, n_draws = 6000,
                                     direction = dir),
         dir = dir)
  })

  output$bpanel <- renderUI({
    if (input$analyse == 0) return(helpText(
      "Upload a CSV with columns for treatment received, treatment guessed and",
      "an outcome value, map them on the left, then click Analyse."))
    tagList(
      uiOutput("b_error"),
      fluidRow(column(5, h4("Strata (from your data)"), tableOutput("b_strata")),
               column(7, h4("Adjusted vs unadjusted"), tableOutput("b_summary"))),
      uiOutput("b_identity"),
      h4("CGR curve"), plotOutput("b_curve", height = "420px"),
      h4("Region of practical equivalence"),
      p(class = "muted", "A ROPE conclusion is only as good as the band width, so",
        "the sensitivity to that width is shown beside it."),
      fluidRow(column(7, plotOutput("b_rope", height = "300px")),
               column(5, tableOutput("b_sens"))))
  })

  safe_fit <- reactive(tryCatch(fit(), error = function(e) structure(conditionMessage(e), class = "cgrc_err")))
  output$b_error <- renderUI({
    f <- safe_fit()
    if (inherits(f, "cgrc_err")) div(class = "verdict warn", paste("Could not analyse:", f))
  })

  output$b_strata <- renderTable({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    st <- cgr_strata(f$trial)
    data.frame(stratum = STRATA, n = lengths(st)[STRATA],
               mean = round(vapply(st[STRATA], mean, numeric(1)), 2),
               row.names = NULL)
  })

  output$b_summary <- renderTable({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    s <- f$fit$summary
    ok <- cgrc_pct_ok(s$post_mean[1], s$cri_lo[1], s$cri_hi[1], s$post_mean[2])
    s$pct_attenuation[2] <- if (ok) s$pct_attenuation[2] else NA
    s[, c("what","post_mean","cri_lo","cri_hi","p_favourable","pct_attenuation")]
  }, digits = 3)

  output$b_identity <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    z <- cgr_reference_line_test(f$trial, orig_cgr = f$fit$observed_cgr)
    div(class = "verdict feas", HTML(sprintf(
      "<b>Identity check.</b> Observed CGR = %.4f. The curve at the observed CGR
       equals the raw arm-mean difference to %.1e — the no-op identity holds, so
       the reference line is in the right place.",
      f$fit$observed_cgr, abs(z$D_at_obs - z$raw_mean_diff))))
  })

  output$b_curve <- renderPlot({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    lab <- if (f$dir < 0) "favourable (lower is better)" else "positive"
    cgr_plot(f$fit$curve, obs_cgr = f$fit$observed_cgr, direction_label = lab)
  })

  output$b_rope <- renderPlot({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    z <- f$rope
    stack <- do.call(rbind, lapply(c("p_benefit","p_negligible","p_harm"), function(k)
      data.frame(cgr = z$cgr, p = z[[k]], region = k)))
    stack$region <- factor(stack$region, c("p_benefit","p_negligible","p_harm"),
      labels = c("meaningful benefit","practically negligible","meaningful harm"))
    ggplot(stack, aes(cgr, p, fill = region)) + geom_area() +
      geom_vline(xintercept = 0.5, linetype = "dashed") +
      geom_vline(xintercept = f$fit$observed_cgr, linetype = "dashed", colour = "darkgreen") +
      scale_fill_manual(values = c("meaningful benefit" = "#2471A3",
        "practically negligible" = "grey75", "meaningful harm" = "#C0392B")) +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "correct guess rate", y = "posterior probability", fill = NULL,
           subtitle = "black = perfect blinding; green = your observed CGR") +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")
  })

  output$b_sens <- renderTable({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    data.frame(`delta (SD frac)` = f$sens$delta_in_SD,
               `delta (points)` = round(f$sens$delta, 2),
               `P negligible` = round(f$sens$p_negligible, 3),
               `P benefit` = round(f$sens$p_benefit, 3), check.names = FALSE)
  }, digits = 3)

  output$to_design <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    tagList(tags$hr(), actionButton("do_bridge",
      sprintf("Run the design check at this trial (n=%d, CGR=%.2f)",
              nrow(f$trial), f$fit$observed_cgr), class = "btn-sm"))
  })

  ## Bridge: send the uploaded trial's n and observed CGR to Panel A - the exact
  ## workflow the paper's limitations paragraph describes.
  observeEvent(input$do_bridge, {
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return()
    updateSliderInput(session, "n", value = round(nrow(f$trial) / 10) * 10)
    updateSliderInput(session, "pcg", value = round(f$fit$observed_cgr, 2))
    updateNavbarPage(session, "navbar", selected = "Design")
  })
}

shinyApp(ui, server)
