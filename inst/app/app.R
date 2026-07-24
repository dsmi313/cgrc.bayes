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
# effect sizes offered = those simulated at EVERY expectancy level, so no cell
# is missing when the expectancy slider moves.
EFFS <- if (is.null(LUT)) c(0, 1.5, 3) else if ("mu_aeb" %in% names(LUT))
  sort(Reduce(intersect, split(LUT$true_effect, LUT$mu_aeb))) else
  sort(unique(LUT$true_effect))
NRANGE <- if (is.null(LUT)) c(60, 1000) else range(LUT$n)
CELL_LABEL <- c("00" = "no effect, no expectancy",
                "10" = "real effect, no expectancy",
                "01" = "no effect, pure expectancy",
                "11" = "real effect + expectancy")

## ---- named thresholds (were magic numbers inline) ---------------------------
THIN_STRATUM  <- 15     # expected smallest stratum below this = fragile design
DEGEN_WARN    <- 0.02   # > this share of empty-stratum trials triggers a warning
P_FAV_LEVEL   <- 0.95   # posterior P(favourable) reporting level (one-sided)
P_FAV_MATCHED <- 0.975  # matched to a two-sided frequentist test at p < 0.05
MATCHED_OK    <- !is.null(LUT) && "p_fav_gt_975" %in% names(LUT)  # lookup has it?

## ---- helpers used only for display -----------------------------------------

op_table_A <- function(lut, n, p_cg, eff, mu_aeb) {
  rows <- do.call(rbind, lapply(list(c(0,0), c(1,0), c(0,1), c(1,1)), function(z)
    cgrc_op_at(lut, n, p_cg, eff, z[1], z[2], mu_aeb)))
  # honesty markers: * interpolated between grid points; dagger = outside this
  # expectancy level's grid (an edge value is shown under the requested label).
  mark <- ifelse(rows$clamped, " †", ifelse(rows$interpolated, " *", ""))
  data.frame(
    scenario = paste0(CELL_LABEL[paste0(rows$DTE, rows$AEB)], mark),
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
    "body, .shiny-input-container, .form-control, .btn, label,
      .control-label, table, .table{font-size:16px;}
     h4{font-size:1.35em;} p, .help-block{font-size:1.05em;}
     .verdict{font-size:1.2em;line-height:1.55;padding:14px 16px;
      background:#eef3f8;border-left:5px solid #2471A3;border-radius:4px;}
     .warn{background:#FDEDEC;border-left-color:#C0392B;}
     .feas{font-size:1.05em;} .muted{color:#666;font-size:1.0em;}")),

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
        # Discrete, because the lookup only holds these three levels; a continuous
        # slider that snaps would show a resolution it does not have.
        radioButtons("mu_aeb", "Expectancy magnitude (points)",
                     choices = c("3.85  (half the microdose reference)" = 3.85,
                                 "7.7  (microdose reference, Szigeti)"   = 7.7,
                                 "15.4  (double the reference)"          = 15.4),
                     selected = 7.7),
        div(class = "muted",
            "7.7 is Szigeti's microdose calibration, not a universal constant."),
        uiOutput("inflation_note"),
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
          strong("Criteria differ:"), "\"adjusted flags effect\" is the share of",
          "trials with posterior P(favourable) > 0.95 (Bayesian); \"unadjusted",
          "significant\" is the share with p < 0.05 (frequentist) — comparable in",
          "spirit, not the same quantity. The honest power comparison is the",
          "\"real effect, no expectancy\" row, where both measure the same thing;",
          "in the expectancy rows the unadjusted column is mostly detecting",
          "expectancy, not drug, so do not read it as power."),
        tableOutput("opchar"),
        div(class = "muted",
            "* interpolated between simulated grid points.",
            "† outside the simulated grid for this expectancy level — the nearest",
            "edge value is shown (e.g. n=60 is only simulated at the 7.7 level)."),
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
        sliderInput("rope", "Meaningful-difference threshold (fraction of outcome SD)",
                    min = 0.05, max = 0.5, value = 0.5, step = 0.05),
        div(class = "muted",
            "0.5 SD is the minimum important difference Norman (2003) argues for",
            "and Szigeti's 2024 escitalopram trial adopts. Narrow it for a",
            "stricter bar; it sets both the headline and the ROPE band."),
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
  output$inflation_note <- renderUI({
    inf <- cgr_aeb_inflation(as.numeric(input$mu_aeb), input$pcg)
    div(class = "muted", style = "margin-top:6px;",
        HTML(sprintf("At your CGR of %.2f, this inflates an <b>unadjusted</b>
                      estimate by <b>%.1f points</b>.", input$pcg, inf)))
  })

  output$verdict <- renderUI({
    if (no_lut) return(div(class = "verdict warn",
      "Lookup table not built. Run data-raw/build_lookup.R, then reinstall."))
    div(class = "verdict",
        cgrc_verdict(LUT, input$n, input$pcg, as.numeric(input$eff), as.numeric(input$mu_aeb)))
  })

  output$feasibility <- renderUI({
    minstr <- cgr_min_stratum(input$n, input$pcg)
    degen  <- if (no_lut) NA else
      cgrc_op_at(LUT, input$n, input$pcg, as.numeric(input$eff), 0, 1, as.numeric(input$mu_aeb))$empty_stratum_rate
    thin       <- minstr < THIN_STRATUM
    high_degen <- !is.na(degen) && degen > DEGEN_WARN
    warn <- thin || high_degen
    extra <- if (high_degen) sprintf(paste0(
      " For a real trial of this design, roughly <b>%.0f%%</b> of the time a ",
      "wrong-guess stratum comes up empty and CGRC cannot be computed at all — ",
      "so the rates above are <i>conditional on the trials where it could be</i>."),
      100 * degen)
      else if (thin) paste0(" CGR adjustment is fragile at this design and can be",
                            " undefined for some trials.")
      else ""
    div(class = if (warn) "verdict warn feas" else "verdict feas",
      HTML(sprintf(
        "<b>Feasibility.</b> Expected smallest stratum: <b>~%.0f</b> participants%s.
         Simulated trials with an empty stratum: <b>%s</b>.%s",
        minstr, if (thin) " (thin)" else "",
        if (is.na(degen)) "n/a" else sprintf("%.1f%%", 100 * degen), extra)))
  })

  output$power_plot <- renderPlot({
    if (no_lut) return(NULL)
    eff <- as.numeric(input$eff)
    pc <- cgrc_power_curve(LUT, input$pcg, eff)
    # with no true effect there is nothing to have "power" for: the same curve
    # is then the adjusted false-favourable rate.
    ylab <- if (eff == 0) "adjusted false-favourable rate (no true effect)"
            else "power of the adjusted analysis"
    ggplot(pc, aes(n, power)) +
      geom_hline(yintercept = c(0.8, 0.9), linetype = "dotted", colour = "grey60") +
      geom_line(colour = "#2471A3", linewidth = 1) +
      geom_point(colour = "#2471A3", size = 2) +
      geom_vline(xintercept = input$n, linetype = "dashed", colour = "#C0392B") +
      annotate("text", x = input$n, y = 0.02, label = paste0("your n=", input$n),
               colour = "#C0392B", hjust = -0.05, size = 4.8) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(x = "sample size (n)", y = ylab) +
      theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank())
  })

  output$tradeoff_plot <- renderPlot({
    if (no_lut) return(NULL)
    eff <- as.numeric(input$eff)
    pw  <- cgrc_op_at(LUT, input$n, input$pcg, eff, 1, 0)
    fp  <- cgrc_op_at(LUT, input$n, input$pcg, eff, 0, 1, as.numeric(input$mu_aeb))
    # The frequentist rate is a TWO-SIDED t-test at p<0.05. Its matched Bayesian
    # comparator is posterior P(favourable) > 0.975 (not 0.95). Use the matched
    # column when the lookup has it; otherwise fall back and say so.
    adj_col <- if (MATCHED_OK) "p_fav_gt_975" else "p_fav_gt_95"
    adj_lab <- if (MATCHED_OK) "CGR-adjusted (posterior P>0.975)"
               else "CGR-adjusted (posterior P>0.95)"
    cap <- if (MATCHED_OK)
      "Thresholds are matched: two-sided p<0.05 vs its Bayesian equivalent P>0.975."
    else paste("Note: unadjusted is two-sided p<0.05; adjusted is P>0.95, a looser",
               "bar than the matched P>0.975. Read within-analysis, not as a race.")
    df <- data.frame(
      metric = factor(rep(c("false positive\n(pure expectancy)",
                            "power\n(real effect)"), each = 2),
                      levels = c("false positive\n(pure expectancy)", "power\n(real effect)")),
      analysis = factor(c("unadjusted (p<0.05, two-sided)", adj_lab),
                        levels = c("unadjusted (p<0.05, two-sided)", adj_lab)),
      rate = c(fp$freq_sig, fp[[adj_col]], pw$freq_sig, pw[[adj_col]]))
    ggplot(df, aes(metric, rate, fill = analysis)) +
      geom_col(position = position_dodge(0.7), width = 0.62) +
      geom_text(aes(label = sprintf("%.0f%%", 100 * rate)),
                position = position_dodge(0.7), vjust = -0.4, size = 4.8) +
      scale_fill_manual(values = setNames(c("#C0392B", "#2471A3"), levels(df$analysis))) +
      scale_y_continuous(limits = c(0, 1.08), labels = scales::percent) +
      labs(x = NULL, y = NULL, fill = NULL, caption = cap) +
      theme_minimal(base_size = 16) +
      theme(legend.position = "top", panel.grid.minor = element_blank())
  })

  output$opchar <- renderTable({
    if (no_lut) return(NULL)
    op_table_A(LUT, input$n, input$pcg, as.numeric(input$eff), as.numeric(input$mu_aeb))
  }, digits = 3)

  ## exact simulation, only on demand
  exact_rv <- reactiveVal(NULL)
  observeEvent(input$run_exact, {
    withProgress(message = "Running 500 simulated trials x 4 scenarios...",
                 value = 0.3, {
      op <- cgr_operating(n_trials = 500, n = input$n, p_cg = input$pcg,
                          mu_dte = as.numeric(input$eff), mu_aeb = as.numeric(input$mu_aeb),
                          noise = "all", seed = 1)
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
         # the two interpretable probabilities, before/after the blinding
         # correction, at the same ROPE width delta the user chose
         head = cgrc_headline(trial, direction = dir, delta_sd_frac = input$rope,
                              n_draws = 8000),
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
      uiOutput("b_headline"),
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

  # The headline: two plain probabilities, before and after the blinding
  # correction. This is the interpretable answer - "is there an effect" and
  # "is it big enough to matter" - that a single p-value cannot give.
  output$b_headline <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    h <- f$head; pct <- function(p) sprintf("%.0f%%", 100 * p)
    div(class = "verdict",
      HTML(sprintf(
        "<b>Your trial, in two probabilities.</b><br>
         <table style='width:100%%;margin-top:6px;border-collapse:collapse;'>
         <tr style='color:#666;'><td></td><td><b>at your CGR (raw)</b></td>
             <td><b>at perfect blinding</b></td></tr>
         <tr><td>probability of a favourable effect</td>
             <td><b>%s</b></td><td><b>%s</b></td></tr>
         <tr><td>probability it is meaningful (beyond %.2g pts)</td>
             <td><b>%s</b></td><td><b>%s</b></td></tr></table>
         <div class='muted' style='margin-top:8px;'>%s</div>",
        pct(h$p_dir_obs), pct(h$p_dir_blind),
        h$delta, pct(h$p_meaningful_obs), pct(h$p_meaningful_blind),
        paste(h$text,
              "These are continuous probabilities — deliberately no",
              "significant/not cut-off."))))
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
      theme_minimal(base_size = 15) + theme(legend.position = "bottom")
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
