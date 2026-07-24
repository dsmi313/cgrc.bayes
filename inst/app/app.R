# Shiny app for cgrc.bayes: "How reliable is CGR adjustment for my trial design?"
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
# The UNKNOWN-preserving design lookup (six-stratum). NULL until built with
# data-raw/build_unknown_lookup.R; the app then falls back to binary-only.
LUT_U <- tryCatch(cgrc_unknown_lookup(), error = function(e) NULL)
U_MAX <- if (is.null(LUT_U)) 0 else max(LUT_U$u)     # slider cap = highest built u
# effect sizes offered = those simulated at EVERY expectancy level, so no cell
# is missing when the expectancy slider moves.
EFFS <- if (is.null(LUT)) c(0, 1.5, 3) else if ("mu_aeb" %in% names(LUT))
  # unique(): with a single mu_aeb level (e.g. a partial build) Reduce() returns
  # that one group verbatim, duplicates and all - dedupe so `eff` stays a short list.
  sort(unique(Reduce(intersect, split(LUT$true_effect, LUT$mu_aeb)))) else
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

# The operating-characteristics table shows BOTH Bayesian flag levels, clearly
# labelled: P>0.95 (the standard one-sided reporting flag) and P>0.975 (the rough
# Bayesian-tail comparator to a direction-filtered two-sided p<0.05). The
# frequentist column is a two-sided p<0.05 filtered to the favourable direction
# (~0.025 in that tail under the null), so it is compared like for like against
# the one-tail Bayesian flags. `matched` guards older lookups without P>0.975.
op_table_A <- function(lut, n, p_cg, eff, mu_aeb) {
  rows <- do.call(rbind, lapply(list(c(0,0), c(1,0), c(0,1), c(1,1)), function(z)
    cgrc_op_at(lut, n, p_cg, eff, z[1], z[2], mu_aeb)))
  # honesty markers: * interpolated between grid points; dagger = outside this
  # expectancy level's grid (an edge value is shown under the requested label).
  mark <- ifelse(rows$clamped, " †", ifelse(rows$interpolated, " *", ""))
  matched <- "p_fav_gt_975" %in% names(rows)
  out <- data.frame(
    scenario = paste0(CELL_LABEL[paste0(rows$DTE, rows$AEB)], mark),
    `true effect` = round(ifelse(rows$DTE == 1, eff, 0), 2),
    `adjusted bias` = round(rows$adj_bias, 2),
    `95% coverage` = round(rows$coverage95, 3),
    `adjusted flags (P>0.95)` = round(rows$p_fav_gt_95, 3),
    check.names = FALSE)
  if (matched) out[["matched flag (P>0.975)"]] <- round(rows$p_fav_gt_975, 3)
  out[["unadjusted significant (direction-filtered p<0.05)"]] <- round(rows$freq_sig, 3)
  out
}

# Same four-scenario operating-characteristics table for the UNKNOWN estimator,
# read from the UNKNOWN lookup (u snapped to its nearest grid level).
op_table_A_unknown <- function(lut, n, p_cg, eff, u) {
  rows <- do.call(rbind, lapply(list(c(0,0), c(1,0), c(0,1), c(1,1)), function(z)
    cgrc_unknown_op_at(lut, n, p_cg, eff, u, z[1], z[2])))
  mark <- ifelse(rows$clamped, " †", ifelse(rows$interpolated, " *", ""))
  matched <- "p_fav_gt_975" %in% names(rows)
  out <- data.frame(
    scenario = paste0(CELL_LABEL[paste0(rows$DTE, rows$AEB)], mark),
    `true effect` = round(ifelse(rows$DTE == 1, eff, 0), 2),
    `adjusted bias` = round(rows$adj_bias, 2),
    `95% coverage` = round(rows$coverage95, 3),
    `adjusted flags (P>0.95)` = round(rows$p_fav_gt_95, 3),
    check.names = FALSE)
  if (matched) out[["matched flag (P>0.975)"]] <- round(rows$p_fav_gt_975, 3)
  out[["unadjusted significant (direction-filtered p<0.05)"]] <- round(rows$freq_sig, 3)
  out
}

## ---- UI ---------------------------------------------------------------------

ui <- navbarPage(
  "CGRC — how reliable is adjustment for my trial design?",
  id = "navbar",
  header = tagList(
    tags$style(HTML(
    "body{font-size:18px;}
     .shiny-input-container, label, .control-label, .selectize-input,
       .selectize-dropdown, .irs, .radio label, input, select, textarea{font-size:18px;}
     .form-control, .selectize-input{font-size:17px;padding:8px 10px;}
     .btn{font-size:17px;padding:8px 14px;margin-top:4px;}
     .radio, .form-group{margin-bottom:14px;}
     h3{font-size:26px;} h4{font-size:23px;margin-top:0.4em;}
     p, .help-block{font-size:16px;} table, .table{font-size:17px;}
     .table td, .table th{padding:6px 8px;}
     .verdict{font-size:20px;line-height:1.55;padding:16px 18px;
      background:#eef3f8;border-left:5px solid #2471A3;border-radius:4px;}
     .warn{background:#FDEDEC;border-left-color:#C0392B;}
     .caution{background:#FEF9E7;border-left-color:#B9770E;}
     .feas{font-size:17px;}
     .reliability{font-size:20px;font-weight:600;padding:10px 16px;border-radius:4px;
      display:inline-block;margin-bottom:6px;}
     .rel-ok{background:#E9F7EF;border-left:6px solid #1E8449;color:#145A32;}
     .rel-caution{background:#FEF9E7;border-left:6px solid #B9770E;color:#7E5109;}
     .rel-warn{background:#FDEDEC;border-left:6px solid #C0392B;color:#922B21;}
     .done{background:#E9F7EF;border-left:5px solid #1E8449;padding:10px 14px;
      border-radius:4px;font-size:17px;}
     .muted{color:#555;font-size:16px;}
     .sidebar-wide .well{min-width:300px;}
     img.shiny-plot-output, .shiny-plot-output{max-width:100%;}")),
    # Disable the exact-simulation button the instant it is clicked (client-side,
    # no server round-trip) and re-enable it when the server signals completion.
    tags$script(HTML(
    "$(document).on('click','#run_exact,#run_unknown_design',function(){
       $(this).prop('disabled',true).addClass('disabled').data('label',$(this).html());
       $(this).html('Running...');});
     Shiny.addCustomMessageHandler('btn_reenable',function(id){
       var b=$('#'+id); if(b.data('label')) b.html(b.data('label'));
       b.prop('disabled',false).removeClass('disabled');});"))),

  ## ===== Panel A: Design ====================================================
  tabPanel(
    "Design",
    sidebarLayout(
      sidebarPanel(
        width = 4,
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
        if (!is.null(LUT_U)) tagList(
          tags$hr(),
          sliderInput("u_rate", "Expected UNKNOWN-response rate", min = 0,
                      max = U_MAX, value = 0, step = 0.05),
          div(class = "muted",
              "0 = the original binary design tool. Above 0 switches to the",
              "UNKNOWN-preserving six-stratum design: the correct-guess rate above",
              "is then the DIRECTIONAL rate (among AC/PL responders) and this",
              "UNKNOWN rate is held fixed. Snapped to the nearest simulated level.")),
        tags$hr(),
        helpText(class = "muted",
          "Curves and tables are read from a precomputed simulation grid and",
          "interpolated. For the exact numbers at these settings, run the",
          "simulation below. More trials = less Monte Carlo noise but a longer",
          "wait (each trial fits 4 posteriors; ~10 s per 500 trials)."),
        sliderInput("n_trials", "Simulated trials", min = 10, max = 1000,
                    value = 500, step = 10),
        actionButton("run_exact", "Run exact simulation",
                     class = "btn-primary btn-sm")
      ),
      mainPanel(
        width = 8,
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
        h4("Trade-off using the stricter matched threshold, P(favourable) > 0.975"),
        p(class = "muted",
          strong("The values differ from the summary above because this plot uses",
                 "P>0.975 rather than the standard P>0.95 reporting threshold."),
          "\"False treatment attribution\" is a real observed arm difference driven",
          "by expectancy when the direct treatment effect is zero — the unadjusted",
          "analysis attributes it to the drug. This is the one plot that uses the",
          "stricter, approximately direction-matched comparator (posterior P>0.975",
          "vs a direction-filtered two-sided p<0.05); see the note under the plot."),
        plotOutput("tradeoff_plot", height = "260px"),
        br(),
        h4("Operating characteristics at your settings"),
        p(class = "muted",
          strong("Two Bayesian flags are shown, both labelled."),
          "\"adjusted flags (P>0.95)\" is the standard one-sided Bayesian flag;",
          "\"matched flag (P>0.975)\" is the stricter level shown only as a rough",
          "comparator to the \"unadjusted significant\" column — a two-sided p<0.05",
          "filtered to the prespecified favourable direction (≈0.025 in that tail",
          "under the null). Neither threshold is universally correct — 0.95 is the reporting",
          "default, 0.975 is the comparator. The Bayesian and frequentist columns use",
          "different estimators and are not inferentially identical. The honest power",
          "comparison is the \"real effect, no expectancy\" row; in the expectancy",
          "rows the unadjusted column is mostly detecting expectancy, not drug."),
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
        width = 4,
        fileInput("csv", "Upload trial CSV", accept = ".csv"),
        helpText(class = "muted", "One row per participant."),
        selectInput("col_cond", "Column: treatment received", choices = NULL),
        selectInput("col_guess", "Column: treatment guessed", choices = NULL),
        selectInput("col_value", "Column: outcome value", choices = NULL),
        textInput("unknown_level", "Extra label meaning \"I do not know\"",
                  value = "UNKNOWN"),
        div(class = "muted",
            "UNKNOWN, unsure, uncertain, \"don't know\", DK/IDK are recognised",
            "automatically; add your own token here if it differs."),
        uiOutput("unknown_mode_ui"),
        radioButtons("direction", "Favourable direction",
                     c("higher is better" = "1", "lower is better" = "-1"), "1"),
        radioButtons("threshold_mode", "Meaningful-difference threshold",
                     c("fraction of outcome SD" = "sd", "outcome units" = "units"), "sd"),
        conditionalPanel("input.threshold_mode == 'sd'",
          numericInput("rope", "Threshold (fraction of outcome SD)",
                       value = 0.5, min = 0.01, step = 0.05)),
        conditionalPanel("input.threshold_mode == 'units'",
          numericInput("rope_units", "Threshold (outcome units)",
                       value = NA, min = 0, step = 0.1)),
        div(class = "muted",
            "0.5 SD is a common minimum important difference (Norman 2003;",
            "Szigeti 2024). No upper cap — widen it if your field's meaningful",
            "difference is larger. It sets both the headline and the ROPE band."),
        numericInput("seed_b", "Random seed", value = 1, min = 0, step = 1),
        actionButton("analyse", "Analyse", class = "btn-primary"),
        uiOutput("download_ui"),
        uiOutput("to_design")
      ),
      mainPanel(
        width = 8,
        uiOutput("bpanel")
      )
    )
  )
)

## ---- server -----------------------------------------------------------------

server <- function(input, output, session) {

  no_lut <- is.null(LUT)

  ## UNKNOWN-design mode is on when the UNKNOWN lookup exists and u > 0. At u = 0
  ## every Panel A output is exactly the original binary design tool.
  u_on <- reactive(!is.null(LUT_U) && !is.null(input$u_rate) && isTRUE(input$u_rate > 0))

  ## ---- Panel A ----
  output$inflation_note <- renderUI({
    inf <- cgr_aeb_inflation(as.numeric(input$mu_aeb), input$pcg)
    div(class = "muted", style = "margin-top:6px;",
        HTML(sprintf("At your CGR of %.2f, this inflates an <b>unadjusted</b>
                      estimate by <b>%.1f points</b>.", input$pcg, inf)))
  })

  # Shared feasibility numbers, so the reliability badge and the feasibility
  # readout can never disagree.
  feas_numbers <- reactive({
    if (u_on()) {
      minstr <- cgr_unknown_min_stratum(input$n, input$pcg, input$u_rate)
      degen  <- cgrc_unknown_op_at(LUT_U, input$n, input$pcg, as.numeric(input$eff),
                                   input$u_rate, 0, 1)$empty_stratum_rate
    } else {
      minstr <- cgr_min_stratum(input$n, input$pcg)
      degen  <- if (no_lut) NA_real_ else
        cgrc_op_at(LUT, input$n, input$pcg, as.numeric(input$eff), 0, 1,
                   as.numeric(input$mu_aeb))$empty_stratum_rate
    }
    list(minstr = minstr, degen = degen)
  })

  output$verdict <- renderUI({
    if (no_lut) return(div(class = "verdict warn",
      "Lookup table not built. Run data-raw/build_lookup.R, then reinstall."))
    fn  <- feas_numbers()
    rel <- cgrc_reliability(fn$minstr, fn$degen, thin = THIN_STRATUM,
                            degen_warn = DEGEN_WARN)
    txt <- if (u_on())
      cgrc_unknown_verdict(LUT_U, input$n, input$pcg, as.numeric(input$eff), input$u_rate)
    else
      cgrc_verdict(LUT, input$n, input$pcg, as.numeric(input$eff), as.numeric(input$mu_aeb))
    tagList(
      div(class = paste0("reliability rel-", rel$class), rel$category),
      div(class = "muted", style = "margin:-2px 0 6px;",
          "Feasibility only (smallest expected stratum and empty-stratum rate).",
          "Bias, coverage and power are in the operating-characteristics table below."),
      div(class = "verdict", HTML(txt)))
  })

  output$feasibility <- renderUI({
    fn <- feas_numbers(); minstr <- fn$minstr; degen <- fn$degen
    # all four (or six) EXPECTED stratum sizes, not just the smallest, so an arm
    # imbalance is never hidden behind a single number.
    if (u_on()) {
      u <- input$u_rate; n <- input$n; p <- input$pcg
      es <- c("received active / guessed active"   = 0.5 * n * (1 - u) * p,
              "received active / guessed placebo"   = 0.5 * n * (1 - u) * (1 - p),
              "received active / UNKNOWN"            = 0.5 * n * u,
              "received placebo / guessed active"   = 0.5 * n * (1 - u) * (1 - p),
              "received placebo / guessed placebo"  = 0.5 * n * (1 - u) * p,
              "received placebo / UNKNOWN"           = 0.5 * n * u)
      strata_word <- "smallest of six strata"
    } else {
      es4 <- cgr_expected_strata(input$n, input$pcg)
      es <- c("received active / guessed active"  = es4[["ACAC"]],
              "received active / guessed placebo" = es4[["ACPL"]],
              "received placebo / guessed active" = es4[["PLAC"]],
              "received placebo / guessed placebo"= es4[["PLPL"]])
      strata_word <- "smallest stratum"
    }
    strata_rows <- paste(sprintf(
      "<tr><td>%s</td><td style='text-align:right'><b>~%.0f</b></td></tr>",
      names(es), es), collapse = "")
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
    tagList(
      div(class = if (warn) "verdict warn feas" else "verdict feas",
        HTML(sprintf(
          "<b>Feasibility.</b> Expected %s: <b>~%.0f</b> participants%s.
           Simulated trials with an empty stratum: <b>%s</b>.%s
           <table style='margin-top:8px;border-collapse:collapse;'>
           <tr style='color:#555;'><td><b>expected stratum</b></td>
               <td style='text-align:right'><b>size</b></td></tr>%s</table>",
          strata_word, minstr, if (thin) " (thin)" else "",
          if (is.na(degen)) "n/a" else sprintf("%.1f%%", 100 * degen), extra,
          strata_rows))),
      tags$details(style = "margin-top:8px;",
        tags$summary(class = "muted",
          "Assumptions behind this design simulation (click to expand)"),
        tags$ul(class = "muted", if (u_on()) tagList(
          tags$li("Balanced random assignment (each arm with probability 0.5)."),
          tags$li("The selected DIRECTIONAL correct-guess rate applies symmetrically among directional (AC/PL) responders; the UNKNOWN rate is equal across arms and held fixed."),
          tags$li("The UNKNOWN-aware generative model is used: UNKNOWN responders carry no expectancy (B_TE = 0)."),
          tags$li("The target is the direct treatment effect (the reweighted estimate at a 50% directional guess rate, UNKNOWN rate held fixed)."),
          tags$li("Both directional guess classes must be present, and the UNKNOWN class when the UNKNOWN rate is above zero; structural-zero cells are permitted (the estimand errors only when a required class is absent)."),
          tags$li("Results are conditional on the selected true-effect and expectancy magnitudes.")
        ) else tagList(
          tags$li("Balanced random assignment (each arm with probability 0.5)."),
          tags$li("The selected correct-guess rate applies symmetrically to both arms unless otherwise modelled."),
          tags$li("The activated-expectancy-bias (AEB) generative model is used."),
          tags$li("The target is the direct treatment effect (the CGR-adjusted estimate at a 50% guess rate)."),
          tags$li("All four treatment-by-guess strata must be nonempty for the estimand to be defined."),
          tags$li("Results are conditional on the selected true-effect and expectancy magnitudes.")))))
  })

  output$power_plot <- renderPlot({
    if (no_lut) return(NULL)
    eff <- as.numeric(input$eff)
    pc <- if (u_on()) cgrc_unknown_power_curve(LUT_U, input$pcg, eff, input$u_rate)
          else cgrc_power_curve(LUT, input$pcg, eff)
    # with no true effect there is nothing to have "power" for: the same curve
    # is then the adjusted false-favourable rate.
    ylab <- if (eff == 0) "adjusted false-favourable rate (no true effect)"
            else "power of the adjusted analysis"
    g <- ggplot(pc, aes(n, power)) +
      geom_hline(yintercept = c(0.8, 0.9), linetype = "dotted", colour = "grey60") +
      geom_line(colour = "#2471A3", linewidth = 1) +
      geom_point(colour = "#2471A3", size = 2) +
      geom_vline(xintercept = input$n, linetype = "dashed", colour = "#C0392B") +
      annotate("text", x = input$n, y = 0.02, label = paste0("your n=", input$n),
               colour = "#C0392B", hjust = -0.05, size = 4.8) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(x = "sample size (n)", y = ylab) +
      theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank())
    # overlay the exact-simulation point at your n, when the run still matches
    em <- exact_match()
    if (!is.null(em)) {
      opp <- em$op[em$op$DTE == 1 & em$op$AEB == 0, ]
      g <- g + geom_point(data = data.frame(n = em$n, power = opp$p_fav_gt_95),
                          aes(n, power), colour = "#E67E22", size = 4.5, shape = 18) +
        annotate("text", x = em$n, y = opp$p_fav_gt_95, label = "exact",
                 colour = "#E67E22", vjust = -1, hjust = -0.1, size = 4.6) +
        labs(subtitle = "orange ◆ = exact simulation at your n")
    }
    g
  })

  output$tradeoff_plot <- renderPlot({
    if (no_lut) return(NULL)
    eff <- as.numeric(input$eff)
    if (u_on()) {
      pw <- cgrc_unknown_op_at(LUT_U, input$n, input$pcg, eff, input$u_rate, 1, 0)
      fp <- cgrc_unknown_op_at(LUT_U, input$n, input$pcg, eff, input$u_rate, 0, 1)
    } else {
      pw <- cgrc_op_at(LUT, input$n, input$pcg, eff, 1, 0)
      fp <- cgrc_op_at(LUT, input$n, input$pcg, eff, 0, 1, as.numeric(input$mu_aeb))
    }
    # The frequentist rate is now the FAVOURABLE-TAIL t-test rate (p<0.05 in the
    # favourable direction). Its approximately matched Bayesian comparator is
    # posterior P(favourable) > 0.975. They use different estimators and are not
    # inferentially identical - this is the one plot that uses the 0.975 level.
    matched <- "p_fav_gt_975" %in% names(pw)
    adj_col <- if (matched) "p_fav_gt_975" else "p_fav_gt_95"
    adj_lab <- if (matched) "CGR-adjusted (posterior P>0.975)"
               else "CGR-adjusted (posterior P>0.95)"
    cap <- if (matched)
      paste("Rough tail comparator: two-sided p<0.05 filtered to the favourable",
            "direction (≈0.025 in that tail under the null) vs posterior",
            "P(favourable)>0.975. Different estimators — not inferentially identical.")
    else paste("Note: unadjusted is a direction-filtered two-sided p<0.05; adjusted",
               "is P>0.95, a looser bar than the matched P>0.975.")
    mlev <- c("false treatment attribution\n(pure expectancy)", "power\n(real effect)")
    df <- data.frame(
      metric = factor(rep(mlev, each = 2), levels = mlev),
      analysis = factor(c("unadjusted (direction-filtered p<0.05)", adj_lab),
                        levels = c("unadjusted (direction-filtered p<0.05)", adj_lab)),
      rate = c(fp$freq_sig, fp[[adj_col]], pw$freq_sig, pw[[adj_col]]))
    g <- ggplot(df, aes(metric, rate, fill = analysis)) +
      geom_col(position = position_dodge(0.7), width = 0.62) +
      geom_text(aes(label = sprintf("%.0f%%", 100 * rate)),
                position = position_dodge(0.7), vjust = -0.4, size = 4.8) +
      scale_fill_manual(values = setNames(c("#C0392B", "#2471A3"), levels(df$analysis))) +
      scale_y_continuous(limits = c(0, 1.12), labels = scales::percent) +
      labs(x = NULL, y = NULL, fill = NULL, caption = cap) +
      theme_minimal(base_size = 16) +
      theme(legend.position = "top", panel.grid.minor = element_blank())
    em <- exact_match()
    if (!is.null(em)) {
      fpe <- em$op[em$op$DTE == 0 & em$op$AEB == 1, ]
      pwe <- em$op[em$op$DTE == 1 & em$op$AEB == 0, ]
      ecol <- if (matched && "p_fav_gt_975" %in% names(em$op)) "p_fav_gt_975" else "p_fav_gt_95"
      edf <- data.frame(metric = df$metric, analysis = df$analysis,
        rate = c(fpe$freq_sig, fpe[[ecol]], pwe$freq_sig, pwe[[ecol]]))
      g <- g + geom_point(data = edf, aes(metric, rate, group = analysis),
                          position = position_dodge(0.7), shape = 18, size = 4.5,
                          colour = "#E67E22") +
        labs(subtitle = "orange ◆ = exact simulation at these settings")
    }
    g
  })

  output$opchar <- renderTable({
    if (no_lut) return(NULL)
    if (u_on())
      op_table_A_unknown(LUT_U, input$n, input$pcg, as.numeric(input$eff), input$u_rate)
    else
      op_table_A(LUT, input$n, input$pcg, as.numeric(input$eff), as.numeric(input$mu_aeb))
  }, digits = 3)

  ## exact simulation, only on demand. The full settings are captured at click
  ## time so the readout, the plot markers and the lookup-vs-exact comparison all
  ## describe the SAME run. Seed is fixed so the run is reproducible.
  EXACT_SEED <- 1L
  exact_rv <- reactiveVal(NULL)
  observeEvent(input$run_exact, {
    nt <- input$n_trials; unk <- u_on()
    n <- input$n; pcg <- input$pcg; eff <- as.numeric(input$eff)
    mu_aeb <- as.numeric(input$mu_aeb); u <- if (unk) input$u_rate else 0
    withProgress(message = sprintf("Running %d %ssimulated trials x 4 scenarios...",
                                   nt, if (unk) "UNKNOWN " else ""), value = 0.3, {
      op <- if (unk)
        cgr_unknown_operating(n_trials = nt, n = n, p_cg = pcg, u = u, mu_dte = eff,
                              mu_aeb = mu_aeb, noise = "all", seed = EXACT_SEED)
      else
        cgr_operating(n_trials = nt, n = n, p_cg = pcg, mu_dte = eff,
                      mu_aeb = mu_aeb, noise = "all", seed = EXACT_SEED)
      incProgress(0.7)
      exact_rv(list(op = op, n = n, pcg = pcg, eff = eff, mu_aeb = mu_aeb, u = u,
                    mode = if (unk) "unknown" else "binary",
                    n_trials = nt, seed = EXACT_SEED, time = Sys.time()))
    })
    session$sendCustomMessage("btn_reenable", "run_exact")  # visible completion
  })

  ## The stored exact run only decorates the plots while its settings still match
  ## the current sliders - otherwise the orange markers would misrepresent a run
  ## computed elsewhere. NULL means "do not draw markers"; the comparison table
  ## below always shows, stamped with the settings it was computed at.
  exact_match <- reactive({
    e <- exact_rv(); if (is.null(e)) return(NULL)
    same <- isTRUE(all.equal(e$n, input$n)) && isTRUE(all.equal(e$pcg, input$pcg)) &&
      isTRUE(all.equal(e$eff, as.numeric(input$eff))) &&
      e$mode == (if (u_on()) "unknown" else "binary") &&
      (e$mode == "binary" || isTRUE(all.equal(e$u, input$u_rate))) &&
      (e$mode == "unknown" || isTRUE(all.equal(e$mu_aeb, as.numeric(input$mu_aeb))))
    if (same) e else NULL
  })

  output$exact_out <- renderUI({
    e <- exact_rv(); if (is.null(e)) return(NULL)
    op <- e$op
    # interpolated-lookup values at the SAME settings, per scenario
    look_pfav <- vapply(seq_len(nrow(op)), function(i) {
      r <- if (e$mode == "unknown")
        cgrc_unknown_op_at(LUT_U, e$n, e$pcg, e$eff, e$u, op$DTE[i], op$AEB[i])
      else cgrc_op_at(LUT, e$n, e$pcg, e$eff, op$DTE[i], op$AEB[i], e$mu_aeb)
      r$p_fav_gt_95
    }, numeric(1))
    cmp <- data.frame(
      scenario = CELL_LABEL[paste0(op$DTE, op$AEB)],
      `interpolated lookup` = round(look_pfav, 3),
      `exact simulation`    = round(op$p_fav_gt_95, 3),
      `abs. difference`     = round(abs(op$p_fav_gt_95 - look_pfav), 3),
      `valid sims`          = op$n_valid,
      `empty-stratum rate`  = paste0(round(100 * op$empty_stratum_rate, 1), "%"),
      check.names = FALSE)
    settings_txt <- sprintf("n=%d, CGR=%.2f, effect=%.1f%s", e$n, e$pcg, e$eff,
      if (e$mode == "unknown") sprintf(", UNKNOWN=%.0f%%", 100 * e$u)
      else sprintf(", expectancy=%.1f", e$mu_aeb))
    tagList(
      br(),
      div(class = "done", HTML(sprintf(
        "✓ Exact simulation completed at <b>%s</b> — <b>%d</b> trials per scenario, seed <b>%d</b>.",
        format(e$time, "%H:%M:%S"), e$n_trials, e$seed))),
      h4("Interpolated lookup vs exact simulation"),
      div(class = "muted", HTML(sprintf(paste(
        "Both columns are the adjusted-flags rate P(favourable)>0.95. The curves and",
        "the trade-off plot above are read from the <b>interpolated lookup</b>; the",
        "orange markers (shown while your settings still match this run) are the",
        "<b>exact simulation</b> points. Run settings: %s. Small differences are Monte",
        "Carlo noise, not scientific change."), settings_txt))),
      renderTable(cmp, digits = 3),
      h4(sprintf("Exact simulation detail (%d trials, seed %d)", e$n_trials, e$seed)),
      renderTable({
        data.frame(scenario = CELL_LABEL[paste0(op$DTE, op$AEB)],
                   `adj bias` = round(op$adj_bias, 3),
                   `95% coverage` = round(op$coverage95, 3),
                   `adjusted flags (P>0.95)` = round(op$p_fav_gt_95, 3),
                   `matched flag (P>0.975)` =
                     if ("p_fav_gt_975" %in% names(op)) round(op$p_fav_gt_975, 3) else NA_real_,
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

  # Does the chosen guess column contain any observed UNKNOWN response? Peeked
  # reactively so the analysis-mode control can appear before Analyse is clicked.
  unknown_present <- reactive({
    req(input$csv, input$col_guess)
    g <- tryCatch(cgrc_normalise_guess(raw_csv()[[input$col_guess]],
                    allow_unknown = TRUE, unknown_labels = input$unknown_level),
                  error = function(e) NULL)
    !is.null(g) && any(g == "UNKNOWN", na.rm = TRUE)
  })

  # The analysis-mode control only appears when UNKNOWN responses are present.
  # Default is to PRESERVE UNKNOWN as a third category.
  output$unknown_mode_ui <- renderUI({
    if (!isTRUE(unknown_present())) return(NULL)
    div(class = "verdict feas", style = "padding:8px 10px;",
      radioButtons("unknown_mode",
        HTML("<b>UNKNOWN guesses detected.</b> How should they be handled?"),
        c("Preserve UNKNOWN as a third response category" = "preserve",
          "Binary complete-case CGRC (exclude UNKNOWN responses)" = "completecase"),
        selected = "preserve"))
  })

  # The meaningful-difference threshold in outcome units, from either control.
  delta_units <- function(sdy) {
    if (identical(input$threshold_mode, "units") &&
        !is.null(input$rope_units) && is.finite(input$rope_units))
      as.numeric(input$rope_units)
    else as.numeric(input$rope) * sdy
  }

  fit <- eventReactive(input$analyse, {
    d <- raw_csv()
    dir <- as.numeric(input$direction)
    seed <- if (is.null(input$seed_b) || !is.finite(input$seed_b)) NULL else as.integer(input$seed_b)
    ulevel <- input$unknown_level

    # 1. Audit every row instead of silently dropping with complete.cases().
    aud <- cgrc_input_audit(d[[input$col_cond]], d[[input$col_guess]],
                            d[[input$col_value]], unknown_level = ulevel)
    clean <- aud$clean

    # 2. Decide the mode. Preserve UNKNOWN by default when present.
    present <- aud$has_unknown
    mode <- if (!present) "binary"
            else if (identical(input$unknown_mode, "completecase")) "binary"
            else "unknown"

    # In binary complete-case with UNKNOWN present, drop the UNKNOWN rows here and
    # count them, so the exclusion is explicit rather than silent.
    n_excl_unknown <- 0L
    if (mode == "binary" && present) {
      is_unk <- clean$guess == "UNKNOWN"
      n_excl_unknown <- sum(is_unk)
      clean <- clean[!is_unk, , drop = FALSE]
    }
    sdy <- if (nrow(clean)) stats::sd(clean$value) else NA_real_
    delta <- delta_units(sdy)

    base <- list(mode = mode, trial = clean, dir = dir, audit = aud,
                 delta = delta, seed = seed, n_excl_unknown = n_excl_unknown)

    if (mode == "unknown") {
      cur_cgr <- NULL
      ufit  <- cgrc_unknown(clean, unknown_level = ulevel, n_draws = 8000,
                            direction = dir, seed = seed)
      grid  <- sort(unique(c(seq(0, 1, length.out = 101),
                             ufit$observed_directional_cgr, 0.5)))
      c(base, list(
        ufit  = ufit,
        uhead = cgrc_unknown_headline(clean, unknown_level = ulevel, direction = dir,
                                      delta = delta, n_draws = 8000, seed = seed),
        urope = cgr_unknown_rope(clean, grid = grid, n_draws = 8000,
                                 delta = delta, direction = dir),
        usens = cgr_unknown_rope_sensitivity(clean, at_cgr = 0.5, n_draws = 6000,
                                             direction = dir)))
    } else {
      grid <- sort(unique(c(seq(0, 1, length.out = 101),
                            cgr_observed(cgr_strata(clean)))))
      c(base, list(
        fit  = cgrc(clean, n_draws = 8000, direction = dir, seed = seed),
        head = cgrc_headline(clean, direction = dir, delta = delta,
                             n_draws = 8000, seed = seed),
        rope = cgr_rope(clean, grid = grid, n_draws = 8000, delta = delta,
                        direction = dir),
        sens = cgr_rope_sensitivity(clean, at_cgr = 0.5, n_draws = 6000,
                                    direction = dir)))
    }
  })

  output$bpanel <- renderUI({
    if (input$analyse == 0) return(helpText(
      "Upload a CSV with columns for treatment received, treatment guessed and",
      "an outcome value, map them on the left, then click Analyse."))
    tagList(
      uiOutput("b_error"),
      uiOutput("b_audit"),
      uiOutput("b_top_summary"),
      uiOutput("b_guess_rates"),
      h4("Mapped data (first rows)"),
      p(class = "muted", "Check the mapping is right before reading the results."),
      tableOutput("b_preview"),
      uiOutput("b_headline"),
      fluidRow(column(5, h4("Strata (from your data)"), tableOutput("b_strata")),
               column(7, h4("Adjusted vs unadjusted"), tableOutput("b_summary"))),
      uiOutput("b_counts"),
      uiOutput("b_identity"),
      h4("CGR curve"), plotOutput("b_curve", height = "420px"),
      h4("Region of practical equivalence"),
      p(class = "muted", "A ROPE conclusion is only as good as the band width, so",
        "the sensitivity to that width is shown beside it."),
      fluidRow(column(7, plotOutput("b_rope", height = "300px")),
               column(5, tableOutput("b_sens"))),
      uiOutput("unknown_design_out"))
  })

  safe_fit <- reactive(tryCatch(fit(), error = function(e) structure(conditionMessage(e), class = "cgrc_err")))
  output$b_error <- renderUI({
    f <- safe_fit()
    if (inherits(f, "cgrc_err")) div(class = "verdict warn", paste("Could not analyse:", f))
  })

  # Missing-data audit: what was excluded and why, and the exact analysis n. No
  # row is dropped silently; the cleaned data and the exclusion log are
  # downloadable from the sidebar.
  output$b_audit <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    s <- f$audit$summary
    modetxt <- if (f$mode == "unknown")
      "UNKNOWN responses are <b>preserved</b> as a third response category."
    else if (f$n_excl_unknown > 0)
      sprintf(paste0("Binary <b>complete-case</b> analysis: %d UNKNOWN response(s) ",
                     "(%.1f%%) were <b>excluded</b> (they were observed, not missing)."),
              f$n_excl_unknown, 100 * f$n_excl_unknown / s[["n_input"]])
    else "No UNKNOWN responses were present; standard binary CGRC."
    seedtxt <- if (is.null(f$seed)) "not set (results may vary run to run)"
               else sprintf("%d", f$seed)
    div(class = "verdict feas", HTML(sprintf(
      "<b>Input audit.</b> %d rows uploaded → <b>%d analysed</b>. %s
       <div class='muted' style='margin-top:6px;'>Excluded — missing condition: %d;
       missing guess: %d; missing outcome: %d; non-numeric outcome: %d.
       Observed UNKNOWN responses: %d. &nbsp;|&nbsp; Random seed: <b>%s</b>.</div>",
      s[["n_input"]], nrow(f$trial), modetxt,
      s[["missing_condition"]], s[["missing_guess"]], s[["missing_outcome"]],
      s[["nonnumeric_outcome"]], s[["observed_unknown"]], seedtxt)))
  })

  # Compact top-of-panel summary (item 9): the numbers a trialist wants first, in
  # plain language, without turning posterior probabilities into significant/not.
  output$b_top_summary <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    summ <- cgrc_analysis_summary(f)
    rows <- paste(sprintf(
      "<tr><td>%s</td><td style='text-align:right'><b>%s</b></td></tr>",
      summ$quantity, ifelse(is.na(summ$value), "—",
        formatC(summ$value, format = "g", digits = 6))), collapse = "")
    div(class = "verdict",
      HTML(sprintf(
        "<b>At a glance.</b>
         <table style='width:100%%;margin-top:6px;border-collapse:collapse;'>%s</table>
         <div class='muted' style='margin-top:8px;'>The <b>adjusted effect</b> is a
         counterfactual under the CGRC assumptions — what the effect would have been
         at a 50%% correct-guess rate — <i>not</i> automatically the true
         pharmacological effect. Probabilities are continuous; no significant/not
         cut-off is imposed.</div>", rows)))
  })

  # Arm-specific correct-guess rates, allocation and guess-response counts, so a
  # single overall rate can never hide severe arm asymmetry (item 7). Works in
  # both modes. Prominently warns when the smallest observed stratum is thin.
  output$b_guess_rates <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    gr <- cgr_guess_rates(f$trial)
    if (f$mode == "unknown") { st <- cgr_unknown_strata(f$trial); nn <- lengths(st) }
    else { st <- cgr_strata(f$trial); nn <- lengths(st)[STRATA] }
    smallest <- min(nn[nn > 0])
    thin <- smallest < THIN_STRATUM
    asym <- is.finite(gr$active) && is.finite(gr$placebo) &&
            abs(gr$active - gr$placebo) > 0.2
    pct <- function(p) if (is.na(p)) "n/a" else sprintf("%.0f%%", 100 * p)
    warn_txt <- if (thin) sprintf(paste0(
      " <b>Warning:</b> the smallest observed stratum has <b>%d</b> participants ",
      "(&lt; %d) — the reweighted estimate is fragile here."), smallest, THIN_STRATUM)
      else ""
    asym_txt <- if (asym) sprintf(paste0(
      " <b>Arm asymmetry:</b> active-arm and placebo-arm correct-guess rates differ ",
      "by %.0f points; a single overall rate would hide this."),
      100 * abs(gr$active - gr$placebo)) else ""
    ug <- if (gr$guess_unknown > 0) sprintf(" / UNKNOWN %d", gr$guess_unknown) else ""
    div(class = if (thin || asym) "verdict warn feas" else "verdict feas",
      HTML(sprintf(
        "<b>Blinding & allocation.</b>
         Overall correct-guess rate <b>%s</b>;
         active-arm <b>%s</b>; placebo-arm <b>%s</b>.
         Allocation — active <b>%d</b>, placebo <b>%d</b>.
         Guesses — active %d, placebo %d%s.
         Smallest occupied stratum <b>%d</b>.%s%s",
        pct(gr$overall), pct(gr$active), pct(gr$placebo),
        gr$n_active, gr$n_placebo, gr$guess_active, gr$guess_placebo, ug,
        smallest, warn_txt, asym_txt)))
  })

  # A short preview of the mapped analysis frame, so the column mapping can be
  # eyeballed before trusting any downstream number.
  output$b_preview <- renderTable({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    utils::head(f$trial, 8)
  })

  # Downloads: cleaned analysis data and the exclusion log.
  output$download_ui <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    tagList(tags$hr(),
      div(class = "muted", "Downloads"),
      downloadButton("dl_summary","Analysis summary (CSV)",     class = "btn-sm"),
      downloadButton("dl_curve",  "CGR curve (CSV)",            class = "btn-sm"),
      downloadButton("dl_png",    "Primary plot (PNG)",         class = "btn-sm"),
      downloadButton("dl_html",   "Analysis report (HTML)",     class = "btn-sm"),
      downloadButton("dl_clean",  "Cleaned analysis data (CSV)", class = "btn-sm"),
      downloadButton("dl_log",    "Exclusion log (CSV)",        class = "btn-sm"),
      downloadButton("dl_report", "Analysis report (Markdown)", class = "btn-sm"))
  })
  # The curve underlying the CGR figure (one row per grid CGR).
  curve_df <- function(f) if (f$mode == "unknown") f$ufit$curve else f$fit$curve
  output$dl_summary <- downloadHandler(
    filename = function() "cgrc_summary.csv",
    content = function(file) {
      f <- safe_fit(); if (inherits(f, "cgrc_err")) return()
      utils::write.csv(cgrc_analysis_summary(f), file, row.names = FALSE)
    })
  output$dl_curve <- downloadHandler(
    filename = function() "cgrc_curve.csv",
    content = function(file) {
      f <- safe_fit(); if (inherits(f, "cgrc_err")) return()
      utils::write.csv(curve_df(f), file, row.names = FALSE)
    })
  output$dl_png <- downloadHandler(
    filename = function() "cgrc_curve.png",
    content = function(file) {
      f <- safe_fit(); if (inherits(f, "cgrc_err")) return()
      lab <- if (f$dir < 0) "favourable" else "positive"
      p <- if (f$mode == "unknown")
        cgr_unknown_plot(f$ufit$curve, obs_cgr = f$ufit$observed_directional_cgr,
                         u = f$ufit$target_unknown_rate, direction_label = lab)
      else cgr_plot(f$fit$curve, obs_cgr = f$fit$observed_cgr, direction_label = lab)
      ggplot2::ggsave(file, p, width = 8, height = 6, dpi = 150)
    })
  output$dl_html <- downloadHandler(
    filename = function() "cgrc_report.html",
    content = function(file) {
      f <- safe_fit(); if (inherits(f, "cgrc_err")) return()
      writeLines(cgrc_build_html_report(f), file)
    })
  output$dl_clean <- downloadHandler(
    filename = function() "cgrc_analysis_data.csv",
    content = function(file) {
      f <- safe_fit(); if (inherits(f, "cgrc_err")) return()
      utils::write.csv(f$trial, file, row.names = FALSE)
    })
  output$dl_log <- downloadHandler(
    filename = function() "cgrc_exclusion_log.csv",
    content = function(file) {
      f <- safe_fit(); if (inherits(f, "cgrc_err")) return()
      utils::write.csv(f$audit$log, file, row.names = FALSE)
    })
  output$dl_report <- downloadHandler(
    filename = function() "cgrc_report.md",
    content = function(file) {
      f <- safe_fit(); if (inherits(f, "cgrc_err")) return()
      writeLines(cgrc_build_report(f), file)
    })

  # The headline: two plain probabilities, before and after reweighting the guess
  # rate to 0.50. This is the interpretable answer - "is there an effect" and "is
  # it big enough to matter" - that a single p-value cannot give.
  output$b_headline <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    pct <- function(p) sprintf("%.0f%%", 100 * p)
    if (f$mode == "unknown") {
      h <- f$uhead
      rhead <- "at directional CGR 0.50<br><span class='muted'>(UNKNOWN rate held fixed)</span>"
      title <- "Your trial, in two probabilities (UNKNOWN preserved)."
      note  <- h$text
    } else {
      h <- f$head
      rhead <- "at guessing-at-chance (CGR 0.50)"
      title <- "Your trial, in two probabilities."
      note  <- h$text
    }
    div(class = "verdict",
      HTML(sprintf(
        "<b>%s</b><br>
         <table style='width:100%%;margin-top:6px;border-collapse:collapse;'>
         <tr style='color:#666;'><td></td><td><b>at your CGR (raw)</b></td>
             <td><b>%s</b></td></tr>
         <tr><td>probability of a favourable effect</td>
             <td><b>%s</b></td><td><b>%s</b></td></tr>
         <tr><td>probability it is meaningful (beyond %.2g pts)</td>
             <td><b>%s</b></td><td><b>%s</b></td></tr></table>
         <div class='muted' style='margin-top:8px;'>%s</div>",
        title, rhead,
        pct(h$p_dir_obs), pct(h$p_dir_blind),
        h$delta, pct(h$p_meaningful_obs), pct(h$p_meaningful_blind),
        paste(note, "These are continuous probabilities — deliberately no",
              "significant/not cut-off."))))
  })

  output$b_strata <- renderTable({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    if (f$mode == "unknown") {
      st <- cgr_unknown_strata(f$trial)
      lab <- c(ACAC = "received active / guessed active",
               ACPL = "received active / guessed placebo",
               ACU  = "received active / UNKNOWN",
               PLAC = "received placebo / guessed active",
               PLPL = "received placebo / guessed placebo",
               PLU  = "received placebo / UNKNOWN")
      data.frame(stratum = lab[UNKNOWN_STRATA], n = lengths(st)[UNKNOWN_STRATA],
                 mean = round(vapply(UNKNOWN_STRATA, function(k)
                   if (length(st[[k]])) mean(st[[k]]) else NA_real_, numeric(1)), 2),
                 row.names = NULL)
    } else {
      st <- cgr_strata(f$trial)
      data.frame(stratum = STRATA, n = lengths(st)[STRATA],
                 mean = round(vapply(st[STRATA], mean, numeric(1)), 2),
                 row.names = NULL)
    }
  })

  output$b_summary <- renderTable({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    s <- if (f$mode == "unknown") f$ufit$summary else f$fit$summary
    ok <- cgrc_pct_ok(s$post_mean[1], s$cri_lo[1], s$cri_hi[1], s$post_mean[2])
    s$pct_attenuation[2] <- if (ok) s$pct_attenuation[2] else NA
    keep <- intersect(c("what","directional_cgr","unknown_rate","post_mean",
                        "cri_lo","cri_hi","p_favourable","pct_attenuation"), names(s))
    s[, keep]
  }, digits = 3)

  # For UNKNOWN mode, the extra design counts the brief asks to show prominently.
  output$b_counts <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err") || f$mode != "unknown") return(NULL)
    st <- cgr_unknown_strata(f$trial); n <- lengths(st)
    minstr <- min(n[n > 0])
    warn <- minstr < THIN_STRATUM
    div(class = if (warn) "verdict warn feas" else "verdict feas", HTML(sprintf(
      "<b>UNKNOWN-preserving design.</b> n total <b>%d</b>; directional <b>%d</b>;
       UNKNOWN <b>%d</b> (<b>%.1f%%</b>). Directional CGR <b>%.3f</b>.
       Active-arm guesses AC/PL/U: %d/%d/%d. Placebo-arm guesses AC/PL/U: %d/%d/%d.
       Smallest occupied stratum: <b>%d</b>%s.",
      f$ufit$n_total, f$ufit$n_directional, f$ufit$n_unknown,
      100 * f$ufit$observed_unknown_rate, f$ufit$observed_directional_cgr,
      n[["ACAC"]], n[["ACPL"]], n[["ACU"]], n[["PLAC"]], n[["PLPL"]], n[["PLU"]],
      minstr, if (warn) " (thin — the reweighted estimate is fragile)" else "")))
  })

  output$b_identity <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    if (f$mode == "unknown") {
      z <- cgr_unknown_reference_line_test(f$trial)
      div(class = "verdict feas", HTML(sprintf(
        "<b>Identity check.</b> At the observed directional CGR (%.4f) and observed
         UNKNOWN rate (%.1f%%), the reweighted curve equals the raw arm-mean
         difference to %.1e — the no-op identity holds.",
        z$computed_obs_directional_cgr, 100 * z$observed_unknown_rate,
        abs(z$D_at_obs - z$raw_mean_diff))))
    } else {
      z <- cgr_reference_line_test(f$trial, orig_cgr = f$fit$observed_cgr)
      div(class = "verdict feas", HTML(sprintf(
        "<b>Identity check.</b> Observed CGR = %.4f. The curve at the observed CGR
         equals the raw arm-mean difference to %.1e — the no-op identity holds, so
         the reference line is in the right place.",
        f$fit$observed_cgr, abs(z$D_at_obs - z$raw_mean_diff))))
    }
  })

  output$b_curve <- renderPlot({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    lab <- if (f$dir < 0) "favourable" else "positive"
    if (f$mode == "unknown")
      cgr_unknown_plot(f$ufit$curve, obs_cgr = f$ufit$observed_directional_cgr,
                       u = f$ufit$target_unknown_rate, direction_label = lab)
    else
      cgr_plot(f$fit$curve, obs_cgr = f$fit$observed_cgr, direction_label = lab)
  })

  output$b_rope <- renderPlot({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    if (f$mode == "unknown") {
      z <- f$urope; obs_cgr <- f$ufit$observed_directional_cgr
      xlab <- "directional correct-guess rate"
      sub  <- sprintf("black = directional 0.50; green = observed directional CGR (UNKNOWN held at %.1f%%)",
                      100 * f$ufit$target_unknown_rate)
    } else {
      z <- f$rope; obs_cgr <- f$fit$observed_cgr
      xlab <- "correct guess rate"
      sub  <- "black = guessing at chance (0.50); green = your observed CGR"
    }
    stack <- do.call(rbind, lapply(c("p_benefit","p_negligible","p_harm"), function(k)
      data.frame(cgr = z$cgr, p = z[[k]], region = k)))
    stack$region <- factor(stack$region, c("p_benefit","p_negligible","p_harm"),
      labels = c("meaningful benefit","practically negligible","meaningful harm"))
    ggplot(stack, aes(cgr, p, fill = region)) + geom_area() +
      geom_vline(xintercept = 0.5, linetype = "dashed") +
      geom_vline(xintercept = obs_cgr, linetype = "dashed", colour = "darkgreen") +
      scale_fill_manual(values = c("meaningful benefit" = "#2471A3",
        "practically negligible" = "grey75", "meaningful harm" = "#C0392B")) +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = xlab, y = "posterior probability", fill = NULL, subtitle = sub) +
      theme_minimal(base_size = 15) + theme(legend.position = "bottom")
  })

  output$b_sens <- renderTable({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    z <- if (f$mode == "unknown") f$usens else f$sens
    data.frame(`delta (SD frac)` = z$delta_in_SD,
               `delta (points)` = round(z$delta, 2),
               `P negligible` = round(z$p_negligible, 3),
               `P benefit` = round(z$p_benefit, 3), check.names = FALSE)
  }, digits = 3)

  output$to_design <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err")) return(NULL)
    if (f$mode == "unknown") {
      # The precomputed Panel A lookup models BINARY guessing only, so it cannot
      # be reused. Instead offer an on-demand UNKNOWN operating-characteristics
      # run using the six-stratum generative model (cgr_unknown_operating).
      return(tagList(tags$hr(),
        actionButton("run_unknown_design",
          sprintf("Run UNKNOWN design check (n=%d, dir CGR=%.2f, UNKNOWN=%.0f%%)",
                  nrow(f$trial), f$ufit$observed_directional_cgr,
                  100 * f$ufit$observed_unknown_rate), class = "btn-sm"),
        div(class = "muted",
            "A fresh six-stratum simulation (not the binary lookup). Assumes",
            "UNKNOWN responders carry no expectancy and the UNKNOWN rate is equal",
            "across arms; uses the effect, expectancy and trial count set on the",
            "Design tab.")))
    }
    tagList(tags$hr(), actionButton("do_bridge",
      sprintf("Run the design check at this trial (n=%d, CGR=%.2f)",
              nrow(f$trial), f$fit$observed_cgr), class = "btn-sm"))
  })

  ## Bridge: send the uploaded trial's n and observed CGR to Panel A - the exact
  ## workflow the paper's limitations paragraph describes. Binary mode only.
  observeEvent(input$do_bridge, {
    f <- safe_fit(); if (inherits(f, "cgrc_err") || f$mode != "binary") return()
    updateSliderInput(session, "n", value = round(nrow(f$trial) / 10) * 10)
    updateSliderInput(session, "pcg", value = round(f$fit$observed_cgr, 2))
    updateNavbarPage(session, "navbar", selected = "Design")
  })

  ## On-demand UNKNOWN design check: operating characteristics of the six-stratum
  ## estimator at THIS trial's n, directional CGR and observed UNKNOWN rate, from
  ## the UNKNOWN-aware generative model. Not a lookup - simulated on the spot.
  udesign_rv <- reactiveVal(NULL)
  observeEvent(input$run_unknown_design, {
    f <- safe_fit(); if (inherits(f, "cgrc_err") || f$mode != "unknown") return()
    nt <- input$n_trials
    withProgress(message = sprintf("Simulating %d UNKNOWN trials x 4 scenarios...", nt),
                 value = 0.3, {
      op <- cgr_unknown_operating(
        n_trials = nt, n = nrow(f$trial),
        p_cg = round(f$ufit$observed_directional_cgr, 2),
        u    = round(f$ufit$observed_unknown_rate, 2),
        mu_dte = as.numeric(input$eff), mu_aeb = as.numeric(input$mu_aeb), seed = 1)
      incProgress(0.7)
      attr(op, "n_trials") <- nt
      udesign_rv(op)
    })
    session$sendCustomMessage("btn_reenable", "run_unknown_design")
  })
  # reset the readout whenever a new analysis is run, so a stale table never shows
  observeEvent(input$analyse, udesign_rv(NULL))
  output$unknown_design_out <- renderUI({
    f <- safe_fit(); if (inherits(f, "cgrc_err") || f$mode != "unknown") return(NULL)
    op <- udesign_rv(); if (is.null(op)) return(NULL)
    minstr <- cgr_unknown_min_stratum(nrow(f$trial),
                round(f$ufit$observed_directional_cgr, 2),
                round(f$ufit$observed_unknown_rate, 2))
    tagList(br(), h4(sprintf("UNKNOWN design check (%d simulated trials)", attr(op, "n_trials"))),
      p(class = "muted", sprintf(paste(
        "Six-stratum operating characteristics at n=%d, directional CGR=%.2f,",
        "UNKNOWN rate=%.0f%%, true effect=%.1f, expectancy=%.1f. Smallest expected",
        "stratum ~%.0f. 'true' is the direct effect the adjusted estimate targets."),
        nrow(f$trial), round(f$ufit$observed_directional_cgr, 2),
        100 * round(f$ufit$observed_unknown_rate, 2),
        as.numeric(input$eff), as.numeric(input$mu_aeb), minstr)),
      renderTable({
        data.frame(scenario = CELL_LABEL[paste0(op$DTE, op$AEB)],
                   `true` = round(op$true, 2),
                   `adj bias` = round(op$adj_bias, 3),
                   `95% coverage` = round(op$coverage95, 3),
                   `adjusted flags` = round(op$p_fav_gt_95, 3),
                   `unadj significant` = round(op$freq_sig, 3),
                   `empty-stratum` = paste0(round(100 * op$empty_stratum_rate, 1), "%"),
                   check.names = FALSE)
      }, digits = 3))
  })
}

shinyApp(ui, server)
