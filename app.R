library(shiny)
library(bslib)
library(plotly)

# ==============================================================================
# 1. Model Defaults & Validation (SPEC Section 2 & 3.4)
# ==============================================================================

STATE_NAMES <- c("Normal", "Adenoma", "Cancer")

DEFAULT_PI <- c(0.90, 0.08, 0.02)

DEFAULT_A <- matrix(
  c(0.92, 0.07, 0.01,
    0.00, 0.88, 0.12,
    0.00, 0.00, 1.00),
  nrow = 3, byrow = TRUE,
  dimnames = list(From = STATE_NAMES, To = STATE_NAMES)
)

DEFAULT_B <- matrix(
  c(0.95, 0.05,
    0.65, 0.35,
    0.15, 0.85),
  nrow = 3, byrow = TRUE,
  dimnames = list(State = STATE_NAMES, Result = c("P(negative)", "P(positive)"))
)

# Color-blind-safe Okabe-Ito palette
PALETTE <- list(
  normal  = "#009E73",
  adenoma = "#E69F00",
  cancer  = "#D55E00",
  viterbi = "#0072B2",
  neg     = "#56B4E9",
  pos     = "#CC79A7"
)

validate_hmm <- function(pi_vec, A_mat, B_mat, tol = 1e-9) {
  stopifnot(
    length(pi_vec) == 3,
    all(pi_vec >= 0 & pi_vec <= 1),
    abs(sum(pi_vec) - 1) < tol,
    all(dim(A_mat) == c(3, 3)),
    all(A_mat >= 0 & A_mat <= 1),
    all(abs(rowSums(A_mat) - 1) < tol),
    all(dim(B_mat) == c(3, 2)),
    all(B_mat >= 0 & B_mat <= 1),
    all(abs(rowSums(B_mat) - 1) < tol)
  )
  invisible(TRUE)
}

# ==============================================================================
# 2. Core Algorithms in Log Space (SPEC Section 3)
# ==============================================================================

logsumexp <- function(x) {
  if (all(x == -Inf)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

safe_log <- function(x) {
  out <- rep(-Inf, length(x))
  pos <- x > 0
  out[pos] <- log(x[pos])
  if (!is.null(dim(x))) dim(out) <- dim(x)
  if (!is.null(dimnames(x))) dimnames(out) <- dimnames(x)
  out
}

#' Seeded Synthetic Patient Simulation (SPEC Section 4.2)
simulate_patient <- function(T_len, seed, pi_vec = DEFAULT_PI, A_mat = DEFAULT_A, B_mat = DEFAULT_B) {
  validate_hmm(pi_vec, A_mat, B_mat)
  stopifnot(T_len >= 1)
  set.seed(as.integer(seed))
  
  states <- integer(T_len)
  obs    <- integer(T_len)
  
  states[1] <- sample.int(3, size = 1, prob = pi_vec)
  obs[1]    <- sample.int(2, size = 1, prob = B_mat[states[1], ]) - 1L
  
  if (T_len >= 2) {
    for (t in 2:T_len) {
      states[t] <- sample.int(3, size = 1, prob = A_mat[states[t - 1], ])
      obs[t]    <- sample.int(2, size = 1, prob = B_mat[states[t], ]) - 1L
    }
  }
  
  list(states = states, obs = obs)
}

#' Log-Space Viterbi Algorithm (SPEC Section 3.1)
#' Tie-breaking: which.max() deterministically selects the lowest state index.
viterbi_log <- function(obs, pi_vec = DEFAULT_PI, A_mat = DEFAULT_A, B_mat = DEFAULT_B) {
  validate_hmm(pi_vec, A_mat, B_mat)
  T_len <- length(obs)
  stopifnot(T_len >= 1)
  
  log_pi <- safe_log(pi_vec)
  log_A  <- safe_log(A_mat)
  log_B  <- safe_log(B_mat)
  
  delta <- matrix(-Inf, nrow = T_len, ncol = 3, dimnames = list(paste0("t=", 1:T_len), STATE_NAMES))
  psi   <- matrix(NA_integer_, nrow = T_len, ncol = 3, dimnames = list(paste0("t=", 1:T_len), STATE_NAMES))
  
  delta[1, ] <- log_pi + log_B[, obs[1] + 1L]
  
  if (T_len >= 2) {
    for (t in 2:T_len) {
      for (j in 1:3) {
        candidates <- delta[t - 1, ] + log_A[, j]
        best_i <- which.max(candidates)
        psi[t, j]   <- best_i
        delta[t, j] <- candidates[best_i] + log_B[j, obs[t] + 1L]
      }
    }
  }
  
  best_log_prob <- max(delta[T_len, ])
  if (!is.finite(best_log_prob)) {
    return(list(path = NULL, delta = delta, psi = psi, log_prob = -Inf))
  }
  
  path <- integer(T_len)
  path[T_len] <- which.max(delta[T_len, ])
  if (T_len >= 2) {
    for (t in (T_len - 1):1) {
      path[t] <- psi[t + 1, path[t + 1]]
    }
  }
  
  list(path = path, delta = delta, psi = psi, log_prob = best_log_prob)
}

#' Log-Space Forward Algorithm (SPEC Section 3.2)
forward_log <- function(obs, pi_vec = DEFAULT_PI, A_mat = DEFAULT_A, B_mat = DEFAULT_B) {
  validate_hmm(pi_vec, A_mat, B_mat)
  T_len <- length(obs)
  stopifnot(T_len >= 1)
  
  log_pi <- safe_log(pi_vec)
  log_A  <- safe_log(A_mat)
  log_B  <- safe_log(B_mat)
  
  log_alpha <- matrix(-Inf, nrow = T_len, ncol = 3, dimnames = list(paste0("t=", 1:T_len), STATE_NAMES))
  log_alpha[1, ] <- log_pi + log_B[, obs[1] + 1L]
  
  if (T_len >= 2) {
    for (t in 2:T_len) {
      for (j in 1:3) {
        log_alpha[t, j] <- log_B[j, obs[t] + 1L] + logsumexp(log_alpha[t - 1, ] + log_A[, j])
      }
    }
  }
  
  ll <- logsumexp(log_alpha[T_len, ])
  list(
    log_alpha      = log_alpha,
    log_likelihood = ll,
    likelihood     = if (is.finite(ll)) exp(ll) else 0
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.infinite(x), "-Inf", formatC(x, format = "f", digits = digits))
}

# ==============================================================================
# 3. User Interface (SPEC Sections 4 & 5)
# ==============================================================================

ui <- page_sidebar(
  title = "HMM Visualizer for Colorectal Cancer Transition",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  withMathJax(),
  
  tags$head(tags$style(HTML("
    html { scroll-behavior: smooth; }
    .section-card { margin-bottom: 1.5rem; }
    .card-header { font-weight: 600; font-size: 1.05rem; }
    .tab-content { padding-top: 1rem; }
    .obs-btn { width: 56px; line-height: 1.2; padding: 0.35rem 0.25rem; }
    .metric-box { border-radius: 0.5rem; padding: 1rem 1.25rem; color: #fff; }
    .metric-label { font-size: 0.85rem; opacity: 0.9; margin-bottom: 0.25rem; }
    .metric-value { font-size: 1.75rem; font-weight: 700; line-height: 1.2; }
    .jump-links a { display: block; padding: 0.2rem 0; text-decoration: none; font-size: 0.9rem; }
    .jump-links a:hover { text-decoration: underline; }
    table.table { margin-bottom: 0; font-size: 0.9rem; }
  "))),
  
  sidebar = sidebar(
    width = 280,
    h5("Simulation Controls", class = "mb-3"),
    numericInput("seed", "Random Seed:", value = 42, min = 1, max = 999999, step = 1),
    sliderInput("seq_len", "Sequence Length (T):", min = 1, max = 20, value = 10, step = 1, ticks = FALSE),
    actionButton("sim_new", "Simulate new patient", class = "btn-primary w-100 mb-2"),
    actionButton("reset_btn", "Reset to defaults", class = "btn-outline-secondary w-100 mb-3"),
    checkboxInput("reveal_true", "Reveal true hidden states", value = TRUE),
    p(class = "small text-muted mb-3",
      "Teaching toggle: real clinicians never observe the true hidden state directly."),
    hr(),
    h6("Jump to Section", class = "text-uppercase text-muted small fw-bold"),
    div(
      class = "jump-links",
      tags$a(href = "#sec-overview", "1. Model Overview"),
      tags$a(href = "#sec-timeline", "2. Timeline & Observations"),
      tags$a(href = "#sec-viterbi", "3. Viterbi Algorithm"),
      tags$a(href = "#sec-forward", "4. Forward Algorithm"),
      tags$a(href = "#sec-compare", "5. Viterbi vs. Forward"),
      tags$a(href = "#sec-assumptions", "6. Assumptions & Limitations")
    )
  ),
  
  # Main content wrapper with generous spacing
  div(
    class = "container-fluid px-2 py-2",
    
    # Subtitle & Non-clinical Disclaimer Banner
    div(
      class = "mb-4",
      p(class = "lead mb-2",
        "An Interactive Hidden Markov Model for Simulated Colorectal Cancer Screening and Disease Progression"),
      div(
        class = "alert alert-warning mb-0 py-2 px-3",
        strong("Educational tool only. "),
        "All data are synthetic and all probabilities are hypothetical. Nothing in this app is a clinical conclusion, diagnosis, or medical advice."
      )
    ),
    
    # --------------------------------------------------------------------------
    # Section 1: Model Overview
    # --------------------------------------------------------------------------
    card(
      id = "sec-overview",
      class = "section-card shadow-sm",
      fill = FALSE,
      card_header("1. Model Overview"),
      card_body(
        fill = FALSE,
        div(
          class = "text-center p-3 bg-light rounded border mb-4",
          div(
            class = "mb-1",
            HTML("<strong>Progression-Only State Transitions:</strong> &nbsp;
                  <span style='color:#009E73;font-weight:700;'>Normal</span> &mdash;(0.07)&rarr;
                  <span style='color:#E69F00;font-weight:700;'>Adenoma</span> &mdash;(0.12)&rarr;
                  <span style='color:#D55E00;font-weight:700;'>Cancer (absorbing)</span>")
          ),
          div(
            class = "small text-muted",
            HTML("Direct Normal &rarr; Cancer jump: 0.01 &nbsp;|&nbsp; Self-loops: Normal 0.92, Adenoma 0.88, Cancer 1.00")
          )
        ),
        layout_columns(
          col_widths = c(4, 4, 4),
          fill = FALSE,
          div(
            h6("Initial Distribution (π)", class = "fw-semibold mb-2"),
            div(class = "table-responsive", tableOutput("tbl_pi"))
          ),
          div(
            h6("Transition Matrix (A)", class = "fw-semibold mb-2"),
            div(class = "table-responsive", tableOutput("tbl_A"))
          ),
          div(
            h6("Emission Matrix (B)", class = "fw-semibold mb-2"),
            div(class = "table-responsive", tableOutput("tbl_B"))
          )
        ),
        hr(class = "my-4"),
        navset_tab(
          nav_panel(
            "Intuition",
            p("A patient's true biological state (Normal, Adenoma, or Cancer) is ", strong("hidden"),
              " from the clinician. At each time step, a screening test produces an ", strong("observed"),
              " result (negative or positive). Because screening tests are imperfect—false positives and false negatives both happen—we use an HMM to infer the hidden states from the pattern of test results over time."),
            p("The ", strong("transition matrix A"), " governs how the hidden disease progresses from one step to the next, while the ",
              strong("emission matrix B"), " governs how likely a positive or negative test result is given the current hidden state.")
          ),
          nav_panel(
            "Formal Definition",
            p("Let \\(z_t \\in \\{0, 1, 2\\}\\) denote the hidden state (Normal, Adenoma, Cancer) and \\(x_t \\in \\{0, 1\\}\\) denote the screening result (negative, positive) at time \\(t\\). With parameters \\(\\theta = (\\pi, A, B)\\):"),
            p("$$\\pi_i = P(z_1 = i), \\qquad A_{ij} = P(z_{t+1} = j \\mid z_t = i), \\qquad B_j(x) = P(x_t = x \\mid z_t = j)$$"),
            p("The joint probability of a full hidden path and observation sequence is:"),
            p("$$P(z_{1:T}, x_{1:T}) = \\pi_{z_1} B_{z_1}(x_1) \\prod_{t=2}^{T} A_{z_{t-1} z_t} B_{z_t}(x_t)$$")
          ),
          nav_panel(
            "Worked Example (Live)",
            uiOutput("model_worked")
          )
        )
      )
    ),
    
    # --------------------------------------------------------------------------
    # Section 2: Synthetic Patient Timeline & Editable Observations
    # --------------------------------------------------------------------------
    card(
      id = "sec-timeline",
      class = "section-card shadow-sm",
      fill = FALSE,
      card_header("2. Synthetic Patient Timeline & Editable Observations"),
      card_body(
        fill = FALSE,
        plotlyOutput("timeline_plot", height = "360px"),
        hr(class = "my-4"),
        h6("Edit Observations", class = "fw-semibold mb-1"),
        p(class = "small text-muted mb-3",
          "Click any time step below to toggle between Negative (-) and Positive (+). All downstream inference updates automatically."),
        uiOutput("obs_controls")
      )
    ),
    
    # --------------------------------------------------------------------------
    # Section 3: Viterbi Algorithm
    # --------------------------------------------------------------------------
    card(
      id = "sec-viterbi",
      class = "section-card shadow-sm",
      fill = FALSE,
      card_header("3. Viterbi Algorithm: Most Likely Hidden Path"),
      card_body(
        fill = FALSE,
        uiOutput("viterbi_alert"),
        layout_columns(
          col_widths = c(7, 5),
          fill = FALSE,
          div(
            h6("Viterbi Trellis (δ values & winning back-pointers)", class = "fw-semibold mb-2"),
            plotlyOutput("trellis_plot", height = "360px"),
            p(class = "small text-muted mt-2 mb-0",
              "Blue highlighted nodes and edges trace the optimal Viterbi path; gray edges show non-winning back-pointers.")
          ),
          div(
            h6("Step-by-Step Calculation Table", class = "fw-semibold mb-2"),
            layout_columns(
              col_widths = c(6, 6),
              fill = FALSE,
              uiOutput("calc_t_ui"),
              selectInput("calc_j", "Target State (j):", choices = setNames(1:3, STATE_NAMES), selected = 2)
            ),
            div(class = "table-responsive mb-3", tableOutput("viterbi_calc_tbl")),
            uiOutput("viterbi_calc_summary")
          )
        ),
        hr(class = "my-4"),
        navset_tab(
          nav_panel(
            "Intuition",
            p("The Viterbi algorithm finds the single most probable sequence of hidden disease states that could have produced the observed screening results. Instead of enumerating all \\(3^T\\) possible paths, it steps forward one time point at a time, recording only the best preceding state for each current state, and then backtracks from the end to recover the winning path.")
          ),
          nav_panel(
            "Formal Definition",
            p("All inference runs in log space to prevent underflow, with \\(\\log 0 = -\\infty\\):"),
            p("$$\\delta_1(j) = \\log \\pi_j + \\log B_j(x_1)$$"),
            p("$$\\delta_t(j) = \\max_i \\left[ \\delta_{t-1}(i) + \\log A_{ij} \\right] + \\log B_j(x_t)$$"),
            p("$$\\psi_t(j) = \\arg\\max_i \\left[ \\delta_{t-1}(i) + \\log A_{ij} \\right]$$"),
            p("Backtrack from \\(\\hat{z}_T = \\arg\\max_j \\delta_T(j)\\) via \\(\\hat{z}_t = \\psi_{t+1}(\\hat{z}_{t+1})\\). Ties are broken deterministically by choosing the lowest state index.")
          ),
          nav_panel(
            "Worked Example (Live)",
            uiOutput("viterbi_worked")
          )
        )
      )
    ),
    
    # --------------------------------------------------------------------------
    # Section 4: Forward Algorithm
    # --------------------------------------------------------------------------
    card(
      id = "sec-forward",
      class = "section-card shadow-sm",
      fill = FALSE,
      card_header("4. Forward Algorithm: Sequence Likelihood"),
      card_body(
        fill = FALSE,
        layout_columns(
          col_widths = c(6, 6),
          fill = FALSE,
          class = "mb-4",
          div(
            class = "metric-box bg-primary",
            div(class = "metric-label", "Log-Likelihood: log P(x)"),
            div(class = "metric-value", textOutput("fwd_ll_val", inline = TRUE))
          ),
          div(
            class = "metric-box bg-info",
            div(class = "metric-label", "Likelihood: P(x)"),
            div(class = "metric-value", textOutput("fwd_lik_val", inline = TRUE))
          )
        ),
        h6("Forward Log-Alpha Values (first 3 time steps)", class = "fw-semibold mb-2"),
        div(class = "table-responsive", tableOutput("fwd_alpha_tbl")),
        hr(class = "my-4"),
        navset_tab(
          nav_panel(
            "Intuition",
            p("The Forward algorithm computes the total probability of observing this exact sequence of screening results under the model, summing across every possible hidden disease trajectory. It measures how plausible the observed sequence is overall.")
          ),
          nav_panel(
            "Formal Definition",
            p("$$\\log \\alpha_1(j) = \\log \\pi_j + \\log B_j(x_1)$$"),
            p("$$\\log \\alpha_t(j) = \\log B_j(x_t) + \\operatorname{logsumexp}_i \\left[ \\log \\alpha_{t-1}(i) + \\log A_{ij} \\right]$$"),
            p("$$\\log P(x_{1:T}) = \\operatorname{logsumexp}_j \\log \\alpha_T(j)$$")
          ),
          nav_panel(
            "Worked Example (Live)",
            uiOutput("fwd_worked")
          )
        )
      )
    ),
    
    # --------------------------------------------------------------------------
    # Section 5: Viterbi vs Forward
    # --------------------------------------------------------------------------
    card(
      id = "sec-compare",
      class = "section-card shadow-sm",
      fill = FALSE,
      card_header("5. Viterbi vs. Forward Algorithm"),
      card_body(
        fill = FALSE,
        layout_columns(
          col_widths = c(6, 6),
          fill = FALSE,
          div(
            class = "p-3 bg-light rounded border h-100",
            h6("Viterbi: Best Single Path", class = "fw-bold text-primary"),
            p(class = "mb-0",
              "Uses ", strong("max"), " at each step to identify the single most likely hidden trajectory. Its log-probability reflects the joint probability of that one winning path and the observations.")
          ),
          div(
            class = "p-3 bg-light rounded border h-100",
            h6("Forward: Total Probability Over All Paths", class = "fw-bold text-info"),
            p(class = "mb-0",
              "Uses ", strong("sum (logsumexp)"), " at each step to combine probability mass across all possible hidden trajectories. Because it sums over every path, the Forward likelihood is always greater than or equal to the Viterbi path probability.")
          )
        )
      )
    ),
    
    # --------------------------------------------------------------------------
    # Section 6: Model Assumptions and Limitations
    # --------------------------------------------------------------------------
    card(
      id = "sec-assumptions",
      class = "section-card shadow-sm",
      fill = FALSE,
      card_header("6. Model Assumptions and Limitations"),
      card_body(
        fill = FALSE,
        tags$ul(
          class = "mb-0",
          tags$li(class = "mb-2", strong("Progression-only:"), " The model only permits Normal → Adenoma → Cancer transitions; no regression or spontaneous recovery is allowed."),
          tags$li(class = "mb-2", strong("No treatment modeled:"), " Clinical interventions such as polypectomy or surgery are not represented."),
          tags$li(class = "mb-2", strong("Abstract time steps:"), " Time steps (t = 1, 2, ...) are abstract units, not calibrated to months or years."),
          tags$li(class = "mb-2", strong("Hypothetical probabilities:"), " Transition and emission values are illustrative defaults, not calibrated to any real screening test."),
          tags$li(class = "mb-2", strong("Markov & conditional independence:"), " The next state depends only on the current state, and each screening result depends only on the current hidden state."),
          tags$li(class = "mb-2", strong("Homogeneous transitions:"), " The same transition matrix A applies at every step; age, sex, and individual risk factors are excluded."),
          tags$li(strong("Not clinical conclusions:"), " Viterbi paths and likelihoods are mathematical outputs of a toy model, not medical diagnoses or advice.")
        )
      )
    ),
    
    # Footer
    div(
      class = "text-center text-muted small py-3",
      "HMM Visualizer for Colorectal Cancer Transition  |  Educational use only  |  Synthetic data"
    )
  )
)

# ==============================================================================
# 4. Server Logic
# ==============================================================================

server <- function(input, output, session) {
  
  init_sim <- simulate_patient(10, 42)
  sim_data <- reactiveVal(init_sim)
  obs_vec  <- reactiveVal(init_sim$obs)
  
  # Re-simulate whenever seed or sequence length changes
  observeEvent(list(input$seed, input$seq_len), {
    req(is.numeric(input$seed), !is.na(input$seed), is.numeric(input$seq_len), input$seq_len >= 1)
    res <- simulate_patient(input$seq_len, input$seed)
    sim_data(res)
    obs_vec(res$obs)
  }, ignoreInit = TRUE)
  
  # Generate a new random seed
  observeEvent(input$sim_new, {
    new_seed <- sample.int(100000, 1)
    if (identical(new_seed, as.integer(input$seed))) new_seed <- new_seed + 1L
    updateNumericInput(session, "seed", value = new_seed)
  })
  
  # Reset all controls and edited observations back to defaults
  observeEvent(input$reset_btn, {
    updateNumericInput(session, "seed", value = 42)
    updateSliderInput(session, "seq_len", value = 10)
    updateCheckboxInput(session, "reveal_true", value = TRUE)
    res <- simulate_patient(10, 42)
    sim_data(res)
    obs_vec(res$obs)
  })
  
  # Render editable observation toggle buttons using a single JS event callback
  output$obs_controls <- renderUI({
    obs <- obs_vec()
    T_len <- length(obs)
    btns <- lapply(seq_len(T_len), function(t) {
      is_pos <- obs[t] == 1L
      tags$button(
        type = "button",
        class = paste("btn btn-sm obs-btn", ifelse(is_pos, "btn-danger", "btn-outline-secondary")),
        onclick = sprintf("Shiny.setInputValue('toggle_obs', %d, {priority: 'event'})", t),
        HTML(paste0("<div style='font-size:0.75rem;'>t", t, "</div><div style='font-size:1rem;font-weight:700;'>", ifelse(is_pos, "+", "&minus;"), "</div>"))
      )
    })
    div(class = "d-flex flex-wrap gap-2", btns)
  })
  
  # Single observer for all observation toggles (no duplicate listeners)
  observeEvent(input$toggle_obs, {
    t_idx <- as.integer(input$toggle_obs)
    cur <- obs_vec()
    if (!is.na(t_idx) && t_idx >= 1L && t_idx <= length(cur)) {
      cur[t_idx] <- 1L - cur[t_idx]
      obs_vec(cur)
    }
  })
  
  # Reactive inference
  vit_res <- reactive({
    viterbi_log(obs_vec(), DEFAULT_PI, DEFAULT_A, DEFAULT_B)
  })
  
  fwd_res <- reactive({
    forward_log(obs_vec(), DEFAULT_PI, DEFAULT_A, DEFAULT_B)
  })
  
  # --- Section 1: Parameter Tables & Live Model Example ---
  output$tbl_pi <- renderTable({
    data.frame(State = STATE_NAMES, Probability = DEFAULT_PI)
  }, digits = 2, bordered = TRUE, hover = TRUE, width = "100%", spacing = "s")
  
  output$tbl_A <- renderTable({
    as.data.frame(DEFAULT_A)
  }, rownames = TRUE, digits = 2, bordered = TRUE, hover = TRUE, width = "100%", spacing = "s")
  
  output$tbl_B <- renderTable({
    as.data.frame(DEFAULT_B)
  }, rownames = TRUE, digits = 2, bordered = TRUE, hover = TRUE, width = "100%", spacing = "s")
  
  output$model_worked <- renderUI({
    obs <- obs_vec()
    st  <- sim_data()$states
    s1  <- st[1]
    x1  <- obs[1]
    p_init <- DEFAULT_PI[s1]
    p_emit <- DEFAULT_B[s1, x1 + 1L]
    tagList(
      p(strong("At t = 1 for this patient: "),
        sprintf("The simulated true state is %s (initial probability π = %.2f) and the observed test result is %s (emission probability B = %.2f).",
                STATE_NAMES[s1], p_init, ifelse(x1 == 1L, "Positive (+)", "Negative (-)"), p_emit)),
      if (length(obs) >= 2) {
        s2 <- st[2]
        p_trans <- DEFAULT_A[s1, s2]
        p(strong("Transition t = 1 → 2: "),
          sprintf("%s → %s has transition probability A = %.2f.",
                  STATE_NAMES[s1], STATE_NAMES[s2], p_trans))
      }
    )
  })
  
  # --- Section 2: Timeline Chart ---
  output$timeline_plot <- renderPlotly({
    obs <- obs_vec()
    T_len <- length(obs)
    t_seq <- seq_len(T_len)
    true_states <- sim_data()$states
    vit <- vit_res()
    
    p <- plot_ly()
    
    if (isTRUE(input$reveal_true)) {
      p <- p %>% add_trace(
        x = t_seq, y = true_states,
        type = "scatter", mode = "lines+markers",
        name = "True Hidden State",
        line = list(color = PALETTE$normal, width = 2.5, shape = "hv"),
        marker = list(size = 13, symbol = "circle-open", line = list(width = 2.5, color = PALETTE$normal))
      )
    }
    
    if (!is.null(vit$path)) {
      p <- p %>% add_trace(
        x = t_seq, y = vit$path,
        type = "scatter", mode = "lines+markers",
        name = "Viterbi Inferred Path",
        line = list(color = PALETTE$viterbi, width = 2, dash = "dot"),
        marker = list(size = 10, symbol = "x", color = PALETTE$viterbi)
      )
    }
    
    obs_colors <- ifelse(obs == 1L, PALETTE$pos, PALETTE$neg)
    obs_text   <- ifelse(obs == 1L, "+", "-")
    p <- p %>% add_trace(
      x = t_seq, y = rep(0, T_len),
      type = "scatter", mode = "markers+text",
      name = "Screening Result (+/-)",
      text = obs_text,
      textposition = "middle center",
      textfont = list(color = "white", size = 13, family = "Arial Black"),
      marker = list(size = 22, symbol = "square", color = obs_colors),
      hovertext = paste0("t = ", t_seq, ": ", ifelse(obs == 1L, "Positive (+)", "Negative (-)")),
      hoverinfo = "text"
    )
    
    p %>% layout(
      margin = list(t = 55, r = 20, b = 50, l = 90),
      xaxis = list(title = "Time Step (t)", tickmode = "linear", dtick = 1, range = c(0.4, T_len + 0.6)),
      yaxis = list(
        title = "",
        tickvals = c(0, 1, 2, 3),
        ticktext = c("Screening", "Normal", "Adenoma", "Cancer"),
        range = c(-0.6, 3.6)
      ),
      legend = list(orientation = "h", x = 0, y = 1.15)
    ) %>% config(displayModeBar = FALSE)
  })
  
  # --- Section 3: Viterbi Alert, Trellis, Calculation Table ---
  output$viterbi_alert <- renderUI({
    if (is.null(vit_res()$path)) {
      div(class = "alert alert-danger mb-3",
          "No hidden-state path can produce this sequence under the current model.")
    }
  })
  
  output$trellis_plot <- renderPlotly({
    vit <- vit_res()
    T_len <- nrow(vit$delta)
    p <- plot_ly()
    
    if (T_len >= 2) {
      for (t in 2:T_len) {
        for (j in 1:3) {
          prev_i <- vit$psi[t, j]
          if (!is.na(prev_i) && is.finite(vit$delta[t, j])) {
            is_opt <- !is.null(vit$path) && vit$path[t] == j && vit$path[t - 1] == prev_i
            p <- p %>% add_segments(
              x = t - 1, xend = t, y = prev_i, yend = j,
              line = list(
                color = ifelse(is_opt, PALETTE$viterbi, "#CCCCCC"),
                width = ifelse(is_opt, 3.5, 1.2)
              ),
              showlegend = FALSE, hoverinfo = "none"
            )
          }
        }
      }
    }
    
    grid_df <- expand.grid(t = 1:T_len, s = 1:3)
    grid_df$val <- mapply(function(t, s) vit$delta[t, s], grid_df$t, grid_df$s)
    grid_df$lbl <- fmt_num(grid_df$val, 2)
    grid_df$is_opt <- if (!is.null(vit$path)) mapply(function(t, s) vit$path[t] == s, grid_df$t, grid_df$s) else FALSE
    
    p <- p %>% add_trace(
      data = grid_df,
      x = ~t, y = ~s,
      type = "scatter", mode = "markers+text",
      text = ~lbl, textposition = "top center",
      textfont = list(size = 10),
      marker = list(
        size = ifelse(grid_df$is_opt, 13, 8),
        color = ifelse(grid_df$is_opt, PALETTE$viterbi, "#888888")
      ),
      hovertext = ~paste0("t = ", t, ", ", STATE_NAMES[s], "<br>δ = ", lbl),
      hoverinfo = "text",
      showlegend = FALSE
    )
    
    p %>% layout(
      margin = list(t = 25, r = 20, b = 50, l = 85),
      xaxis = list(title = "Time Step (t)", tickmode = "linear", dtick = 1, range = c(0.4, T_len + 0.6)),
      yaxis = list(title = "", tickvals = 1:3, ticktext = STATE_NAMES, range = c(0.5, 3.7))
    ) %>% config(displayModeBar = FALSE)
  })
  
  output$calc_t_ui <- renderUI({
    T_len <- length(obs_vec())
    if (T_len == 1) {
      selectInput("calc_t", "Time Step (t):", choices = 1, selected = 1)
    } else {
      selectInput("calc_t", "Time Step (t):", choices = 2:T_len, selected = min(4, T_len))
    }
  })
  
  output$viterbi_calc_tbl <- renderTable({
    req(input$calc_t, input$calc_j)
    t_val <- as.integer(input$calc_t)
    j_val <- as.integer(input$calc_j)
    vit <- vit_res()
    req(t_val >= 1, t_val <= nrow(vit$delta))
    
    if (t_val == 1) {
      data.frame(
        Term = c("log π_j", "log B_j(x_1)"),
        Value = c(fmt_num(safe_log(DEFAULT_PI[j_val])), fmt_num(safe_log(DEFAULT_B[j_val, obs_vec()[1] + 1L])))
      )
    } else {
      prev_d <- vit$delta[t_val - 1, ]
      log_a  <- safe_log(DEFAULT_A[, j_val])
      sums   <- prev_d + log_a
      winner <- which.max(sums)
      data.frame(
        `Prev State` = STATE_NAMES,
        `δ(t-1)`     = fmt_num(prev_d),
        `log A_ij`   = fmt_num(log_a),
        Sum          = fmt_num(sums),
        Winner       = ifelse(seq_len(3) == winner & is.finite(sums), "◀ max", ""),
        check.names  = FALSE
      )
    }
  }, bordered = TRUE, hover = TRUE, width = "100%", spacing = "s")
  
  output$viterbi_calc_summary <- renderUI({
    req(input$calc_t, input$calc_j)
    t_val <- as.integer(input$calc_t)
    j_val <- as.integer(input$calc_j)
    vit <- vit_res()
    req(t_val >= 1, t_val <= nrow(vit$delta))
    obs_t <- obs_vec()[t_val]
    emit_val <- safe_log(DEFAULT_B[j_val, obs_t + 1L])
    res_val  <- vit$delta[t_val, j_val]
    
    div(
      class = "small bg-light p-2 px-3 rounded border",
      strong(sprintf("δ(%d, %s) = %s", t_val, STATE_NAMES[j_val], fmt_num(res_val))),
      sprintf("  (includes emission log B(%s | %s) = %s)",
              ifelse(obs_t == 1L, "+", "-"), STATE_NAMES[j_val], fmt_num(emit_val))
    )
  })
  
  output$viterbi_worked <- renderUI({
    vit <- vit_res()
    obs <- obs_vec()
    if (is.null(vit$path)) {
      return(p("No valid path exists for the current observations."))
    }
    path_str <- paste(STATE_NAMES[vit$path], collapse = " → ")
    tagList(
      p(strong("Current observation sequence: "), paste(ifelse(obs == 1L, "+", "-"), collapse = " ")),
      p(strong("Initialization (t = 1): "),
        sprintf("δ(1, Normal) = %s,  δ(1, Adenoma) = %s,  δ(1, Cancer) = %s",
                fmt_num(vit$delta[1, 1]), fmt_num(vit$delta[1, 2]), fmt_num(vit$delta[1, 3]))),
      p(strong("Winning Viterbi path: "), path_str),
      p(class = "mb-0", strong("Best joint log-probability: "), fmt_num(vit$log_prob, 4))
    )
  })
  
  # --- Section 4: Forward Algorithm ---
  output$fwd_ll_val <- renderText({
    fmt_num(fwd_res()$log_likelihood, 4)
  })
  
  output$fwd_lik_val <- renderText({
    lik <- fwd_res()$likelihood
    if (lik == 0) "0" else formatC(lik, format = "e", digits = 3)
  })
  
  output$fwd_alpha_tbl <- renderTable({
    la <- fwd_res()$log_alpha
    n_show <- min(3, nrow(la))
    mat <- la[seq_len(n_show), , drop = FALSE]
    df <- as.data.frame(matrix(fmt_num(mat, 3), nrow = n_show, ncol = 3))
    colnames(df) <- STATE_NAMES
    cbind(`Time Step` = paste0("t = ", seq_len(n_show)), df)
  }, bordered = TRUE, hover = TRUE, width = "100%", spacing = "s")
  
  output$fwd_worked <- renderUI({
    fwd <- fwd_res()
    obs <- obs_vec()
    la <- fwd$log_alpha
    tagList(
      p(strong(sprintf("Step t = 1 (x₁ = %s): ", ifelse(obs[1] == 1L, "+", "-"))),
        sprintf("log α₁ = (%s, %s, %s)", fmt_num(la[1, 1]), fmt_num(la[1, 2]), fmt_num(la[1, 3]))),
      if (length(obs) >= 2) {
        p(strong(sprintf("Step t = 2 (x₂ = %s): ", ifelse(obs[2] == 1L, "+", "-"))),
          sprintf("log α₂ = (%s, %s, %s)", fmt_num(la[2, 1]), fmt_num(la[2, 2]), fmt_num(la[2, 3])))
      },
      p(class = "mb-0",
        strong(sprintf("Total log-likelihood at t = %d: ", length(obs))),
        sprintf("logsumexp(log α_%d) = %s  (P(x) = %s)",
                length(obs), fmt_num(fwd$log_likelihood, 4),
                if (fwd$likelihood == 0) "0" else formatC(fwd$likelihood, format = "e", digits = 3)))
    )
  })
}

shinyApp(ui = ui, server = server)