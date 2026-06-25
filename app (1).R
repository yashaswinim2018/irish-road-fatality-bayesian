# ============================================================
# SHINY APP: Bayesian Ranking of Traffic Fatality Rates
# per Distance Travelled in Ireland (2015-2024)
#
# Authors: Sara Raees & Yashaswini Malleshwarppa
#
# HOW TO RUN:
# 1. Place this file (app.R) in your HDS Project folder
#    -- the SAME folder that has road_accidents.csv,
#    THA17.....csv, and panel.rds / posterior_M3.rds if
#    you already saved them from your Rmd.
# 2. Install any missing packages once:
#    install.packages(c("shiny","shinydashboard","DT","plotly",
#                        "scales","leaflet","sf","rnaturalearth"))
# 3. Open app.R in RStudio
# 4. Click "Run App" (top right of the script editor)
#
# The app will:
# - Load and clean the data exactly as in your report
# - Fit the M3 Bayesian model (or load a cached version if
#   posterior_M3.rds already exists, so it opens instantly)
# - Build the Ireland county map once and cache it to
#   ie_counties_sf.rds (only downloads boundaries the first time)
# - Let your professor explore EVERY result interactively,
#   including hovering/clicking counties on the map for
#   live Bayesian results
# ============================================================

# ---- Packages ----
library(shiny)
library(shinydashboard)
library(rjags)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(viridis)
library(DT)
library(plotly)
library(scales)
library(leaflet)
library(sf)
library(rnaturalearth)

# ---- Detect whether JAGS is actually usable in this environment ----
# (shinyapps.io has no JAGS binary installed; live model fitting in the
# DIC tab is only offered when running locally where JAGS is present.)
jags_available <- tryCatch({
  requireNamespace("rjags", quietly = TRUE) &&
    !inherits(try(rjags::jags.model(textConnection(
      "model{x~dnorm(0,1)}"), quiet = TRUE), silent = TRUE), "try-error")
}, error = function(e) FALSE)


load_panel <- function() {

  if (file.exists("panel.rds") && file.exists("county_lookup.rds")) {
    panel         <- readRDS("panel.rds")
    county_lookup <- readRDS("county_lookup.rds")
    return(list(panel = panel, county_lookup = county_lookup))
  }

  df1 <- read.csv("road_accidents.csv")
  df2 <- read.csv("THA17.20251105T151118.csv")

  acc_wide <- df1 %>%
    filter(Statistic.Label %in% c("Persons Killed", "Persons Injured"),
           Year >= 2015, Year <= 2024) %>%
    mutate(County = str_replace(County, "^Co\\.\\s*", "")) %>%
    select(Year, County, Statistic.Label, VALUE) %>%
    pivot_wider(names_from = Statistic.Label, values_from = VALUE) %>%
    rename(deaths = `Persons Killed`, injuries = `Persons Injured`)

  vmt <- df2 %>%
    mutate(County = str_replace(County.of.Ownership, "^Co\\.\\s*", "")) %>%
    filter(Year >= 2015, Year <= 2024,
           County.of.Ownership != "Ireland",
           Type.of.Vehicle == "All vehicle types",
           Year.of.Registration == "All years") %>%
    transmute(Year, County, km = VALUE * 1e6) %>%
    mutate(VM = km / 1e9)

  panel <- acc_wide %>%
    left_join(vmt, by = c("Year", "County")) %>%
    filter(!is.na(VM), VM > 0) %>%
    mutate(
      covid        = ifelse(Year %in% c(2020, 2021), 1L, 0L),
      Year_c       = Year - mean(Year, na.rm = TRUE),
      county_id    = as.integer(factor(County)),
      rate_billion = deaths / VM,
      rate_million = deaths / (km / 1e6)
    ) %>%
    arrange(County, Year)

  county_lookup <- panel %>%
    select(County, county_id) %>% distinct() %>% arrange(county_id)

  saveRDS(panel,         "panel.rds")
  saveRDS(county_lookup, "county_lookup.rds")

  list(panel = panel, county_lookup = county_lookup)
}

dat           <- load_panel()
panel         <- dat$panel
county_lookup <- dat$county_lookup

X      <- panel$deaths
VM_vec <- panel$VM
Year_c <- panel$Year_c
covid  <- panel$covid
county <- panel$county_id
nObs   <- nrow(panel)
N      <- length(unique(county))

# ---- Road infrastructure data (TII 2015) ----
road_data <- tribble(
  ~County,~motorway_km,~dual_km,~single_km,~total_km,~area_km2,
  "Carlow",23.533,0,54.329,77.862,897,
  "Cavan",0,0,123.272,123.272,1932,
  "Clare",31.631,19.546,181.546,232.723,3449,
  "Cork",48.915,62.651,402.680,514.246,7500,
  "Donegal",0,5.714,298.185,303.899,4861,
  "Dublin",81.560,41.397,15.003,138.960,921,
  "Galway",61.538,12.312,379.938,453.788,6148,
  "Kerry",0,11.056,413.165,424.221,4807,
  "Kildare",108.233,8.006,17.888,134.127,1693,
  "Kilkenny",67.602,17.032,111.946,196.580,2073,
  "Laois",66.979,0,101.508,168.487,1719,
  "Leitrim",0,7.282,49.095,56.377,1590,
  "Limerick",27.250,8.927,156.361,192.538,2756,
  "Longford",0,0.156,96.618,96.774,1091,
  "Louth",39.682,9.350,48.754,97.786,826,
  "Mayo",0,0,397.316,397.316,5586,
  "Meath",88.063,12.087,103.632,203.782,2342,
  "Monaghan",0,15.222,88.478,103.700,1295,
  "Offaly",14.684,0,101.247,115.931,2000,
  "Roscommon",21.143,1.151,224.793,247.087,2544,
  "Sligo",0,12.800,140.628,153.428,1836,
  "Tipperary",121.973,3.294,210.023,335.290,4305,
  "Waterford",0,9.690,97.061,106.751,1837,
  "Westmeath",56.229,17.632,101.538,175.399,1840,
  "Wexford",21.869,0,142.148,164.017,2353,
  "Wicklow",35.327,18.332,38.146,91.805,2024
) %>%
  mutate(
    pct_dual     = (motorway_km + dual_km) / total_km * 100,
    pct_single   = single_km / total_km * 100,
    pct_motorway = motorway_km / total_km * 100,
    road_density = total_km / area_km2
  )


fit_m3 <- function() {

  if (file.exists("posterior_M3.rds")) {
    return(readRDS("posterior_M3.rds"))
  }

  stop(
    "posterior_M3.rds not found.\n",
    "This app does not fit the JAGS model live (shinyapps.io has no JAGS\n",
    "installed and a 60-second startup limit). Run the model-fitting code\n",
    "locally first to create posterior_M3.rds, then redeploy with that\n",
    "file included in the app bundle."
  )
}

samp_m3 <- fit_m3()
pm3     <- as.matrix(samp_m3)
n_iter  <- nrow(pm3) / 4

b0_cols  <- grep("^b0\\[", colnames(pm3), value = TRUE)
b0_samps <- pm3[, b0_cols]

rank_matrix <- t(apply(b0_samps, 1, function(x) rank(-x)))

rr <- data.frame(
  County      = county_lookup$County,
  mean_b0     = colMeans(b0_samps),
  lower_b0    = apply(b0_samps, 2, quantile, 0.025),
  upper_b0    = apply(b0_samps, 2, quantile, 0.975),
  mean_rank   = colMeans(rank_matrix),
  lower_rank  = apply(rank_matrix, 2, quantile, 0.025),
  upper_rank  = apply(rank_matrix, 2, quantile, 0.975),
  prob_higher = colMeans(b0_samps > 0),
  prob_lower  = colMeans(b0_samps < 0),
  RR          = exp(colMeans(b0_samps)) / mean(exp(colMeans(b0_samps)))
) %>%
  mutate(
    Evidence = case_when(
      pmax(prob_higher, prob_lower) > 0.999 ~ "Decisive",
      pmax(prob_higher, prob_lower) > 0.975 ~ "Strong",
      pmax(prob_higher, prob_lower) > 0.950 ~ "Moderate",
      pmax(prob_higher, prob_lower) > 0.750 ~ "Weak",
      TRUE ~ "None"
    ),
    Conclusion = case_when(
      prob_higher > 0.975 ~ "Significantly HIGHER RISK",
      prob_lower  > 0.975 ~ "Significantly LOWER RISK",
      TRUE ~ "No significant difference"
    )
  ) %>%
  left_join(road_data %>% select(County, pct_single, pct_motorway, road_density),
            by = "County") %>%
  arrange(mean_b0)

# Average traffic volume per county for bubble sizes
avg_km <- panel %>% group_by(County) %>%
  summarise(avg_km = mean(km, na.rm = TRUE), .groups = "drop")
rr <- rr %>% left_join(avg_km, by = "County")


build_county_geometry <- function() {

  if (file.exists("ie_counties_sf.rds")) {
    return(readRDS("ie_counties_sf.rds"))
  }

  ie <- ne_states(country = "Ireland", returnclass = "sf")

  ie_cleaned <- ie %>%
    mutate(County = case_when(
      name %in% c("Dublin", "D\u00fan Laoghaire\u2013Rathdown", "Fingal", "South Dublin") ~ "Dublin",
      name %in% c("North Tipperary", "South Tipperary") ~ "Tipperary",
      name == "Laoighis" ~ "Laois",
      TRUE ~ name
    )) %>%
    group_by(County) %>%
    summarise(geometry = st_union(geometry), .groups = "drop") %>%
    st_transform(4326)

  saveRDS(ie_cleaned, "ie_counties_sf.rds")
  ie_cleaned
}

ie_geom <- tryCatch(build_county_geometry(), error = function(e) NULL)

if (!is.null(ie_geom)) {
  rr_map <- ie_geom %>% left_join(rr, by = "County")
} else {
  rr_map <- NULL
}


ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "Irish Road Fatality Risk — Bayesian Model (M3)",
    titleWidth = 420
  ),

  dashboardSidebar(
    width = 240,
    sidebarMenu(
      menuItem("Overview",            tabName = "overview",  icon = icon("home")),
      menuItem("Ireland Map",         tabName = "map",        icon = icon("map-marked-alt")),
      menuItem("Exploratory Data",    tabName = "eda",        icon = icon("chart-line")),
      menuItem("Model & Diagnostics", tabName = "diagnostics",icon = icon("flask")),
      menuItem("County Rankings",     tabName = "rankings",   icon = icon("list-ol")),
      menuItem("Road Infrastructure", tabName = "roads",      icon = icon("road")),
      menuItem("Model Selection (DIC)",tabName = "dic",       icon = icon("balance-scale")),
      menuItem("Explore a County",    tabName = "explore",    icon = icon("search")),
      menuItem("About / Methods",     tabName = "about",      icon = icon("info-circle"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-top-color: #2c3e6b; }
    "))),

    tabItems(

      tabItem(tabName = "overview",
        fluidRow(
          valueBoxOutput("vb_obs"),
          valueBoxOutput("vb_counties"),
          valueBoxOutput("vb_draws")
        ),
        fluidRow(
          valueBoxOutput("vb_safest"),
          valueBoxOutput("vb_riskiest"),
          valueBoxOutput("vb_sig")
        ),
        fluidRow(
          box(title = "Project Summary", width = 12, status = "primary",
              solidHeader = TRUE,
              p("This dashboard accompanies the report ", strong("\"Bayesian Ranking of Traffic Fatality Rates per Distance Travelled in Ireland (2015-2024)\""), " by Sara Raees and Yashaswini Malleshwarappa."),
              p("The analysis uses a Bayesian hierarchical Poisson regression (Model M3) fitted in JAGS, with vehicle kilometres travelled (VKT) as the exposure offset and a county-level random intercept capturing baseline risk differences. The national time trend and COVID-19 effect are modelled as global fixed effects, since model selection (DIC) showed county-specific deviations for these were not supported by the data."),
              p("Use the tabs on the left to explore the raw data, the model diagnostics, the full county risk rankings with credible intervals, the road infrastructure validation, and the DIC model comparison.")
          )
        )
      ),

      tabItem(tabName = "map",
        fluidRow(
          box(title = "Interactive county risk map", width = 9, status = "primary",
              solidHeader = TRUE,
              p("Hover over a county to see its Bayesian results in a tooltip. Click a county to pin the full result panel on the right and keep it open."),
              radioButtons("map_metric", "Colour counties by:",
                           choices = c("Relative Risk (RR)" = "RR",
                                       "Posterior mean b0" = "mean_b0",
                                       "Posterior mean rank" = "mean_rank",
                                       "% Single carriageway" = "pct_single"),
                           selected = "RR", inline = TRUE),
              leafletOutput("ireland_map", height = 620)
          ),
          box(title = "Selected county", width = 3, status = "primary",
              solidHeader = TRUE,
              uiOutput("map_county_panel")
          )
        )
      ),


      tabItem(tabName = "eda",
        fluidRow(
          box(title = "Controls", width = 3, status = "primary",
              selectInput("eda_counties", "Select counties to highlight:",
                          choices = sort(unique(panel$County)),
                          selected = c("Dublin","Louth","Donegal"),
                          multiple = TRUE),
              radioButtons("eda_metric", "Metric:",
                           choices = c("Raw Deaths" = "deaths",
                                       "Raw Injuries" = "injuries",
                                       "Rate per billion km" = "rate_billion",
                                       "Rate per million km" = "rate_million"),
                           selected = "deaths")
          ),
          box(title = "County trends over time", width = 9, status = "primary",
              plotlyOutput("eda_trend_plot", height = 480))
        ),
        fluidRow(
          box(title = "National fatality rate (deaths per billion km)", width = 6,
              status = "primary", plotlyOutput("eda_national_plot", height = 380)),
          box(title = "Fatality rate heatmap — county x year", width = 6,
              status = "primary", plotOutput("eda_heatmap", height = 380))
        )
      ),

      tabItem(tabName = "diagnostics",
        fluidRow(
          box(title = "Model M3 specification", width = 12, status = "primary",
              solidHeader = TRUE,
              withMathJax(),
              p("$$X_{it} \\sim \\text{Poisson}(\\mu_{it})$$"),
              p("$$\\log(\\mu_{it}) = \\log(\\text{VM}_{it}) + (\\beta_0 + b_{0c}) + \\beta_1 \\cdot \\text{Year\\_c}_t + \\beta_2 \\cdot \\text{COVID}_t$$"),
              p("$$b_{0c} \\sim \\text{Normal}(0,\\tau_0), \\quad \\tau_0 \\sim \\text{Gamma}(0.001,0.001)$$"),
              p("4 chains, 1000 adaptation, 5000 burn-in, 20000 sampling iterations per chain — 80,000 total posterior draws.")
          )
        ),
        fluidRow(
          box(title = "Posterior summary — global parameters", width = 5,
              status = "primary", DTOutput("posterior_table")),
          box(title = "Gelman-Rubin R-hat and Effective Sample Size", width = 7,
              status = "primary", DTOutput("convergence_table"))
        ),
        fluidRow(
          box(title = "Trace plots (select parameter)", width = 6, status = "primary",
              selectInput("trace_param", "Parameter:",
                          choices = c("beta0","beta1","beta2","sd0"), selected = "beta0"),
              plotOutput("trace_plot", height = 350)),
          box(title = "Density plot — 4 chains overlaid", width = 6, status = "primary",
              plotOutput("density_plot", height = 350))
        ),
        fluidRow(
          box(title = "Prior vs posterior", width = 12, status = "primary",
              plotOutput("prior_post_plot", height = 320))
        )
      ),

      tabItem(tabName = "rankings",
        fluidRow(
          box(title = "Bayesian county baseline risk rankings (b0)", width = 6,
              status = "primary", plotlyOutput("rank_caterpillar", height = 600)),
          box(title = "Credible intervals on the RANK itself", width = 6,
              status = "primary", plotlyOutput("rank_ci_plot", height = 600))
        ),
        fluidRow(
          box(title = "Full results table", width = 12, status = "primary",
              DTOutput("rankings_table"))
        )
      ),

      tabItem(tabName = "roads",
        fluidRow(
          box(title = "Single carriageway % vs fatality rate", width = 6,
              status = "primary", plotlyOutput("road_scatter", height = 480)),
          box(title = "Bayesian risk vs road type and traffic volume", width = 6,
              status = "primary", plotlyOutput("road_combined", height = 480))
        ),
        fluidRow(
          box(title = "Road infrastructure data (TII 2015)", width = 12,
              status = "primary", DTOutput("road_table"))
        )
      ),

      tabItem(tabName = "dic",
        fluidRow(
          box(title = "Model selection (DIC) -- 5 nested models", width = 12,
              status = "warning", solidHeader = TRUE,
              uiOutput("dic_intro"),
              br(),
              DTOutput("dic_table")
          )
        )
      ),

      tabItem(tabName = "explore",
        fluidRow(
          box(title = "Choose a county", width = 3, status = "primary",
              selectInput("county_pick", "County:",
                          choices = sort(unique(panel$County)), selected = "Dublin")),
          box(title = "County summary", width = 9, status = "primary",
              tableOutput("county_summary"))
        ),
        fluidRow(
          box(title = "Raw deaths and injuries over time", width = 6,
              status = "primary", plotlyOutput("county_raw_plot", height = 380)),
          box(title = "Exposure-adjusted rate over time", width = 6,
              status = "primary", plotlyOutput("county_rate_plot", height = 380))
        ),
        fluidRow(
          box(title = "Pairwise comparison", width = 12, status = "primary",
              selectInput("county_pick2", "Compare against:",
                          choices = sort(unique(panel$County)), selected = "Louth"),
              verbatimTextOutput("pairwise_result"))
        )
      ),

      tabItem(tabName = "about",
        fluidRow(
          box(title = "About this app", width = 12, status = "primary",
              solidHeader = TRUE,
              h4("Data sources"),
              tags$ul(
                tags$li("RSA — annual county-level fatality and injury counts (2015-2024)."),
                tags$li("CSO Table THA17 — vehicle kilometres travelled (VKT) by county and year."),
                tags$li("TII National Road Lengths 2015 — motorway, dual, and single carriageway km by county."),
                tags$li("CSO intercensal population estimates and Census 2016/2022.")
              ),
              h4("Methodology"),
              p("A Bayesian hierarchical Poisson regression model was fitted in JAGS with VKT as a fixed-coefficient exposure offset. County-level random intercepts capture baseline risk heterogeneity; the national time trend and COVID-19 effect are modelled as global fixed effects, following model selection via the Deviance Information Criterion (DIC)."),
              h4("Authors"),
              p("Sara Raees & Yashaswini Malleshwarppa, MSc Health Data Science, University of Galway. Supervisors: Carl Scarrott & John Ferguson.")
          )
        )
      )
    )
  )
)


server <- function(input, output, session) {

  # ---------------- Overview value boxes ----------------
  output$vb_obs <- renderValueBox({
    valueBox(nObs, "Observations (county-years)", icon = icon("database"), color = "blue")
  })
  output$vb_counties <- renderValueBox({
    valueBox(N, "Counties", icon = icon("map"), color = "purple")
  })
  output$vb_draws <- renderValueBox({
    valueBox(format(nrow(pm3), big.mark = ","), "Posterior draws",
             icon = icon("dice"), color = "green")
  })
  output$vb_safest <- renderValueBox({
    safest <- rr$County[which.min(rr$mean_b0)]
    valueBox(safest, "Lowest risk county", icon = icon("shield-alt"), color = "green")
  })
  output$vb_riskiest <- renderValueBox({
    riskiest <- rr$County[which.max(rr$mean_b0)]
    valueBox(riskiest, "Highest risk county", icon = icon("exclamation-triangle"), color = "red")
  })
  output$vb_sig <- renderValueBox({
    n_sig <- sum(rr$Conclusion != "No significant difference")
    valueBox(n_sig, "Counties significantly different from average",
             icon = icon("check-circle"), color = "yellow")
  })


  pinned_county <- reactiveVal(NULL)

  output$ireland_map <- renderLeaflet({

    if (is.null(rr_map)) {
      return(
        leaflet() %>% addTiles() %>%
          setView(lng = -8, lat = 53.4, zoom = 6)
      )
    }

    metric <- "RR"   # default colouring on first draw; updated via observer below
    pal <- colorNumeric("viridis", domain = rr_map[[metric]])

    labels <- sprintf(
      "<strong>%s</strong><br/>RR: %.2f<br/>Mean b0: %.3f<br/>Mean rank: %.1f<br/>%s",
      rr_map$County, rr_map$RR, rr_map$mean_b0, rr_map$mean_rank, rr_map$Conclusion
    ) %>% lapply(htmltools::HTML)

    leaflet(rr_map) %>%
      addProviderTiles("CartoDB.Positron") %>%
      setView(lng = -8, lat = 53.4, zoom = 6.4) %>%
      addPolygons(
        layerId = ~County,
        fillColor = ~pal(get(metric)),
        color = "white", weight = 1.2, fillOpacity = 0.78,
        highlightOptions = highlightOptions(
          weight = 3, color = "#222", fillOpacity = 0.9, bringToFront = TRUE
        ),
        label = labels,
        labelOptions = labelOptions(
          style = list("font-weight" = "normal", padding = "6px 10px"),
          textsize = "13px", direction = "auto"
        )
      ) %>%
      addLegend(pal = pal, values = rr_map[[metric]],
                title = metric, position = "bottomright")
  })

  observeEvent(input$map_metric, {
    if (is.null(rr_map)) return(NULL)
    metric <- input$map_metric
    pal <- colorNumeric("viridis", domain = rr_map[[metric]])

    labels <- sprintf(
      "<strong>%s</strong><br/>RR: %.2f<br/>Mean b0: %.3f<br/>Mean rank: %.1f<br/>%% single: %.0f%%<br/>%s",
      rr_map$County, rr_map$RR, rr_map$mean_b0, rr_map$mean_rank,
      rr_map$pct_single, rr_map$Conclusion
    ) %>% lapply(htmltools::HTML)

    leafletProxy("ireland_map", data = rr_map) %>%
      clearShapes() %>%
      clearControls() %>%
      addPolygons(
        layerId = ~County,
        fillColor = ~pal(get(metric)),
        color = "white", weight = 1.2, fillOpacity = 0.78,
        highlightOptions = highlightOptions(
          weight = 3, color = "#222", fillOpacity = 0.9, bringToFront = TRUE
        ),
        label = labels,
        labelOptions = labelOptions(
          style = list("font-weight" = "normal", padding = "6px 10px"),
          textsize = "13px", direction = "auto"
        )
      ) %>%
      addLegend(pal = pal, values = rr_map[[metric]],
                title = metric, position = "bottomright")
  })

  observeEvent(input$ireland_map_shape_click, {
    click <- input$ireland_map_shape_click
    if (!is.null(click$id)) pinned_county(click$id)
  })

  output$map_county_panel <- renderUI({
    if (is.null(pinned_county())) {
      return(p(em("Click a county on the map to pin its full results here.")))
    }
    c1 <- pinned_county()
    rd <- rr %>% filter(County == c1)
    if (nrow(rd) == 0) return(p("No data for this county."))

    tagList(
      h4(c1),
      tags$table(class = "table table-condensed",
        tags$tr(tags$td(strong("Posterior mean b0")), tags$td(round(rd$mean_b0, 3))),
        tags$tr(tags$td(strong("95% CrI on b0")),
                tags$td(paste0("[", round(rd$lower_b0,3), ", ", round(rd$upper_b0,3), "]"))),
        tags$tr(tags$td(strong("Mean rank")), tags$td(round(rd$mean_rank, 1))),
        tags$tr(tags$td(strong("95% CrI on rank")),
                tags$td(paste0("[", round(rd$lower_rank,1), ", ", round(rd$upper_rank,1), "]"))),
        tags$tr(tags$td(strong("Relative Risk (RR)")), tags$td(round(rd$RR, 2))),
        tags$tr(tags$td(strong("P(higher risk)")), tags$td(round(rd$prob_higher, 3))),
        tags$tr(tags$td(strong("P(lower risk)")), tags$td(round(rd$prob_lower, 3))),
        tags$tr(tags$td(strong("Evidence level")), tags$td(rd$Evidence)),
        tags$tr(tags$td(strong("% Single carriageway")), tags$td(round(rd$pct_single,0))),
        tags$tr(tags$td(strong("% Motorway")), tags$td(round(rd$pct_motorway,0)))
      ),
      tags$div(style = paste0(
        "padding:6px 10px;border-radius:4px;margin-top:8px;color:white;background-color:",
        ifelse(rd$Conclusion == "Significantly HIGHER RISK", "#d9534f",
               ifelse(rd$Conclusion == "Significantly LOWER RISK", "#2c7a4b", "#888"))
      ), rd$Conclusion)
    )
  })

  output$eda_trend_plot <- renderPlotly({
    df <- panel
    df$highlight <- ifelse(df$County %in% input$eda_counties, df$County, "Other")
    p <- ggplot(df, aes(x = Year, y = .data[[input$eda_metric]],
                         color = highlight, group = County,
                         text = paste0(County, ", ", Year, ": ",
                                       round(.data[[input$eda_metric]], 2)))) +
      geom_line(data = df %>% filter(highlight == "Other"),
                color = "gray80", linewidth = 0.4) +
      geom_line(data = df %>% filter(highlight != "Other"), linewidth = 1.1) +
      geom_point(data = df %>% filter(highlight != "Other"), size = 1.6) +
      scale_x_continuous(breaks = 2015:2024) +
      labs(x = "Year", y = input$eda_metric, color = "County") +
      theme_minimal(base_size = 12)
    ggplotly(p, tooltip = "text")
  })

  output$eda_national_plot <- renderPlotly({
    nat <- panel %>% group_by(Year) %>%
      summarise(rate = sum(deaths, na.rm = TRUE) / (sum(km, na.rm = TRUE) / 1e9),
                .groups = "drop")
    p <- ggplot(nat, aes(x = Year, y = rate,
                         text = paste0("Year: ", Year, "<br>Rate: ", round(rate,2)))) +
      geom_line(color = "#d73027", linewidth = 1.2) +
      geom_point(color = "#d73027", size = 2.5) +
      geom_vline(xintercept = 2020, linetype = "dashed", color = "gray50") +
      scale_x_continuous(breaks = 2015:2024) +
      labs(x = "Year", y = "Deaths per billion km") +
      theme_minimal(base_size = 12)
    ggplotly(p, tooltip = "text")
  })

  output$eda_heatmap <- renderPlot({
    ggplot(panel, aes(x = Year, y = County, fill = rate_billion)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_viridis(option = "C", name = "Deaths per\nbillion km") +
      scale_x_continuous(breaks = 2015:2024) +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$posterior_table <- renderDT({
    ps <- c("beta0","beta1","beta2","sd0")
    pt <- data.frame(
      Parameter = c("\u03b2\u2080 (baseline)","\u03b2\u2081 (trend)",
                    "\u03b2\u2082 (COVID)","\u03c3\u2080 (county SD)"),
      Mean   = sapply(ps, function(p) round(mean(pm3[,p]), 4)),
      L95    = sapply(ps, function(p) round(quantile(pm3[,p],0.025), 4)),
      U95    = sapply(ps, function(p) round(quantile(pm3[,p],0.975), 4)),
      Pgt0   = sapply(ps, function(p) round(mean(pm3[,p] > 0), 4))
    )
    datatable(pt, rownames = FALSE, options = list(dom = 't'))
  })

  output$convergence_table <- renderDT({
    ps <- c("beta0","beta1","beta2","sd0")
    gd <- gelman.diag(samp_m3)
    es <- effectiveSize(samp_m3)
    ct <- data.frame(
      Parameter = ps,
      Rhat      = round(gd$psrf[ps,1], 4),
      ESS       = round(es[ps]),
      Passed    = ifelse(gd$psrf[ps,1] < 1.05 & es[ps] > 1000, "YES", "CHECK")
    )
    datatable(ct, rownames = FALSE, options = list(dom = 't'))
  })

  output$trace_plot <- renderPlot({
    p <- input$trace_param
    ad <- pm3[, p]
    chains <- lapply(1:4, function(i) ad[((i-1)*n_iter+1):(i*n_iter)])
    cols4 <- c("blue","red","darkgreen","purple")
    plot(chains[[1]], type = "l", col = cols4[1],
         main = paste("Trace:", p), xlab = "Iteration", ylab = p,
         ylim = range(ad))
    for (i in 2:4) lines(chains[[i]], col = cols4[i])
    legend("topright", legend = paste("Chain", 1:4), col = cols4, lwd = 1.5, cex = 0.8, bty = "n")
  })

  output$density_plot <- renderPlot({
    p <- input$trace_param
    ad <- pm3[, p]
    chains <- lapply(1:4, function(i) ad[((i-1)*n_iter+1):(i*n_iter)])
    dl <- lapply(chains, density)
    xlim <- range(sapply(dl, function(d) range(d$x)))
    ylim <- range(sapply(dl, function(d) range(d$y)))
    cols4 <- c("blue","red","darkgreen","purple")
    plot(dl[[1]], col = cols4[1], lwd = 2, main = paste("Density:", p),
         xlab = p, xlim = xlim, ylim = ylim)
    for (i in 2:4) lines(dl[[i]], col = cols4[i], lwd = 2)
    abline(v = 0, lty = 3)
    legend("topright", legend = paste("Chain", 1:4), col = cols4, lwd = 1.5, cex = 0.8, bty = "n")
  })

  output$prior_post_plot <- renderPlot({
    set.seed(123); ns <- 80000
    b0p <- rnorm(ns, 15.015, 1/sqrt(9.905))
    b1p <- rnorm(ns, -0.002, 1/sqrt(84.665))
    b2p <- rnorm(ns, 0.053,  1/sqrt(1.989))
    par(mfrow = c(1,3), mar = c(4,3,3,1))
    for (i in 1:3) {
      param <- c("beta0","beta1","beta2")[i]
      pv  <- list(b0p, b1p, b2p)[[i]]
      pov <- pm3[, param]
      xl  <- quantile(c(pv, pov), c(0.001, 0.999))
      plot(density(pv), col = "red", lwd = 2, main = param, xlab = param, xlim = xl)
      lines(density(pov), col = "blue", lwd = 2)
      abline(v = 0, lty = 3)
      legend("topright", legend = c("Prior","Posterior"), col = c("red","blue"), lwd = 2, cex = 0.8)
    }
  })

  output$rank_caterpillar <- renderPlotly({
    p <- ggplot(rr, aes(x = mean_b0, y = reorder(County, mean_b0), color = Conclusion,
                        text = paste0(County, "<br>b0: ", round(mean_b0,3),
                                      "<br>RR: ", round(RR,2)))) +
      geom_errorbarh(aes(xmin = lower_b0, xmax = upper_b0), height = 0.3) +
      geom_point(size = 2.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("Significantly HIGHER RISK" = "red",
                                     "Significantly LOWER RISK" = "darkgreen",
                                     "No significant difference" = "gray50")) +
      labs(x = "Posterior mean b0 with 95% CrI", y = "County") +
      theme_minimal(base_size = 10)
    ggplotly(p, tooltip = "text")
  })

  output$rank_ci_plot <- renderPlotly({
    p <- ggplot(rr, aes(x = mean_rank, y = reorder(County, -mean_rank),
                        text = paste0(County, "<br>Mean rank: ", round(mean_rank,1),
                                      "<br>95% CrI: [", round(lower_rank,1), ", ",
                                      round(upper_rank,1), "]"))) +
      geom_errorbarh(aes(xmin = lower_rank, xmax = upper_rank), height = 0.3, color = "steelblue") +
      geom_point(size = 2.2, color = "darkblue") +
      geom_vline(xintercept = 13.5, linetype = "dashed", color = "gray50") +
      labs(x = "Posterior mean rank (1=riskiest, 26=safest)", y = "County") +
      theme_minimal(base_size = 10)
    ggplotly(p, tooltip = "text")
  })

  output$rankings_table <- renderDT({
    tab <- rr %>%
      select(County, mean_b0, lower_b0, upper_b0, mean_rank, lower_rank, upper_rank,
             RR, Evidence, Conclusion) %>%
      mutate(across(where(is.numeric), round, 3))
    datatable(tab, rownames = FALSE,
              colnames = c("County","Mean b0","Lower 95%","Upper 95%",
                           "Mean rank","Lower rank","Upper rank","RR","Evidence","Conclusion"),
              options = list(pageLength = 26))
  })

  output$road_scatter <- renderPlotly({
    rate_road <- panel %>% group_by(County) %>%
      summarise(adr = mean(deaths/(km/1e9), na.rm = TRUE)/1e6, .groups = "drop") %>%
      left_join(road_data, by = "County")
    p <- ggplot(rate_road, aes(x = pct_single, y = adr, text = County, color = pct_motorway)) +
      geom_point(size = 3) +
      geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
      scale_color_gradient(low = "red", high = "darkgreen", name = "% Motorway") +
      labs(x = "% single carriageway", y = "Avg deaths per million km") +
      theme_minimal(base_size = 11)
    ggplotly(p, tooltip = c("text","x","y"))
  })

  output$road_combined <- renderPlotly({
    p <- ggplot(rr, aes(x = mean_b0, y = reorder(County, mean_b0), color = pct_single,
                        size = avg_km, text = paste0(County, "<br>b0: ", round(mean_b0,3)))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(alpha = 0.85) +
      scale_color_gradient(low = "darkgreen", high = "red", name = "% Single\ncarriageway") +
      labs(x = "Posterior mean b0", y = "County") +
      theme_minimal(base_size = 10)
    ggplotly(p, tooltip = "text")
  })

  output$road_table <- renderDT({
    datatable(road_data %>%
                select(County, motorway_km, dual_km, single_km, total_km,
                       pct_motorway, pct_single) %>%
                mutate(across(where(is.numeric), round, 1)),
              rownames = FALSE,
              colnames = c("County","Motorway km","Dual km","Single km",
                           "Total km","% Motorway","% Single"))
  })


  # Pre-fitted DIC table, saved locally with saveRDS(dic_table, "dic_results.rds")
  # after running the comparison once on your own machine. The deployed app
  # reads this directly -- it never fits JAGS on the server.
  cached_dic <- if (file.exists("dic_results.rds")) readRDS("dic_results.rds") else NULL

  output$dic_intro <- renderUI({
    if (!is.null(cached_dic)) {
      p("The table below shows the DIC comparison across 5 nested models (M1-M5), pre-computed locally and bundled with this app. Model M3 -- county baseline effects only -- has the lowest DIC and is therefore the model used throughout the rest of this dashboard.")
    } else if (jags_available) {
      tagList(
        p("This refits 5 nested JAGS models (M1 to M5) and can take several minutes."),
        actionButton("run_dic", "Run DIC comparison", icon = icon("play"), class = "btn-warning")
      )
    } else {
      p(em("DIC results are not available: no cached dic_results.rds was found, and JAGS is not available in this deployment environment to fit the comparison live. Run the DIC comparison locally and redeploy with dic_results.rds included to enable this tab."))
    }
  })

  dic_result <- eventReactive(input$run_dic, {

    jd <- list(X = X, VM = VM_vec, Year_c = Year_c, covid = covid,
               county = county, nObs = nObs, N = N)

    rdic <- function(ms, nm) {
      m <- jags.model(textConnection(ms), data = jd, n.chains = 4,
                       n.adapt = 1000, quiet = TRUE)
      update(m, 5000, progress.bar = "none")
      d <- dic.samples(m, n.iter = 10000, type = "pD", progress.bar = "none")
      data.frame(Model = nm, DIC = sum(d$deviance) + sum(d$penalty),
                 pD = sum(d$penalty), Deviance = sum(d$deviance))
    }

    withProgress(message = "Fitting 5 models for DIC comparison...", value = 0, {
      incProgress(0.2, detail = "M1: Full model")
      r1 <- rdic("model{for(i in 1:nObs){X[i]~dpois(mu[i])
        log(mu[i])<-log(VM[i])+(beta0+b0[county[i]])+
        (beta1+b1[county[i]])*Year_c[i]+(beta2+b2[county[i]])*covid[i]}
        for(c in 1:N){b0[c]~dnorm(0,tau0);b1[c]~dnorm(0,tau1);b2[c]~dnorm(0,tau2)}
        beta0~dnorm(15.015,9.905);beta1~dnorm(-0.002,84.665);beta2~dnorm(0.053,1.989)
        tau0~dgamma(.001,.001);tau1~dgamma(.001,.001);tau2~dgamma(.001,.001)
        sd0<-1/sqrt(tau0);sd1<-1/sqrt(tau1);sd2<-1/sqrt(tau2)}", "M1: Full (b0+b1+b2)")

      incProgress(0.2, detail = "M2: No county COVID")
      r2 <- rdic("model{for(i in 1:nObs){X[i]~dpois(mu[i])
        log(mu[i])<-log(VM[i])+(beta0+b0[county[i]])+
        (beta1+b1[county[i]])*Year_c[i]+beta2*covid[i]}
        for(c in 1:N){b0[c]~dnorm(0,tau0);b1[c]~dnorm(0,tau1)}
        beta0~dnorm(15.015,9.905);beta1~dnorm(-0.002,84.665);beta2~dnorm(0.053,1.989)
        tau0~dgamma(.001,.001);tau1~dgamma(.001,.001)
        sd0<-1/sqrt(tau0);sd1<-1/sqrt(tau1)}", "M2: No county COVID (b0+b1)")

      incProgress(0.2, detail = "M3: No county trend (best)")
      r3 <- rdic("model{for(i in 1:nObs){X[i]~dpois(mu[i])
        log(mu[i])<-log(VM[i])+(beta0+b0[county[i]])+
        beta1*Year_c[i]+beta2*covid[i]}
        for(c in 1:N){b0[c]~dnorm(0,tau0)}
        beta0~dnorm(15.015,9.905);beta1~dnorm(-0.002,84.665);beta2~dnorm(0.053,1.989)
        tau0~dgamma(.001,.001);sd0<-1/sqrt(tau0)}", "M3: No county trend (b0) BEST")

      incProgress(0.2, detail = "M4: No county effects")
      r4 <- rdic("model{for(i in 1:nObs){X[i]~dpois(mu[i])
        log(mu[i])<-log(VM[i])+beta0+beta1*Year_c[i]+beta2*covid[i]}
        beta0~dnorm(15.015,9.905);beta1~dnorm(-0.002,84.665);beta2~dnorm(0.053,1.989)}",
        "M4: No county effects NULL")

      incProgress(0.2, detail = "M5: No COVID at all")
      r5 <- rdic("model{for(i in 1:nObs){X[i]~dpois(mu[i])
        log(mu[i])<-log(VM[i])+(beta0+b0[county[i]])+
        (beta1+b1[county[i]])*Year_c[i]}
        for(c in 1:N){b0[c]~dnorm(0,tau0);b1[c]~dnorm(0,tau1)}
        beta0~dnorm(15.015,9.905);beta1~dnorm(-0.002,84.665)
        tau0~dgamma(.001,.001);tau1~dgamma(.001,.001)
        sd0<-1/sqrt(tau0);sd1<-1/sqrt(tau1)}", "M5: No COVID at all")
    })

    result <- bind_rows(r1, r2, r3, r4, r5) %>%
      mutate(across(where(is.numeric), round, 2),
             Delta_DIC = round(DIC - min(DIC), 2),
             Verdict = case_when(
               Delta_DIC == 0 ~ "BEST MODEL",
               Delta_DIC < 2  ~ "Equivalent",
               Delta_DIC < 10 ~ "Some evidence against",
               TRUE ~ "REJECTED")) %>%
      arrange(DIC)

    saveRDS(result, "dic_results.rds")
    result
  })

  output$dic_table <- renderDT({
    tab <- if (!is.null(cached_dic)) {
      cached_dic
    } else if (jags_available && !is.null(tryCatch(dic_result(), error = function(e) NULL))) {
      dic_result()
    } else {
      NULL
    }

    if (is.null(tab)) {
      return(datatable(data.frame(Message = "No DIC results available yet.")))
    }
    datatable(tab, rownames = FALSE,
              colnames = c("Model","DIC","pD","Deviance","Delta DIC","Verdict"))
  })

  output$county_summary <- renderTable({
    c1 <- input$county_pick
    cd <- panel %>% filter(County == c1)
    rd <- rr %>% filter(County == c1)
    data.frame(
      Metric = c("Total deaths (2015-2024)","Total injuries (2015-2024)",
                 "Mean rate (deaths/billion km)","Posterior mean b0",
                 "Posterior mean rank","Relative risk (RR)","Conclusion"),
      Value = c(sum(cd$deaths, na.rm=TRUE), sum(cd$injuries, na.rm=TRUE),
                round(mean(cd$rate_billion, na.rm=TRUE),2),
                round(rd$mean_b0,3), round(rd$mean_rank,1),
                round(rd$RR,2), rd$Conclusion)
    )
  })

  output$county_raw_plot <- renderPlotly({
    c1 <- input$county_pick
    cd <- panel %>% filter(County == c1) %>%
      pivot_longer(cols = c(deaths, injuries), names_to = "type", values_to = "count")
    p <- ggplot(cd, aes(x = Year, y = count, color = type)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      scale_x_continuous(breaks = 2015:2024) +
      labs(x = "Year", y = "Count", color = "") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$county_rate_plot <- renderPlotly({
    c1 <- input$county_pick
    cd <- panel %>% filter(County == c1)
    p <- ggplot(cd, aes(x = Year, y = rate_billion)) +
      geom_line(color = "#d73027", linewidth = 1) +
      geom_point(color = "#d73027", size = 2) +
      scale_x_continuous(breaks = 2015:2024) +
      labs(x = "Year", y = "Deaths per billion km") +
      theme_minimal(base_size = 11)
    ggplotly(p)
  })

  output$pairwise_result <- renderPrint({
    a <- input$county_pick; b <- input$county_pick2
    if (a == b) {
      cat("Please select two different counties.")
      return(invisible(NULL))
    }
    ia <- which(county_lookup$County == a)
    ib <- which(county_lookup$County == b)
    p_a <- mean(b0_samps[, ia] > b0_samps[, ib])
    cat(sprintf("P(%s riskier than %s) = %.4f\n", a, b, p_a))
    cat(sprintf("P(%s riskier than %s) = %.4f\n", b, a, 1 - p_a))
    if (p_a > 0.975) {
      cat(sprintf("\n--> %s is significantly riskier than %s.\n", a, b))
    } else if ((1 - p_a) > 0.975) {
      cat(sprintf("\n--> %s is significantly riskier than %s.\n", b, a))
    } else {
      cat("\n--> No significant difference between these two counties.\n")
    }
  })
}

shinyApp(ui = ui, server = server)
