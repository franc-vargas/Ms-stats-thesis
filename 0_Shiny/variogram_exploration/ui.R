#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinythemes)
library(shinydashboard)
library(leaflet)
# Define UI for application that draws a histogram
navbarPage("Dashboard", theme = shinytheme("flatly"),
           tabPanel("Los Lagos",
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("lagos.breaks",
                        "Number of breaks:",
                        min = 1,
                        max = 1000,
                        value = 30),
            sliderInput("lagos.distmax",
                        "Maximum Distance:",
                        min = 1,
                        max = 300,
                        value = 100),
            selectInput("lagos.Response", 
                        label = "Select Response Variable",
                        choices = c(
                          "Residuals" = "residuos", 
                          "Residuals of Log Total Adults LM" = "residuos.log",
                          "Avg. Max. Total Load per Cycle" = "max.carga_total",
                          "Avg. Max. Total Adults per Cycle"= "max.adultos_totales",
                          "Log-Transformed Avg. Max. Total Load per Cycle" = "log_max.carga_total",
                          "Log-Transformed Avg. Max. Total Adults per Cycle"= "log_max.adultos_totales") )
        ),
        mainPanel(
          fluidRow(
            splitLayout(
              plotOutput("LC_Semivariogram.Lagos"), 
              plotOutput("LC_Robust_Semivariogram.Lagos")
              )),
          fluidRow(
            splitLayout(plotOutput("Geodesic.Semivariogram.Lagos"),
                        plotOutput("Geodesic.Robust.Lagos")))
          )
        )
    ),
    tabPanel("Aysen",
             # Sidebar with a slider input for number of bins
             sidebarLayout(
               sidebarPanel(
                 sliderInput("aysen.breaks",
                             "Number of breaks:",
                             min = 1,
                             max = 1000,
                             value = 30),
                 sliderInput("aysen.distmax",
                             "Maximum Distance:",
                             min = 1,
                             max = 300,
                             value = 100),
                 selectInput("aysen.Response", 
                             label = "Select Response Variable",
                             choices = c(
                               "Residuals of Total Adults LM" = "residuos",
                               "Residuals of Log Total Adults LM" = "residuos.log",
                               "Avg. Max. Total Load per Cycle" = "max.carga_total",
                               "Avg. Max. Total Adults per Cycle"= "max.adultos_totales",
                               "Log-Transformed Avg. Max. Total Load per Cycle" = "log_max.carga_total",
                               "Log-Transformed Avg. Max. Total Adults per Cycle"= "log_max.adultos_totales"))               ),
               mainPanel(
                 fluidRow(
                   splitLayout(
                     plotOutput("LC_Semivariogram.Aysen"), 
                     plotOutput("LC_Robust_Semivariogram.Aysen")
                   )),
                 fluidRow(
                   splitLayout(plotOutput("Geodesic.Semivariogram.Aysen"),
                               plotOutput("Geodesic.Robust.Aysen")))
               )
             )
    ),
    tabPanel("Map",
             leafletOutput("leaflet.map", height = 800),
             )
             
             
    )
