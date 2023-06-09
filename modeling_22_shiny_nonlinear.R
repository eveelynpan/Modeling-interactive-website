# Demonstration of building curves to fit data using parametric and
# non-parametric fitting algorithms.
# Author: Greg Hamerly (hamerly@cs.baylor.edu)
#
# Disclaimer: this code is not robust against errors and is for demonstration
# only.

library(shiny)
library(glmnet)

########################
# Pre-Load Data Sets (try from libraries, fall back to disk)
########################

tryCatch({
  library("mowateR")
  data("eml")
  data("tds_rivers")
  data("boulder_ammonia")
  print("Loaded mowater datasets")
}, error=function(e) {
  print(e)
  print("Could not load mowater datasets from library... trying loading from disk")
  load('dat/eml.rda')
  load('dat/tds_rivers.rda')
  load('dat/boulder_ammonia.rda')
})

tryCatch({
  library("fields")
  data("WorldBankCO2")
}, error=function(e) {
  print(e)
  print("Could not load fields dataset from library... trying loading from disk")
  load('dat/worldbankco2.rda')
})

# Global constants
MAX_POLY_ORDER = 20 # largest polynomial degree we will try
MAX_SAMPLE_SIZE = 1000 # largest data sample (to keep the app responsive)

# Metadata about some datasets that can be loaded -- used to construct the UI
# and load data
datasets <- list(
     eml = c(name="Eagle Mountain Lake", dataset="eml", xname="pH", yname="DO"),
     co2 = c(name="World Bank CO2", dataset="WorldBankCO2", xname="GDP.cap", yname="CO2.cap"),
     tds_rivers = c(name="TDS Rivers", dataset="tds_rivers", xname="lees", yname="hoover"),
     boulder_ammonia = c(name="Boulder Ammonia", dataset="boulder_ammonia", xname="AB3.Z7.DO.mg.L", yname="AB3.Z7.Ammonia.mg.N.L")
)

# set up a vector of datasets to select for the UI
dataset_select_keys <- NULL
for (key in names(datasets)) {
    ds <- datasets[[key]]
    dataset_select_keys[ds['name']] <- key
}

# Global variables indicating which dataset is currently loaded
dataset <- NULL
datasetkey <- "none"

########################
# Helper Functions
########################

# Load a dataset from a library
loadDataset <- function(name, key) {
    # load the dataset by name
    data(list=name)
    d <- get(name)
    d <- as.data.frame(d) # treat it as a data frame

    # only keep numeric variables...
    numeric_vars <- sapply(d, is.numeric)
    d <- d[,numeric_vars]

    # update our global variables
    assign("dataset", d, inherits=T)
    assign("datasetkey", key, inherits=T)
}

# for the currently-selected dataset, choose the x and y variables
updateXandY <- function(xname, yname) {
    # take a sample if the amount of data is too large
    n <- nrow(dataset)
    s <- 1:n
    if (n > MAX_SAMPLE_SIZE) { s <- sample(n, MAX_SAMPLE_SIZE) }
    x <- dataset[s,xname]
    y <- dataset[s,yname]

    # this "if" is a kludge to avoid errors that were popping up
    if (!(is.null(x) || is.null(y))) {
        # sort the data first, to avoid having to sort it every time we plot it
        o <- order(x)
        x <- x[o]
        y <- y[o]

        assign("x", x, envir=.GlobalEnv)
        assign("y", y, envir=.GlobalEnv)

        assign("xname", xname, envir=.GlobalEnv)
        assign("yname", yname, envir=.GlobalEnv)
    }
}

# Load a dataset specified by the given key (which is a lookup into the global
# "datasets"). If the dataset is already loaded, do nothing. Return true if
# a new dataset was loaded; false otherwise.
useDataset <- function(key) {
    if (datasetkey == key) {
        F # nothing changed
    } else {
        print(sprintf("changing to dataset '%s'", key))
        ds <- datasets[[key]]
        loadDataset(ds["dataset"], key)
        # choose an initial x variable from the default
        updateXandY(ds["xname"], ds["yname"])
        T # something changed
    }
}

# Compute the R^2 from two vectors of observations and predictions
rsquared <- function(y, p) {
  # y: observations
  # p: predictions
  ssr <- sum((y - p) ^ 2)
  sst <- sum((y - mean(y)) ^ 2)
  1 - ssr / sst
}

########################
# Shiny Code
########################

ui <- fluidPage(
    headerPanel('Fitting Data with Curves'),
    fluidRow(
        column(4,
            selectInput("datasetkey", "Dataset:", dataset_select_keys),
            selectInput("xname", "X Variable:", choices=NULL),
            checkboxGroupInput("log_vars", "Log transform:", c("X" = "x", "Y" = "y"), inline=T),
            sliderInput(inputId="loess_span", label="Loess span", 
                        min=0.02, max=2, value=0.7, step=0.01),
            sliderInput(inputId="lowess_f", label="Lowess f", 
                        min=0.01, max=1, value=0.2, step=0.01),
            sliderInput(inputId="polynomial_degree", label="Polynomial degree", 
                        min=0, max=MAX_POLY_ORDER, value=1),
            sliderInput(inputId="lasso_log_lambda", label="LASSO log(lambda)", 
                        min=-50, max=30, value=0, step=0.01)
        ),
        column(8,
            plotOutput('plot1'),
            h3("R^2 values:"),
            div(tableOutput('r2_table'), style="font-size: 200%")
        )
    )
)

server <- function(input, output, session) {
    output$plot1 <- renderPlot({
        # dataset loading, wrangling, and sampling
        if (useDataset(input$datasetkey) || is.null(input$xname) || input$xname == "") {
            updateSelectInput(session=session, inputId="xname", choices=colnames(dataset))
        }

        if ((!is.null(input$xname)) && (input$xname != xname)) {
            updateXandY(input$xname, yname)
        }

        # prepare labels for the plot, and r2 vector
        xlab <- xname
        ylab <- yname
        r2 <- NULL

        # handle log scale if the user selected it
        log_x <- "x" %in% input$log_vars
        log_y <- "y" %in% input$log_vars

        # remove any observations <= 0 if we are using log scaling
        sx <- sy <- rep(T, length(x))
        if (log_x) { sx <- x > 0 }
        if (log_y) { sy <- y > 0 }
        s <- which(sx & sy)

        xx <- x[s]
        if (log_x) {
            xx <- log(xx)
            xlab <- paste(xlab, "(Log)")
        }
        yy <- y[s]
        if (log_y) {
            yy <- log(yy)
            ylab <- paste(ylab, "(Log)")
        }
        XX <- poly(xx, MAX_POLY_ORDER, raw=T)

        ##############################################################
        # evenly-spaced range of x values for making predictions and plotting
        # fit lines (at more than just the original observations)
        predict_x <- seq(min(xx), max(xx), length.out=200)
        predict_X <- data.frame(poly(predict_x, MAX_POLY_ORDER, raw=T))
        colnames(predict_X) <- paste("x", 1:ncol(predict_X), sep='')

        # plot the points and fit lines
        plot(yy ~ xx, xlab=xlab, ylab=ylab,
             cex.lab=1.5, cex=2, pch=19, col=rgb(0, 0, 0, 0.3))

        # LOESS
        loess_fit <- loess(yy ~ xx, span=input$loess_span)
        lines(predict_x, predict(loess_fit, predict_x), col="green", lwd=5)
        r2["loess"] <- rsquared(yy, predict(loess_fit))

        # LOWESS
        lowess_fitted <- lowess(xx, yy, f=input$lowess_f)
        lines(lowess_fitted, col="cyan", lwd=5)
        r2["lowess"] <- rsquared(yy, lowess_fitted$y)

        # Polynomial
        if (input$polynomial_degree == 0) {
            # special case for 0-degree polynomial
            poly_fit <- lm(yy ~ 1)
            p <- predict(poly_fit, newdata=data.frame(x1=predict_x))
        } else {
            lm_data <- data.frame(XX[,1:input$polynomial_degree])
            colnames(lm_data) <- paste("x", 1:input$polynomial_degree, sep='')
            lm_data <- cbind(lm_data, list(yy=yy))
            poly_fit <- lm(yy ~ ., data=lm_data)
            p <- predict(poly_fit, newdata=predict_X)
        }
        lines(predict_x, p, col="pink", lwd=5)
        r2["polynomial"] <- rsquared(yy, poly_fit$fitted.values)

        # LASSO (of the full polynomial)
        lasso_fit <- glmnet(XX, yy, lambda=exp(input$lasso_log_lambda))
        p <- predict(lasso_fit, as.matrix(predict_X))
        lines(predict_x, p, col="purple", lwd=5)
        r2["lasso"] <- rsquared(yy, predict(lasso_fit, XX))

        legend("topleft", 
               legend=c("loess", "lowess", "polynomial", "lasso"),
               col=c("green", "cyan", "pink", "purple"),
               lty=1, lwd=5)

        # render the R^2 values
        output$r2_table <- renderTable({ t(r2) }, colnames=T, digits=3)
    })
}

shinyApp(ui = ui, server = server)
