test_that("HolzingerSwineford works", {
  testthat::skip_on_cran()
  library(lessTREE)

  set.seed(123)

  ## The famous Holzinger and Swineford (1939) example
  HS.model <- ' visual  =~ 1*x1 + l21*x2 + l31*x3
              textual =~ 1*x4 + l52*x5 + l62*x6
              speed   =~ 1*x7 + l83*x8 + l93*x9
              visual ~~ v1*visual + v12*textual + v13*speed
textual ~~ v2*textual + v23*speed
speed ~~ v3*speed
x1 ~~ d1*x1
x2 ~~ d2*x2
x3 ~~ d3*x3
x4 ~~ d4*x4
x5 ~~ d5*x5
x6 ~~ d6*x6
x7 ~~ d7*x7
x8 ~~ d8*x8
x9 ~~ d9*x9
'

  fit <- cfa(HS.model, data = HolzingerSwineford1939[complete.cases(HolzingerSwineford1939),])

  ctrl_tree <- semtree.control(method = "score", min.bucket = 10)

  tree <- semtree(model = fit,
                  data = HolzingerSwineford1939[complete.cases(HolzingerSwineford1939),],
                  control = ctrl_tree,
                  predictors = c("sex", "school"))

  rtree <- regularize_semtree(tree = tree) |>
    semtree_lasso(lambdas = seq(0.1,1,.1))

  plot.lessTREE(rtree)

  out <- select_final(rtree, criterion = "BIC")

  # let's test a model, where some parameters have names that are not
  # allowed:
  ## The famous Holzinger and Swineford (1939) example
  HS.model <- ' visual  =~ 1*x1 + l21*x2 + l31*x3
              textual =~ 1*x4 + l52*x5 + l62*x6
              speed   =~ 1*x7 + l83*x8 + l93*x9
              visual ~~ v1*visual + v12*textual + speed
textual ~~ v2*textual + v23*speed
speed ~~ v3*speed
x1 ~~ x1
x2 ~~ x2
x3 ~~ d3*x3
x4 ~~ d4*x4
x5 ~~ d5*x5
x6 ~~ d6*x6
x7 ~~ d7*x7
x8 ~~ x8
x9 ~~ d9*x9
'
  fit <- cfa(HS.model, data = HolzingerSwineford1939[complete.cases(HolzingerSwineford1939),])

  ctrl_tree <- semtree.control(method = "score", min.bucket = 10)

  tree <- semtree(model = fit,
                  data = HolzingerSwineford1939[complete.cases(HolzingerSwineford1939),],
                  control = ctrl_tree,
                  predictors = c("sex", "school"))

  rtree <- regularize_semtree(tree = tree) |>
    semtree_lasso(lambdas = seq(0.1,1,.1))

  out <- select_final(rtree, criterion = "BIC")

  plot.lessTREE(rtree)

  })
