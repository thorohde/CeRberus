test_that("position-agnostic aggregate slots have typed defaults", {
  screen <- methods::new("PosAgnMultiplexScreen")

  expect_s4_class(aggregatedGuideGIs(screen), "gRNA_GI")
  expect_identical(aggregatedLimmaModels(screen), list())
})

test_that("position-agnostic aggregate accessors round-trip values", {
  screen <- methods::new("PosAgnMultiplexScreen")
  guide_gis <- methods::new("gRNA_GI")
  models <- list(A = list(fit = "aggregate"))

  aggregatedGuideGIs(screen) <- guide_gis
  aggregatedLimmaModels(screen) <- models

  expect_identical(aggregatedGuideGIs(screen), guide_gis)
  expect_identical(aggregatedLimmaModels(screen), models)
})

test_that("aggregate accessors reject non-position-agnostic screens", {
  screen <- methods::new("MultiplexScreen")

  expect_error(
    aggregatedGuideGIs(screen),
    "only available for PosAgnMultiplexScreen"
  )
  expect_error(
    aggregatedLimmaModels(screen),
    "only available for PosAgnMultiplexScreen"
  )
})
