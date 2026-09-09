test_that("label_ntc_pairs adds an exact binary NTC flag", {
  results <- data.table::data.table(
    query_gene = c("NTC", "TARGET", "NTC_like", "TARGET"),
    library_gene = c("TARGET", "NTC", "TARGET", "TARGET")
  )

  expect_identical(
    label_ntc_pairs(results, "NTC")$has_NTC,
    c(TRUE, TRUE, FALSE, FALSE)
  )
  expect_false("has_NTC" %in% names(results))
})

test_that("label_ntc_pairs returns false when no NTCs are configured", {
  results <- data.frame(query_gene = "A", library_gene = "B")

  expect_false(label_ntc_pairs(results)$has_NTC)
})
