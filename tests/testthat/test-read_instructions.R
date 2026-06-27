write_instruction_file <- function(path, scores_file, output_directory, ...) {
  yaml::write_yaml(
    c(
      list(
        scores_file = scores_file,
        output_directory = output_directory
      ),
      list(...)
    ),
    file = path
  )
}

test_that("read_instructions reads a valid YAML instruction file", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempdir()
  yaml_fpath <- tempfile(fileext = ".yaml")
  file.create(scores_file)

  write_instruction_file(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    FDR = "bonferroni",
    overwrite_output = FALSE,
    verbose = TRUE
  )

  result <- read_instructions(yaml_fpath)

  expect_type(result, "list")
  expect_equal(result$scores_file, scores_file)
  expect_equal(
    result$output_directory,
    normalizePath(
      output_directory,
      winslash = "/",
      mustWork = FALSE
    )
  )
  expect_equal(result$FDR, "bonferroni")
  expect_false(result$overwrite_output)
  expect_true(result$verbose)
})

test_that("read_instructions validates the YAML file path", {
  expect_error(
    read_instructions(character()),
    "yaml_fpath must be a single path string"
  )
  expect_error(
    read_instructions(NA_character_),
    "yaml_fpath must be a single path string"
  )
  expect_error(
    read_instructions(""),
    "yaml_fpath must be a single path string"
  )
  expect_error(
    read_instructions(tempfile(fileext = ".yaml")),
    "Instruction file does not exist"
  )
})

test_that("read_instructions requires a YAML mapping/object", {
  yaml_fpath <- tempfile(fileext = ".yaml")
  writeLines("- one\n- two", yaml_fpath)

  expect_error(
    read_instructions(yaml_fpath),
    "Instruction file must contain a YAML mapping/object"
  )
})

test_that("read_instructions validates required scores_file", {
  yaml_fpath <- tempfile(fileext = ".yaml")
  yaml::write_yaml(list(output_directory = tempdir()), file = yaml_fpath)

  expect_error(
    read_instructions(yaml_fpath),
    "Please provide a scores file using the 'scores_file' argument"
  )

  yaml::write_yaml(
    list(
      scores_file = tempfile(fileext = ".csv"),
      output_directory = tempdir()
    ),
    file = yaml_fpath
  )

  expect_error(
    read_instructions(yaml_fpath),
    "Given scores file does not exist"
  )
})

test_that("read_instructions validates required output_directory", {
  scores_file <- tempfile(fileext = ".csv")
  yaml_fpath <- tempfile(fileext = ".yaml")
  file.create(scores_file)

  yaml::write_yaml(list(scores_file = scores_file), file = yaml_fpath)

  expect_error(
    read_instructions(yaml_fpath),
    "Output directory required"
  )

  yaml::write_yaml(
    list(scores_file = scores_file, output_directory = ""),
    file = yaml_fpath
  )

  expect_error(
    read_instructions(yaml_fpath),
    "Output directory must not be empty"
  )
})

test_that("read_instructions defaults invalid or missing FDR to BH", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempdir()
  yaml_fpath <- tempfile(fileext = ".yaml")
  file.create(scores_file)

  write_instruction_file(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory
  )

  expect_equal(read_instructions(yaml_fpath)$FDR, "BH")

  write_instruction_file(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    FDR = "invalid"
  )

  expect_equal(read_instructions(yaml_fpath)$FDR, "BH")
})

test_that("read_instructions parses boolean-like values", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- tempdir()
  yaml_fpath <- tempfile(fileext = ".yaml")
  file.create(scores_file)

  write_instruction_file(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = "no",
    verbose = "yes"
  )

  result <- read_instructions(yaml_fpath)

  expect_false(result$overwrite_output)
  expect_true(result$verbose)

  write_instruction_file(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory,
    overwrite_output = "not-a-bool",
    verbose = "not-a-bool"
  )

  result <- read_instructions(yaml_fpath)

  expect_true(result$overwrite_output)
  expect_false(result$verbose)
})

test_that("read_instructions returns output_directory without trailing separators", {
  scores_file <- tempfile(fileext = ".csv")
  output_directory <- file.path(tempdir(), "")
  yaml_fpath <- tempfile(fileext = ".yaml")
  file.create(scores_file)

  write_instruction_file(
    yaml_fpath,
    scores_file = scores_file,
    output_directory = output_directory
  )

  result <- read_instructions(yaml_fpath)

  expect_false(grepl("[/\\\\]$", result$output_directory))
})
