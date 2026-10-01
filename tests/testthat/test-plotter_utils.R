test_that("signif_to_zero truncates to the requested number of significant digits", {
  expect_equal(signif_to_zero(123.456, digits = 4), 123.4)
  expect_equal(signif_to_zero(0.012345, digits = 2), 0.012)
})

test_that("signif_to_zero preserves sign", {
  expect_true(signif_to_zero(-5.6789, digits = 3) < 0)
  expect_true(signif_to_zero(5.6789, digits = 3) > 0)
  expect_equal(
    abs(signif_to_zero(-5.6789, digits = 3)),
    signif_to_zero(5.6789, digits = 3)
  )
})

test_that("scale_color_de_gradient returns a ggplot2 continuous colour scale", {
  sc <- scale_color_de_gradient(abs_max = 2)
  expect_s3_class(sc, "ScaleContinuous")
  expect_equal(sc$limits, c(-2, 2))
})

test_that("scale_color_de_gradient respects custom limits and breaks", {
  sc <- scale_color_de_gradient(abs_max = 5, limits = c(-1, 1), breaks = c(-1, 0, 1))
  expect_equal(sc$limits, c(-1, 1))
  expect_equal(sc$breaks, c(-1, 0, 1))
})

test_that("small_axis returns coord/theme/annotation components by default", {
  components <- small_axis(label = "PC1")
  expect_type(components, "list")
  expect_length(components, 4)
  expect_false(is.null(components[[1]])) # coord_fixed
  expect_false(is.null(components[[2]])) # axis theme
  expect_false(is.null(components[[3]])) # arrow annotation
  expect_false(is.null(components[[4]])) # label annotation
})

test_that("small_axis omits optional components when disabled", {
  components <- small_axis(label = NULL, fix_coord = FALSE, remove_axes = FALSE)
  expect_null(components[[1]]) # no coord_fixed requested
  expect_null(components[[2]]) # axis theme not removed
  expect_false(is.null(components[[3]])) # the axis lines are always drawn
  expect_null(components[[4]]) # no label supplied
})
