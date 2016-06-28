context("qseaSet")
qs <- getExampleQseaSet()
test_that("getExampleQseaSet works properly", {
    expect_is(qs, "qseaSet")
})

