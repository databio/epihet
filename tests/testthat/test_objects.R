context("object handling")

test_that("accept allowed objects for IM calculations", {

    # only first line of message is necessary in this case
    expect_error(RPIM("foo"), NULL)

})

test_that("PIM and RPIM warn based on detected number of samples", {

    skip_on_bioc()

    data("exampleBSDT", package = "epihet")
    data("BSDTlist", package = "epihet")

    # PIM errors when given multiple samples
    expect_error(PIM(BSDTlist, cacheDir = "~"), NULL)

    # RPIM errors when given single sample
    expect_error(RPIM(exampleBSDT, cacheDir = "~"), NULL)

})
