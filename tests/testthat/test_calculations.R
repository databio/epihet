context("calculations")

test_that("calculations return vectors of correct length", {

    skip_on_bioc()

    data("exampleBSDT", package = "RPIM")
    data("BSDTlist", package = "RPIM")

    # only first line of message is necessary in this case
    pimres = PIM(exampleBSDT, cacheDir = "~")
    rpimres = RPIM(BSDTlist, cacheDir = "~")

    # PIM
    expect_is(pimres, "numeric")
    expect_length(pimres, 1)

    # RPIM
    expect_is(rpimres, "numeric")
    expect_length(rpimres, length(BSDTlist))

})
