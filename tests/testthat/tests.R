library(RPIM)

context("object handling")

test_that("accept allowed objects for IM calculations", {

    # only first line of message is necessary in this case
    expect_error(RPIM("foo"), "the following are allowed:\n")

})

