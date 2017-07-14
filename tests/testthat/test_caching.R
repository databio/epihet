context("caching")

test_that("cached binomial returns object as expected", {

    confdt = cacheBinomConfIntervals(200,200,.95)

    expect_is(confdt, "data.table")

})
