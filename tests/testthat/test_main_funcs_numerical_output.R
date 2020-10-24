
# testthat::test_dir("tests/")
# NB devtools::test() crashes help index & need to reload package

context("check numerical output of functions in main.R")

# NB test includes variation numberical values &
# in how d1,d2 with capitals or accents
# guess should check separately with specific tests

test_that("av_long_infected works", {

  t01 <- av_long_infected(
    a1 = 3.0, b1 = 0.6, d1 = "Weibull",
    a2 = 3.0, b2 = 0.6, d2 = "weibull")

  t02 <- av_long_infected(
    a1 = 2.0, b1 = 1.0, d1 = "Fréchet",
    a2 = 2.5, b2 = 1.0, d2 = "frechet")

  t03 <- av_long_infected(
    a1 = 23.0, b1 = 4.0, d1 = "Gumbel",
    a2 = 22.0, b2 = 5.0, d2 = "gumbel")

  expect_equal(t01$value, 11.840, tolerance = 0.01)
  expect_equal(t02$value, 12.973, tolerance = 0.01)
  expect_equal(t03$value, 16.795, tolerance = 0.01)
}
)

test_that("av_long_uninfected works", {

  t01 <- av_long_uninfected(
    a1 = 3.0, b1 = 0.6, d1 = "Weibull")

  t02 <- av_long_uninfected(
    a1 = 1.0, b1 = 0.5, d1 = "Fréchet")

  t03 <- av_long_uninfected(
    a1 = 23.0, b1 = 4.0, d1 = "Gumbel")

  expect_equal(t01$value, 17.947, tolerance = 0.01)
  expect_equal(t02$value,  4.818, tolerance = 0.01)
  expect_equal(t03$value, 20.704, tolerance = 0.01)
}
)

test_that("etd_infected works", {

  t01 <- etd_infected(a1 = 2, b1 = 0.5, a2 = 30, b2 = 5,
                      d1 = "Weibull", d2 = "Gumbel", tmax = 100)

  expect_equal(t01, 7.34, tolerance = 0.1)

}
)

test_that("etd_uninfected works", {

  t01 <- etd_uninfected(a1 = 20, b1 = 5,
                        d1 = "Gumbel", tmax = 100)

  expect_equal(t01, 20.0, tolerance = 0.1)

}
)


# checks output is a matrix &
# first & last values of t, h, sd[h], ci95-,ci95+ are > 0
# values should be > 0
# but can get 'Inf' with some values
# so not sure what testing here, as choose values to test

test_that("conf_inf_virulence works", {

  tmax = 15

  t01 <- conf_ints_virulence(
    a2 = 2, b2 = 0.5,
    var_a2 = 0.001, var_b2 = 0.001,
    cov_a2b2 = -0.001,
    d2 = "Weibull", tmax = tmax)

  expect_is(t01, 'matrix')

  expect_gt(t01[1,1], 0)
  expect_gt(t01[1,2], 0)
  expect_gt(t01[1,5], 0)
  expect_gt(t01[1,6], 0)
  expect_gt(t01[1,7], 0)

  expect_gt(t01[tmax,2], 0)
  expect_gt(t01[tmax,5], 0)
  expect_gt(t01[tmax,6], 0)
  expect_gt(t01[tmax,7], 0)

}
)

