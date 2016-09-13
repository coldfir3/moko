test_that('mkm can build properly the objects',{

n <- 10
d <- 2
doe <- replicate(d,sample(0:n,n))/n
res <- t(apply(doe, 1, nowacki_beam))
model <- mkm(doe, res, modelcontrol = list(objective = 1:2))

expect_is(model, 'mkm')
expect_equal(model@d, d)
expect_equal(model@n, n)
expect_equal(model@j, 5)
expect_equal(model@m, 2)

})
