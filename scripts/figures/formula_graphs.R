x <- 1:100
y <- exp(-(.5 * (4*(x/100))))

write_tsv(tibble(x, y), path = 'exp_fun_distance.tsv')

plot(x, y,
     xlab = 'x = 1:100',
     ylab = 'y = exp(-(0.5 * (4 * (x/100))))')

set.seed(123)
g <- sample(1:20, 100, replace = TRUE)

z <- aggregate(y, by = list(g = g), sum)


write_tsv(tibble(x = as.numeric(table(g)), y = z$x), 'agg_sum.tsv')

plot(as.numeric(table(g)), z$x,
     xlab = 'N_g of y (g = 1:20)',
     ylab = 'aggregate(y, by = g, sum)')

set.seed(123)
t <- rnorm(40, 0, 2)

write_tsv(tibble(x = c(-z$x, z$x), y = t), 'distance_stats.tsv')

plot(c(-z$x, z$x), t,
     xlab = 'c(-y, y)',
     ylab = 't = rnorm(40, 0, 2)')

r <- (rank(-abs(c(-z$x, z$x))) * rank(-abs(t)))/40^2
min(r)

points(c(-z$x, z$x)[which(r == min(r))],
       t[which(r == min(r))],
       col = 'red')

points(c(-z$x, z$x)[which(r == max(r))],
       t[which(r == max(r))],
       col = 'blue')
