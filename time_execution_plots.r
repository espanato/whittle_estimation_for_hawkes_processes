library(hawkesbow)
library(Rcpp)


## The code below was used to plot Figure 3.3 of our paper.

# bin = 0.25
l3 <- c()
l2 <- c()
l1 <- c()
somme1 <- 0
somme2 <- 0
somme3 <- 0

# bin = 1
l_11 <- c()
l_21 <- c()
l_31 <- c()
s_11 <- 0
s_21 <- 0
s_31 <- 0

# bin = 2
l_12 <- c()
l_22 <- c()
l_32 <- c()
s_12 <- 0
s_22 <- 0
s_32 <- 0

# mle
mle1 <- c()
mle2 <- c()
mle3 <- c()
sm1 <- 0
sm2 <- 0
sm3 <- 0

# T=1000
w3 <- c()
w2 <- c()
w1 <- c()

w_11 <- c()
w_21 <- c()
w_31 <- c()

w_12 <- c()
w_22 <- c()
w_32 <- c()

me1 <- c()
me2 <- c()
me3 <- c()

# T=100
w1b <- c()
w2b <- c()
w3b <- c()

w_11b <- c()
w_21b <- c()
w_31b <- c()

w_12b <- c()
w_22b <- c()
w_32b <- c()

me1b <- c()
me2b <- c()
me3b <- c()


iter <- 100

for (end in c(100, 200, 400, 800, 1600, 3200)) {
    for (j in seq(1, iter)) {
        x <- hawkes(end, fun = 1, repr = .5, family = "exp", rate = 1)


        # mle
        mle <- mle(x$p, "Exponential", x$end)$par
        sm1 <- ((1 - mle[1])**2) + sm1
        sm2 <- ((0.5 - mle[2])**2) + sm2
        sm3 <- ((1 - mle[3])**2) + sm3

        # binsize = 0.25
        y <- discrete(x, binsize = 0.25)
        z <- whittle(y, "Exponential", 0.25)$par
        somme3 <- ((1 - z[3])**2) + somme3
        somme2 <- ((0.5 - z[2])**2) + somme2
        somme1 <- ((1 - z[1])**2) + somme1


        # binsize = 1
        y <- discrete(x, binsize = 1)
        z_1 <- whittle(y, "Exponential", 1)$par
        s_11 <- ((1 - z_1[1])**2) + s_11
        s_21 <- ((0.5 - z_1[2])**2) + s_21
        s_31 <- ((1 - z_1[3])**2) + s_31


        # binsize = 2
        y <- discrete(x, binsize = 2)
        z_2 <- whittle(y, "Exponential", 2)$par
        s_12 <- ((1 - z_2[1])**2) + s_12
        s_22 <- ((0.5 - z_2[2])**2) + s_22
        s_32 <- ((1 - z_2[3])**2) + s_32

        if (end == 3200) {
            w3 <- c(w3, (z[3]))
            w2 <- c(w2, (z[2]))
            w1 <- c(w1, (z[1]))

            w_11 <- c(w_11, (z_1[1]))
            w_21 <- c(w_21, (z_1[2]))
            w_31 <- c(w_31, (z_1[3]))

            w_12 <- c(w_12, (z_2[1]))
            w_22 <- c(w_22, (z_2[2]))
            w_32 <- c(w_32, (z_2[3]))

            me1 <- c(me1, (mle[1]))
            me2 <- c(me2, (mle[2]))
            me3 <- c(me3, (mle[3]))
        }


        if (end == 100) {
            w3b <- c(w3b, (z[3]))
            w2b <- c(w2b, (z[2]))
            w1b <- c(w1b, (z[1]))

            w_11b <- c(w_11b, (z_1[1]))
            w_21b <- c(w_21b, (z_1[2]))
            w_31b <- c(w_31b, (z_1[3]))

            w_12b <- c(w_12b, (z_2[1]))
            w_22b <- c(w_22b, (z_2[2]))
            w_32b <- c(w_32b, (z_2[3]))

            me1b <- c(me1b, (mle[1]))
            me2b <- c(me2b, (mle[2]))
            me3b <- c(me3b, (mle[3]))
        }
    }

    mle1 <- c(mle1, sm1 / iter)
    mle2 <- c(mle2, sm2 / iter)
    mle3 <- c(mle3, sm3 / iter)
    sm1 <- 0
    sm2 <- 0
    sm3 <- 0

    l1 <- c(l1, somme1 / iter)
    l2 <- c(l2, somme2 / iter)
    l3 <- c(l3, somme3 / iter)
    somme3 <- 0
    somme2 <- 0
    somme1 <- 0

    l_11 <- c(l_11, s_11 / iter)
    l_21 <- c(l_21, s_21 / iter)
    l_31 <- c(l_31, s_31 / iter)
    s_11 <- 0
    s_21 <- 0
    s_31 <- 0

    l_12 <- c(l_12, s_12 / iter)
    l_22 <- c(l_22, s_22 / iter)
    l_32 <- c(l_32, s_32 / iter)
    s_12 <- 0
    s_22 <- 0
    s_32 <- 0
}


# We export all the data in a csv file to make the plots in Python
my_data <- data.frame(l1, l2, l3, l_11, l_21, l_31, l_12, l_22, l_32, mle1, mle2, mle3)
boxplot <- data.frame(w1, w2, w3, w_11, w_21, w_31, w_12, w_22, w_32, me1, me2, me3, w1b, w2b, w3b, w_11b, w_21b, w_31b, w_12b, w_22b, w_32b, me1b, me2b, me3b)
write.csv(my_data, "C:\\Users\\etoma\\OneDrive\\Bureau\\Hawkes\\data_mse.csv", , row.names = FALSE)
write.csv(boxplot, "C:\\Users\\etoma\\OneDrive\\Bureau\\Hawkes\\data_boxplots.csv", , row.names = FALSE)




###




## The programs below were used to plot all the execution times of the section "Python vs R: execution time" of our Jupyter notebook.


# Generate Hawkes processes
temps_hawkes <- function(end) {
    for (j in seq(1, 100)) {
        hawkes(end, fun = 1, repr = .5, family = "exp", rate = 1)
    }
}

l <- c()
for (end in c(100, 200, 400, 800, 1600, 3200)) {
    l <- c(l, system.time(temps_hawkes(end))[1])
}

print(l)



# Generate the count series
temps_discrete <- function(hawkes) {
    for (j in seq(1, 100)) {
        discrete(hawkes, binsize = 2)
    }
}

l <- c()
for (end in c(100, 200, 400, 800, 1600, 3200)) {
    x <- hawkes(end, fun = 1, repr = .5, family = "exp", rate = 1)
    l <- c(l, system.time(temps_discrete(x))[1])
}

print(l)



# Generate the periodogram
temps_periodogram <- function(counts) {
    for (j in seq(1, 1000)) {
        n <- length(counts)
        dft <- fft(counts - mean(counts))
        I <- Mod(dft)^2 / n
        I <- I[-1]
    }
}

l <- c()
for (end in c(100, 200, 400, 800, 1600, 3200)) {
    x <- hawkes(end, fun = 1, repr = .5, family = "exp", rate = 1)
    y <- discrete(x, binsize = 2)
    l <- c(l, system.time(temps_periodogram(y))[1])
}

print(l)



# The Whittle estimation for exponential kernels
temps_whittle <- function(counts) {
    for (j in seq(1, 100)) {
        whittle(counts, "Exponential", binsize = 2)
    }
}

l <- c()
for (end in c(100, 200, 400, 800, 1600, 3200)) {
    x <- hawkes(end, fun = 1, repr = .5, family = "exp", rate = 1)
    y <- discrete(x, binsize = 2)
    l <- c(l, system.time(temps_whittle(y)))
}

print(l)



# The Whittle estimation for powerlaw kernels
temps_whittle <- function(counts) {
    for (j in seq(1, 10)) {
        whittle(counts, "powerlaw")
    }
}

l <- c()
for (end in c(100, 200, 400, 800, 1600, 3200)) {
    x <- hawkes(end, fun = 1, repr = .3, family = "powerlaw", shape = 3.5, scale = 1.0)
    y <- discrete(x, binsize = 1)
    l <- c(l, system.time(temps_whittle(y)))
}

print(l)


# MLE estimation for exponential kernels
temps_whittle <- function(x) {
    for (j in seq(1, 10)) {
        mle(x$p, "Exponential", x$end)
    }
}

l <- c()
for (end in c(100, 200, 400, 800, 1600, 3200)) {
    x <- hawkes(end, fun = 1, repr = .5, family = "exp", rate = 1)
    l <- c(l, system.time(temps_whittle(x)))
}

print(l)
