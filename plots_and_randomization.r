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




########################################## SIMULATIONS for RANDOMIZATION ##########################################


randomization <- function(y, binsize) {
  time = numeric()
  t = 0
  L = length(y)
  for (i in 1:L) {
    if (y[i] == 0) {
        # dont do anything
    } else {
      for (j in 1:y[i]) {
        time = c(time, t+j*binsize/(y[i]+1))
      }
    }
    t = t + binsize
  }
  return(time)
}

N = 500
end = 1000
#bins is a list of binsizes from 0.125 to 2 in 0.125 increments
bins = cbind(0.125, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)

# results for whittle 
muw = numeric(length(bins))
alphaw = numeric(length(bins))
betaw = numeric(length(bins))

# results for mle_randomized
mum = numeric(length(bins))
alpham = numeric(length(bins))
betam = numeric(length(bins))

# results for mle not randomized
mumle = numeric(length(bins))
alphamle = numeric(length(bins))
betamle = numeric(length(bins))

for (i in 1:length(bins)) {
  sw = c(0 , 0 , 0)
  sm = c(0 , 0 , 0)
  mlem = c(0 , 0 , 0)

  for (k in 1:N) {
    x <- hawkes(end, fun = 1, repr = .5, family = "exp", rate = 1)
    mle_original <- mle(x$p, "Exponential", x$end)$par

    y = discrete(x, binsize = bins[i])
    whi = whittle(y, "exp", binsize = bins[i])
    whi = whi$par

    x_new = randomization(y, bins[i])
    ml <- mle(x_new, "Exponential", x$end)$par
    # compute the squared error (absolute value of the difference between the two values squared)
    for (j in 1:3) {
      a = c(1,0.5,1)
      sw[j] = sw[j] + (a[j]- whi[j])^2
      sm[j] = sm[j] + (a[j]- ml[j])^2
      mlem[j] = mlem[j] + (a[j]- mle_original[j])^2
    }

  }

  for (j in 1:3) {
    sw[j] = sw[j] / N
    sm[j] = sm[j] / N
    mlem[j] = mlem[j] / N
  }
  # store the results
  muw[i] = sw[1]
  alphaw[i] = sw[2]
  betaw[i] = sw[3]

  mum[i] = sm[1]
  alpham[i] = sm[2]
  betam[i] = sm[3]

  mumle[i] = mlem[1]
  alphamle[i] = mlem[2]
  betamle[i] = mlem[3]

  print(paste("binsize: ", bins[i]))
}

# store the results
write.table(cbind(muw, mum, mumle), file = "results_randomization_mu.csv", row.names = FALSE, col.names = FALSE)
write.table(cbind(alphaw, alpham, alphamle), file = "results_randomization_alpha.csv", row.names = FALSE, col.names = FALSE)
write.table(cbind(betaw, betam, betamle), file = "results_randomization_beta.csv", row.names = FALSE, col.names = FALSE)

""" PYTHON CODE FOR RENDERING
import pandas as pd
import matplotlib.pyplot as plt

binlist = [0.125, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
# load the data from results_randomization_alpha.csv, results_randomization_beta.csv and results_randomization_mu.csv which are parsed using space as delimiter
# N.B. each row corresponds to simulation results with binsize store in binlist
# N.B. each column corresponds to whittle, mle_randomized, mle_not_randomized
data_mu = pd.read_csv('results_randomization_mu.csv',  header=None, names=['whittle', 'mle_randomized', 'mle_not_randomized'], sep=' ')
data_alpha = pd.read_csv('results_randomization_alpha.csv', header=None, names=['whittle', 'mle_randomized', 'mle_not_randomized'] , sep=' ')
data_beta = pd.read_csv('results_randomization_beta.csv', header=None, names=['whittle', 'mle_randomized', 'mle_not_randomized'] , sep=' ')

# 3 subplots horizontally
fig, ax = plt.subplots(3,1, figsize=(10,20))
# plot results against binsize for mu 
ax[0].plot(binlist, data_mu['whittle'],"o-", label='whittle',  color='midnightblue')
ax[0].plot(binlist, data_mu['mle_randomized'], label='mle_randomized', marker='.', ls='-', color='darkorange')
ax[0].plot(binlist, data_mu['mle_not_randomized'], label='mle_not_randomized', marker='.', ls='-', color='darkgreen')
ax[0].set_xlabel('Binsize')
ax[0].set_ylabel('MSE')
ax[0].set_title('MSE for mu', fontsize=15)

# plot results against binsize for alpha
ax[1].plot(binlist, data_alpha['whittle'], "o-", label='whittle', color='midnightblue')
ax[1].plot(binlist, data_alpha['mle_randomized'], label='mle_randomized', marker='.', ls='-', color='darkorange')
ax[1].plot(binlist, data_alpha['mle_not_randomized'], label='mle_not_randomized', marker='.', ls='-', color='darkgreen')
ax[1].set_xlabel('Binsize')
ax[1].set_ylabel('MSE')
ax[1].set_title('MSE for alpha', fontsize=15)

# plot results against binsize for beta
ax[2].plot(binlist, data_beta['whittle'], "o-", label='whittle', marker='.', ls='-', color='midnightblue')
ax[2].plot(binlist, data_beta['mle_randomized'], label='mle_randomized', marker='.', ls='-', color='darkorange')
ax[2].plot(binlist, data_beta['mle_not_randomized'], label='mle_not_randomized', marker='.', ls='-', color='darkgreen')
ax[2].set_xlabel('Binsize')
ax[2].set_ylabel('MSE')
ax[2].set_title('MSE for beta', fontsize=15)

#set y log scale
ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[2].set_yscale('log')

fig.suptitle('MSE for against binsize for T = 1000 over 50 simulations', fontsize=16)
# show legend 
ax[0].legend(loc='upper left')
ax[1].legend(loc='upper left')
ax[2].legend(loc='upper left')

plt.show()
"""



