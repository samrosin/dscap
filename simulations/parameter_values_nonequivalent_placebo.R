# each theta_at has values for a=0 and a=t, for t=1,2,3,4,5

theta_a1 <- list(matrix(c(-6, 1.5, .1), ncol = 1), # theta_01 vector (theta_010, theta_011, theta_012)
                matrix(c(-6.69, 1.5, .1), ncol = 1)) # theta_11 vector (theta_110, theta_111, theta_112)

theta_a2 <- list(matrix(c(-5.65, 1.5, .1), ncol = 1), # theta_02 vector
                matrix(c(-6.42, 1.5, .081), ncol = 1)) # theta_22 vecctor

theta_a3 <- list(matrix(c(-5.8, 1.5, .1), ncol = 1), # theta_03
                matrix(c(-9.7, 1, .088), ncol = 1)) # theta_33

theta_a4 <- list(matrix(c(-6.2, 1.5, .1), ncol = 1), # theta_04
                matrix(c(-8.1, 1.15, .1), ncol = 1)) # theta_44

theta_a5 <- list(matrix(c(-6.35, 1.5, .1), ncol = 1), # theta_05 
                matrix(c(-9.52, 1.5, .095), ncol = 1)) # theta_55

theta <- list(theta_a1, theta_a2, theta_a3, theta_a4, theta_a5)

gamma_a1 <- list(matrix(c(-5.5, 1.5, .1), ncol = 1),
                matrix(c(-2.971, 1, .08), ncol = 1))
gamma_a2 <- list(matrix(c(-5.425, 1.5, .1), ncol = 1),
                matrix(c(-3.775, 1.5, .081), ncol = 1))
gamma_a3 <- list(matrix(c(-5.45, 1.5, .1), ncol = 1),
                 matrix(c(-0.718, 1, .09), ncol = 1))
gamma_a4 <- list(matrix(c(-5.55, 1.5, .1), ncol = 1),
                 matrix(c(-2.818, 1.15, .1), ncol = 1))
gamma_a5 <- list(matrix(c(-5.575, 1.5, .1), ncol = 1),
                 matrix(c(-1.668, 1.5, .1), ncol = 1))
gamma = list(gamma_a1, gamma_a2, gamma_a3, gamma_a4, gamma_a5)
